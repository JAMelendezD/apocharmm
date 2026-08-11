// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: Samarjeet Prasad, James E. Gonzales II
//
// ENDLICENSE

#include "CharmmContext.h"

#include "ApoCharmmError.h"
#include "Constants.h"
#include "PBC.h"
#include "cuda_utils.h"
#include "gpu_utils.h"

#include <cmath>
#include <cuda.h>
#include <cuda_runtime.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>

static bool hasValidBoxDimensions(const std::vector<double> &boxDimensions) {
  if (boxDimensions.size() != 3)
    return false;
  return ((boxDimensions[0] > 0.0) && (boxDimensions[1] > 0.0) &&
          (boxDimensions[2] > 0.0));
}

CharmmContext::CharmmContext(void)
    : m_RandomSeed(0), m_Psf(nullptr), m_Prm(nullptr),
      m_BoxDimensions({-9999.9999, -9999.9999, -9999.9999}),
      m_HasBoxDimensions(false), m_ForceManager(nullptr), m_NumAtoms(-1),
      m_NumDegreesOfFreedom(-1), m_Pbc(PBC::P1), m_HasPbc(false),
      m_CoordinatesChargesSP(), m_CoordinatesChargesDP(),
      m_HasCoordinates(false), m_VelocitiesInverseMasses(), m_KineticEnergy(1),
      m_Pressure(9), m_VirialKineticEnergyTensor(9),
      m_EnergyTableEvaluationCount(0), m_PreviousPrintedPotentialEnergy(0.0),
      m_HasPreviousPrintedPotentialEnergy(false), m_Temperature(0.0f),
      m_UsingHolonomicConstraints(false) {}

CharmmContext::CharmmContext(std::shared_ptr<CharmmPSF> psf,
                             std::shared_ptr<CharmmParameters> prm)
    : CharmmContext() {
  APOCHARMM_REQUIRE(psf != nullptr, ApoCharmmErrorCode::InvalidArgument,
                    "CharmmPSF must not be null");

  APOCHARMM_REQUIRE(prm != nullptr, ApoCharmmErrorCode::InvalidArgument,
                    "CharmmParameters must not be null");

  std::vector<int> devices = {0, 1, 2, 3};
  start_gpu(1, 1, 0, devices);

  this->setPsf(psf);
  this->setPrm(prm);
  this->setForceManager(std::make_shared<ForceManager>());

  this->useHolonomicConstraints(true);

  std::random_device rd{};
  m_RandomSeed = rd();
}

CharmmContext::CharmmContext(std::shared_ptr<ForceManager> fm)
    : CharmmContext() {
  std::vector<int> devices = {0, 1, 2, 3};
  start_gpu(1, 1, 0, devices);

  this->setForceManager(fm);
  this->requireForceManager();

  this->syncStateFromForceManager();

  if (!m_ForceManager->isInitialized() &&
      this->hasCompleteForceManagerState()) {
    m_ForceManager->initialize();
    this->syncStateFromForceManager();
    this->useHolonomicConstraints(true);
  }

  std::random_device rd{};
  m_RandomSeed = rd();
}

CharmmContext::CharmmContext(const CharmmContext &other)
    : std::enable_shared_from_this<CharmmContext>(other) {
  m_Psf = other.m_Psf;
  m_Prm = other.m_Prm;
  m_BoxDimensions = other.m_BoxDimensions;
  m_HasBoxDimensions = other.m_HasBoxDimensions;
  m_Pbc = other.m_Pbc;
  m_HasPbc = other.m_HasPbc;

  m_NumAtoms = other.m_NumAtoms;
  m_CoordinatesChargesSP = other.m_CoordinatesChargesSP;
  m_CoordinatesChargesDP = other.m_CoordinatesChargesDP;
  m_HasCoordinates = other.m_HasCoordinates;
  m_VelocitiesInverseMasses = other.m_VelocitiesInverseMasses;
  m_KineticEnergy = other.m_KineticEnergy;
  m_Pressure = other.m_Pressure;
  m_EnergyTableEvaluationCount = 0;
  m_PreviousPrintedPotentialEnergy = 0.0;
  m_HasPreviousPrintedPotentialEnergy = false;
  m_Temperature = other.m_Temperature;
  m_NumDegreesOfFreedom = other.m_NumDegreesOfFreedom;
  m_UsingHolonomicConstraints = other.m_UsingHolonomicConstraints;
  m_RandomSeed = other.m_RandomSeed;
}

void CharmmContext::setPrm(std::shared_ptr<CharmmParameters> prm) {
  APOCHARMM_REQUIRE(prm != nullptr, ApoCharmmErrorCode::InvalidArgument,
                    "CharmmParameters must not be null");

  m_Prm = prm;

  if (m_ForceManager != nullptr)
    m_ForceManager->setPrm(m_Prm);

  this->finalizeSetupIfReady();

  return;
}

void CharmmContext::setPsf(std::shared_ptr<CharmmPSF> psf) {
  APOCHARMM_REQUIRE(psf != nullptr, ApoCharmmErrorCode::InvalidArgument,
                    "CharmmPSF must not be null");

  APOCHARMM_REQUIRE(
      (m_NumAtoms <= 0) || (m_NumAtoms == psf->getNumAtoms()),
      ApoCharmmErrorCode::InvalidArgument,
      "Atom count mismatch between CharmmContext and CharmmPSF; expected " +
          std::to_string(m_NumAtoms) + ", observed " +
          std::to_string(psf->getNumAtoms()));

  m_Psf = psf;

  this->setNumAtoms(psf->getNumAtoms());
  this->setCharges(psf->getCharges());
  this->setMasses(psf->getMasses());

  if (m_ForceManager != nullptr)
    m_ForceManager->setPsf(m_Psf);

  this->finalizeSetupIfReady();

  return;
}

void CharmmContext::setForceManager(std::shared_ptr<ForceManager> fm) {
  APOCHARMM_REQUIRE(fm != nullptr, ApoCharmmErrorCode::InvalidArgument,
                    "ForceManager must not be null");

  m_ForceManager = fm;

  this->syncForceManagerFromState();

  this->finalizeSetupIfReady();

  return;
}

void CharmmContext::setNumAtoms(const int numAtoms) {
  APOCHARMM_REQUIRE(numAtoms > 0, ApoCharmmErrorCode::InvalidArgument,
                    "Atom count must be positive; observed " +
                        std::to_string(numAtoms));

  const bool coordinatesSizeChanged =
      (m_CoordinatesChargesSP.size() != static_cast<std::size_t>(numAtoms)) ||
      (m_CoordinatesChargesDP.size() != static_cast<std::size_t>(numAtoms));

  m_NumAtoms = numAtoms;

  // If array sizes change, reallocate and set to initial values

  if (coordinatesSizeChanged)
    m_HasCoordinates = false;

  if (m_CoordinatesChargesSP.size() != static_cast<std::size_t>(numAtoms)) {
    m_CoordinatesChargesSP.resize(numAtoms);
    m_CoordinatesChargesSP.setToValue(make_float4(0.0f, 0.0f, 0.0f, 0.0f));
  }

  if (m_CoordinatesChargesDP.size() != static_cast<std::size_t>(numAtoms)) {
    m_CoordinatesChargesDP.resize(numAtoms);
    m_CoordinatesChargesDP.setToValue(make_double4(0.0, 0.0, 0.0, 0.0));
  }

  if (m_VelocitiesInverseMasses.size() != static_cast<std::size_t>(numAtoms)) {
    m_VelocitiesInverseMasses.resize(numAtoms);
    m_VelocitiesInverseMasses.setToValue(make_double4(0.0, 0.0, 0.0, 0.0));
  }

  return;
}

void CharmmContext::setCoordinatesCharges(
    const std::vector<double4> &coordinatesCharges) {
  APOCHARMM_REQUIRE(
      m_NumAtoms > 0, ApoCharmmErrorCode::NotInitialized,
      "Atom count must be set before assigning coordinates and charges");

  APOCHARMM_REQUIRE(coordinatesCharges.size() ==
                        static_cast<std::size_t>(m_NumAtoms),
                    ApoCharmmErrorCode::InvalidArgument,
                    "Coordinate and charge count mismatch; expected " +
                        std::to_string(m_NumAtoms) + ", observed " +
                        std::to_string(coordinatesCharges.size()));

  if ((m_ForceManager != nullptr) && m_ForceManager->isInitialized() &&
      !m_ForceManager->isComposite()) {
    APOCHARMM_REQUIRE(
        coordinatesCharges.size() ==
            static_cast<std::size_t>(m_ForceManager->getNumAtoms()),
        ApoCharmmErrorCode::InvalidArgument,
        "Coordinate and charge count mismatch with ForceManager; expected " +
            std::to_string(m_ForceManager->getNumAtoms()) + ", observed " +
            std::to_string(coordinatesCharges.size()));
  }

  for (int i = 0; i < m_NumAtoms; i++) {
    const double4 &xyzq = coordinatesCharges[i];
    m_CoordinatesChargesSP[i] =
        make_float4(static_cast<float>(xyzq.x), static_cast<float>(xyzq.y),
                    static_cast<float>(xyzq.z), static_cast<float>(xyzq.w));
    m_CoordinatesChargesDP[i] = xyzq;
  }

  m_CoordinatesChargesSP.transferToDevice();
  m_CoordinatesChargesDP.transferToDevice();

  // JEG260707: This is done to prevent building the neighbor list twice.
  const bool forceManagerWasInitialized =
      ((m_ForceManager != nullptr) && m_ForceManager->isInitialized());

  m_HasCoordinates = true;

  this->finalizeSetupIfReady();

  if (forceManagerWasInitialized)
    this->resetNeighborList();

  return;
}

void CharmmContext::setCoordinatesCharges(
    const std::vector<std::vector<double>> &coordinatesCharges) {
  std::vector<double4> coordsCharges;
  coordsCharges.reserve(coordinatesCharges.size());

  for (const std::vector<double> &xyzq : coordinatesCharges) {
    APOCHARMM_REQUIRE(
        xyzq.size() == 4, ApoCharmmErrorCode::InvalidArgument,
        "Each coordinate and charge entry must contain exactly 4 values");

    coordsCharges.push_back(make_double4(xyzq[0], xyzq[1], xyzq[2], xyzq[3]));
  }

  this->setCoordinatesCharges(coordsCharges);

  return;
}

void CharmmContext::setCoordinatesCharges(
    const std::vector<double> &coordinatesCharges) {
  APOCHARMM_REQUIRE(
      coordinatesCharges.size() % 4 == 0, ApoCharmmErrorCode::InvalidArgument,
      "Coordinate and charge array length must be a multiple of 4; observed " +
          std::to_string(coordinatesCharges.size()));

  const std::size_t natom = coordinatesCharges.size() / 4;
  std::vector<double4> coordsCharges(natom);

  for (std::size_t i = 0; i < natom; i++) {
    coordsCharges[i] = make_double4(
        coordinatesCharges[i * 4 + 0], coordinatesCharges[i * 4 + 1],
        coordinatesCharges[i * 4 + 2], coordinatesCharges[i * 4 + 3]);
  }

  this->setCoordinatesCharges(coordsCharges);

  return;
}

void CharmmContext::setCoordinates(const std::vector<double3> &coordinates) {
  APOCHARMM_REQUIRE(m_NumAtoms > 0, ApoCharmmErrorCode::NotInitialized,
                    "Atom count must be set before assigning coordinates");

  APOCHARMM_REQUIRE(coordinates.size() == static_cast<std::size_t>(m_NumAtoms),
                    ApoCharmmErrorCode::InvalidArgument,
                    "Coordinate count mismatch; expected " +
                        std::to_string(m_NumAtoms) + ", observed " +
                        std::to_string(coordinates.size()));

  for (int i = 0; i < m_NumAtoms; i++) {
    m_CoordinatesChargesSP[i].x = static_cast<float>(coordinates[i].x);
    m_CoordinatesChargesSP[i].y = static_cast<float>(coordinates[i].y);
    m_CoordinatesChargesSP[i].z = static_cast<float>(coordinates[i].z);

    m_CoordinatesChargesDP[i].x = coordinates[i].x;
    m_CoordinatesChargesDP[i].y = coordinates[i].y;
    m_CoordinatesChargesDP[i].z = coordinates[i].z;
  }

  m_CoordinatesChargesSP.transferToDevice();
  m_CoordinatesChargesDP.transferToDevice();

  // JEG260707: This is done to prevent building the neighbor list twice.
  const bool forceManagerWasInitialized =
      ((m_ForceManager != nullptr) && m_ForceManager->isInitialized());

  m_HasCoordinates = true;

  this->finalizeSetupIfReady();

  if (forceManagerWasInitialized)
    this->resetNeighborList();

  return;
}

void CharmmContext::setCoordinates(
    const std::vector<std::vector<double>> &coordinates) {
  std::vector<double3> coords;
  coords.reserve(coordinates.size());

  for (const std::vector<double> &xyz : coordinates) {
    APOCHARMM_REQUIRE(xyz.size() == 3, ApoCharmmErrorCode::InvalidArgument,
                      "Each coordinate entry must contain exactly 3 values");

    coords.push_back(make_double3(xyz[0], xyz[1], xyz[2]));
  }

  this->setCoordinates(coords);

  return;
}

void CharmmContext::setCoordinates(const std::vector<double> &coordinates) {
  APOCHARMM_REQUIRE(
      coordinates.size() % 3 == 0, ApoCharmmErrorCode::InvalidArgument,
      "Coordinate array length must be a multiple of 3; observed " +
          std::to_string(coordinates.size()));

  const std::size_t natom = coordinates.size() / 3;
  std::vector<double3> coords(natom);

  for (std::size_t i = 0; i < natom; i++) {
    coords[i] = make_double3(coordinates[i * 3 + 0], coordinates[i * 3 + 1],
                             coordinates[i * 3 + 2]);
  }

  this->setCoordinates(coords);

  return;
}

void CharmmContext::setCoordinates(const std::shared_ptr<Coordinates> crd) {
  APOCHARMM_REQUIRE(crd != nullptr, ApoCharmmErrorCode::InvalidArgument,
                    "Coordinates must not be null");

  this->setCoordinates(crd->getCoordinatesDP());

  return;
}

void CharmmContext::setCharges(const std::vector<double> &charges) {
  APOCHARMM_REQUIRE(m_NumAtoms > 0, ApoCharmmErrorCode::NotInitialized,
                    "Atom count must be set before assigning charges");

  APOCHARMM_REQUIRE(charges.size() == static_cast<std::size_t>(m_NumAtoms),
                    ApoCharmmErrorCode::InvalidArgument,
                    "Charge count mismatch; expected " +
                        std::to_string(m_NumAtoms) + ", observed " +
                        std::to_string(charges.size()));

  for (int i = 0; i < m_NumAtoms; i++) {
    m_CoordinatesChargesSP[i].w = static_cast<float>(charges[i]);
    m_CoordinatesChargesDP[i].w = charges[i];
  }

  m_CoordinatesChargesSP.transferToDevice();
  m_CoordinatesChargesDP.transferToDevice();

  return;
}

void CharmmContext::setVelocitiesInverseMasses(
    const std::vector<double4> &velocitiesInverseMasses) {
  APOCHARMM_REQUIRE(
      m_NumAtoms > 0, ApoCharmmErrorCode::NotInitialized,
      "Atom count must be set before assigning velocities and inverse masses");

  APOCHARMM_REQUIRE(velocitiesInverseMasses.size() ==
                        static_cast<std::size_t>(m_NumAtoms),
                    ApoCharmmErrorCode::InvalidArgument,
                    "Velocity and inverse-mass count mismatch; expected " +
                        std::to_string(m_NumAtoms) + ", observed " +
                        std::to_string(velocitiesInverseMasses.size()));

  for (int i = 0; i < m_NumAtoms; i++)
    m_VelocitiesInverseMasses[i] = velocitiesInverseMasses[i];

  m_VelocitiesInverseMasses.transferToDevice();

  return;
}

void CharmmContext::setVelocitiesInverseMasses(
    const std::vector<std::vector<double>> &velocitiesInverseMasses) {
  std::vector<double4> velMass;
  velMass.reserve(velocitiesInverseMasses.size());

  for (const std::vector<double> &xyzm : velocitiesInverseMasses) {
    APOCHARMM_REQUIRE(
        xyzm.size() == 4, ApoCharmmErrorCode::InvalidArgument,
        "Each velocity and inverse-mass entry must contain exactly 4 values");

    velMass.push_back(make_double4(xyzm[0], xyzm[1], xyzm[2], xyzm[3]));
  }

  this->setVelocitiesInverseMasses(velMass);

  return;
}

void CharmmContext::setVelocitiesInverseMasses(
    const std::vector<double> &velocitiesInverseMasses) {
  APOCHARMM_REQUIRE(velocitiesInverseMasses.size() % 4 == 0,
                    ApoCharmmErrorCode::InvalidArgument,
                    "Velocity and inverse-mass array length must be a multiple "
                    "of 4; observed " +
                        std::to_string(velocitiesInverseMasses.size()));

  const std::size_t natom = velocitiesInverseMasses.size() / 4;
  std::vector<double4> velMass(natom);

  for (std::size_t i = 0; i < natom; i++) {
    velMass[i] = make_double4(
        velocitiesInverseMasses[i * 4 + 0], velocitiesInverseMasses[i * 4 + 1],
        velocitiesInverseMasses[i * 4 + 2], velocitiesInverseMasses[i * 4 + 3]);
  }

  this->setVelocitiesInverseMasses(velMass);

  return;
}

void CharmmContext::setVelocities(const std::vector<double3> &velocities) {
  APOCHARMM_REQUIRE(m_NumAtoms > 0, ApoCharmmErrorCode::NotInitialized,
                    "Atom count must be set before assigning velocities");

  APOCHARMM_REQUIRE(velocities.size() == static_cast<std::size_t>(m_NumAtoms),
                    ApoCharmmErrorCode::InvalidArgument,
                    "Velocity count mismatch; expected " +
                        std::to_string(m_NumAtoms) + ", observed " +
                        std::to_string(velocities.size()));

  for (int i = 0; i < m_NumAtoms; i++) {
    m_VelocitiesInverseMasses[i].x = velocities[i].x;
    m_VelocitiesInverseMasses[i].y = velocities[i].y;
    m_VelocitiesInverseMasses[i].z = velocities[i].z;
  }

  m_VelocitiesInverseMasses.transferToDevice();

  return;
}

void CharmmContext::setVelocities(
    const std::vector<std::vector<double>> &velocities) {
  std::vector<double3> vels;
  vels.reserve(velocities.size());

  for (const std::vector<double> &xyz : velocities) {
    APOCHARMM_REQUIRE(xyz.size() == 3, ApoCharmmErrorCode::InvalidArgument,
                      "Each velocity entry must contain exactly 3 values");

    vels.push_back(make_double3(xyz[0], xyz[1], xyz[2]));
  }

  this->setVelocities(vels);

  return;
}

void CharmmContext::setVelocities(const std::vector<double> &velocities) {
  APOCHARMM_REQUIRE(velocities.size() % 3 == 0,
                    ApoCharmmErrorCode::InvalidArgument,
                    "Velocity array length must be a multiple of 3; observed " +
                        std::to_string(velocities.size()));

  const std::size_t natom = velocities.size() / 3;
  std::vector<double3> vels(natom);

  for (std::size_t i = 0; i < natom; i++) {
    vels[i] = make_double3(velocities[i * 3 + 0], velocities[i * 3 + 1],
                           velocities[i * 3 + 2]);
  }

  this->setVelocities(vels);

  return;
}

void CharmmContext::setVelocitiesFromCHARMMVelocityFile(
    const std::string &fileName) {
  APOCHARMM_REQUIRE(m_NumAtoms > 0, ApoCharmmErrorCode::NotInitialized,
                    "Atom count must be set before loading velocities from a "
                    "CHARMM velocity file");
  APOCHARMM_REQUIRE(!fileName.empty(), ApoCharmmErrorCode::InvalidArgument,
                    "CHARMM velocity file name must not be empty");

  std::ifstream fin(fileName);
  APOCHARMM_REQUIRE(fin.is_open(), ApoCharmmErrorCode::Runtime,
                    "Could not open CHARMM velocity file \"" + fileName + "\"");

  std::string line;
  bool foundAtomCount = false;
  while (std::getline(fin, line)) {
    if (!line.empty() && line[0] != '*') {
      foundAtomCount = true;
      break;
    }
  }

  APOCHARMM_REQUIRE(foundAtomCount, ApoCharmmErrorCode::Runtime,
                    "Could not read atom count from CHARMM velocity file \"" +
                        fileName + "\"");

  int fileAtomCount = 0;
  std::stringstream atomCountStream(line);
  APOCHARMM_REQUIRE(static_cast<bool>(atomCountStream >> fileAtomCount),
                    ApoCharmmErrorCode::Runtime,
                    "Could not parse atom count from CHARMM velocity file \"" +
                        fileName + "\"");

  APOCHARMM_REQUIRE(fileAtomCount == m_NumAtoms,
                    ApoCharmmErrorCode::InvalidArgument,
                    "CHARMM velocity atom count mismatch; expected " +
                        std::to_string(m_NumAtoms) + ", observed " +
                        std::to_string(fileAtomCount));

  for (int i = 0; i < m_NumAtoms; i++) {
    APOCHARMM_REQUIRE(
        static_cast<bool>(std::getline(fin, line)), ApoCharmmErrorCode::Runtime,
        "Could not read velocity record " + std::to_string(i + 1) +
            " from CHARMM velocity file \"" + fileName + "\"");

    int atomId = 0;
    int residueId = 0;
    int residueIdInSegment = 0;
    std::string residueName = "";
    std::string atomName = "";
    std::string segmentName = "";
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double bFactor = 0.0;

    std::stringstream recordStream(line);
    APOCHARMM_REQUIRE(
        static_cast<bool>(recordStream >> atomId >> residueId >> residueName >>
                          atomName >> x >> y >> z >> segmentName >>
                          residueIdInSegment >> bFactor),
        ApoCharmmErrorCode::Runtime,
        "Could not parse velocity record " + std::to_string(i + 1) +
            " from CHARMM velocity file \"" + fileName + "\"");

    m_VelocitiesInverseMasses[i].x = x;
    m_VelocitiesInverseMasses[i].y = y;
    m_VelocitiesInverseMasses[i].z = z;
  }

  m_VelocitiesInverseMasses.transferToDevice();

  return;
}

void CharmmContext::setMasses(const std::vector<double> &masses) {
  APOCHARMM_REQUIRE(m_NumAtoms > 0, ApoCharmmErrorCode::NotInitialized,
                    "Atom count must be set before assigning masses");

  APOCHARMM_REQUIRE(masses.size() == static_cast<std::size_t>(m_NumAtoms),
                    ApoCharmmErrorCode::InvalidArgument,
                    "Mass count mismatch; expected " +
                        std::to_string(m_NumAtoms) + ", observed " +
                        std::to_string(masses.size()));

  for (int i = 0; i < m_NumAtoms; i++)
    m_VelocitiesInverseMasses[i].w = 1.0 / masses[i];

  m_VelocitiesInverseMasses.transferToDevice();

  return;
}

void CharmmContext::setTemperature(const double temperature) {
  APOCHARMM_REQUIRE(std::isfinite(temperature) && (temperature >= 0.0),
                    ApoCharmmErrorCode::InvalidArgument,
                    "Temperature must be finite and non-negative");

  m_Temperature = temperature;

  return;
}

void CharmmContext::setPeriodicBoundaryCondition(const PBC pbc) {
  m_Pbc = pbc;
  m_HasPbc = true;

  if (m_ForceManager != nullptr) {
    m_ForceManager->setPeriodicBoundaryCondition(m_Pbc);
    m_Pbc = m_ForceManager->getPeriodicBoundaryCondition();
  }

  // JEG260707: This is done to prevent building the neighbor list twice.
  const bool forceManagerWasInitialized =
      ((m_ForceManager != nullptr) && m_ForceManager->isInitialized());

  this->finalizeSetupIfReady();

  if (forceManagerWasInitialized && m_HasCoordinates)
    this->resetNeighborList();

  return;
}

void CharmmContext::setBoxDimensions(const std::vector<double> &boxDimensions) {
  APOCHARMM_REQUIRE(hasValidBoxDimensions(boxDimensions),
                    ApoCharmmErrorCode::InvalidArgument,
                    "Box dimensions must contain exactly 3 positive values");

  m_BoxDimensions = boxDimensions;
  m_HasBoxDimensions = true;

  if (m_ForceManager != nullptr) {
    m_ForceManager->setBoxDimensions(m_BoxDimensions);
    m_BoxDimensions = m_ForceManager->getBoxDimensions();
    m_HasBoxDimensions = hasValidBoxDimensions(m_BoxDimensions);
  }

  // JEG260707: This is done to prevent building the neighbor list twice.
  const bool forceManagerWasInitialized =
      ((m_ForceManager != nullptr) && m_ForceManager->isInitialized());

  this->finalizeSetupIfReady();

  if (forceManagerWasInitialized && m_HasCoordinates)
    this->resetNeighborList();

  return;
}

void CharmmContext::setRandomSeed(const std::uint64_t randomSeed) {
  m_RandomSeed = randomSeed;
  return;
}

void CharmmContext::useHolonomicConstraints(
    const bool usingHolonomicConstraints) {
  this->requirePsf();
  this->requireForceManager();

  m_UsingHolonomicConstraints = usingHolonomicConstraints;

  int ndegf = m_NumAtoms * 3;

  if (m_Pbc == PBC::P1)
    ndegf -= 3;
  else if (m_Pbc == PBC::P21)
    ndegf -= 1;

  if (m_UsingHolonomicConstraints) {
    ndegf -= this->getWaterMolecules().size() * 3;
    int numShakeConstraints = 0;
    const std::vector<int4> &shakeAtoms =
        m_ForceManager->getShakeAtoms().getHostArray();
    for (std::size_t i = 0; i < shakeAtoms.size(); i++) {
      numShakeConstraints++;
      if (shakeAtoms[i].z != -1)
        numShakeConstraints++;
      if (shakeAtoms[i].w != -1)
        numShakeConstraints++;
    }
    ndegf -= numShakeConstraints;
  }

  m_NumDegreesOfFreedom = ndegf;

  return;
}

void CharmmContext::setKappa(const float kappa) {
  this->requireForceManager();
  m_ForceManager->setKappa(kappa);
  return;
}

void CharmmContext::setCutoff(const float cutoff) {
  this->requireForceManager();
  m_ForceManager->setCutoff(cutoff);
  return;
}

void CharmmContext::setCtonnb(const float ctonnb) {
  this->requireForceManager();
  m_ForceManager->setCtonnb(ctonnb);
  return;
}

void CharmmContext::setCtofnb(const float ctofnb) {
  this->requireForceManager();
  m_ForceManager->setCtofnb(ctofnb);
  return;
}

void CharmmContext::setFFTGrid(const int nfftx, const int nffty,
                               const int nfftz) {
  this->requireForceManager();
  m_ForceManager->setFFTGrid(nfftx, nffty, nfftz);
  return;
}

void CharmmContext::setPmeSplineOrder(const int pmeSplineOrder) {
  this->requireForceManager();
  m_ForceManager->setPmeSplineOrder(pmeSplineOrder);
  return;
}

void CharmmContext::setVdwType(const int vdwType) {
  this->requireForceManager();
  m_ForceManager->setVdwType(vdwType);
  return;
}

int CharmmContext::getNumAtoms(void) const { return m_NumAtoms; }

int CharmmContext::getNumDegreesOfFreedom(void) const {
  return m_NumDegreesOfFreedom;
}

const CudaContainer<float4> &
CharmmContext::getCoordinatesChargesSP(void) const {
  return m_CoordinatesChargesSP;
}

const CudaContainer<double4> &
CharmmContext::getCoordinatesChargesDP(void) const {
  return m_CoordinatesChargesDP;
}

const CudaContainer<double4> &
CharmmContext::getVelocitiesInverseMasses(void) const {
  return m_VelocitiesInverseMasses;
}

double CharmmContext::getTemperature(void) const { return m_Temperature; }

PBC CharmmContext::getPeriodicBoundaryCondition(void) const { return m_Pbc; }

const std::vector<double> &CharmmContext::getBoxDimensions(void) const {
  return m_BoxDimensions;
}

std::uint64_t CharmmContext::getRandomSeed(void) const { return m_RandomSeed; }

bool CharmmContext::usingHolonomicConstraints(void) const {
  return m_UsingHolonomicConstraints;
}

double CharmmContext::getVolume(void) const {
  APOCHARMM_REQUIRE(
      m_HasBoxDimensions && hasValidBoxDimensions(m_BoxDimensions),
      ApoCharmmErrorCode::NotInitialized, "Box dimensions are not set");

  return (m_BoxDimensions[0] * m_BoxDimensions[1] * m_BoxDimensions[2]);
}

const CudaContainer<double> &CharmmContext::getPressure(void) const {
  return m_Pressure;
}

const std::vector<Bond> &CharmmContext::getBonds(void) const {
  this->requirePsf();
  return m_Psf->getBonds();
}

const CudaContainer<int4> &CharmmContext::getWaterMolecules(void) const {
  this->requirePsf();
  return m_Psf->getWaterMolecules();
}

const CudaContainer<int4> &CharmmContext::getShakeAtoms(void) const {
  this->requireInitializedForceManager();
  return m_ForceManager->getShakeAtoms();
}

const CudaContainer<float4> &CharmmContext::getShakeParams(void) const {
  this->requireInitializedForceManager();
  return m_ForceManager->getShakeParams();
}

float CharmmContext::getKappa(void) const {
  this->requireForceManager();
  return m_ForceManager->getKappa();
}

float CharmmContext::getCutoff(void) const {
  this->requireForceManager();
  return m_ForceManager->getCutoff();
}

float CharmmContext::getCtonnb(void) const {
  this->requireForceManager();
  return m_ForceManager->getCtonnb();
}

float CharmmContext::getCtofnb(void) const {
  this->requireForceManager();
  return m_ForceManager->getCtofnb();
}

std::vector<int> CharmmContext::getFFTGrid(void) const {
  this->requireForceManager();
  return m_ForceManager->getFFTGrid();
}

int CharmmContext::getPmeSplineOrder(void) const {
  this->requireForceManager();
  return m_ForceManager->getPmeSplineOrder();
}

int CharmmContext::getVdwType(void) const {
  this->requireForceManager();
  return m_ForceManager->getVdwType();
}

int CharmmContext::getForceStride(void) const {
  this->requireInitializedForceManager();
  return m_ForceManager->getForceStride();
}

CudaContainer<float4> &CharmmContext::getCoordinatesChargesSP(void) {
  return m_CoordinatesChargesSP;
}

CudaContainer<double4> &CharmmContext::getCoordinatesChargesDP(void) {
  return m_CoordinatesChargesDP;
}

CudaContainer<double4> &CharmmContext::getVelocitiesInverseMasses(void) {
  return m_VelocitiesInverseMasses;
}

std::vector<double> &CharmmContext::getBoxDimensions(void) {
  return m_BoxDimensions;
}

double CharmmContext::getKineticEnergy(void) {
  this->calculateKineticEnergy();
  m_KineticEnergy.transferFromDevice();
  return m_KineticEnergy.getHostArray()[0];
}

CudaContainer<double> &CharmmContext::getPressure(void) { return m_Pressure; }

std::shared_ptr<CharmmPSF> CharmmContext::getPsf(void) { return m_Psf; }

std::shared_ptr<CharmmParameters> CharmmContext::getPrm(void) { return m_Prm; }

std::shared_ptr<ForceManager> CharmmContext::getForceManager(void) {
  return m_ForceManager;
}

CudaContainer<double> &CharmmContext::getPotentialEnergy(void) {
  this->requireInitializedForceManager();
  return m_ForceManager->getPotentialEnergy();
}

float CharmmContext::getPotentialEnergies(void) {
  this->requireInitializedForceManager();
  return m_ForceManager->getPotentialEnergies();
}

std::shared_ptr<Force<double>> CharmmContext::getForces(void) {
  this->requireInitializedForceManager();
  return m_ForceManager->getForces();
}

CudaContainer<double> &CharmmContext::getVirial(void) {
  this->requireInitializedForceManager();
  return m_ForceManager->getVirial();
}

std::vector<Bond> &CharmmContext::getBonds(void) {
  this->requirePsf();
  return m_Psf->getBonds();
}

CudaContainer<int4> &CharmmContext::getWaterMolecules(void) {
  this->requirePsf();
  return m_Psf->getWaterMolecules();
}

CudaContainer<int4> &CharmmContext::getShakeAtoms(void) {
  this->requireInitializedForceManager();
  return m_ForceManager->getShakeAtoms();
}

CudaContainer<float4> &CharmmContext::getShakeParams(void) {
  this->requireInitializedForceManager();
  return m_ForceManager->getShakeParams();
}

void CharmmContext::assignVelocitiesAtTemperature(const double temperature) {
  APOCHARMM_REQUIRE(
      m_NumAtoms > 0, ApoCharmmErrorCode::NotInitialized,
      "Atom count must be set before assigning velocities at a temperature");

  this->setTemperature(temperature);
  const double kbt = charmm::constants::kBoltz * temperature;

  std::mt19937 rng{m_RandomSeed};

  for (int i = 0; i < m_NumAtoms; i++) {
    const double sd = std::sqrt(kbt * m_VelocitiesInverseMasses[i].w);
    std::normal_distribution<double> dist(0, sd);

    m_VelocitiesInverseMasses[i].x = dist(rng);
    m_VelocitiesInverseMasses[i].y = dist(rng);
    m_VelocitiesInverseMasses[i].z = dist(rng);
  }

  m_VelocitiesInverseMasses.transferToDevice();

  return;
}

double CharmmContext::computeTemperature(void) {
  APOCHARMM_REQUIRE(m_NumAtoms > 0 && m_VelocitiesInverseMasses.size() ==
                                          static_cast<std::size_t>(m_NumAtoms),
                    ApoCharmmErrorCode::NotInitialized,
                    "Atom count and velocity storage must be initialized "
                    "before computing temperature");

  APOCHARMM_REQUIRE(
      m_NumDegreesOfFreedom > 0, ApoCharmmErrorCode::NotInitialized,
      "Degrees of freedom must be initialized before computing temperature");

  const double kineticEnergy = this->getKineticEnergy();
  return kineticEnergy / (0.5 * static_cast<double>(m_NumDegreesOfFreedom) *
                          charmm::constants::kBoltz);
}

void CharmmContext::computePressure(void) {
  APOCHARMM_THROW(ApoCharmmErrorCode::NotImplemented,
                  "Pressure computation is not implemented");
}

__global__ static void ImageCenteringKernel(
    double4 *__restrict__ coordinatesChargesDP,
    float4 *__restrict__ coordinatesChargesSP,
    double4 *__restrict__ velocitiesInverseMasses, double *__restrict__ forces,
    const int forceStride, const int2 *__restrict__ groups, const int numGroups,
    const double boxX, const double boxY, const double boxZ, const PBC pbc) {
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  const int stride = gridDim.x * blockDim.x;

  for (int i = index; i < numGroups; i += stride) {
    const int2 group = groups[i];

    double gx = 0.0;
    double gy = 0.0;
    double gz = 0.0;
    for (int j = group.x; j <= group.y; j++) {
      gx += coordinatesChargesDP[j].x;
      gy += coordinatesChargesDP[j].y;
      gz += coordinatesChargesDP[j].z;
    }

    const double ng = static_cast<double>(group.y - group.x + 1);
    gx /= ng;
    gy /= ng;
    gz /= ng;

    if (gx < -0.5 * boxX) {
      for (int j = group.x; j <= group.y; j++) {
        coordinatesChargesDP[j].x += boxX;
        coordinatesChargesSP[j].x += static_cast<float>(boxX);
        if (pbc == PBC::P21) {
          coordinatesChargesDP[j].y = -coordinatesChargesDP[j].y;
          coordinatesChargesDP[j].z = -coordinatesChargesDP[j].z;
          coordinatesChargesSP[j].y = -coordinatesChargesSP[j].y;
          coordinatesChargesSP[j].z = -coordinatesChargesSP[j].z;
          velocitiesInverseMasses[j].y = -velocitiesInverseMasses[j].y;
          velocitiesInverseMasses[j].z = -velocitiesInverseMasses[j].z;
          forces[1 * forceStride + j] = -forces[1 * forceStride + j];
          forces[2 * forceStride + j] = -forces[2 * forceStride + j];
          gy = -gy;
          gz = -gz;
        }
      }
    }

    if (gx > 0.5 * boxX) {
      for (int j = group.x; j <= group.y; j++) {
        coordinatesChargesDP[j].x -= boxX;
        coordinatesChargesSP[j].x -= static_cast<float>(boxX);
        if (pbc == PBC::P21) {
          coordinatesChargesDP[j].y = -coordinatesChargesDP[j].y;
          coordinatesChargesDP[j].z = -coordinatesChargesDP[j].z;
          coordinatesChargesSP[j].y = -coordinatesChargesSP[j].y;
          coordinatesChargesSP[j].z = -coordinatesChargesSP[j].z;
          velocitiesInverseMasses[j].y = -velocitiesInverseMasses[j].y;
          velocitiesInverseMasses[j].z = -velocitiesInverseMasses[j].z;
          forces[1 * forceStride + j] = -forces[1 * forceStride + j];
          forces[2 * forceStride + j] = -forces[2 * forceStride + j];
          gy = -gy;
          gz = -gz;
        }
      }
    }

    if (gy < -0.5 * boxY) {
      for (int j = group.x; j <= group.y; j++) {
        coordinatesChargesDP[j].y += boxY;
        coordinatesChargesSP[j].y += static_cast<float>(boxY);
      }
    }

    if (gy > 0.5 * boxY) {
      for (int j = group.x; j <= group.y; j++) {
        coordinatesChargesDP[j].y -= boxY;
        coordinatesChargesSP[j].y -= static_cast<float>(boxY);
      }
    }

    if (gz < -0.5 * boxZ) {
      for (int j = group.x; j <= group.y; j++) {
        coordinatesChargesDP[j].z += boxZ;
        coordinatesChargesSP[j].z += static_cast<float>(boxZ);
      }
    }

    if (gz > 0.5 * boxZ) {
      for (int j = group.x; j <= group.y; j++) {
        coordinatesChargesDP[j].z -= boxZ;
        coordinatesChargesSP[j].z -= static_cast<float>(boxZ);
      }
    }
  }

  return;
}

void CharmmContext::imageCentering(void) {
  this->requireInitializedForceManager();
  this->requirePsf();
  APOCHARMM_REQUIRE(
      m_HasBoxDimensions && hasValidBoxDimensions(m_BoxDimensions),
      ApoCharmmErrorCode::NotInitialized, "Box dimensions are not set");

  double *forces = this->getForces()->xyz();
  const int forceStride = this->getForceStride();

  const CudaContainer<int2> &groups = m_Psf->getGroups();
  const int numGroups = groups.size();

  if (numGroups == 0)
    return;

  constexpr int numThreads = 128;
  const int numBlocks = (numGroups + numThreads - 1) / numThreads;

  cudaCheckLaunch(ImageCenteringKernel<<<numBlocks, numThreads>>>(
      m_CoordinatesChargesDP.getDeviceArray().data(),
      m_CoordinatesChargesSP.getDeviceArray().data(),
      m_VelocitiesInverseMasses.getDeviceArray().data(), forces, forceStride,
      groups.getDeviceArray().data(), numGroups, m_BoxDimensions[0],
      m_BoxDimensions[1], m_BoxDimensions[2], m_Pbc));

  cudaCheck(cudaDeviceSynchronize());

  return;
}

void CharmmContext::resetNeighborList(void) {
  this->requireInitializedForceManager();

  this->imageCentering();

  m_ForceManager->resetNeighborList(
      m_CoordinatesChargesSP.getDeviceArray().data());

  return;
}

__global__ static void CalculateKineticEnergyKernel(
    double *__restrict__ kineticEnergy,
    const double4 *__restrict__ velocitiesInverseMasses, const int numAtoms) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int stride = gridDim.x * blockDim.x;

  double ke = 0.0;
  for (int i = idx; i < numAtoms; i += stride) {
    ke += 0.5 *
          ((velocitiesInverseMasses[i].x * velocitiesInverseMasses[i].x) +
           (velocitiesInverseMasses[i].y * velocitiesInverseMasses[i].y) +
           (velocitiesInverseMasses[i].z * velocitiesInverseMasses[i].z)) /
          velocitiesInverseMasses[i].w;
  }

  ke = BlockReduceSum<double>(ke);

  if (threadIdx.x == 0)
    atomicAdd(kineticEnergy, ke);

  return;
}

void CharmmContext::calculateKineticEnergy(void) {
  APOCHARMM_REQUIRE(m_NumAtoms > 0 && m_VelocitiesInverseMasses.size() ==
                                          static_cast<std::size_t>(m_NumAtoms),
                    ApoCharmmErrorCode::NotInitialized,
                    "Atom count and velocity storage must be initialized "
                    "before calculating kinetic energy");

  cudaCheck(
      cudaMemset(static_cast<void *>(m_KineticEnergy.getDeviceArray().data()),
                 0, sizeof(double)));

  cudaCheckLaunch(CalculateKineticEnergyKernel<<<1, 1024>>>(
      m_KineticEnergy.getDeviceArray().data(),
      m_VelocitiesInverseMasses.getDeviceArray().data(), m_NumAtoms));

  cudaCheck(cudaDeviceSynchronize());

  return;
}

void CharmmContext::calculatePotentialEnergy(const bool reset,
                                             const bool print) {
  this->requireInitializedForceManager();

  APOCHARMM_REQUIRE(
      !print || !m_ForceManager->isComposite(),
      ApoCharmmErrorCode::NotImplemented,
      "Energy-table printing is not implemented for composite ForceManagers");

  m_ForceManager->calcForce(m_CoordinatesChargesSP.getDeviceArray().data(),
                            reset, true, true);

  if (print)
    this->printEnergyTable();

  return;
}

void CharmmContext::calculateForces(const bool reset, const bool calcEnergy,
                                    const bool calcVirial) {
  this->requireInitializedForceManager();

  m_ForceManager->calcForce(m_CoordinatesChargesSP.getDeviceArray().data(),
                            reset, calcEnergy, calcVirial);

  return;
}

void CharmmContext::linkBackForceManager(void) {
  this->requireForceManager();

  std::shared_ptr<CharmmContext> self = this->weak_from_this().lock();
  if (self == nullptr)
    return;

  m_ForceManager->setContext(self);

  return;
}

void CharmmContext::printEnergyTable(void) {
  const std::map<std::string, double> energyComponents =
      m_ForceManager->getEnergyComponents();

  const auto component =
      [&energyComponents](const std::string &name) -> double {
    return energyComponents.at(name);
  };

  CudaContainer<double> &potentialEnergy = m_ForceManager->getPotentialEnergy();
  potentialEnergy.transferToHost();

  const double totalEnergy = potentialEnergy[0];
  const double deltaEnergy =
      m_HasPreviousPrintedPotentialEnergy
          ? totalEnergy - m_PreviousPrintedPotentialEnergy
          : 0.0;

  const std::shared_ptr<Force<double>> forces = m_ForceManager->getForces();
  const std::size_t numAtoms = static_cast<std::size_t>(m_NumAtoms);

  std::vector<double> forceX(numAtoms);
  std::vector<double> forceY(numAtoms);
  std::vector<double> forceZ(numAtoms);

  forces->getXYZ(forceX.data(), forceY.data(), forceZ.data());

  double squaredGradientNorm = 0.0;
  for (std::size_t i = 0; i < numAtoms; i++) {
    squaredGradientNorm +=
        forceX[i] * forceX[i] + forceY[i] * forceY[i] + forceZ[i] * forceZ[i];
  }

  const double gradientRms =
      std::sqrt(squaredGradientNorm / static_cast<double>(numAtoms));

  std::ostringstream output;

  output << "ENER ENR:  Eval#     ENERgy      Delta-E         GRMS\n"
         << "ENER INTERN:          BONDs       ANGLes       UREY-b"
            "    DIHEdrals    IMPRopers\n"
         << "ENER CROSS:           CMAPs        PMF1D        PMF2D"
            "        PRIMO\n"
         << "ENER EXTERN:        VDWaals         ELEC       HBONds"
            "          ASP         USER\n"
         << "ENER EWALD:          EWKSum       EWSElf       EWEXcl"
            "       EWQCor       EWUTil\n"
         << " ----------       ---------    ---------    ---------"
            "    ---------    ---------\n";

  output << std::fixed << std::setprecision(5) << std::right;

  output << "ENER>" << std::setw(9) << m_EnergyTableEvaluationCount
         << std::setw(13) << totalEnergy << std::setw(13) << deltaEnergy
         << std::setw(13) << gradientRms << '\n';

  output << std::left << std::setw(14) << "ENER INTERN>" << std::right
         << std::setw(13) << component("bond") << std::setw(13)
         << component("angle") << std::setw(13) << component("ureyb")
         << std::setw(13) << component("dihe") << std::setw(13)
         << component("imdihe") << '\n';

  output << std::left << std::setw(14) << "ENER CROSS>" << std::right
         << std::setw(13) << component("cmap") << std::setw(13) << 0.0
         << std::setw(13) << 0.0 << std::setw(13) << 0.0 << '\n';

  output << std::left << std::setw(14) << "ENER EXTERN>" << std::right
         << std::setw(13) << component("vdw") << std::setw(13)
         << component("elec") << std::setw(13) << 0.0 << std::setw(13) << 0.0
         << std::setw(13) << component("user") << '\n';

  output << std::left << std::setw(14) << "ENER EWALD>" << std::right
         << std::setw(13) << component("ewks") << std::setw(13)
         << component("ewse") << std::setw(13) << component("ewex")
         << std::setw(13) << 0.0 << std::setw(13) << 0.0 << '\n';

  output << " ----------       ---------    ---------    ---------"
            "    ---------    ---------\n";

  std::cout << output.str() << std::flush;

  m_PreviousPrintedPotentialEnergy = totalEnergy;
  m_HasPreviousPrintedPotentialEnergy = true;
  m_EnergyTableEvaluationCount++;

  return;
}

void CharmmContext::syncStateFromForceManager(void) {
  this->requireForceManager();

  std::shared_ptr<CharmmPSF> psf = m_ForceManager->getPsf();
  if (psf != nullptr) {
    const bool needToInitializeAtoms = (m_NumAtoms != psf->getNumAtoms());

    m_Psf = psf;
    this->setNumAtoms(psf->getNumAtoms());

    if (needToInitializeAtoms) {
      this->setCharges(psf->getCharges());
      this->setMasses(psf->getMasses());
    }
  }

  std::shared_ptr<CharmmParameters> prm = m_ForceManager->getPrm();
  if (prm != nullptr)
    m_Prm = prm;

  m_BoxDimensions = m_ForceManager->getBoxDimensions();
  m_HasBoxDimensions = hasValidBoxDimensions(m_BoxDimensions);

  m_Pbc = m_ForceManager->getPeriodicBoundaryCondition();
  m_HasPbc = true;

  return;
}

void CharmmContext::syncForceManagerFromState(void) {
  this->requireForceManager();

  if (m_Psf != nullptr)
    m_ForceManager->setPsf(m_Psf);
  else {
    std::shared_ptr<CharmmPSF> psf = m_ForceManager->getPsf();
    if (psf != nullptr) {
      const bool needToInitializeAtoms = (m_NumAtoms != psf->getNumAtoms());

      m_Psf = psf;
      this->setNumAtoms(psf->getNumAtoms());

      if (needToInitializeAtoms) {
        this->setCharges(psf->getCharges());
        this->setMasses(psf->getMasses());
      }
    }
  }

  if (m_Prm != nullptr)
    m_ForceManager->setPrm(m_Prm);
  else {
    std::shared_ptr<CharmmParameters> prm = m_ForceManager->getPrm();
    if (prm != nullptr)
      m_Prm = prm;
  }

  if (m_HasBoxDimensions)
    m_ForceManager->setBoxDimensions(m_BoxDimensions);
  else {
    m_BoxDimensions = m_ForceManager->getBoxDimensions();
    m_HasBoxDimensions = hasValidBoxDimensions(m_BoxDimensions);
  }

  if (m_HasPbc)
    m_ForceManager->setPeriodicBoundaryCondition(m_Pbc);
  else {
    m_Pbc = m_ForceManager->getPeriodicBoundaryCondition();
    m_HasPbc = true;
  }

  return;
}

bool CharmmContext::hasCompleteForceManagerState(void) const {
  return ((m_ForceManager != nullptr) && (m_Psf != nullptr) &&
          (m_Prm != nullptr) && (m_HasBoxDimensions) &&
          (hasValidBoxDimensions(m_BoxDimensions)));
}

void CharmmContext::initializeForceManagerIfReady(void) {
  if (m_ForceManager == nullptr)
    return;

  if (m_ForceManager->isInitialized())
    return;

  if (!this->hasCompleteForceManagerState())
    return;

  m_ForceManager->initialize();
  this->syncStateFromForceManager();

  return;
}

void CharmmContext::requirePsf(void) const {
  APOCHARMM_REQUIRE(m_Psf != nullptr, ApoCharmmErrorCode::NotInitialized,
                    "CharmmPSF is not set");
  return;
}

void CharmmContext::requireForceManager(void) const {
  APOCHARMM_REQUIRE(m_ForceManager != nullptr,
                    ApoCharmmErrorCode::NotInitialized,
                    "ForceManager is not set");
  return;
}

void CharmmContext::requireInitializedForceManager(void) const {
  this->requireForceManager();

  APOCHARMM_REQUIRE(m_ForceManager->isInitialized(),
                    ApoCharmmErrorCode::NotInitialized,
                    "ForceManager is not initialized");

  return;
}

void CharmmContext::finalizeSetupIfReady(void) {
  if (m_ForceManager == nullptr)
    return;

  const bool forceManagerWasInitialized = m_ForceManager->isInitialized();

  if (!m_ForceManager->hasCharmmContext())
    this->linkBackForceManager();

  this->initializeForceManagerIfReady();

  const bool forceManagerWasInitializedHere =
      (!forceManagerWasInitialized && m_ForceManager->isInitialized());

  if (m_ForceManager->isInitialized())
    this->useHolonomicConstraints(m_UsingHolonomicConstraints);

  if (forceManagerWasInitializedHere && m_HasCoordinates)
    this->resetNeighborList();

  return;
}

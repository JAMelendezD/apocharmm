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

#include "Constants.h"
#include "PBC.h"
#include "cuda_utils.h"
#include "gpu_utils.h"

#include <cmath>
#include <cstdlib>
#include <cuda.h>
#include <cuda_runtime.h>
#include <fstream>
#include <iostream>
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
      m_Pressure(9), m_VirialKineticEnergyTensor(9), m_Temperature(0.0f),
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
  constexpr std::string_view functionName =
      "CharmmContext::CharmmContext(std::shared_ptr<ForceManager>)";

  std::vector<int> devices = {0, 1, 2, 3};
  start_gpu(1, 1, 0, devices);

  this->setForceManager(fm);
  this->requireForceManager(functionName);

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
  m_Temperature = other.m_Temperature;
  m_NumDegreesOfFreedom = other.m_NumDegreesOfFreedom;
  m_UsingHolonomicConstraints = other.m_UsingHolonomicConstraints;
  m_RandomSeed = other.m_RandomSeed;
}

void CharmmContext::setPrm(std::shared_ptr<CharmmParameters> prm) {
  constexpr std::string_view functionName =
      "void CharmmContext::setPrm(std::shared_ptr<CharmmParameters>)";

  if (prm == nullptr) {
    std::string msg(functionName);
    msg += ": CharmmParameters == nullptr";
    throw std::invalid_argument(msg);
  }

  m_Prm = prm;

  if (m_ForceManager != nullptr)
    m_ForceManager->setPrm(m_Prm);

  this->finalizeSetupIfReady();

  return;
}

void CharmmContext::setPsf(std::shared_ptr<CharmmPSF> psf) {
  constexpr std::string_view functionName =
      "void CharmmContext::setPsf(std::shared_ptr<CharmmPSF>)";

  if (psf == nullptr) {
    std::string msg(functionName);
    msg += ": CharmmPSF == nullptr";
    throw std::invalid_argument(msg);
  }

  if ((m_NumAtoms > 0) && (m_NumAtoms != psf->getNumAtoms())) {
    std::string msg(functionName);
    msg += ": Number of atoms in CharmmContext and CharmmPSF do not match\n";
    msg += " CharmmContext NATOM: " + std::to_string(m_NumAtoms) + "\n";
    msg += " CharmmPSF NATOM: " + std::to_string(psf->getNumAtoms()) + "\n";
    throw std::invalid_argument(msg);
  }

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
  constexpr std::string_view functionName =
      "void CharmmContext::setForceManager(std::shared_ptr<ForceManager>)";

  if (fm == nullptr) {
    std::string msg(functionName);
    msg += ": ForceManager == nullptr";
    throw std::invalid_argument(msg);
  }

  m_ForceManager = fm;

  this->syncForceManagerFromState();

  this->finalizeSetupIfReady();

  return;
}

void CharmmContext::setNumAtoms(const int numAtoms) {
  constexpr std::string_view functionName =
      "void CharmmContext::setNumAtoms(const int)";

  if (numAtoms <= 0) {
    std::string msg(functionName);
    msg += ": Number of atoms must be positive (numAtoms == " +
           std::to_string(numAtoms) + ")";
    throw std::invalid_argument(msg);
  }

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
  constexpr std::string_view functionName =
      "void CharmmContext::setCoordinatesCharges(const std::vector<double4> &)";

  if (m_NumAtoms <= 0) {
    std::string msg(functionName);
    msg += ": Number of atoms must be set before assigning coordinates";
    throw std::runtime_error(msg);
  }

  if (coordinatesCharges.size() != static_cast<std::size_t>(m_NumAtoms)) {
    std::string msg(functionName);
    msg += ": Array length and number of atoms in CharmmContext do not match\n";
    msg += " NATOM: " + std::to_string(m_NumAtoms) + "\n";
    msg += " Array length: " + std::to_string(coordinatesCharges.size()) + "\n";
    throw std::invalid_argument(msg);
  }

  if ((m_ForceManager != nullptr) && m_ForceManager->isInitialized() &&
      !m_ForceManager->isComposite()) {
    if (coordinatesCharges.size() !=
        static_cast<std::size_t>(m_ForceManager->getNumAtoms())) {
      std::string msg(functionName);
      msg +=
          ": Array length and number of atoms in ForceManager do not match\n";
      msg +=
          " FM NATOM: " + std::to_string(m_ForceManager->getNumAtoms()) + "\n";
      msg +=
          " Array length: " + std::to_string(coordinatesCharges.size()) + "\n";
      throw std::invalid_argument(msg);
    }
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
  constexpr std::string_view functionName =
      "void CharmmContext::setCoordinatesCharges(const "
      "std::vector<std::vector<double>> &)";

  std::vector<double4> coordsCharges;
  coordsCharges.reserve(coordinatesCharges.size());

  for (const std::vector<double> &xyzq : coordinatesCharges) {
    if (xyzq.size() != 4) {
      std::string msg(functionName);
      msg += ": Entry must have 4 values";
      throw std::invalid_argument(msg);
    }

    coordsCharges.push_back(make_double4(xyzq[0], xyzq[1], xyzq[2], xyzq[3]));
  }

  this->setCoordinatesCharges(coordsCharges);

  return;
}

void CharmmContext::setCoordinatesCharges(
    const std::vector<double> &coordinatesCharges) {
  constexpr std::string_view functionName =
      "void CharmmContext::setCoordinatesCharges(const std::vector<double> &)";

  if (coordinatesCharges.size() % 4 != 0) {
    std::string msg(functionName);
    msg += ": Array length must be a multiple of 4";
    throw std::invalid_argument(msg);
  }

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
  constexpr std::string_view functionName =
      "void CharmmContext::setCoordinates(const std::vector<double3> &)";

  if (m_NumAtoms <= 0) {
    std::string msg(functionName);
    msg += ": Number of atoms must be set before assigning coordinates";
    throw std::runtime_error(msg);
  }

  if (coordinates.size() != static_cast<std::size_t>(m_NumAtoms)) {
    std::string msg(functionName);
    msg += ": Array length and number of atoms do not match\n";
    msg += " NATOM: " + std::to_string(m_NumAtoms) + "\n";
    msg += " Array length: " + std::to_string(coordinates.size()) + "\n";
    throw std::invalid_argument(msg);
  }

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
  constexpr std::string_view functionName =
      "void CharmmContext::setCoordinates(const "
      "std::vector<std::vector<double>> &)";

  std::vector<double3> coords;
  coords.reserve(coordinates.size());

  for (const std::vector<double> &xyz : coordinates) {
    if (xyz.size() != 3) {
      std::string msg(functionName);
      msg += ": Entry must have 3 values";
      throw std::invalid_argument(msg);
    }

    coords.push_back(make_double3(xyz[0], xyz[1], xyz[2]));
  }

  this->setCoordinates(coords);

  return;
}

void CharmmContext::setCoordinates(const std::vector<double> &coordinates) {
  constexpr std::string_view functionName =
      "void CharmmContext::setCoordinates(const std::vector<double> &)";

  if (coordinates.size() % 3 != 0) {
    std::string msg(functionName);
    msg += ": Array length must be a multiple of 3";
    throw std::invalid_argument(msg);
  }

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
  constexpr std::string_view functionName =
      "void CharmmContext::setCoordinates(const std::shared_ptr<Coordinates>)";

  if (crd == nullptr) {
    std::string msg(functionName);
    msg += ": Coordinates == nullptr";
    throw std::invalid_argument(msg);
  }

  this->setCoordinates(crd->getCoordinatesDP());

  return;
}

void CharmmContext::setCharges(const std::vector<double> &charges) {
  constexpr std::string_view functionName =
      "void CharmmContext::setCharges(const std::vector<double> &)";

  if (charges.size() != static_cast<std::size_t>(m_NumAtoms)) {
    std::string msg(functionName);
    msg += ": Array length and number of atoms do not match\n";
    msg += " NATOM: " + std::to_string(m_NumAtoms) + "\n";
    msg += " Array length: " + std::to_string(charges.size()) + "\n";
    throw std::invalid_argument(msg);
  }

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
  constexpr std::string_view functionName =
      "CharmmContext::setVelocitiesInverseMasses";

  if (m_NumAtoms == -1) {
    std::string msg(functionName);
    msg += ": Number of atoms must be set before assigning velocities.";
    throw std::runtime_error(msg);
  }

  if (velocitiesInverseMasses.size() != static_cast<std::size_t>(m_NumAtoms)) {
    std::string msg(functionName);
    msg += ": Array length and number of atoms do not match\n";
    msg += " NATOM: " + std::to_string(m_NumAtoms) + "\n";
    msg += " Array length: " + std::to_string(velocitiesInverseMasses.size()) +
           "\n";
    throw std::invalid_argument(msg);
  }

  for (int i = 0; i < m_NumAtoms; i++)
    m_VelocitiesInverseMasses[i] = velocitiesInverseMasses[i];

  m_VelocitiesInverseMasses.transferToDevice();

  return;
}

void CharmmContext::setVelocitiesInverseMasses(
    const std::vector<std::vector<double>> &velocitiesInverseMasses) {
  constexpr std::string_view functionName =
      "void CharmmContext::setVelocitiesInverseMasses(const "
      "std::vector<std::vector<double>> &)";

  std::vector<double4> velMass;
  velMass.reserve(velocitiesInverseMasses.size());

  for (const std::vector<double> &xyzm : velocitiesInverseMasses) {
    if (xyzm.size() != 4) {
      std::string msg(functionName);
      msg += ": Entry must have 4 values";
      throw std::invalid_argument(msg);
    }

    velMass.push_back(make_double4(xyzm[0], xyzm[1], xyzm[2], xyzm[3]));
  }

  this->setVelocitiesInverseMasses(velMass);

  return;
}

void CharmmContext::setVelocitiesInverseMasses(
    const std::vector<double> &velocitiesInverseMasses) {
  constexpr std::string_view functionName =
      "void CharmmContext::setVelocitiesInverseMasses(const "
      "std::vector<double> &)";

  if (velocitiesInverseMasses.size() % 4 != 0) {
    std::string msg(functionName);
    msg += ": Array length must be a multiple of 4";
    throw std::invalid_argument(msg);
  }

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
  constexpr std::string_view functionName =
      "void CharmmContext::setVelocities(const std::vector<double3> &)";

  if (velocities.size() != static_cast<std::size_t>(m_NumAtoms)) {
    std::string msg(functionName);
    msg += ": Array length and number of atoms do not match\n";
    msg += " NATOM: " + std::to_string(m_NumAtoms) + "\n";
    msg += " Array length: " + std::to_string(velocities.size()) + "\n";
    throw std::invalid_argument(msg);
  }

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
  constexpr std::string_view functionName =
      "void CharmmContext::setVelocities(const "
      "std::vector<std::vector<double>> &)";

  std::vector<double3> vels;
  vels.reserve(velocities.size());

  for (const std::vector<double> &xyz : velocities) {
    if (xyz.size() != 3) {
      std::string msg(functionName);
      msg += ": Entry must have 3 values";
      throw std::invalid_argument(msg);
    }

    vels.push_back(make_double3(xyz[0], xyz[1], xyz[2]));
  }

  this->setVelocities(vels);

  return;
}

void CharmmContext::setVelocities(const std::vector<double> &velocities) {
  constexpr std::string_view functionName =
      "void CharmmContext::setVelocities(const std::vector<double> &)";

  if (velocities.size() % 3 != 0) {
    std::string msg(functionName);
    msg += ": Array length must be a multiple of 3";
    throw std::invalid_argument(msg);
  }

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
  constexpr std::string_view functionName =
      "void CharmmContext::setVelocitiesFromCHARMMVelocityFile(const "
      "std::string &)";

  std::ifstream fin(fileName);

  if (!fin.is_open()) {
    std::string msg(functionName);
    msg += ": Could not open CHARMM velocity file \"" + fileName + "\"";
    throw std::runtime_error(msg);
  }

  std::string line;
  while (std::getline(fin, line)) {
    if (!line.empty() && line[0] != '*')
      break;
  }

  if (line.empty()) {
    std::string msg(functionName);
    msg += ": Could not read atom count from CHARMM velocity file \"" +
           fileName + "\"";
    throw std::runtime_error(msg);
  }

  const int nAtoms = std::stoul(line);
  if (nAtoms != m_NumAtoms) {
    std::string msg(functionName);
    msg += ": Number of atom mismatch in CHARMM velocities file\n";
    msg += " NATOM: " + std::to_string(m_NumAtoms) + "\n";
    msg += " MATOM: " + std::to_string(nAtoms) + "\n";
    throw std::invalid_argument(msg);
  }

  int i = 0;

  std::getline(fin, line);
  while (i < m_NumAtoms) {
    if (line.size() == 0) {
      std::string msg(functionName);
      msg += ": Blank line read in \"" + fileName + "\"";
      throw std::runtime_error(msg);
    }

    // auto content = split(line);
    // static_cast<void>(content);

    int atomId, resId, resIdInSeg;
    std::string resName, atomName, segName;
    float x, y, z, bFactor;

    std::stringstream ss(line);
    ss >> atomId >> resId >> resName >> atomName >> x >> y >> z >> segName >>
        resIdInSeg >> bFactor;

    m_VelocitiesInverseMasses[i].x = x;
    m_VelocitiesInverseMasses[i].y = y;
    m_VelocitiesInverseMasses[i].z = z;

    std::getline(fin, line);
    ++i;
  }

  m_VelocitiesInverseMasses.transferToDevice();

  float kineticEnergy = 0.0;
  for (int i = 0; i < m_NumAtoms; ++i) {
    kineticEnergy += 1.0 / m_VelocitiesInverseMasses[i].w *
                     (pow(m_VelocitiesInverseMasses[i].x, 2) +
                      pow(m_VelocitiesInverseMasses[i].y, 2) +
                      pow(m_VelocitiesInverseMasses[i].z, 2));
  }
  kineticEnergy *= 0.5;

  // JEG260702: Debug
  /* *
  const int ndegf = this->getDegreesOfFreedom();

  std::cout << "dof : " << ndegf << "\n";
  const float backTemp =
      kineticEnergy / (1 / 2.0 * ndegf * charmm::constants::kBoltz);
  std::cout << "calculated temp from ke (host) : " << backTemp << "\n";

  m_VelocitiesInverseMasses.transferToDevice();
  std::cout << "calculated temp from ke : " << this->computeTemperature()
            << "\n";
  * */

  return;
}

void CharmmContext::setMasses(const std::vector<double> &masses) {
  constexpr std::string_view functionName =
      "void CharmmContext::setMasses(const std::vector<double> &)";

  if (masses.size() != static_cast<std::size_t>(m_NumAtoms)) {
    std::string msg(functionName);
    msg += ": Array length and number of atoms do not match\n";
    msg += " NATOM: " + std::to_string(m_NumAtoms) + "\n";
    msg += " Array length: " + std::to_string(masses.size()) + "\n";
    throw std::invalid_argument(msg);
  }

  for (int i = 0; i < m_NumAtoms; i++)
    m_VelocitiesInverseMasses[i].w = 1.0 / masses[i];

  m_VelocitiesInverseMasses.transferToDevice();

  return;
}

void CharmmContext::setTemperature(const double temperature) {
  constexpr std::string_view functionName =
      "void CharmmContext::setTemperature(const double)";

  if (!std::isfinite(temperature) || (temperature < 0.0)) {
    std::string msg(functionName);
    msg += ": Temperature must be finite and non-negative";
    throw std::invalid_argument(msg);
  }

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
  constexpr std::string_view functionName =
      "void CharmmContext::setBoxDimensions(const std::vector<double> &)";

  if (!hasValidBoxDimensions(boxDimensions)) {
    std::string msg(functionName);
    msg += ": Box dimensions must be exactly 3 positive values";
    throw std::invalid_argument(msg);
  }

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
  constexpr std::string_view functionName =
      "void CharmmContext::useHolonomicConstraints(const bool)";

  this->requirePsf(functionName);
  this->requireForceManager(functionName);

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
  this->requireForceManager("void CharmmContext::setKappa(const float)");
  m_ForceManager->setKappa(kappa);
  return;
}

void CharmmContext::setCutoff(const float cutoff) {
  this->requireForceManager("void CharmmContext::setCutoff(const float)");
  m_ForceManager->setCutoff(cutoff);
  return;
}

void CharmmContext::setCtonnb(const float ctonnb) {
  this->requireForceManager("void CharmmContext::setCtonnb(const float)");
  m_ForceManager->setCtonnb(ctonnb);
  return;
}

void CharmmContext::setCtofnb(const float ctofnb) {
  this->requireForceManager("void CharmmContext::setCtofnb(const float)");
  m_ForceManager->setCtofnb(ctofnb);
  return;
}

void CharmmContext::setFFTGrid(const int nfftx, const int nffty,
                               const int nfftz) {
  this->requireForceManager(
      "void CharmmContext::setFFTGrid(const int, const int, const int)");
  m_ForceManager->setFFTGrid(nfftx, nffty, nfftz);
  return;
}

void CharmmContext::setPmeSplineOrder(const int pmeSplineOrder) {
  this->requireForceManager("void CharmmContext::setPmeSplineOrder(const int)");
  m_ForceManager->setPmeSplineOrder(pmeSplineOrder);
  return;
}

void CharmmContext::setVdwType(const int vdwType) {
  this->requireForceManager("void CharmmContext::setVdwType(const int)");
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
  if (!m_HasBoxDimensions || !hasValidBoxDimensions(m_BoxDimensions)) {
    throw std::runtime_error("double CharmmContext::getVolume(void) const: Box "
                             "dimensions are not set");
  }

  return (m_BoxDimensions[0] * m_BoxDimensions[1] * m_BoxDimensions[2]);
}

const CudaContainer<double> &CharmmContext::getPressure(void) const {
  return m_Pressure;
}

const std::vector<Bond> &CharmmContext::getBonds(void) const {
  this->requirePsf("CharmmContext::getBonds");
  return m_Psf->getBonds();
}

const CudaContainer<int4> &CharmmContext::getWaterMolecules(void) const {
  this->requirePsf("CharmmContext::getWaterMolecules");
  return m_Psf->getWaterMolecules();
}

const CudaContainer<int4> &CharmmContext::getShakeAtoms(void) const {
  this->requireInitializedForceManager(
      "CudaContainer<int4>& CharmmContext::getShakeAtoms(void)");
  return m_ForceManager->getShakeAtoms();
}

const CudaContainer<float4> &CharmmContext::getShakeParams(void) const {
  this->requireInitializedForceManager(
      "CudaContainer<float4> &CharmmContext::getShakeParams(void)");
  return m_ForceManager->getShakeParams();
}

float CharmmContext::getKappa(void) const {
  this->requireForceManager("CharmmContext::getKappa");
  return m_ForceManager->getKappa();
}

float CharmmContext::getCutoff(void) const {
  this->requireForceManager("CharmmContext::getCutoff");
  return m_ForceManager->getCutoff();
}

float CharmmContext::getCtonnb(void) const {
  this->requireForceManager("CharmmContext::getCtonnb");
  return m_ForceManager->getCtonnb();
}

float CharmmContext::getCtofnb(void) const {
  this->requireForceManager("CharmmContext::getCtofnb");
  return m_ForceManager->getCtofnb();
}

std::vector<int> CharmmContext::getFFTGrid(void) const {
  this->requireForceManager("CharmmContext::getFFTGrid");
  return m_ForceManager->getFFTGrid();
}

int CharmmContext::getPmeSplineOrder(void) const {
  this->requireForceManager("CharmmContext::getPmeSplineOrder");
  return m_ForceManager->getPmeSplineOrder();
}

int CharmmContext::getVdwType(void) const {
  this->requireForceManager("CharmmContext::getVdwType");
  return m_ForceManager->getVdwType();
}

int CharmmContext::getForceStride(void) const {
  this->requireInitializedForceManager("CharmmContext::getForceStride");
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
  this->requireInitializedForceManager("CharmmContext::getPotentialEnergy");
  return m_ForceManager->getPotentialEnergy();
}

float CharmmContext::getPotentialEnergies(void) {
  this->requireInitializedForceManager(
      "float CharmmContext::getPotentialEnergies(void)");
  return m_ForceManager->getPotentialEnergies();
}

std::shared_ptr<Force<double>> CharmmContext::getForces(void) {
  this->requireInitializedForceManager(
      "std::shared_ptr<Force<double>> CharmmContext::getForces(void)");
  return m_ForceManager->getForces();
}

CudaContainer<double> &CharmmContext::getVirial(void) {
  this->requireInitializedForceManager("CharmmContext::getVirial");
  return m_ForceManager->getVirial();
}

std::vector<Bond> &CharmmContext::getBonds(void) {
  this->requirePsf("CharmmContext::getBonds");
  return m_Psf->getBonds();
}

CudaContainer<int4> &CharmmContext::getWaterMolecules(void) {
  this->requirePsf("CharmmContext::getWaterMolecules");
  return m_Psf->getWaterMolecules();
}

CudaContainer<int4> &CharmmContext::getShakeAtoms(void) {
  this->requireInitializedForceManager(
      "CudaContainer<int4>& CharmmContext::getShakeAtoms(void)");
  return m_ForceManager->getShakeAtoms();
}

CudaContainer<float4> &CharmmContext::getShakeParams(void) {
  this->requireInitializedForceManager(
      "CudaContainer<float4> &CharmmContext::getShakeParams(void)");
  return m_ForceManager->getShakeParams();
}

void CharmmContext::assignVelocitiesAtTemperature(const double temperature) {
  constexpr std::string_view functionName =
      "void CharmmContext::assignVelocitiesAtTemperature(const double)";

  if (m_NumAtoms == -1) {
    std::string msg(functionName);
    msg += ": The number of atoms has not been set.\n";
    msg += "HINT: Trying setting coordinates or number of atoms before "
           "assigning velocities.";
    throw std::invalid_argument(msg);
  }

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

  // this->removeCenterOfMassMotion();

  // JEG260702: This code is not used for anything.
  /* *
  double kineticEnergy = 0.0;
  for (int i = 0; i < m_NumAtoms; ++i) {
    kineticEnergy +=
        1.0 / m_VelocitiesInverseMasses[i].w *
        (m_VelocitiesInverseMasses[i].x * m_VelocityInverseMasses[i].x +
         m_VelocitiesInverseMasses[i].y * m_VelocityInverseMasses[i].y +
         m_VelocitiesInverseMasses[i].z * m_VelocityInverseMasses[i].z);
  }
  kineticEnergy *= 0.5;

  const double ndegf = static_cast<double>(this->getDegreesOfFreedom());
  const double backTemp =
      kineticEnergy / (1 / 2.0 * ndegf * charmm::constants::kBoltz);
  static_cast<void>(backTemp);
  * */

  return;
}

double CharmmContext::computeTemperature(void) {
  if ((m_NumAtoms == -1) || (m_VelocitiesInverseMasses.size() == 0)) {
    throw std::invalid_argument(
        "No atoms in the system -- coordinates have not been loaded and/or "
        "velocities not assigned.");
  }

  const double kineticEnergy = this->getKineticEnergy();
  return kineticEnergy / (0.5 * static_cast<double>(m_NumDegreesOfFreedom) *
                          charmm::constants::kBoltz);
}

/* *
__global__ static void CalculateKineticPressureKernel(
    double *accumulant, const double4 *__restrict__ velocitiesInverseMasses,
    const int numAtoms) {
  constexpr int blockSize = 128 * 9;
  __shared__ double sdata[blockSize];
  int threadId = threadIdx.x;

  for (int i = 0; i < 9; ++i)
    sdata[9 * threadId + i] = 0.0;

  int index = threadIdx.x + blockIdx.x * blockDim.x;

  for (int i = index; i < numAtoms; i += blockDim.x)
    if (index < numAtoms) {
      const double rvc = 0.5 / velocitiesInverseMasses[i].w;
      sdata[threadId + 0] =
          rvc * velocitiesInverseMasses[i].x * velocitiesInverseMasses[i].x;
      sdata[threadId + 1] =
          rvc * velocitiesInverseMasses[i].x * velocitiesInverseMasses[i].y;
      sdata[threadId + 2] =
          rvc * velocitiesInverseMasses[i].x * velocitiesInverseMasses[i].z;
      sdata[threadId + 3] =
          rvc * velocitiesInverseMasses[i].y * velocitiesInverseMasses[i].x;
      sdata[threadId + 4] =
          rvc * velocitiesInverseMasses[i].y * velocitiesInverseMasses[i].y;
      sdata[threadId + 5] =
          rvc * velocitiesInverseMasses[i].y * velocitiesInverseMasses[i].z;
      sdata[threadId + 6] =
          rvc * velocitiesInverseMasses[i].z * velocitiesInverseMasses[i].x;
      sdata[threadId + 7] =
          rvc * velocitiesInverseMasses[i].z * velocitiesInverseMasses[i].y;
      sdata[threadId + 8] =
          rvc * velocitiesInverseMasses[i].z * velocitiesInverseMasses[i].z;
    }

  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    __syncthreads();
    if (threadId < s) {
      sdata[threadId + 0] += sdata[threadId + s * 9 + 0];
      sdata[threadId + 1] += sdata[threadId + s * 9 + 1];
      sdata[threadId + 2] += sdata[threadId + s * 9 + 2];
      sdata[threadId + 3] += sdata[threadId + s * 9 + 3];
      sdata[threadId + 4] += sdata[threadId + s * 9 + 4];
      sdata[threadId + 5] += sdata[threadId + s * 9 + 5];
      sdata[threadId + 6] += sdata[threadId + s * 9 + 6];
      sdata[threadId + 7] += sdata[threadId + s * 9 + 7];
      sdata[threadId + 8] += sdata[threadId + s * 9 + 8];
    }
  }

  if (threadId == 0) {
#pragma unroll
    for (int i = 0; i < 9; i++)
      atomicAdd(accumulant + i, sdata[i]);
  }

  return;
}
* */

void CharmmContext::computePressure(void) {
  throw std::runtime_error("void CharmmContext::computePressure(void): Has not "
                           "been implemented yet");

  /* *
  this->requireInitializedForceManager(
      "void CharmmContext::computePressure(void)");

  CalculateKineticPressureKernel<<<1, 1024>>>(
      m_VirialKineticEnergyTensor.getDeviceArray().data(),
      m_VelocitiesInverseMasses.getDeviceArray().data(), m_NumAtoms);

  cudaCheck(cudaDeviceSynchronize());

  const double vcell = charmm::constants::patmos / this->getVolume();

  CudaContainer<double> &virial = m_ForceManager->getVirial();
  virial.transferFromDevice();

  m_VirialKineticEnergyTensor.transferFromDevice();

  for (int i = 0; i < 9; i++) {
    m_Pressure[i] = (2.0 * m_VirialKineticEnergyTensor.getHostArray()[i] +
                     virial.getHostArray()[i]) /
                    vcell;
  }

  m_Pressure.transferToDevice();
  * */

  return;
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
  this->requireInitializedForceManager("CharmmContext::imageCentering");
  this->requirePsf("CharmmContext::imageCentering");

  double *forces = this->getForces()->xyz();
  const int forceStride = this->getForceStride();

  const CudaContainer<int2> &groups = m_Psf->getGroups();
  const int numGroups = groups.size();

  constexpr int numThreads = 128;
  const int numBlocks = (numGroups + numThreads - 1) / numThreads;

  ImageCenteringKernel<<<numBlocks, numThreads>>>(
      m_CoordinatesChargesDP.getDeviceArray().data(),
      m_CoordinatesChargesSP.getDeviceArray().data(),
      m_VelocitiesInverseMasses.getDeviceArray().data(), forces, forceStride,
      groups.getDeviceArray().data(), numGroups, m_BoxDimensions[0],
      m_BoxDimensions[1], m_BoxDimensions[2], m_Pbc);

  cudaCheck(cudaDeviceSynchronize());

  return;
}

void CharmmContext::resetNeighborList(void) {
  this->requireInitializedForceManager("CharmmContext::resetNeighborList");

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
  cudaCheck(
      cudaMemset(static_cast<void *>(m_KineticEnergy.getDeviceArray().data()),
                 0, sizeof(double)));

  CalculateKineticEnergyKernel<<<1, 1024>>>(
      m_KineticEnergy.getDeviceArray().data(),
      m_VelocitiesInverseMasses.getDeviceArray().data(), m_NumAtoms);

  cudaCheck(cudaDeviceSynchronize());

  return;
}

void CharmmContext::calculatePotentialEnergy(const bool reset,
                                             const bool print) {
  this->requireInitializedForceManager(
      "CharmmContext::calculatePotentialEnergy");

  m_ForceManager->calcForce(m_CoordinatesChargesSP.getDeviceArray().data(),
                            reset, true, true);

  if (print == true)
    std::cout << "Implement some nice printing" << std::endl;

  return;
}

void CharmmContext::calculateForces(bool reset, bool calcEnergy,
                                    bool calcVirial) {
  this->requireInitializedForceManager("CharmmContext::calculateForces");

  m_ForceManager->calcForce(m_CoordinatesChargesSP.getDeviceArray().data(),
                            reset, calcEnergy, calcVirial);

  return;
}

void CharmmContext::linkBackForceManager(void) {
  this->requireForceManager("void CharmmContext::linkBackForceManager(void)");

  std::shared_ptr<CharmmContext> self = this->weak_from_this().lock();
  if (self == nullptr)
    return;

  m_ForceManager->setContext(self);

  return;
}

void CharmmContext::syncStateFromForceManager(void) {
  this->requireForceManager(
      "void CharmmContext::syncStateFromForceManager(void)");

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
  this->requireForceManager(
      "void CharmmContext::syncForceManagerFromState(void)");

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

void CharmmContext::requirePsf(const std::string_view functionName) const {
  if (m_Psf == nullptr) {
    std::string msg(functionName);
    msg += ": CharmmPSF is not set";
    throw std::invalid_argument(msg);
  }
  return;
}

void CharmmContext::requireForceManager(
    const std::string_view functionName) const {
  if (m_ForceManager == nullptr) {
    std::string msg(functionName);
    msg += ": ForceManager is not set";
    throw std::invalid_argument(msg);
  }
  return;
}

void CharmmContext::requireInitializedForceManager(
    const std::string_view functionName) const {
  this->requireForceManager(functionName);

  if (!m_ForceManager->isInitialized()) {
    std::string msg(functionName);
    msg += ": ForceManager is not initialized";
    throw std::runtime_error(msg);
  }

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

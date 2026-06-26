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

#include "CharmmCrd.h"
// #include "Checkpoint.h"
#include "Constants.h"
#include "PBC.h"
#include "cuda_utils.h"
#include "gpu_utils.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cuda.h>
#include <cuda_runtime.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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
      m_VelocitiesInverseMasses(), m_KineticEnergy(1), m_Pressure(9),
      m_VirialKineticEnergyTensor(9), m_Temperature(0.0f),
      m_UsingHolonomicConstraints(false) {}

CharmmContext::CharmmContext(std::shared_ptr<ForceManager> forceManager)
    : CharmmContext() {
  std::vector<int> devices = {0, 1, 2, 3};
  start_gpu(1, 1, 0, devices);

  m_ForceManager = forceManager;

  this->requireForceManager(
      "CharmmContext::CharmmContext(std::shared_ptr<ForceManager>)");

  this->syncStateFromForceManager();

  this->useHolonomicConstraints(true);

  if (!m_ForceManager->isInitialized())
    m_ForceManager->initialize();

  this->syncStateFromForceManager();

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
  m_VelocitiesInverseMasses = other.m_VelocitiesInverseMasses;
  m_KineticEnergy = other.m_KineticEnergy;
  m_Pressure = other.m_Pressure;
  m_Temperature = other.m_Temperature;
  m_NumDegreesOfFreedom = other.m_NumDegreesOfFreedom;
  m_UsingHolonomicConstraints = other.m_UsingHolonomicConstraints;
  m_RandomSeed = other.m_RandomSeed;
}

// void CharmmContext::setupFromCheckpoint(
//     std::shared_ptr<Checkpoint> checkpoint) {
//   m_NumAtoms = checkpoint->get<int>("numAtoms");
//   return;
// }

void CharmmContext::setMasses(const std::vector<double> &masses) {
  if (masses.size() != static_cast<std::size_t>(m_NumAtoms)) {
    std::stringstream tmpexc;
    tmpexc << "Masses vector size does not match numAtoms (" << masses.size()
           << " != " << m_NumAtoms << ")\n";
    throw std::invalid_argument(tmpexc.str());
  }

  for (int i = 0; i < m_NumAtoms; i++)
    m_VelocitiesInverseMasses[i].w = 1.0 / masses[i];

  m_VelocitiesInverseMasses.transferToDevice();

  return;
}

void CharmmContext::setNumAtoms(const int numAtoms) {
  m_NumAtoms = numAtoms;
  return;
}

void CharmmContext::setCoordinates(const std::shared_ptr<Coordinates> crd) {
  this->requireForceManager("CharmmContext::setCoordinates");
  this->setCoordinates(crd->getCoordinatesD());
  return;
}

void CharmmContext::setCoordinates(const std::vector<double4> &coordinates) {
  this->requirePsf("CharmmContext::setCoordinates");

  if (!m_ForceManager->isComposite())
    assert(coordinates.size() == m_ForceManager->getNumAtoms());

  if (!m_ForceManager->hasCharmmContext())
    this->linkBackForceManager();

  this->setNumAtoms(coordinates.size());

  m_VelocitiesInverseMasses.resize(m_NumAtoms);
  this->setMasses(m_Psf->getMasses());

  this->useHolonomicConstraints(m_UsingHolonomicConstraints);
  std::vector<double> charges = m_Psf->getCharges();

  m_CoordinatesChargesSP.resize(m_NumAtoms);
  m_CoordinatesChargesDP.resize(m_NumAtoms);
  for (int i = 0; i < m_NumAtoms; i++) {
    m_CoordinatesChargesSP[i] = make_float4(
        static_cast<float>(coordinates[i].x),
        static_cast<float>(coordinates[i].y),
        static_cast<float>(coordinates[i].z), static_cast<float>(charges[i]));
    m_CoordinatesChargesDP[i] = make_double4(coordinates[i].x, coordinates[i].y,
                                             coordinates[i].z, charges[i]);
  }

  m_CoordinatesChargesSP.transferToDevice();
  m_CoordinatesChargesDP.transferToDevice();

  this->resetNeighborList();

  return;
}

std::vector<std::vector<double>> CharmmContext::getCoordinates(void) {
  m_CoordinatesChargesDP.transferFromDevice();

  std::vector<std::vector<double>> coordinates;
  for (int i = 0; i < m_NumAtoms; i++) {
    std::vector<double> coordinate = {m_CoordinatesChargesDP[i].x,
                                      m_CoordinatesChargesDP[i].y,
                                      m_CoordinatesChargesDP[i].z};
    coordinates.push_back(coordinate);
  }

  return coordinates;
}

void CharmmContext::setCoords(const std::vector<float> &coordinates) {
  assert(coordinates.size() == static_cast<std::size_t>(m_NumAtoms * 3));

  for (int i = 0; i < m_NumAtoms; i++) {
    m_CoordinatesChargesSP[i].x = coordinates[i * 3 + 0];
    m_CoordinatesChargesSP[i].y = coordinates[i * 3 + 1];
    m_CoordinatesChargesSP[i].z = coordinates[i * 3 + 2];
  }

  m_CoordinatesChargesSP.transferToDevice();
  this->resetNeighborList();

  return;
}

int CharmmContext::getNumAtoms(void) const { return m_NumAtoms; }

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
  this->requireForceManager("CharmmContext::imageCentering");
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
  this->requireForceManager("CharmmContext::resetNeighborList");

  this->imageCentering();

  m_ForceManager->resetNeighborList(
      m_CoordinatesChargesSP.getDeviceArray().data());

  return;
}

void CharmmContext::calculateForces(bool reset, bool calcEnergy,
                                    bool calcVirial) {
  this->requireForceManager("CharmmContext::calculateForces");

  m_ForceManager->calcForce(m_CoordinatesChargesSP.getDeviceArray().data(),
                            reset, calcEnergy, calcVirial);

  return;
}

float CharmmContext::getPotentialEnergies(void) {
  return m_ForceManager->getPotentialEnergies();
}

std::shared_ptr<Force<double>> CharmmContext::getForces(void) {
  this->requireForceManager("CharmmContext::getForces");
  return m_ForceManager->getForces();
}

float CharmmContext::getTemperature(void) const { return m_Temperature; }

void CharmmContext::setTemperature(const float temperature) {
  m_Temperature = temperature;
  return;
}

__global__ void calculateCenterOfMassMomemtumKernel(
    int numAtoms, double4 *__restrict__ velocitiesInverseMasses) {
  return;
}

void CharmmContext::removeCenterOfMassMotion(void) { return; }

void CharmmContext::assignVelocitiesAtTemperature(const float temperature) {
  if (m_NumAtoms == -1) {
    throw std::invalid_argument(
        "numAtoms = -1 in CharmmContext::assignVelocitiesAtTemperature -- This "
        "Context object was not initialized properly (no Coordinate given?).\n "
        "Make sure you used setCoordinates before trying to assignVelocities.");
  }

  this->setTemperature(temperature);
  const double boltz = charmm::constants::kBoltz * m_Temperature;

  std::mt19937 gen{m_RandomSeed};

  for (int i = 0; i < m_NumAtoms; i++) {
    double sd = boltz * m_VelocitiesInverseMasses[i].w;
    sd = sqrt(sd);

    std::normal_distribution<> d(0, sd);

    m_VelocitiesInverseMasses[i].x = d(gen);
    m_VelocitiesInverseMasses[i].y = d(gen);
    m_VelocitiesInverseMasses[i].z = d(gen);
  }

  m_VelocitiesInverseMasses.transferToDevice();
  this->removeCenterOfMassMotion();

  float kineticEnergy = 0.0;
  for (int i = 0; i < m_NumAtoms; ++i) {
    kineticEnergy += 1.0 / m_VelocitiesInverseMasses[i].w *
                     (pow(m_VelocitiesInverseMasses[i].x, 2) +
                      pow(m_VelocitiesInverseMasses[i].y, 2) +
                      pow(m_VelocitiesInverseMasses[i].z, 2));
  }
  kineticEnergy *= 0.5;
  const int ndegf = this->getDegreesOfFreedom();

  const float backTemp =
      kineticEnergy / (1 / 2.0 * ndegf * charmm::constants::kBoltz);
  static_cast<void>(backTemp);

  m_VelocitiesInverseMasses.transferToDevice();

  return;
}

static std::vector<std::string> split(std::string line) {
  std::stringstream ss(line);
  std::string atomId, resId, resName, atom, x, y, z;
  ss >> atomId >> resId >> resName >> atom >> x >> y >> z;
  std::vector<std::string> content = {atomId, resId, resName, atom, x, y, z};

  return content;
}

void CharmmContext::assignVelocitiesFromCHARMMVelocityFile(
    const std::string &fileName) {
  std::ifstream fin(fileName);

  if (!fin.is_open()) {
    throw std::invalid_argument("Could not open CHARMM velocity file ");
    exit(0);
  }

  std::string line;
  while (1) {
    std::getline(fin, line);
    if (line[0] != '*')
      break;
  }

  const int nAtoms = std::stoul(line);
  assert(nAtoms == m_NumAtoms);

  int i = 0;

  std::getline(fin, line);
  while (i < m_NumAtoms) {
    if (line.size() == 0) {
      throw std::invalid_argument("ERROR: Blank line read in " + fileName +
                                  "\n. Exiting\n");
      exit(0);
    }

    auto content = split(line);
    static_cast<void>(content);

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

  float kineticEnergy = 0.0;
  for (int i = 0; i < m_NumAtoms; ++i) {
    kineticEnergy += 1.0 / m_VelocitiesInverseMasses[i].w *
                     (pow(m_VelocitiesInverseMasses[i].x, 2) +
                      pow(m_VelocitiesInverseMasses[i].y, 2) +
                      pow(m_VelocitiesInverseMasses[i].z, 2));
  }
  kineticEnergy *= 0.5;

  const int ndegf = this->getDegreesOfFreedom();

  std::cout << "dof : " << ndegf << "\n";
  const float backTemp =
      kineticEnergy / (1 / 2.0 * ndegf * charmm::constants::kBoltz);
  std::cout << "calculated temp from ke (host) : " << backTemp << "\n";

  m_VelocitiesInverseMasses.transferToDevice();
  std::cout << "calculated temp from ke : " << this->computeTemperature()
            << "\n";
  return;
}

void CharmmContext::assignVelocities(const std::vector<double> &velocities) {
  assert(velocities.size() == static_cast<std::size_t>(m_NumAtoms * 3));

  for (int i = 0; i < m_NumAtoms; ++i) {
    m_VelocitiesInverseMasses[i].x = velocities[i * 3];
    m_VelocitiesInverseMasses[i].y = velocities[i * 3 + 1];
    m_VelocitiesInverseMasses[i].z = velocities[i * 3 + 2];
  }

  m_VelocitiesInverseMasses.transferToDevice();
  return;
}

void CharmmContext::assignVelocities(
    const std::vector<std::vector<double>> &velocities) {
  assert(velocities.size() == static_cast<std::size_t>(m_NumAtoms));

  for (int i = 0; i < m_NumAtoms; ++i) {
    assert(velocities[i].size() == 3);
    m_VelocitiesInverseMasses[i].x = velocities[i][0];
    m_VelocitiesInverseMasses[i].y = velocities[i][1];
    m_VelocitiesInverseMasses[i].z = velocities[i][2];
  }

  m_VelocitiesInverseMasses.transferToDevice();
  return;
}

CudaContainer<double4> &CharmmContext::getVelocityMass(void) {
  return m_VelocitiesInverseMasses;
}

CudaContainer<double4> &CharmmContext::getCoordinatesCharges(void) {
  return m_CoordinatesChargesDP;
}

const CudaContainer<float4> &CharmmContext::getXYZQ(void) const {
  return m_CoordinatesChargesSP;
}

CudaContainer<float4> &CharmmContext::getXYZQ(void) {
  return m_CoordinatesChargesSP;
}

int CharmmContext::getForceStride(void) const {
  return m_ForceManager->getForceStride();
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

double CharmmContext::getKineticEnergy(void) {
  this->calculateKineticEnergy();
  m_KineticEnergy.transferFromDevice();
  return m_KineticEnergy.getHostArray()[0];
}

float CharmmContext::computeTemperature(void) {
  if (m_NumAtoms == -1 || m_VelocitiesInverseMasses.size() == 0) {
    throw std::invalid_argument(
        "No atoms in the system -- coordinates have not been loaded and/or "
        "velocities not assigned.");
  }

  const double kineticEnergy = this->getKineticEnergy();
  const int numDegreesOfFreedom = this->getDegreesOfFreedom();
  return kineticEnergy /
         (0.5 * numDegreesOfFreedom * charmm::constants::kBoltz);
}

void CharmmContext::setPeriodicBoundaryCondition(const PBC pbc) {
  m_Pbc = pbc;
  m_HasPbc = true;

  if (m_ForceManager != nullptr) {
    m_ForceManager->setPeriodicBoundaryCondition(m_Pbc);
    m_Pbc = m_ForceManager->getPeriodicBoundaryCondition();
    this->resetNeighborList();
  }

  return;
}

PBC CharmmContext::getPeriodicBoundaryCondition(void) const { return m_Pbc; }

const std::vector<double> &CharmmContext::getBoxDimensions(void) const {
  return m_BoxDimensions;
}

std::vector<double> &CharmmContext::getBoxDimensions(void) {
  return m_BoxDimensions;
}

void CharmmContext::setBoxDimensions(const std::vector<double> &boxDimensions) {
  if (!hasValidBoxDimensions(boxDimensions)) {
    throw std::invalid_argument("CharmmContext::setBoxDimensions: Box "
                                "dimensions must be exactly 3 positive values");
  }

  m_BoxDimensions = boxDimensions;
  m_HasBoxDimensions = true;

  if (m_ForceManager != nullptr) {
    m_ForceManager->setBoxDimensions(m_BoxDimensions);
    m_BoxDimensions = m_ForceManager->getBoxDimensions();
    m_HasBoxDimensions = hasValidBoxDimensions(m_BoxDimensions);
    this->initializeForceManagerIfReady();
  }

  return;
}

std::vector<Bond> CharmmContext::getBonds(void) {
  this->requirePsf("CharmmContext::getBonds");
  return m_Psf->getBonds();
}

int CharmmContext::getDegreesOfFreedom(void) const {
  return m_NumDegreesOfFreedom;
}

void CharmmContext::calculatePotentialEnergy(bool reset, bool print) {
  m_ForceManager->calcForce(m_CoordinatesChargesSP.getDeviceArray().data(),
                            reset, true, true);
  if (print == true)
    std::cout << "Implement some nice printing" << std::endl;
  return;
}

CudaContainer<double> &CharmmContext::getPotentialEnergy(void) {
  return m_ForceManager->getPotentialEnergy();
}

double CharmmContext::getVolume(void) const {
  return (m_BoxDimensions[0] * m_BoxDimensions[1] * m_BoxDimensions[2]);
}

static __global__ void
calculateKineticKernel(int numAtoms,
                       const double4 *__restrict__ velocitiesInverseMasses,
                       double *accumulant) {
  constexpr int blockSize = 128 * 9;
  __shared__ double sdata[blockSize];
  int threadId = threadIdx.x;

  for (int i = 0; i < 9; ++i)
    sdata[9 * threadId + i] = 0.0;

  int index = threadIdx.x + blockIdx.x * blockDim.x;

  for (int i = index; i < numAtoms; i += blockDim.x)
    if (index < numAtoms) {
      double rvc = 0.5 / velocitiesInverseMasses[index].w;
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
}

void CharmmContext::computePressure(void) {
  auto kineticEnergy = this->getKineticEnergy();
  static_cast<void>(kineticEnergy);

  calculateKineticKernel<<<1, 1024>>>(
      m_NumAtoms, m_VelocitiesInverseMasses.getDeviceArray().data(),
      m_VirialKineticEnergyTensor.getDeviceArray().data());

  cudaCheck(cudaDeviceSynchronize());

  const double vcell = charmm::constants::patmos / this->getVolume();

  auto virial = m_ForceManager->getVirial();
  virial.transferFromDevice();

  m_VirialKineticEnergyTensor.transferFromDevice();

  for (int i = 0; i < 9; ++i) {
    m_Pressure[i] = (2.0 * m_VirialKineticEnergyTensor.getHostArray()[i] +
                     virial.getHostArray()[i]) /
                    vcell;
  }

  m_Pressure.transferToDevice();
  return;
}

CudaContainer<double> &CharmmContext::getVirial(void) {
  return m_ForceManager->getVirial();
}

CudaContainer<int4> CharmmContext::getWaterMolecules(void) {
  this->requirePsf("CharmmContext::getWaterMolecules");
  return m_Psf->getWaterMolecules();
}

CudaContainer<int4> CharmmContext::getShakeAtoms(void) {
  return m_ForceManager->getShakeAtoms();
}

CudaContainer<float4> CharmmContext::getShakeParams(void) {
  return m_ForceManager->getShakeParams();
}

void CharmmContext::useHolonomicConstraints(const bool useConstraints) {
  m_UsingHolonomicConstraints = useConstraints;
  int ndegf = m_NumAtoms * 3;

  if (m_Pbc == PBC::P1)
    ndegf -= 3;
  else if (m_Pbc == PBC::P21)
    ndegf -= 1;

  if (m_UsingHolonomicConstraints) {
    ndegf -= this->getWaterMolecules().size() * 3;
    int numShakeConstraints = 0;
    auto shakeAtoms = m_ForceManager->getShakeAtoms().getHostArray();
    for (std::size_t i = 0; i < shakeAtoms.size(); ++i) {
      ++numShakeConstraints;
      if (shakeAtoms[i].z != -1)
        ++numShakeConstraints;
      if (shakeAtoms[i].w != -1)
        ++numShakeConstraints;
    }
    ndegf -= numShakeConstraints;
  }

  m_NumDegreesOfFreedom = ndegf;

  return;
}

void CharmmContext::orient(void) { return; }

void CharmmContext::setPsf(std::shared_ptr<CharmmPSF> psf) {
  if (psf == nullptr)
    throw std::invalid_argument("CharmmContext::setPsf: psf == nullptr");

  m_Psf = psf;

  if (m_ForceManager != nullptr) {
    m_ForceManager->setPsf(m_Psf);
    this->initializeForceManagerIfReady();
  }

  return;
}

void CharmmContext::setPrm(std::shared_ptr<CharmmParameters> prm) {
  if (prm == nullptr)
    throw std::invalid_argument("CharmmContext::setPrm: prm == nullptr");

  m_Prm = prm;

  if (m_ForceManager != nullptr) {
    m_ForceManager->setPrm(m_Prm);
    this->initializeForceManagerIfReady();
  }

  return;
}

void CharmmContext::setForceManager(
    std::shared_ptr<ForceManager> forceManager) {
  if (forceManager == nullptr) {
    throw std::invalid_argument(
        "CharmmContext::setForceManager: forceManager == nullptr");
  }

  m_ForceManager = forceManager;
  this->syncForceManagerFromState();
  this->initializeForceManagerIfReady();

  return;
}

std::shared_ptr<ForceManager> CharmmContext::getForceManager(void) {
  return m_ForceManager;
}

std::shared_ptr<CharmmPSF> CharmmContext::getPsf(void) const { return m_Psf; }

std::shared_ptr<CharmmParameters> CharmmContext::getPrm(void) const {
  return m_Prm;
}

void CharmmContext::linkBackForceManager(void) {
  m_ForceManager->setContext(shared_from_this());
  return;
}

void CharmmContext::writeCrd(std::string fileName) {
  std::ofstream fout(fileName);

  if (!fout.is_open())
    throw std::invalid_argument("ERROR! Can't open the crd file to write \n");

  return;
}

int CharmmContext::getNumDegreesOfFreedom(void) const {
  return m_NumDegreesOfFreedom;
}

void CharmmContext::setRandomSeedForVelocities(const std::uint64_t randomSeed) {
  m_RandomSeed = randomSeed;
  return;
}

std::uint64_t CharmmContext::getRandomSeedForVelocities(void) const {
  return m_RandomSeed;
}

CudaContainer<double> CharmmContext::getPressure(void) const {
  return m_Pressure;
}

bool CharmmContext::isUsingHolonomicConstraints(void) const {
  return m_UsingHolonomicConstraints;
}

void CharmmContext::syncStateFromForceManager(void) {
  this->requireForceManager("CharmmContext::syncStateFromForceManager");

  m_Psf = m_ForceManager->getPsf();
  m_Prm = m_ForceManager->getPrm();

  m_BoxDimensions = m_ForceManager->getBoxDimensions();
  m_HasBoxDimensions = hasValidBoxDimensions(m_BoxDimensions);

  m_Pbc = m_ForceManager->getPeriodicBoundaryCondition();
  m_HasPbc = true;

  return;
}

void CharmmContext::syncForceManagerFromState(void) {
  this->requireForceManager("CharmmContext::syncForceManagerFromState");

  if (m_Psf != nullptr)
    m_ForceManager->setPsf(m_Psf);
  else
    m_Psf = m_ForceManager->getPsf();

  if (m_Prm != nullptr)
    m_ForceManager->setPrm(m_Prm);
  else
    m_Prm = m_ForceManager->getPrm();

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

void CharmmContext::requirePsf(const std::string &functionName) const {
  if (m_Psf == nullptr)
    throw std::invalid_argument(functionName + ": CharmmPSF is not set");
  return;
}

void CharmmContext::requireForceManager(const std::string &functionName) const {
  if (m_ForceManager == nullptr)
    throw std::invalid_argument(functionName + ": ForceManager is not set");
  return;
}

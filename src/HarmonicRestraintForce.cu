// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author:  James E. Gonzales II, Samarjeet Prasad
//
// ENDLICENSE

#include "HarmonicRestraintForce.h"

#include "ApoCharmmError.h"
#include "cuda_utils.h"
#include "gpu_utils.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

template <typename AT, typename CT>
HarmonicRestraintForce<AT, CT>::HarmonicRestraintForce(const int numAtoms)
    : m_NumAtoms(numAtoms), m_ForceConstants(), m_ReferenceCoordinates(),
      m_BoxDimensions(3), m_Selection(0, AtomSelection::InitialValue::ALL),
      m_EnergyVirial(nullptr), m_Forces(nullptr), m_Stream(nullptr) {
  APOCHARMM_REQUIRE(numAtoms > 0, ApoCharmmErrorCode::InvalidArgument,
                    "Atom count must be positive; observed " +
                        std::to_string(numAtoms));

  m_ForceConstants.resize(static_cast<std::size_t>(numAtoms));
  m_ReferenceCoordinates.resize(static_cast<std::size_t>(numAtoms));
  m_Selection.setNumAtoms(numAtoms, AtomSelection::InitialValue::ALL);

  m_ForceConstants.set(0.0);

  m_ReferenceCoordinates.set(make_double4(0.0, 0.0, 0.0, 1.0));

  m_BoxDimensions.set(0.0);

  m_EnergyVirial = std::make_shared<CudaEnergyVirial>();
  m_EnergyVirial->insert("harm");

  m_Forces = std::make_shared<Force<AT>>();
  m_Forces->realloc(numAtoms, 1.5f);

  m_Stream = std::make_shared<cudaStream_t>();
  cudaCheck(cudaStreamCreate(m_Stream.get()));
}

template <typename AT, typename CT>
HarmonicRestraintForce<AT, CT>::~HarmonicRestraintForce(void) noexcept {
  if (m_Stream != nullptr) {
    destroy_cuda_stream_noexcept(m_Stream.get());
    m_Stream.reset();
  }
}

template <typename AT, typename CT>
void HarmonicRestraintForce<AT, CT>::setSelection(
    const AtomSelection &selection) {
  APOCHARMM_REQUIRE(selection.getNumAtoms() == m_NumAtoms,
                    ApoCharmmErrorCode::InvalidArgument,
                    "Selection atom count mismatch; expected " +
                        std::to_string(m_NumAtoms) + ", observed " +
                        std::to_string(selection.getNumAtoms()));

  m_Selection = selection;

  for (int i = 0; i < m_NumAtoms; i++) {
    if (!selection.contains(i))
      m_ForceConstants[i] = 0.0;
  }
  m_ForceConstants.transferToDevice();

  return;
}

template <typename AT, typename CT>
void HarmonicRestraintForce<AT, CT>::setForceConstant(
    const double forceConstant) {
  APOCHARMM_REQUIRE(std::isfinite(forceConstant),
                    ApoCharmmErrorCode::InvalidArgument,
                    "Force constant must be finite; observed " +
                        std::to_string(forceConstant));

  APOCHARMM_REQUIRE(forceConstant >= 0.0, ApoCharmmErrorCode::InvalidArgument,
                    "Force constant must be non-negative; observed " +
                        std::to_string(forceConstant));

  m_ForceConstants.set(0.0);

  for (const int i : m_Selection.getAtomIndices())
    m_ForceConstants[i] = forceConstant;

  m_ForceConstants.transferToDevice();

  return;
}

template <typename AT, typename CT>
void HarmonicRestraintForce<AT, CT>::setForceConstants(
    const std::vector<double> &forceConstants) {
  APOCHARMM_REQUIRE(forceConstants.size() ==
                        static_cast<std::size_t>(m_NumAtoms),
                    ApoCharmmErrorCode::InvalidArgument,
                    "Force-constant array size mismatch; expected " +
                        std::to_string(m_NumAtoms) + ", observed " +
                        std::to_string(forceConstants.size()));

  const std::vector<int> selectedAtomIndices = m_Selection.getAtomIndices();

  for (const int i : selectedAtomIndices) {
    APOCHARMM_REQUIRE(
        std::isfinite(forceConstants[i]), ApoCharmmErrorCode::InvalidArgument,
        "Force constant at index " + std::to_string(i) +
            " must be finite; observed " + std::to_string(forceConstants[i]));

    APOCHARMM_REQUIRE(forceConstants[i] >= 0.0,
                      ApoCharmmErrorCode::InvalidArgument,
                      "Force constant at index " + std::to_string(i) +
                          " must be non-negative; observed " +
                          std::to_string(forceConstants[i]));
  }

  m_ForceConstants.set(0.0);

  for (const int i : selectedAtomIndices)
    m_ForceConstants[i] = forceConstants[i];

  m_ForceConstants.transferToDevice();

  return;
}

template <typename AT, typename CT>
void HarmonicRestraintForce<AT, CT>::setReferenceCoordinates(
    const std::vector<double3> &referenceCoordinates) {
  APOCHARMM_REQUIRE(referenceCoordinates.size() ==
                        static_cast<std::size_t>(m_NumAtoms),
                    ApoCharmmErrorCode::InvalidArgument,
                    "Reference-coordinate array size mismatch; expected " +
                        std::to_string(m_NumAtoms) + ", observed " +
                        std::to_string(referenceCoordinates.size()));

  for (int i = 0; i < m_NumAtoms; i++) {
    APOCHARMM_REQUIRE(std::isfinite(referenceCoordinates[i].x),
                      ApoCharmmErrorCode::InvalidArgument,
                      "Reference coordinate at atom index " +
                          std::to_string(i) +
                          ", X component must be finite; observed " +
                          std::to_string(referenceCoordinates[i].x));

    APOCHARMM_REQUIRE(std::isfinite(referenceCoordinates[i].y),
                      ApoCharmmErrorCode::InvalidArgument,
                      "Reference coordinate at atom index " +
                          std::to_string(i) +
                          ", Y component must be finite; observed " +
                          std::to_string(referenceCoordinates[i].y));

    APOCHARMM_REQUIRE(std::isfinite(referenceCoordinates[i].z),
                      ApoCharmmErrorCode::InvalidArgument,
                      "Reference coordinate at atom index " +
                          std::to_string(i) +
                          ", Z component must be finite; observed " +
                          std::to_string(referenceCoordinates[i].z));

    m_ReferenceCoordinates[i].x = referenceCoordinates[i].x;
    m_ReferenceCoordinates[i].y = referenceCoordinates[i].y;
    m_ReferenceCoordinates[i].z = referenceCoordinates[i].z;
  }
  m_ReferenceCoordinates.transferToDevice();

  return;
}

template <typename AT, typename CT>
void HarmonicRestraintForce<AT, CT>::setReferenceCoordinates(
    const std::vector<std::vector<double>> &referenceCoordinates) {
  APOCHARMM_REQUIRE(referenceCoordinates.size() ==
                        static_cast<std::size_t>(m_NumAtoms),
                    ApoCharmmErrorCode::InvalidArgument,
                    "Reference-coordinate array size mismatch; expected " +
                        std::to_string(m_NumAtoms) + ", observed " +
                        std::to_string(referenceCoordinates.size()));

  std::vector<double3> convertedCoordinates(referenceCoordinates.size());

  for (int i = 0; i < m_NumAtoms; i++) {
    APOCHARMM_REQUIRE(referenceCoordinates[i].size() == 3,
                      ApoCharmmErrorCode::InvalidArgument,
                      "Reference coordinate at atom index " +
                          std::to_string(i) +
                          " has invalid size; expected 3, observed " +
                          std::to_string(referenceCoordinates[i].size()));

    convertedCoordinates[i] =
        make_double3(referenceCoordinates[i][0], referenceCoordinates[i][1],
                     referenceCoordinates[i][2]);
  }

  this->setReferenceCoordinates(convertedCoordinates);

  return;
}

template <typename AT, typename CT>
void HarmonicRestraintForce<AT, CT>::setMasses(
    const std::vector<double> &masses) {
  APOCHARMM_REQUIRE(masses.size() == static_cast<std::size_t>(m_NumAtoms),
                    ApoCharmmErrorCode::InvalidArgument,
                    "Mass array size mismatch; expected " +
                        std::to_string(m_NumAtoms) + ", observed " +
                        std::to_string(masses.size()));

  for (int i = 0; i < m_NumAtoms; i++) {
    APOCHARMM_REQUIRE(
        std::isfinite(masses[i]), ApoCharmmErrorCode::InvalidArgument,
        "Mass at index " + std::to_string(i) + " must be finite; observed " +
            std::to_string(masses[i]));

    APOCHARMM_REQUIRE(masses[i] >= 0.0, ApoCharmmErrorCode::InvalidArgument,
                      "Mass at index " + std::to_string(i) +
                          " must be non-negative; observed " +
                          std::to_string(masses[i]));

    m_ReferenceCoordinates[i].w = masses[i];
  }

  m_ReferenceCoordinates.transferToDevice();

  return;
}

template <typename AT, typename CT>
std::shared_ptr<cudaStream_t> HarmonicRestraintForce<AT, CT>::getStream(void) {
  return m_Stream;
}

template <typename AT, typename CT>
void HarmonicRestraintForce<AT, CT>::initialize(
    const int numAtoms, const std::vector<double> &boxDimensions) {
  APOCHARMM_REQUIRE(m_NumAtoms == numAtoms, ApoCharmmErrorCode::InvalidArgument,
                    "Initialization atom count mismatch; expected " +
                        std::to_string(m_NumAtoms) + ", observed " +
                        std::to_string(numAtoms));

  this->setBoxDimensions(boxDimensions);

  return;
}

template <typename AT, typename CT>
void HarmonicRestraintForce<AT, CT>::clear(void) {
  m_EnergyVirial->clear(*m_Stream);
  m_Forces->clear(*m_Stream);
  return;
}

template <typename AT, typename CT, bool calcEnergy, bool calcVirial>
__global__ static void HarmonicRestraintForceKernel(
    AT *__restrict__ forces, const int forceStride, double *__restrict__ energy,
    const double *__restrict__ forceConstants, const float4 *__restrict__ xyzq,
    const double4 *__restrict__ coordsRef, const int numAtoms) {
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  const int stride = gridDim.x * blockDim.x;

  double epot = 0.0;
  for (int i = index; i < numAtoms; i += stride) {
    const CT kf = static_cast<CT>(forceConstants[i]);

    const CT rx = static_cast<CT>(xyzq[i].x);
    const CT ry = static_cast<CT>(xyzq[i].y);
    const CT rz = static_cast<CT>(xyzq[i].z);

    const CT rx0 = static_cast<CT>(coordsRef[i].x);
    const CT ry0 = static_cast<CT>(coordsRef[i].y);
    const CT rz0 = static_cast<CT>(coordsRef[i].z);
    const CT mass = static_cast<CT>(coordsRef[i].w);

    const CT dx = rx - rx0;
    const CT dy = ry - ry0;
    const CT dz = rz - rz0;
    const CT dr2 = dx * dx + dy * dy + dz * dz;

    if (calcEnergy == true)
      epot += static_cast<double>(kf * mass * dr2);

    const CT forceCoef = static_cast<CT>(2.0) * kf * mass;

    AT fx = static_cast<AT>(0), fy = static_cast<AT>(0),
       fz = static_cast<AT>(0);
    calc_component_force<AT, CT>(forceCoef, dx, dy, dz, fx, fy, fz);

    // Store forces
    write_force<AT>(fx, fy, fz, i, forceStride, forces);

    // JEG260501: We can add this calculation later if we implement the external
    // virial. For now, we only have the internal virial, and absolute harmonic
    // restraints only contribute to the external virial.
    //   if (calcVirial) {
    // #ifdef USE_DP_SFORCE
    //     atomicAdd((unsigned long long int *)&virial->sforce_dp[][0],
    //               static_cast<double>(f * dx));
    //     atomicAdd((unsigned long long int *)&virial->sforce_dp[][0],
    //               static_cast<double>(f * dy));
    //     atomicAdd((unsigned long long int *)&virial->sforce_dp[][0],
    //               static_cast<double>(f * dz));
    // #else
    //     fx /= CONVERT_TO_VIR;
    //     fy /= CONVERT_TO_VIR;
    //     fz /= CONVERT_TO_VIR;
    //     atomicAdd((unsigned long long int *)&virial->sforce_fp[][0],
    //     llitoulli(fx)); atomicAdd((unsigned long long int
    //     *)&virial->sforce_fp[][0], llitoulli(fy)); atomicAdd((unsigned long
    //     long int *)&virial->sforce_fp[][0], llitoulli(fz));
    // #endif
    //   }
  }

  if (calcEnergy == true) {
    epot = BlockReduceSum<double>(epot);
    if (threadIdx.x == 0)
      atomicAdd(energy, epot);
  }

  return;
}

template <typename AT, typename CT>
void HarmonicRestraintForce<AT, CT>::calcForce(const float4 *xyzq,
                                               const bool calcEnergy,
                                               const bool calcVirial) {
  constexpr int numThreads = 256;
  const int numBlocks = (m_NumAtoms + numThreads - 1) / numThreads;

  if ((calcEnergy == true) && (calcVirial == true)) {
    cudaCheckLaunch(
        HarmonicRestraintForceKernel<AT, CT, true, true>
        <<<numBlocks, numThreads, 0, *m_Stream>>>(
            m_Forces->xyz(), m_Forces->stride(),
            m_EnergyVirial->getEnergyPointer("harm"),
            m_ForceConstants.getDeviceArray().data(), xyzq,
            m_ReferenceCoordinates.getDeviceArray().data(), m_NumAtoms));
  } else if ((calcEnergy == true) && (calcVirial == false)) {
    cudaCheckLaunch(
        HarmonicRestraintForceKernel<AT, CT, true, false>
        <<<numBlocks, numThreads, 0, *m_Stream>>>(
            m_Forces->xyz(), m_Forces->stride(),
            m_EnergyVirial->getEnergyPointer("harm"),
            m_ForceConstants.getDeviceArray().data(), xyzq,
            m_ReferenceCoordinates.getDeviceArray().data(), m_NumAtoms));
  } else if ((calcEnergy == false) && (calcVirial == true)) {
    cudaCheckLaunch(
        HarmonicRestraintForceKernel<AT, CT, false, true>
        <<<numBlocks, numThreads, 0, *m_Stream>>>(
            m_Forces->xyz(), m_Forces->stride(),
            m_EnergyVirial->getEnergyPointer("harm"),
            m_ForceConstants.getDeviceArray().data(), xyzq,
            m_ReferenceCoordinates.getDeviceArray().data(), m_NumAtoms));
  } else if ((calcEnergy == false) && (calcVirial == false)) {
    cudaCheckLaunch(
        HarmonicRestraintForceKernel<AT, CT, false, false>
        <<<numBlocks, numThreads, 0, *m_Stream>>>(
            m_Forces->xyz(), m_Forces->stride(),
            m_EnergyVirial->getEnergyPointer("harm"),
            m_ForceConstants.getDeviceArray().data(), xyzq,
            m_ReferenceCoordinates.getDeviceArray().data(), m_NumAtoms));
  }

  return;
}

template <typename AT, typename CT>
void HarmonicRestraintForce<AT, CT>::setBoxDimensions(
    const std::vector<double> &boxDimensions) {
  APOCHARMM_REQUIRE(boxDimensions.size() == 3,
                    ApoCharmmErrorCode::InvalidArgument,
                    "Box-dimension array size mismatch; expected 3, observed " +
                        std::to_string(boxDimensions.size()));

  for (std::size_t i = 0; i < 3; i++) {
    APOCHARMM_REQUIRE(
        std::isfinite(boxDimensions[i]), ApoCharmmErrorCode::InvalidArgument,
        "Box dimension at index " + std::to_string(i) +
            " must be finite; observed " + std::to_string(boxDimensions[i]));

    APOCHARMM_REQUIRE(
        boxDimensions[i] > 0.0, ApoCharmmErrorCode::InvalidArgument,
        "Box dimension at index " + std::to_string(i) +
            " must be positive; observed " + std::to_string(boxDimensions[i]));
  }

  if ((m_BoxDimensions[0] == boxDimensions[0]) &&
      (m_BoxDimensions[1] == boxDimensions[1]) &&
      (m_BoxDimensions[2] == boxDimensions[2]))
    return;

  // JEG260615: After some thought, I am not entirely convinced this is correct.
  // So I have commented it out.
  /* *
  if ((m_BoxDimensions[0] != 0.0) && (m_BoxDimensions[1] != 0.0) &&
      (m_BoxDimensions[2] != 0.0)) {
    // JEG260512: This only works for cubic, tetragonal, and orthorhombic boxes.
    // If we ever want to do other crystals (e.g. triclinic) we would need to
    // update all of the code to use the A, B, C, alpha, beta, gamma matrix.
    const double dx = boxDimensions[0] / m_BoxDimensions[0];
    const double dy = boxDimensions[1] / m_BoxDimensions[1];
    const double dz = boxDimensions[2] / m_BoxDimensions[2];

    for (int i = 0; i < m_NumAtoms; i++) {
      m_ReferenceCoordinates[i].x *= dx;
      m_ReferenceCoordinates[i].y *= dy;
      m_ReferenceCoordinates[i].z *= dz;
    }
    m_ReferenceCoordinates.transferToDevice();
  }
  * */

  m_BoxDimensions = boxDimensions;

  return;
}

template <typename AT, typename CT>
std::shared_ptr<Force<AT>> HarmonicRestraintForce<AT, CT>::getForce(void) {
  return m_Forces;
}

template <typename AT, typename CT>
std::shared_ptr<CudaEnergyVirial>
HarmonicRestraintForce<AT, CT>::getEnergyVirial(void) {
  return m_EnergyVirial;
}

//
// Explicit instances of HarmonicRestraintForce
//
template class HarmonicRestraintForce<long long int, float>;
template class HarmonicRestraintForce<long long int, double>;

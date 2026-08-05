// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author:  James E. Gonzales II
//
// ENDLICENSE

#include "HarmonicCenterOfMassRestraintForce.h"

#include "cuda_utils.h"
#include "gpu_utils.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

template <typename AT, typename CT>
HarmonicCenterOfMassRestraintForce<AT, CT>::HarmonicCenterOfMassRestraintForce(
    const int numAtoms)
    : m_NumAtoms(numAtoms), m_NumSelected(numAtoms), m_ForceConstant(0.0),
      m_ReferenceDistance(0.0), m_UseMassWeighting(false),
      m_AtomIndices(numAtoms), m_AtomWeights(numAtoms), m_Masses(numAtoms),
      m_PartialSums(1), m_RestraintState(2), m_BoxDimX(0.0), m_BoxDimY(0.0),
      m_BoxDimZ(0.0), m_ReferencePosition(make_double4(0.0, 0.0, 0.0, 0.0)),
      m_ReferenceMaskX(1), m_ReferenceMaskY(1), m_ReferenceMaskZ(1),
      m_Selection(numAtoms, AtomSelection::InitialValue::ALL),
      m_EnergyVirial(nullptr), m_Forces(nullptr), m_Stream(nullptr) {
  m_Masses.set(1.0);
  m_PartialSums.set(make_double4(0.0, 0.0, 0.0, 0.0));
  m_RestraintState.set(make_double4(0.0, 0.0, 0.0, 0.0));

  this->updateSelectedAtoms();

  m_EnergyVirial = std::make_shared<CudaEnergyVirial>();
  m_EnergyVirial->insert("hmcm");

  m_Forces = std::make_shared<Force<AT>>();
  m_Forces->realloc(numAtoms, 1.5f);

  m_Stream = std::make_shared<cudaStream_t>();
  cudaCheck(cudaStreamCreate(m_Stream.get()));
}

template <typename AT, typename CT>
HarmonicCenterOfMassRestraintForce<AT, CT>::~HarmonicCenterOfMassRestraintForce(
    void) noexcept {
  if (m_Stream != nullptr) {
    destroy_cuda_stream_noexcept(m_Stream.get());
    m_Stream.reset();
  }
}

template <typename AT, typename CT>
void HarmonicCenterOfMassRestraintForce<AT, CT>::setSelection(
    const AtomSelection &selection) {
  if (selection.getNumAtoms() != m_NumAtoms) {
    throw std::invalid_argument("HarmonicCenterOfMassRestraintForce selection "
                                "has incorrect number of atoms");
  }

  if (selection.getNumSelected() <= 0) {
    throw std::invalid_argument("HarmonicCenterOfMassRestraintForce selection "
                                "must contain at least one atom");
  }

  m_Selection = selection;
  this->updateSelectedAtoms();

  return;
}

template <typename AT, typename CT>
void HarmonicCenterOfMassRestraintForce<AT, CT>::setForceConstant(
    const double forceConstant) {
  if (!std::isfinite(forceConstant)) {
    throw std::invalid_argument(
        "HarmonicCenterOfMassRestraintForce force constant must be finite");
  }

  if (forceConstant < 0.0) {
    throw std::invalid_argument("HarmonicCenterOfMassRestraintForce force "
                                "constant must be non-negative");
  }

  m_ForceConstant = forceConstant;

  return;
}

template <typename AT, typename CT>
void HarmonicCenterOfMassRestraintForce<AT, CT>::setReferencePosition(
    const std::vector<double> &referencePosition) {
  this->setReferencePosition(referencePosition, {1, 1, 1});
  return;
}

template <typename AT, typename CT>
void HarmonicCenterOfMassRestraintForce<AT, CT>::setReferencePosition(
    const std::vector<double> &referencePosition,
    const std::vector<int> &referenceMask) {
  if (referencePosition.size() != 3) {
    std::string msg =
        "ERROR: HarmonicCenterOfMassRestraintForce::setReferencePosition: "
        "referencePosition must contain exactly 3 elements\n";
    msg += "referencePosition.size() = " +
           std::to_string(referencePosition.size()) + "\n";
    throw std::invalid_argument(msg);
  }

  if (referenceMask.size() != 3) {
    std::string msg =
        "ERROR: HarmonicCenterOfMassRestraintForce::setReferencePosition: "
        "referenceMask must contain exactly 3 elements\n";
    msg +=
        "referenceMask.size() = " + std::to_string(referenceMask.size()) + "\n";
    throw std::invalid_argument(msg);
  }

  int numActiveComponents = 0;
  for (std::size_t i = 0; i < 3; i++) {
    if (!std::isfinite(referencePosition[i])) {
      throw std::invalid_argument("HarmonicCenterOfMassRestraintForce "
                                  "reference position must be finite");
    }

    if ((referenceMask[i] != 0) && (referenceMask[i] != 1)) {
      throw std::invalid_argument("HarmonicCenterOfMassRestraintForce "
                                  "reference mask values must be 0 or 1");
    }

    numActiveComponents += referenceMask[i];
  }

  if (numActiveComponents == 0) {
    throw std::invalid_argument("HarmonicCenterOfMassRestraintForce reference "
                                "mask must activate at least one coordinate");
  }

  m_ReferencePosition.x = referencePosition[0];
  m_ReferencePosition.y = referencePosition[1];
  m_ReferencePosition.z = referencePosition[2];

  m_ReferenceMaskX = referenceMask[0];
  m_ReferenceMaskY = referenceMask[1];
  m_ReferenceMaskZ = referenceMask[2];

  return;
}

template <typename AT, typename CT>
void HarmonicCenterOfMassRestraintForce<AT, CT>::setReferenceDistance(
    const double referenceDistance) {
  if (!std::isfinite(referenceDistance)) {
    throw std::invalid_argument(
        "HarmonicCenterOfMassRestraintForce reference distance must be finite");
  }

  if (referenceDistance < 0.0) {
    throw std::invalid_argument("HarmonicCenterOfMassRestraintForce reference "
                                "distance must be non-negative");
  }

  m_ReferenceDistance = referenceDistance;

  return;
}

template <typename AT, typename CT>
void HarmonicCenterOfMassRestraintForce<AT, CT>::setMasses(
    const std::vector<double> &masses) {
  if (masses.size() != static_cast<std::size_t>(m_NumAtoms)) {
    std::string msg =
        "ERROR: HarmonicCenterOfMassRestraintForce::setMasses(const "
        "std::vector<double> &): Size of input vector must match the total "
        "number of atoms\n";
    msg += "        NATOM = " + std::to_string(m_NumAtoms) + "\n";
    msg += "masses.size() = " + std::to_string(masses.size()) + "\n";
    throw std::invalid_argument(msg);
  }

  for (std::size_t i = 0; i < masses.size(); i++) {
    if (!std::isfinite(masses[i])) {
      throw std::invalid_argument(
          "HarmonicCenterOfMassRestraintForce masses must be finite");
    }

    if (masses[i] < 0.0) {
      throw std::invalid_argument(
          "HarmonicCenterOfMassRestraintForce masses must be non-negative");
    }
  }

  m_Masses = masses;
  m_UseMassWeighting = true;
  this->updateSelectedAtoms();

  return;
}

template <typename AT, typename CT>
void HarmonicCenterOfMassRestraintForce<AT, CT>::setMassWeighting(
    const bool useMassWeighting) {
  m_UseMassWeighting = useMassWeighting;
  this->updateSelectedAtoms();
  return;
}

template <typename AT, typename CT>
void HarmonicCenterOfMassRestraintForce<AT, CT>::setBoxDimensions(
    const std::vector<double> &boxDimensions) {
  if (boxDimensions.size() != 3) {
    std::string msg =
        "ERROR: HarmonicCenterOfMassRestraintForce::setBoxDimensions(const "
        "std::vector<double> &): boxDimensions must contain exactly 3 "
        "elements\n";
    msg +=
        "boxDimensions.size() = " + std::to_string(boxDimensions.size()) + "\n";
    throw std::invalid_argument(msg);
  }

  for (std::size_t i = 0; i < 3; i++) {
    if (!std::isfinite(boxDimensions[i])) {
      throw std::invalid_argument(
          "HarmonicCenterOfMassRestraintForce box dimensions must be finite");
    }

    if (boxDimensions[i] <= 0.0) {
      throw std::invalid_argument(
          "HarmonicCenterOfMassRestraintForce box dimensions must be positive");
    }
  }

  if ((m_BoxDimX == boxDimensions[0]) && (m_BoxDimY == boxDimensions[1]) &&
      (m_BoxDimZ == boxDimensions[2]))
    return;

  // JEG260615: After some thought, I am not entirely convinced this is correct.
  // So I have commented it out.
  /* *
  if ((m_BoxDimX != 0.0) && (m_BoxDimY != 0.0) && (m_BoxDimZ != 0.0)) {
    // JEG260512: This only works for cubic, tetragonal, and orthorhombic boxes.
    // If we ever want to do other crystals (e.g. triclinic) we would need to
    // update all of the code to use the A, B, C, alpha, beta, gamma matrix.
    const double dx = boxDimensions[0] / m_BoxDimX;
    const double dy = boxDimensions[1] / m_BoxDimY;
    const double dz = boxDimensions[2] / m_BoxDimZ;

    if (m_ReferenceMaskX != 0)
      m_ReferencePosition.x *= dx;
    if (m_ReferenceMaskY != 0)
      m_ReferencePosition.y *= dy;
    if (m_ReferenceMaskZ != 0)
      m_ReferencePosition.z *= dz;
  }
  * */

  m_BoxDimX = boxDimensions[0];
  m_BoxDimY = boxDimensions[1];
  m_BoxDimZ = boxDimensions[2];

  return;
}

template <typename AT, typename CT>
std::shared_ptr<CudaEnergyVirial>
HarmonicCenterOfMassRestraintForce<AT, CT>::getEnergyVirial(void) {
  return m_EnergyVirial;
}

template <typename AT, typename CT>
std::shared_ptr<Force<AT>>
HarmonicCenterOfMassRestraintForce<AT, CT>::getForce(void) {
  return m_Forces;
}

template <typename AT, typename CT>
std::shared_ptr<cudaStream_t>
HarmonicCenterOfMassRestraintForce<AT, CT>::getStream(void) {
  return m_Stream;
}

template <typename AT, typename CT>
void HarmonicCenterOfMassRestraintForce<AT, CT>::initialize(
    const int numAtoms, const std::vector<double> &boxDimensions) {
  if (m_NumAtoms != numAtoms) {
    throw std::runtime_error("HarmonicCenterOfMassRestraintForce: Attempted to "
                             "initialize with different atom count");
  }

  this->setBoxDimensions(boxDimensions);

  return;
}

template <typename AT, typename CT>
void HarmonicCenterOfMassRestraintForce<AT, CT>::clear(void) {
  m_EnergyVirial->clear(*m_Stream);
  m_Forces->clear(*m_Stream);
  return;
}

__global__ static void PartialSumsKernel(double4 *__restrict__ partialSums,
                                         const int *__restrict__ atomIndices,
                                         const double *__restrict__ atomWeights,
                                         const float4 *__restrict__ xyzq,
                                         const int numSelected,
                                         const double boxx, const double boxy,
                                         const double boxz) {
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  const int stride = gridDim.x * blockDim.x;

  const int anchorAtomIndex = atomIndices[0];
  const float4 xyzqAnchor = xyzq[anchorAtomIndex];

  const float3 halfBox = make_float3(static_cast<float>(0.5 * boxx),
                                     static_cast<float>(0.5 * boxy),
                                     static_cast<float>(0.5 * boxz));

  double sumX = 0.0, sumY = 0.0, sumZ = 0.0, sumW = 0.0;
  for (int i = index; i < numSelected; i += stride) {
    const int atomIndex = atomIndices[i];
    const double atomWeight = atomWeights[i];

    const float4 xyzqi = xyzq[atomIndex];

    // Calculate periodic shift for atoms
    const int ish = calc_ishift(xyzqi, xyzqAnchor, halfBox);
    double shx = 0.0, shy = 0.0, shz = 0.0;
    calc_box_shift<double>(ish, boxx, boxy, boxz, shx, shy, shz);

    sumX += atomWeight * (static_cast<double>(xyzqi.x) + shx);
    sumY += atomWeight * (static_cast<double>(xyzqi.y) + shy);
    sumZ += atomWeight * (static_cast<double>(xyzqi.z) + shz);
    sumW += atomWeight;
  }

  sumX = BlockReduceSum<double>(sumX);
  sumY = BlockReduceSum<double>(sumY);
  sumZ = BlockReduceSum<double>(sumZ);
  sumW = BlockReduceSum<double>(sumW);

  if (threadIdx.x == 0)
    partialSums[blockIdx.x] = make_double4(sumX, sumY, sumZ, sumW);

  return;
}

template <bool calcEnergy>
__global__ static void
StateKernel(double4 *__restrict__ restraintState, double *__restrict__ energy,
            const double4 *__restrict__ partialSums, const int numPartialSums,
            const double forceConstant, const double4 referencePosition,
            const int referenceMaskX, const int referenceMaskY,
            const int referenceMaskZ, const double referenceDistance,
            const double boxx, const double boxy, const double boxz) {
  double sumX = 0.0, sumY = 0.0, sumZ = 0.0, sumW = 0.0;
  for (int i = threadIdx.x; i < numPartialSums; i += blockDim.x) {
    const double4 partialSum = partialSums[i];
    sumX += partialSum.x;
    sumY += partialSum.y;
    sumZ += partialSum.z;
    sumW += partialSum.w;
  }
  sumX = BlockReduceSum<double>(sumX);
  sumY = BlockReduceSum<double>(sumY);
  sumZ = BlockReduceSum<double>(sumZ);
  sumW = BlockReduceSum<double>(sumW);

  if (threadIdx.x == 0) {
    const double invTotW = 1.0 / sumW;
    const double centerX = sumX * invTotW;
    const double centerY = sumY * invTotW;
    const double centerZ = sumZ * invTotW;

    int centerShiftCountX = 0, centerShiftCountY = 0, centerShiftCountZ = 0;
    double dx = 0.0, dy = 0.0, dz = 0.0;

    if (referenceMaskX != 0) {
      centerShiftCountX = -static_cast<int>(
          floor((centerX - referencePosition.x) / boxx + 0.5));
      dx = centerX + static_cast<double>(centerShiftCountX) * boxx -
           referencePosition.x;
    }

    if (referenceMaskY != 0) {
      centerShiftCountY = -static_cast<int>(
          floor((centerY - referencePosition.y) / boxy + 0.5));
      dy = centerY + static_cast<double>(centerShiftCountY) * boxy -
           referencePosition.y;
    }

    if (referenceMaskZ != 0) {
      centerShiftCountZ = -static_cast<int>(
          floor((centerZ - referencePosition.z) / boxz + 0.5));
      dz = centerZ + static_cast<double>(centerShiftCountZ) * boxz -
           referencePosition.z;
    }

    const double dr2 = dx * dx + dy * dy + dz * dz;
    double epot = 0.0, gradX = 0.0, gradY = 0.0, gradZ = 0.0;
    double forceCoef = 2.0 * forceConstant;

    if (referenceDistance == 0.0)
      epot = forceConstant * dr2;
    else {
      const double distance = sqrt(dr2);
      const double delta = distance - referenceDistance;

      epot = forceConstant * delta * delta;

      if (distance > 0.0)
        forceCoef *= (delta / distance);
    }
    gradX = forceCoef * dx;
    gradY = forceCoef * dy;
    gradZ = forceCoef * dz;

    restraintState[0] = make_double4(gradX, gradY, gradZ, invTotW);
    restraintState[1] =
        make_double4(static_cast<double>(centerShiftCountX),
                     static_cast<double>(centerShiftCountY),
                     static_cast<double>(centerShiftCountZ), 0.0);

    if (calcEnergy == true)
      atomicAdd(energy, epot);
  }

  return;
}

template <typename AT, typename CT, bool calcVirial>
__global__ static void
ApplyForcesKernel(Virial_t *__restrict__ virial, AT *__restrict__ forces,
                  const int forceStride, const int *__restrict__ atomIndices,
                  const double *__restrict__ atomWeights, const int numSelected,
                  const float4 *__restrict__ xyzq,
                  const double4 *__restrict__ restraintState, const double boxx,
                  const double boxy, const double boxz) {
  const int index = blockIdx.x * blockDim.x + threadIdx.x;
  const int stride = gridDim.x * blockDim.x;

  // Gradient state
  const CT gradX = static_cast<CT>(restraintState[0].x);
  const CT gradY = static_cast<CT>(restraintState[0].y);
  const CT gradZ = static_cast<CT>(restraintState[0].z);
  const CT invTotW = static_cast<CT>(restraintState[0].w);

  // Center periodic image shift relative to the reference
  const int centerShiftCountX = static_cast<int>(restraintState[1].x);
  const int centerShiftCountY = static_cast<int>(restraintState[1].y);
  const int centerShiftCountZ = static_cast<int>(restraintState[1].z);

  const int anchorAtomIndex = atomIndices[0];
  const float4 xyzqAnchor = xyzq[anchorAtomIndex];

  const float3 halfBox = make_float3(static_cast<float>(0.5 * boxx),
                                     static_cast<float>(0.5 * boxy),
                                     static_cast<float>(0.5 * boxz));

  // JEG260617: This assumes the center-reference minimum image shift is within
  // one periodic image in each dimension
  int centerIsh = 13;
  const bool validCenterShift =
      (centerShiftCountX >= -1) && (centerShiftCountX <= 1) &&
      (centerShiftCountY >= -1) && (centerShiftCountY <= 1) &&
      (centerShiftCountZ >= -1) && (centerShiftCountZ <= 1);

  if (validCenterShift == true) {
    centerIsh = (centerShiftCountX + 1) + 3 * (centerShiftCountY + 1) +
                9 * (centerShiftCountZ + 1);
  }

  for (int i = index; i < numSelected; i += stride) {
    const int atomIndex = atomIndices[i];
    const CT weightFraction = static_cast<CT>(atomWeights[i]) * invTotW;

    AT fx = static_cast<AT>(0);
    AT fy = static_cast<AT>(0);
    AT fz = static_cast<AT>(0);

    calc_component_force<AT, CT>(weightFraction, gradX, gradY, gradZ, fx, fy,
                                 fz);

    // Store forces
    write_force<AT>(fx, fy, fz, atomIndex, forceStride, forces);

    if (calcVirial == true) {
      const float4 xyzqi = xyzq[atomIndex];
      const int atomIsh = calc_ishift(xyzqi, xyzqAnchor, halfBox);

      const CT atomGradX = weightFraction * gradX;
      const CT atomGradY = weightFraction * gradY;
      const CT atomGradZ = weightFraction * gradZ;

#ifdef USE_DP_SFORCE
      if (atomIsh != 13) {
        atomicAdd(&virial->sforce_dp[atomIsh][0],
                  static_cast<double>(atomGradX));
        atomicAdd(&virial->sforce_dp[atomIsh][1],
                  static_cast<double>(atomGradY));
        atomicAdd(&virial->sforce_dp[atomIsh][2],
                  static_cast<double>(atomGradZ));
      }

      if ((validCenterShift == true) && (centerIsh != 13)) {
        atomicAdd(&virial->sforce_dp[centerIsh][0],
                  static_cast<double>(atomGradX));
        atomicAdd(&virial->sforce_dp[centerIsh][1],
                  static_cast<double>(atomGradY));
        atomicAdd(&virial->sforce_dp[centerIsh][2],
                  static_cast<double>(atomGradZ));
      }
#else
      const AT virialForceX = fx / CONVERT_TO_VIR;
      const AT virialForceY = fy / CONVERT_TO_VIR;
      const AT virialForceZ = fz / CONVERT_TO_VIR;

      if (atomIsh != 13) {
        atomicAdd((unsigned long long int *)&virial->sforce_fp[atomIsh][0],
                  llitoulli(virialForceX));
        atomicAdd((unsigned long long int *)&virial->sforce_fp[atomIsh][1],
                  llitoulli(virialForceY));
        atomicAdd((unsigned long long int *)&virial->sforce_fp[atomIsh][2],
                  llitoulli(virialForceZ));
      }

      if ((validCenterShift == true) && (centerIsh != 13)) {
        atomicAdd((unsigned long long int *)&virial->sforce_fp[centerIsh][0],
                  llitoulli(virialForceX));
        atomicAdd((unsigned long long int *)&virial->sforce_fp[centerIsh][1],
                  llitoulli(virialForceY));
        atomicAdd((unsigned long long int *)&virial->sforce_fp[centerIsh][2],
                  llitoulli(virialForceZ));
      }
#endif
    }
  }

  return;
}

template <typename AT, typename CT>
void HarmonicCenterOfMassRestraintForce<AT, CT>::calcForce(
    const float4 *xyzq, const bool calcEnergy, const bool calcVirial) {
  if ((m_NumSelected <= 0) || (m_ForceConstant == 0.0))
    return;

  if ((m_BoxDimX <= 0.0) || (m_BoxDimY <= 0.0) || (m_BoxDimZ <= 0.0)) {
    throw std::runtime_error("HarmonicCenterOfMassRestraintForce::calcForce "
                             "requires positive box dimensions");
  }

  constexpr int numThreads = 256;
  const int numBlocks = (m_NumSelected + numThreads - 1) / numThreads;
  const int numPartialSums = static_cast<int>(m_PartialSums.size());

  PartialSumsKernel<<<numBlocks, numThreads, 0, *m_Stream>>>(
      m_PartialSums.getDeviceArray().data(),
      m_AtomIndices.getDeviceArray().data(),
      m_AtomWeights.getDeviceArray().data(), xyzq, m_NumSelected, m_BoxDimX,
      m_BoxDimY, m_BoxDimZ);
  cudaCheck(cudaGetLastError());

  if (calcEnergy == true) {
    StateKernel<true><<<1, numThreads, 0, *m_Stream>>>(
        m_RestraintState.getDeviceArray().data(),
        m_EnergyVirial->getEnergyPointer("hmcm"),
        m_PartialSums.getDeviceArray().data(), numPartialSums, m_ForceConstant,
        m_ReferencePosition, m_ReferenceMaskX, m_ReferenceMaskY,
        m_ReferenceMaskZ, m_ReferenceDistance, m_BoxDimX, m_BoxDimY, m_BoxDimZ);
  } else {
    StateKernel<false><<<1, numThreads, 0, *m_Stream>>>(
        m_RestraintState.getDeviceArray().data(),
        m_EnergyVirial->getEnergyPointer("hmcm"),
        m_PartialSums.getDeviceArray().data(), numPartialSums, m_ForceConstant,
        m_ReferencePosition, m_ReferenceMaskX, m_ReferenceMaskY,
        m_ReferenceMaskZ, m_ReferenceDistance, m_BoxDimX, m_BoxDimY, m_BoxDimZ);
  }
  cudaCheck(cudaGetLastError());

  if (calcVirial == true) {
    ApplyForcesKernel<AT, CT, true><<<numBlocks, numThreads, 0, *m_Stream>>>(
        m_EnergyVirial->getVirialPointer(), m_Forces->xyz(), m_Forces->stride(),
        m_AtomIndices.getDeviceArray().data(),
        m_AtomWeights.getDeviceArray().data(), m_NumSelected, xyzq,
        m_RestraintState.getDeviceArray().data(), m_BoxDimX, m_BoxDimY,
        m_BoxDimZ);
  } else {
    ApplyForcesKernel<AT, CT, false><<<numBlocks, numThreads, 0, *m_Stream>>>(
        nullptr, m_Forces->xyz(), m_Forces->stride(),
        m_AtomIndices.getDeviceArray().data(),
        m_AtomWeights.getDeviceArray().data(), m_NumSelected, xyzq,
        m_RestraintState.getDeviceArray().data(), m_BoxDimX, m_BoxDimY,
        m_BoxDimZ);
  }
  cudaCheck(cudaGetLastError());

  return;
}

template <typename AT, typename CT>
void HarmonicCenterOfMassRestraintForce<AT, CT>::updateSelectedAtoms(void) {
  const std::vector<int> selectedAtomIndices = m_Selection.getAtomIndices();

  if (selectedAtomIndices.empty()) {
    throw std::invalid_argument("HarmonicCenterOfMassRestraintForce selection "
                                "must contain at least one atom");
  }

  m_NumSelected = static_cast<int>(selectedAtomIndices.size());

  std::vector<double> selectedAtomWeights(selectedAtomIndices.size());
  double totalWeight = 0.0;
  for (std::size_t i = 0; i < selectedAtomIndices.size(); i++) {
    const int atomIndex = selectedAtomIndices[i];
    const double weight =
        (m_UseMassWeighting == true) ? m_Masses[atomIndex] : 1.0;
    selectedAtomWeights[i] = weight;
    totalWeight += weight;
  }

  if (totalWeight <= 0.0) {
    throw std::invalid_argument("HarmonicCenterOfMassRestraintForce selected "
                                "atoms must have positive total weight");
  }

  m_AtomIndices = selectedAtomIndices;
  m_AtomWeights = selectedAtomWeights;

  constexpr int numThreads = 256;
  const int numBlocks = (m_NumSelected + numThreads - 1) / numThreads;

  m_PartialSums.resize(numBlocks);
  m_PartialSums.set(make_double4(0.0, 0.0, 0.0, 0.0));

  m_RestraintState.set(make_double4(0.0, 0.0, 0.0, 0.0));

  return;
}

//
// Explicit instances of HarmonicCenterOfMassRestraintForce
//
template class HarmonicCenterOfMassRestraintForce<long long int, float>;
template class HarmonicCenterOfMassRestraintForce<long long int, double>;

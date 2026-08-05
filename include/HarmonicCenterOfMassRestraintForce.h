// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author:  James E. Gonzales II
//
// ENDLICENSE

#pragma once

#include "AtomSelection.h"
#include "CudaContainer.h"
#include "CudaEnergyVirial.h"
#include "Force.h"

#include <cuda_runtime.h>
#include <memory>
#include <vector>

//
// Calculates CONS HMCM-style harmonic restraints on the center of a selection
// of atoms.
//

template <typename AT, typename CT> class HarmonicCenterOfMassRestraintForce {
public:
  static constexpr bool contributesVirial = true;

public:
  HarmonicCenterOfMassRestraintForce(void) = delete;
  HarmonicCenterOfMassRestraintForce(const int numAtoms);
  ~HarmonicCenterOfMassRestraintForce(void) noexcept;

public:
  void setSelection(const AtomSelection &selection);
  void setForceConstant(const double forceConstant);

  void setReferencePosition(const std::vector<double> &referencePosition);

  void setReferencePosition(const std::vector<double> &referencePosition,
                            const std::vector<int> &referenceMask);

  void setReferenceDistance(const double referenceDistance);

  void setMasses(const std::vector<double> &masses);
  void setMassWeighting(const bool useMassWeighting);

  void setBoxDimensions(const std::vector<double> &boxDimensions);

public:
  std::shared_ptr<CudaEnergyVirial> getEnergyVirial(void);
  std::shared_ptr<Force<AT>> getForce(void);
  std::shared_ptr<cudaStream_t> getStream(void);

public:
  void initialize(const int numAtoms, const std::vector<double> &boxDimensions);
  void clear(void);
  void calcForce(const float4 *xyzq, const bool calcEnergy,
                 const bool calcVirial);

private:
  void updateSelectedAtoms(void);

private:
  int m_NumAtoms;
  int m_NumSelected;

  double m_ForceConstant;
  double m_ReferenceDistance;
  bool m_UseMassWeighting;

  CudaContainer<int> m_AtomIndices;
  CudaContainer<double> m_AtomWeights;
  CudaContainer<double> m_Masses;

  CudaContainer<double4> m_PartialSums;

  // m_RestraintState[0] = {dU/dCx, dU/dCy, dU/dCz, 1 / totalWeight}
  // m_RestraintState[1] = {centerImageShiftX, centerImageShiftY,
  //                        centerImageShiftZ, 0}
  CudaContainer<double4> m_RestraintState;

  double m_BoxDimX;
  double m_BoxDimY;
  double m_BoxDimZ;

  double4 m_ReferencePosition;
  int m_ReferenceMaskX;
  int m_ReferenceMaskY;
  int m_ReferenceMaskZ;

  AtomSelection m_Selection;

  std::shared_ptr<CudaEnergyVirial> m_EnergyVirial;
  std::shared_ptr<Force<AT>> m_Forces;
  std::shared_ptr<cudaStream_t> m_Stream;
};

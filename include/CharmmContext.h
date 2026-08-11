// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: Samarjeet Prasad, James E. Gonzales II
//
// ENDLICENSE

#pragma once

#include "CharmmPSF.h"
#include "CharmmParameters.h"
#include "Coordinates.h"
#include "CudaContainer.h"
#include "Force.h"
#include "ForceManager.h"
#include "PBC.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <vector_types.h>

class CharmmContext : public std::enable_shared_from_this<CharmmContext> {
public:
  CharmmContext(void);
  CharmmContext(std::shared_ptr<CharmmPSF> psf,
                std::shared_ptr<CharmmParameters> prm);
  CharmmContext(std::shared_ptr<ForceManager> fm);
  ~CharmmContext(void) noexcept = default;

  CharmmContext(const CharmmContext &other);

public: // Setters
  void setPrm(std::shared_ptr<CharmmParameters> prm);
  void setPsf(std::shared_ptr<CharmmPSF> psf);
  void setForceManager(std::shared_ptr<ForceManager> fm);

  void setNumAtoms(const int numAtoms);

  void setCoordinatesCharges(const std::vector<double4> &coordinatesCharges);
  void setCoordinatesCharges(
      const std::vector<std::vector<double>> &coordinatesCharges);
  void setCoordinatesCharges(const std::vector<double> &coordinatesCharges);
  void setCoordinates(const std::vector<double3> &coordinates);
  void setCoordinates(const std::vector<std::vector<double>> &coordinates);
  void setCoordinates(const std::vector<double> &coordinates);
  void setCoordinates(const std::shared_ptr<Coordinates> crd);
  void setCharges(const std::vector<double> &charges);

  void setVelocitiesInverseMasses(
      const std::vector<double4> &velocitiesInverseMasses);
  void setVelocitiesInverseMasses(
      const std::vector<std::vector<double>> &velocitiesInverseMasses);
  void setVelocitiesInverseMasses(
      const std::vector<double> &velocitiesInverseMasses);
  void setVelocities(const std::vector<double3> &velocities);
  void setVelocities(const std::vector<std::vector<double>> &velocities);
  void setVelocities(const std::vector<double> &velocities);
  void setVelocitiesFromCHARMMVelocityFile(const std::string &fileName);
  void setMasses(const std::vector<double> &masses);

  void setTemperature(const double temperature);
  void setPeriodicBoundaryCondition(const PBC pbc);
  void setBoxDimensions(const std::vector<double> &boxDimensions);
  void setRandomSeed(const std::uint64_t randomSeed);
  void useHolonomicConstraints(const bool usingHolonomicConstraints);

  void setKappa(const float kappa);
  void setCutoff(const float cutoff);
  void setCtonnb(const float ctonnb);
  void setCtofnb(const float ctofnb);
  void setFFTGrid(const int nfftx, const int nffty, const int nfftz);
  void setPmeSplineOrder(const int pmeSplineOrder);
  void setVdwType(const int vdwType);

public: // Getters
  int getNumAtoms(void) const;
  int getNumDegreesOfFreedom(void) const;
  const CudaContainer<float4> &getCoordinatesChargesSP(void) const;
  const CudaContainer<double4> &getCoordinatesChargesDP(void) const;
  const CudaContainer<double4> &getVelocitiesInverseMasses(void) const;
  double getTemperature(void) const;
  PBC getPeriodicBoundaryCondition(void) const;
  const std::vector<double> &getBoxDimensions(void) const;
  std::uint64_t getRandomSeed(void) const;
  bool usingHolonomicConstraints(void) const;
  double getVolume(void) const;
  const CudaContainer<double> &getPressure(void) const;

  const std::vector<Bond> &getBonds(void) const;
  const CudaContainer<int4> &getWaterMolecules(void) const;
  const CudaContainer<int4> &getShakeAtoms(void) const;
  const CudaContainer<float4> &getShakeParams(void) const;

  float getKappa(void) const;
  float getCutoff(void) const;
  float getCtonnb(void) const;
  float getCtofnb(void) const;
  std::vector<int> getFFTGrid(void) const;
  int getPmeSplineOrder(void) const;
  int getVdwType(void) const;
  int getForceStride(void) const;

  CudaContainer<float4> &getCoordinatesChargesSP(void);
  CudaContainer<double4> &getCoordinatesChargesDP(void);
  CudaContainer<double4> &getVelocitiesInverseMasses(void);
  std::vector<double> &getBoxDimensions(void);
  double getKineticEnergy(void);
  CudaContainer<double> &getPressure(void);

  std::shared_ptr<CharmmPSF> getPsf(void);
  std::shared_ptr<CharmmParameters> getPrm(void);
  std::shared_ptr<ForceManager> getForceManager(void);

  CudaContainer<double> &getPotentialEnergy(void);
  float getPotentialEnergies(void);
  std::shared_ptr<Force<double>> getForces(void);
  CudaContainer<double> &getVirial(void);

  std::vector<Bond> &getBonds(void);
  CudaContainer<int4> &getWaterMolecules(void);
  CudaContainer<int4> &getShakeAtoms(void);
  CudaContainer<float4> &getShakeParams(void);

public: // Specialized functions
  void assignVelocitiesAtTemperature(const double temperature);

  double computeTemperature(void);
  void computePressure(void);

  void imageCentering(void);
  void resetNeighborList(void);

  void calculateKineticEnergy(void);
  void calculatePotentialEnergy(const bool reset = false,
                                const bool print = false);
  void calculateForces(const bool reset = false, const bool calcEnergy = false,
                       const bool calcVirial = false);

  void linkBackForceManager(void);

protected:
  void syncStateFromForceManager(void);
  void syncForceManagerFromState(void);

  bool hasCompleteForceManagerState(void) const;
  void initializeForceManagerIfReady(void);

  void requirePsf(void) const;
  void requireForceManager(void) const;
  void requireInitializedForceManager(void) const;

  void finalizeSetupIfReady(void);

protected:
  std::uint64_t m_RandomSeed;

  std::shared_ptr<CharmmPSF> m_Psf;
  std::shared_ptr<CharmmParameters> m_Prm;
  std::vector<double> m_BoxDimensions;
  bool m_HasBoxDimensions;

  std::shared_ptr<ForceManager> m_ForceManager;

  int m_NumAtoms;
  int m_NumDegreesOfFreedom;

  PBC m_Pbc;
  bool m_HasPbc;

  CudaContainer<float4> m_CoordinatesChargesSP;
  CudaContainer<double4> m_CoordinatesChargesDP;
  bool m_HasCoordinates;

  CudaContainer<double4> m_VelocitiesInverseMasses;

  CudaContainer<double> m_KineticEnergy;
  CudaContainer<double> m_Pressure;
  CudaContainer<double> m_VirialKineticEnergyTensor;

  double m_Temperature;
  bool m_UsingHolonomicConstraints;
};

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
  CharmmContext(std::shared_ptr<ForceManager> forceManager);
  ~CharmmContext(void) = default;

  CharmmContext(const CharmmContext &other);

public: // Setters
  void setCoordinates(const std::shared_ptr<Coordinates> crd);
  void setCoordinates(const std::vector<double4> &coordinates);
  void setCoordinates(const std::vector<float> &coordinates);

  void setCoordinatesCharges(const std::vector<float4> &coordinatesCharges);
  void setCoordinatesCharges(
      const std::vector<std::vector<float>> &coordinatesCharges);
  void setCharges(std::vector<float> &charges);

public: // Getters
  int getNumAtoms(void) const;
  const CudaContainer<float4> &getCoordinatesChargesSP(void) const;
  const CudaContainer<double4> &getCoordinatesChargesDP(void) const;

  CudaContainer<float4> &getCoordinatesChargesSP(void);
  CudaContainer<double4> &getCoordinatesChargesDP(void);
  std::vector<std::vector<double>> getCoordinates(void);

public: // Specialized function
public:
  CudaContainer<double4> &getCoordinatesCharges(void);
  const CudaContainer<float4> &getXYZQ(void) const;
  CudaContainer<float4> &getXYZQ(void);

public:
  void setTemperature(const float temperature);
  float getTemperature(void) const;
  float computeTemperature(void);

public:
  void setPeriodicBoundaryCondition(const PBC pbc);
  PBC getPeriodicBoundaryCondition(void) const;
  const std::vector<double> &getBoxDimensions(void) const;
  std::vector<double> &getBoxDimensions(void);
  void setBoxDimensions(const std::vector<double> &boxDimensions);
  double getVolume(void) const;

public: // ForceManager setters and getters
  void setKappa(const float kappa);
  void setCutoff(const float cutoff);
  void setCtonnb(const float ctonnb);
  void setCtofnb(const float ctofnb);
  void setFFTGrid(const int nfftx, const int nffty, const int nfftz);
  void setPmeSplineOrder(const int pmeSplineOrder);
  void setVdwType(const int vdwType);

  float getKappa(void) const;
  float getCutoff(void) const;
  float getCtonnb(void) const;
  float getCtofnb(void) const;
  std::vector<int> getFFTGrid(void) const;
  int getPmeSplineOrder(void) const;
  int getVdwType(void) const;

public:
  void setMasses(const char *fileName);
  void setMasses(const std::vector<double> &masses);
  void setNumAtoms(const int numAtoms);

public:
  void assignVelocitiesAtTemperature(const float temperature);
  void assignVelocitiesFromCHARMMVelocityFile(const std::string &fileName);
  void assignVelocities(const std::vector<double> &velocities);
  void assignVelocities(const std::vector<std::vector<double>> &velocities);
  void setRandomSeedForVelocities(const std::uint64_t randomSeed);
  std::uint64_t getRandomSeedForVelocities(void) const;
  void removeCenterOfMassMotion(void);
  CudaContainer<double4> &getVelocityMass(void);

public:
  void calculateKineticEnergy(void);
  double getKineticEnergy(void);
  CudaContainer<double> getKineticEnergy_(void);

public:
  void resetNeighborList(void);
  void calculatePotentialEnergy(const bool reset = false,
                                const bool print = false);
  void calculateForces(bool reset = false, bool calcEnergy = false,
                       bool calcVirial = false);
  float getPotentialEnergies(void);
  CudaContainer<double> &getPotentialEnergy(void);
  std::shared_ptr<Force<double>> getForces(void);
  CudaContainer<double> &getVirial(void);
  int getForceStride(void) const;
  int *get_loc2glo(void) const;

public:
  void computePressure(void);
  CudaContainer<double> getPressure(void) const;

public:
  std::vector<Bond> getBonds(void);
  int getDegreesOfFreedom(void) const;
  int getNumDegreesOfFreedom(void) const;
  CudaContainer<int4> getWaterMolecules(void);
  CudaContainer<int4> getShakeAtoms(void);
  CudaContainer<float4> getShakeParams(void);
  void useHolonomicConstraints(const bool useConstraints);
  bool isUsingHolonomicConstraints(void) const;

public:
  void imageCentering(void);
  void orient(void);
  void setPsf(std::shared_ptr<CharmmPSF> psf);
  void setPrm(std::shared_ptr<CharmmParameters> prm);
  void setForceManager(std::shared_ptr<ForceManager> forceManager);
  std::shared_ptr<ForceManager> getForceManager(void);
  std::shared_ptr<CharmmPSF> getPsf(void) const;
  std::shared_ptr<CharmmParameters> getPrm(void) const;
  void linkBackForceManager(void);
  void writeCrd(std::string fileName);

protected:
  void syncStateFromForceManager(void);
  void syncForceManagerFromState(void);

  bool hasCompleteForceManagerState(void) const;
  void initializeForceManagerIfReady(void);

  void requirePsf(const std::string &functionName) const;
  void requireForceManager(const std::string &functionName) const;
  void requireInitializedForceManager(const std::string &functionName) const;

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
  CudaContainer<double4> m_VelocitiesInverseMasses;

  CudaContainer<double> m_KineticEnergy;
  CudaContainer<double> m_Pressure;
  CudaContainer<double> m_VirialKineticEnergyTensor;

  float m_Temperature;
  bool m_UsingHolonomicConstraints;
};

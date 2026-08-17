// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II, Samarjeet Prasad
//
// ENDLICENSE

#pragma once

#include "CudaIntegrator.h"

#include "CharmmContext.h"
#include "CudaContainer.h"

#include <cstdint>
#include <curand_kernel.h>
#include <string>
#include <vector>

class CudaLangevinThermostatIntegrator : public CudaIntegrator {
public:
  CudaLangevinThermostatIntegrator(const double timeStep);
  ~CudaLangevinThermostatIntegrator(void) noexcept override;

public:
  void setReferenceTemperature(const double referenceTemperature);
  void setThermostatFriction(const double thermostatFriction);
  void setThermostatRngSeed(const std::uint64_t seed);
  void setRngSequencePos(const unsigned long long int sequencePos);
  void setRngStates(const std::string &rngStateString);
  void resetAverageTemperature(void);

public:
  double getReferenceTemperature(void) const;
  double getThermostatFriction(void) const;
  std::uint64_t getThermostatRngSeed(void) const;
  unsigned long long int getRngSequencePos(void) const;
  std::string getRngStates(void) const;
  int getAverageWindowSize(void) const;
  const CudaContainer<double> &getKineticEnergy(void) const;
  const CudaContainer<double> &getAverageTemperature(void) const;

  CudaContainer<double> &getKineticEnergy(void);
  CudaContainer<double> &getAverageTemperature(void);
  double getInstantaneousTemperature(void);

protected:
  void initializeImpl(void) override;
  void initializeFromRestartFileImpl(const std::string &rstFileName) override;
  void propagateOneStepImpl(void) override;

protected:
  void setRngStateData(const unsigned long long int position,
                       const std::vector<curandStatePhilox4_32_10_t> &states);

  void initializeRng(void);
  void removeCenterOfMassMotion(void);
  void alloc(const int n);
  void dealloc(void);

protected:
  double m_ReferenceTemperature;
  double m_ThermostatFriction;
  double m_ThermostatGamma;

  std::uint64_t m_Seed;
  unsigned long long int m_RngSequencePos;
  curandStatePhilox4_32_10_t *m_RngStates;

  int m_AverageWindowSize;
  CudaContainer<double> m_KineticEnergy;
  CudaContainer<double> m_KineticEnergyPartialSums;
  CudaContainer<double> m_AverageTemperature;
};

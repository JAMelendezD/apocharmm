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

/**
 * @brief Propagates constant-temperature dynamics with Langevin noise.
 *
 * The integrator maintains one Philox RNG state per context atom, applies
 * stochastic position updates, enforces configured holonomic constraints,
 * derives velocities from corrected displacement history, and accumulates two
 * kinetic and temperature estimators.
 *
 * @warning Treat the class as non-copyable. It owns a raw device RNG-state
 * allocation in addition to the base CUDA streams.
 * @see cuda_integrators
 */
class CudaLangevinThermostatIntegrator : public CudaIntegrator {
public:
  /**
   * @brief Constructs an unattached Langevin thermostat integrator.
   *
   * @param[in] timeStep Finite positive time step in picoseconds.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `timeStep` is non-finite or not
   * positive.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if stream or
   * fixed-state allocation fails.
   * @throws std::bad_alloc If host or stream-owner allocation fails.
   * @throws std::system_error If the platform random source used for the
   * initial seed fails.
   *
   * @post The reference temperature is 300 K, friction is zero, a seed has been
   * sampled, the sequence position is zero, and running averages are zero.
   */
  CudaLangevinThermostatIntegrator(const double timeStep);

  /**
   * @brief Destroys the thermostat and its RNG-state allocation.
   *
   * CUDA cleanup failures are discarded and no exception escapes.
   */
  ~CudaLangevinThermostatIntegrator(void) noexcept override;

public:
  /**
   * @brief Sets the reference temperature.
   *
   * @param[in] referenceTemperature Finite non-negative temperature in kelvin.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `referenceTemperature` is
   * non-finite or negative.
   */
  void setReferenceTemperature(const double referenceTemperature);

  /**
   * @brief Sets the thermostat friction coefficient.
   *
   * @param[in] thermostatFriction Finite non-negative friction in inverse
   * picoseconds.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `thermostatFriction` is non-finite
   * or negative.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if changing the
   * value triggers reinitialization and CUDA work fails.
   * @throws std::bad_alloc If triggered reinitialization cannot allocate work
   * storage.
   *
   * @post The dimensionless per-step coefficient is recomputed from the public
   * time step.
   * @note When a context is attached and no propagation steps have yet been
   * requested, this method immediately reruns concrete initialization.
   */
  void setThermostatFriction(const double thermostatFriction);

  /**
   * @brief Sets the thermostat RNG seed.
   *
   * @param[in] seed Dimensionless 64-bit seed.
   *
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if an attached
   * context requires RNG allocation or initialization and CUDA work fails.
   *
   * @post The existing sequence position is preserved. With an attached
   * context, one state per atom is regenerated from the new seed and that
   * position.
   */
  void setThermostatRngSeed(const std::uint64_t seed);

  /**
   * @brief Sets the thermostat RNG sequence position.
   *
   * @param[in] sequencePos Dimensionless Philox sequence position.
   *
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if an attached
   * context requires RNG allocation or initialization and CUDA work fails.
   *
   * @post The current seed is preserved. With an attached context, one state
   * per atom is regenerated and advanced to `sequencePos`.
   */
  void setRngSequencePos(const unsigned long long int sequencePos);

  /**
   * @brief Restores serialized per-atom RNG states.
   *
   * @param[in] rngStateString Borrowed opaque string produced by
   * `getRngStates()`.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::NotInitialized` if no context is attached.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the decoded state count does not
   * equal the context atom count.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Runtime` if the
   * string cannot be decoded.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if device
   * allocation or host-to-device transfer fails.
   *
   * @post Both the encoded position and every decoded state have been copied
   * into the integrator.
   */
  void setRngStates(const std::string &rngStateString);

  /**
   * @brief Resets both running-temperature values and their sample count.
   *
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if the device
   * array cannot be reset.
   */
  void resetAverageTemperature(void);

public:
  /**
   * @brief Returns the reference temperature.
   *
   * @return The stored temperature in kelvin.
   */
  double getReferenceTemperature(void) const;

  /**
   * @brief Returns the thermostat friction coefficient.
   *
   * @return The configured value in inverse picoseconds.
   */
  double getThermostatFriction(void) const;

  /**
   * @brief Returns the thermostat RNG seed.
   *
   * @return The dimensionless 64-bit seed.
   */
  std::uint64_t getThermostatRngSeed(void) const;

  /**
   * @brief Returns the recorded thermostat RNG sequence position.
   *
   * @return The dimensionless Philox position. Successful propagation advances
   * this value by four per step.
   */
  unsigned long long int getRngSequencePos(void) const;

  /**
   * @brief Copies the current per-atom RNG states to a serialized string.
   *
   * @return A newly allocated opaque string containing the recorded sequence
   * position and one Philox state per context atom.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::NotInitialized` if no context is attached or RNG state
   * has not been initialized.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if stream
   * synchronization or device-to-host transfer fails.
   * @throws std::bad_alloc If the host state vector or result string cannot be
   * allocated.
   *
   * @note The method synchronizes the primary integrator stream before copying.
   */
  std::string getRngStates(void) const;

  /**
   * @brief Returns the number of samples in the running average.
   *
   * @return A non-negative dimensionless count.
   */
  int getAverageWindowSize(void) const;

  /**
   * @brief Returns both kinetic-energy estimators.
   *
   * @return Borrowed const alias to a two-element container in kilocalories per
   * mole. Element 0 is the three-point estimator and element 1 is the on-step
   * estimator.
   */
  const CudaContainer<double> &getKineticEnergy(void) const;

  /**
   * @brief Returns both running-temperature estimators.
   *
   * @return Borrowed const alias to a two-element container in kelvin with the
   * same element order as `getKineticEnergy()`.
   */
  const CudaContainer<double> &getAverageTemperature(void) const;

  /**
   * @brief Returns mutable kinetic-energy storage.
   *
   * @return Borrowed mutable alias to the two-element container.
   * @warning Preserve the two-element shape and host/device coherence.
   */
  CudaContainer<double> &getKineticEnergy(void);

  /**
   * @brief Returns mutable running-temperature storage.
   *
   * @return Borrowed mutable alias to the two-element container.
   * @warning Preserve the two-element shape and host/device coherence.
   */
  CudaContainer<double> &getAverageTemperature(void);

  /**
   * @brief Computes the instantaneous three-point temperature.
   *
   * @return Kinetic-energy element 0 converted to kelvin.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::NotInitialized` if no context is attached.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Runtime` if the
   * kinetic-energy container does not contain exactly two elements.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if the kinetic
   * energy cannot be transferred to the host.
   */
  double getInstantaneousTemperature(void);

protected:
  /** @brief Allocates and initializes thermostat state from the context. */
  void initializeImpl(void) override;

  /**
   * @brief Restores thermostat displacement and RNG state from a restart file.
   *
   * @param[in] rstFileName Borrowed restart-file path.
   */
  void initializeFromRestartFileImpl(const std::string &rstFileName) override;

  /** @brief Executes one synchronized stochastic thermostat step. */
  void propagateOneStepImpl(void) override;

protected:
  /**
   * @brief Copies decoded RNG state to device storage.
   *
   * @param[in] position Dimensionless sequence position to store.
   * @param[in] states Borrowed contiguous host states copied during the call.
   *
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if allocation
   * or host-to-device transfer fails.
   */
  void setRngStateData(const unsigned long long int position,
                       const std::vector<curandStatePhilox4_32_10_t> &states);

  /**
   * @brief Allocates and initializes one RNG state per attached-context atom.
   *
   * The method does nothing before context attachment.
   */
  void initializeRng(void);

  /** @brief Removes mass-weighted center-of-mass displacement. */
  void removeCenterOfMassMotion(void);

  /**
   * @brief Allocates raw device storage for `n` RNG states.
   *
   * @param[in] n Number of per-atom states.
   */
  void alloc(const int n);

  /**
   * @brief Releases raw RNG-state storage through the checked native helper.
   */
  void dealloc(void);

protected:
  /** @brief Stores the reference temperature in kelvin. */
  double m_ReferenceTemperature;

  /** @brief Stores the public friction coefficient in inverse picoseconds. */
  double m_ThermostatFriction;

  /** @brief Stores the dimensionless friction accumulated over one time step.
   */
  double m_ThermostatGamma;

  /** @brief Stores the dimensionless 64-bit Philox seed. */
  std::uint64_t m_Seed;

  /** @brief Stores the recorded dimensionless Philox sequence position. */
  unsigned long long int m_RngSequencePos;

  /**
   * @brief Owns one raw device Philox state per context atom.
   *
   * The pointer is null before RNG initialization and is released by the
   * non-throwing destructor.
   */
  curandStatePhilox4_32_10_t *m_RngStates;

  /** @brief Stores the number of samples accumulated into running averages. */
  int m_AverageWindowSize;

  /**
   * @brief Stores two kinetic-energy estimators in kilocalories per mole.
   *
   * Element 0 is the three-point estimator and element 1 is the on-step
   * estimator.
   */
  CudaContainer<double> m_KineticEnergy;

  /** @brief Stores per-block partial sums for both kinetic estimators. */
  CudaContainer<double> m_KineticEnergyPartialSums;

  /** @brief Stores two running-temperature estimators in kelvin. */
  CudaContainer<double> m_AverageTemperature;
};

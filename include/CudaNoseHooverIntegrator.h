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

#include <filesystem>
#include <string>
#include <vector_types.h>

/**
 * @brief Propagates constant-temperature dynamics with Nose-Hoover coupling.
 *
 * The integrator maintains one Nose-Hoover extended variable, two kinetic
 * energy estimators, two running-temperature estimators, and a
 * predictor-corrector displacement buffer. Context attachment initializes the
 * derived state and calculates a default coupling mass when no explicit mass
 * has been supplied.
 *
 * @warning Treat the class as non-copyable. Its base owns CUDA streams whose
 * implicit copy semantics are unsafe.
 * @see cuda_integrators
 */
class CudaNoseHooverIntegrator : public CudaIntegrator {
public:
  /**
   * @brief Constructs an unattached Nose-Hoover integrator.
   *
   * @param[in] timeStep Finite positive time step in picoseconds.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `timeStep` is non-finite or not
   * positive.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if stream or
   * fixed-state allocation fails.
   * @throws std::bad_alloc If host or stream-owner allocation fails.
   *
   * @post The reference temperature is 300 K, the maximum
   * predictor-corrector count is three, the old-temperature selector is false,
   * and the running averages are zero.
   */
  CudaNoseHooverIntegrator(const double timeStep);

public:
  /**
   * @brief Sets the Nose-Hoover reference temperature.
   *
   * @param[in] referenceTemperature Finite non-negative temperature in kelvin.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `referenceTemperature` is
   * non-finite or negative.
   */
  void setReferenceTemperature(const double referenceTemperature);

  /**
   * @brief Sets the scalar Nose-Hoover coupling mass.
   *
   * @param[in] noseHooverPistonMass Finite non-negative native coupling-mass
   * value. The repository does not establish its public dimensional
   * convention.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `noseHooverPistonMass` is
   * non-finite or negative.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if the scalar
   * cannot be copied to device storage.
   *
   * @warning Zero is currently accepted, but active propagation divides by the
   * coupling mass.
   */
  void setNoseHooverPistonMass(const double noseHooverPistonMass);

  /**
   * @brief Sets the current Nose-Hoover extended-variable velocity.
   *
   * @param[in] noseHooverPistonVelocity Finite native extended-state value.
   * The repository does not establish its public dimensional convention.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the value is non-finite.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if the scalar
   * cannot be copied to device storage.
   */
  void setNoseHooverPistonVelocity(const double noseHooverPistonVelocity);

  /**
   * @brief Sets the previous Nose-Hoover extended-variable velocity.
   *
   * @param[in] noseHooverPistonVelocityPrevious Finite native previous-state
   * value. The repository does not establish its public dimensional
   * convention.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the value is non-finite.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if the scalar
   * cannot be copied to device storage.
   */
  void setNoseHooverPistonVelocityPrevious(
      const double noseHooverPistonVelocityPrevious);

  /**
   * @brief Sets the current Nose-Hoover velocity-increment state.
   *
   * Despite the public name, the propagation kernel stores a change in the
   * extended-variable velocity rather than a physical force.
   *
   * @param[in] noseHooverPistonForce Finite native increment value. The
   * repository does not establish its public dimensional convention.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the value is non-finite.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if the scalar
   * cannot be copied to device storage.
   */
  void setNoseHooverPistonForce(const double noseHooverPistonForce);

  /**
   * @brief Sets the previous Nose-Hoover velocity-increment state.
   *
   * @param[in] noseHooverPistonForcePrevious Finite native previous-increment
   * value. The repository does not establish its public dimensional
   * convention.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the value is non-finite.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if the scalar
   * cannot be copied to device storage.
   */
  void
  setNoseHooverPistonForcePrevious(const double noseHooverPistonForcePrevious);

  /**
   * @brief Sets the number of Nose-Hoover predictor-corrector iterations.
   *
   * @param[in] maxPredictorCorrectorIterations Positive, dimensionless
   * iteration count.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if
   * `maxPredictorCorrectorIterations` is not positive.
   */
  void
  setMaxPredictorCorrectorIterations(const int maxPredictorCorrectorIterations);

  /**
   * @brief Selects the on-step temperature estimator.
   *
   * @param[in] usingOldTemperature `true` selects array element 1, the on-step
   * “old” estimator, for native instantaneous temperature and thermostat
   * feedback. `false` selects element 0, the three-point “new” estimator.
   */
  void useOldTemperature(const bool usingOldTemperature);

  /**
   * @brief Resets both running-temperature values and their sample count.
   *
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if the
   * two-element temperature array cannot be reset on the device.
   *
   * @post The sample count and both temperature averages are zero.
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
   * @brief Returns the native temperature-estimator selector.
   *
   * @return `true` when native feedback and instantaneous temperature use
   * element 1; `false` when they use element 0.
   */
  bool usingOldTemperature(void) const;

  /**
   * @brief Returns the scalar coupling-mass container.
   *
   * @return Borrowed const alias to a one-element integrator-owned container.
   * The accessor performs no transfer or synchronization.
   */
  const CudaContainer<double> &getNoseHooverPistonMass(void) const;

  /**
   * @brief Returns the current extended-variable velocity container.
   *
   * @return Borrowed const alias to a one-element integrator-owned container.
   */
  const CudaContainer<double> &getNoseHooverPistonVelocity(void) const;

  /**
   * @brief Returns the previous extended-variable velocity container.
   *
   * @return Borrowed const alias to a one-element integrator-owned container.
   */
  const CudaContainer<double> &getNoseHooverPistonVelocityPrevious(void) const;

  /**
   * @brief Returns the current velocity-increment container.
   *
   * @return Borrowed const alias to a one-element integrator-owned container.
   */
  const CudaContainer<double> &getNoseHooverPistonForce(void) const;

  /**
   * @brief Returns the previous velocity-increment container.
   *
   * @return Borrowed const alias to a one-element integrator-owned container.
   */
  const CudaContainer<double> &getNoseHooverPistonForcePrevious(void) const;

  /**
   * @brief Returns the predictor-corrector iteration count.
   *
   * @return A positive dimensionless iteration count.
   */
  int getMaxPredictorCorrectorIterations(void) const;

  /**
   * @brief Returns both kinetic-energy estimators.
   *
   * @return Borrowed const alias to a two-element container in kilocalories per
   * mole. Element 0 is the three-point estimator and element 1 is the on-step
   * estimator.
   *
   * @note The accessor performs no device-to-host transfer.
   */
  const CudaContainer<double> &getKineticEnergy(void) const;

  /**
   * @brief Returns both running-temperature estimators.
   *
   * @return Borrowed const alias to a two-element container in kelvin. Element
   * 0 corresponds to the three-point estimator and element 1 to the on-step
   * estimator.
   *
   * @note The accessor performs no device-to-host transfer.
   */
  const CudaContainer<double> &getAverageTemperature(void) const;

  /**
   * @brief Returns mutable scalar coupling-mass storage.
   *
   * @return Borrowed mutable alias to the one-element container.
   * @warning Preserve its one-element shape and host/device coherence.
   */
  CudaContainer<double> &getNoseHooverPistonMass(void);

  /**
   * @brief Returns mutable current extended-variable velocity storage.
   *
   * @return Borrowed mutable alias to the one-element container.
   * @warning Preserve its one-element shape and host/device coherence.
   */
  CudaContainer<double> &getNoseHooverPistonVelocity(void);

  /**
   * @brief Returns mutable previous extended-variable velocity storage.
   *
   * @return Borrowed mutable alias to the one-element container.
   * @warning Preserve its one-element shape and host/device coherence.
   */
  CudaContainer<double> &getNoseHooverPistonVelocityPrevious(void);

  /**
   * @brief Returns mutable current velocity-increment storage.
   *
   * @return Borrowed mutable alias to the one-element container.
   * @warning Preserve its one-element shape and host/device coherence.
   */
  CudaContainer<double> &getNoseHooverPistonForce(void);

  /**
   * @brief Returns mutable previous velocity-increment storage.
   *
   * @return Borrowed mutable alias to the one-element container.
   * @warning Preserve its one-element shape and host/device coherence.
   */
  CudaContainer<double> &getNoseHooverPistonForcePrevious(void);

  /**
   * @brief Returns mutable kinetic-energy storage.
   *
   * @return Borrowed mutable alias to the two-element container.
   * @warning Preserve the two-element estimator layout.
   */
  CudaContainer<double> &getKineticEnergy(void);

  /**
   * @brief Returns mutable running-temperature storage.
   *
   * @return Borrowed mutable alias to the two-element container.
   * @warning Preserve the two-element estimator layout.
   */
  CudaContainer<double> &getAverageTemperature(void);

  /**
   * @brief Computes the selected instantaneous temperature.
   *
   * @return The selected kinetic-energy estimator converted to kelvin.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::NotInitialized` if no context is attached.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if kinetic
   * energy cannot be transferred to the host.
   *
   * @note The method selects element 1 when `usingOldTemperature()` is true and
   * element 0 otherwise.
   */
  double getInstantaneousTemperature(void);

protected:
  /** @brief Allocates and initializes Nose-Hoover state from the context. */
  void initializeImpl(void) override;

  /**
   * @brief Restores Nose-Hoover state from a CHARMM-style restart file.
   *
   * @param[in] rstFilePath Borrowed restart-file path. The path is not
   * canonicalized or retained.
   */
  void initializeFromRestartFileImpl(
      const std::filesystem::path &rstFilePath) override;

  /** @brief Executes one synchronized Nose-Hoover propagation step. */
  void propagateOneStepImpl(void) override;

protected:
  /**
   * @brief Computes the default numeric coupling mass.
   *
   * @return `0.2` times the sum of context atom masses. The exact dimensional
   * convention of the returned extended-system value is not established by the
   * repository.
   */
  double computeNoseHooverPistonMass(void);

  /**
   * @brief Removes mass-weighted center-of-mass displacement.
   *
   * The method synchronizes the integrator stream, transfers velocity and
   * previous-displacement storage to the host, applies the PBC-dependent
   * correction, and transfers the displacement storage back to the device.
   */
  void removeCenterOfMassMotion(void);

protected:
  /** @brief Stores the reference temperature in kelvin. */
  double m_ReferenceTemperature;

  /** @brief Stores the one-element Nose-Hoover coupling mass. */
  CudaContainer<double> m_NoseHooverPistonMass;

  /** @brief Stores the current one-element extended-variable velocity. */
  CudaContainer<double> m_NoseHooverPistonVelocity;

  /** @brief Stores the previous extended-variable velocity. */
  CudaContainer<double> m_NoseHooverPistonVelocityPrevious;

  /**
   * @brief Stores the current extended-variable velocity increment.
   *
   * The propagation kernel treats this value as an increment despite the
   * historical member name.
   */
  CudaContainer<double> m_NoseHooverPistonForce;

  /** @brief Stores the previous extended-variable velocity increment. */
  CudaContainer<double> m_NoseHooverPistonForcePrevious;

  /** @brief Stores one predicted `double4` displacement per context atom. */
  CudaContainer<double4> m_CoordsDeltaPredicted;

  /** @brief Stores the positive predictor-corrector iteration count. */
  int m_MaxPredictorCorrectorIterations;

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

  /**
   * @brief Stores two running-temperature estimators in kelvin.
   *
   * The element order matches @ref m_KineticEnergy.
   */
  CudaContainer<double> m_AverageTemperature;

  /** @brief Selects element 1 instead of element 0 for native feedback. */
  bool m_UsingOldTemperature;
};

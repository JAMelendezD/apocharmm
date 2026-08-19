// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author:  James E. Gonzales II, Samarjeet Prasad
//
// ENDLICENSE

#pragma once

#include "Subscriber.h"

#include <string>

/**
 * @brief Replaces a CHARMM-style restart file with the latest simulation state.
 *
 * Each successful @ref update writes version-50 text containing periodic-box
 * state, integrator counters and random state, energy/statistics placeholders,
 * double-precision coordinates, velocities, and the integrator's previous
 * coordinate-delta state. The file is intended for
 * `CudaIntegrator::initializeFromRestartFile()` on supported integrators.
 *
 * Box lengths and coordinates use angstroms. Context velocities use angstroms
 * per AKMA time unit. Average temperature uses kelvin. Arrays are written in
 * atom order with X, Y, and Z components on each line. Several CHARMM energy
 * and statistics fields are currently emitted as zeros.
 *
 * Updating transfers all required CUDA containers to host memory. Each transfer
 * synchronizes all previously requested work on the current CUDA device; no
 * subscriber-specific stream is used.
 *
 * @warning Updates support only `CudaNoseHooverIntegrator`,
 * `CudaLangevinPistonIntegrator`, and
 * `CudaLangevinThermostatIntegrator`.
 * @warning Each update truncates the target before all sections are generated.
 * A later failure can destroy a previously usable restart and leave a partial
 * replacement.
 *
 * @see Subscriber
 * @see subscriber
 */
class RestartSubscriber : public Subscriber {
public:
  /**
   * @brief Constructs an unattached restart writer without an output file.
   *
   * @post The report frequency is `1000`, the logical file name is empty, and
   * the stream is closed.
   */
  RestartSubscriber(void);

  /**
   * @brief Constructs a restart writer with the default reporting interval.
   *
   * Construction creates or truncates the text output immediately. A complete
   * restart is not written until @ref update is invoked.
   *
   * @param[in] fileName Output path copied by the subscriber. The value must be
   * nonempty; when it contains `/`, its nonempty parent path must exist.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `fileName` is empty or its checked
   * parent path does not exist.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Runtime` if the text
   * output file cannot be opened for writing.
   *
   * @post The file has been created or truncated and
   * `getReportFrequency() == 1000`.
   */
  RestartSubscriber(const std::string &fileName);

  /**
   * @brief Constructs a restart writer with an explicit reporting interval.
   *
   * Construction validates the frequency, then creates or truncates the text
   * output. A complete restart is not written until @ref update is invoked.
   *
   * @param[in] fileName Output path copied by the subscriber. The value must be
   * nonempty; when it contains `/`, its nonempty parent path must exist.
   * @param[in] reportFrequency Positive, dimensionless number of propagated
   * steps between restart replacements.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `reportFrequency` is not positive,
   * `fileName` is empty, or its checked parent path does not exist.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Runtime` if the text
   * output file cannot be opened for writing.
   *
   * @post The file has been created or truncated and
   * `getReportFrequency() == reportFrequency`.
   */
  RestartSubscriber(const std::string &fileName, const int reportFrequency);

  /**
   * @brief Destroys the writer and closes its stream.
   *
   * Destruction is non-throwing.
   */
  ~RestartSubscriber(void) noexcept override = default;

  /**
   * @brief Replaces the output with one restart snapshot.
   *
   * The method validates the attached collaborators and box before reopening
   * the path in truncating text mode. It writes one complete restart and closes
   * the stream before checking the final stream state.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::NotInitialized` if the logical file name is empty, no
   * context or integrator is attached, the integrator has no context, or the
   * context does not contain exactly three positive box lengths.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::NotImplemented` if the attached integrator is not a
   * supported Nose-Hoover, Langevin-piston, or Langevin-thermostat integrator.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Runtime` if the
   * subscriber and integrator retain different context objects, the file cannot
   * be reopened, or closing detects failed output.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if any required
   * device-to-host transfer or device-wide synchronization fails.
   *
   * @pre Context per-atom arrays and every integrator-specific state container
   * required by the selected integrator have valid lengths and device storage.
   * @post On success, the stream is closed and the path contains one complete
   * restart representing the state observed during this call.
   * @warning The operation is not transactional. Output and host mirrors
   * updated before a failure remain observable.
   */
  void update(void) override;
};

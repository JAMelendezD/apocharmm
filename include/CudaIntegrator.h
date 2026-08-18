// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author:  Samarjeet Prasad, James E. Gonzales II
//
// ENDLICENSE

#pragma once

#include "CharmmContext.h"
#include "CudaContainer.h"
#include "CudaHolonomicConstraint.h"
#include "Subscriber.h"

#include <cuda_runtime.h>
#include <map>
#include <memory>
#include <string>
#include <vector>

/**
 * @brief Base class for integrators
 * @todo Add CudaHolonomicConstraints as a member here
 * @todo should timeStep be float/double?
 *
 * Base class to build integrators.While propagating, CudaIntegrator will
 * notify the right Subscriber object (through the CharmmContext) using this
 * reportFreqList as an index.
 *
 */
class CudaIntegrator : public std::enable_shared_from_this<CudaIntegrator> {
public:
  CudaIntegrator(void);

  /** @brief Class constructor. Takes timeStep (in ps) as argument.
   *
   * @param[in] timeStep Time step, in ps
   */
  CudaIntegrator(const double timeStep);

  /** @brief Class constructor. Takes timeStep (ps) as argument, as well as
   * debugPrintFrequency (default 0).
   * @param[in] timeStep Time step, in ps
   * @param[in] debugPrintFrequency Frequency (number of timestep) at which to
   * print integrator infos
   */
  CudaIntegrator(const double timeStep, const int debugPrintFrequency);

  virtual ~CudaIntegrator(void) noexcept;

  /** @brief Returns integrator timestep (in ps)
   */
  double getTimeStep(void) const;

  /**
   * @brief Set integrator timestep (in ps)
   *
   * Sets timeStep variable to dt * unit time factor.
   * If the Subscriber list is not empty, then calls setTimeStepFromIntegrator
   * for each member.
   */
  void setTimeStep(const double dt);

  /**
   * @brief Link integrator to  CharmmContext
   *
   * Sets context member variable, allocates coordsRef (?), sets
   * holonomicConstraint. If a CharmmContext was already set, throws an
   * exception.
   *
   * @param[in] ctx CharmmContext to be linked
   */
  virtual void setCharmmContext(std::shared_ptr<CharmmContext> ctx) final;

  const std::shared_ptr<CharmmContext> getCharmmContext(void) const;

  std::shared_ptr<CharmmContext> getCharmmContext(void);

  virtual void initialize(void) final;

  virtual void initializeFromRestartFile(const std::string &rstFileName) final;

  /**
   * @brief Propagate a single time step
   */
  virtual void propagateOneStep(void) final;

  /**
   * @brief Propagates dynamics and invokes scheduled subscriber updates.
   *
   * Steps are numbered locally from one through `numSteps`. After each
   * completed dynamics step, the integrator synchronously calls subscribers
   * whose cached reporting interval divides that local index. A later call to
   * this method restarts the local schedule at one.
   *
   * The total requested-step counter is increased before entering the loop.
   * Subscriber callbacks run after the step, current-step increment, and NaN
   * energy check for that iteration.
   *
   * @param[in] numSteps Positive number of dynamics steps to request.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `numSteps` is not positive.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::NotInitialized` if required integrator, context,
   * force- manager, or subscriber state is missing.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::NotImplemented` if a subscriber update reports an
   * unsupported operation, including restart output for an unsupported
   * integrator type.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if dynamics,
   * force evaluation, a subscriber transfer, a kernel, or CUDA synchronization
   * fails.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Runtime` if a runtime
   * readiness check, NaN check, subscriber file operation, or other native
   * runtime operation fails.
   * @throws std::invalid_argument If a scheduled legacy DynaSubscriber rejects
   * its current integrator.
   * @throws std::runtime_error If a scheduled legacy XYZSubscriber has no
   * attached context.
   *
   * @post On success, exactly `numSteps` dynamics steps and every matching
   * callback have completed.
   * @warning The operation is not transactional. Earlier dynamics steps, state
   * changes, and output bytes remain observable after a later failure.
   */
  void propagate(const int numSteps);

  void setDebugPrintFrequency(const int freq);

  void setNonbondedListUpdateFrequency(const int freq);

  /**
   * @brief Attaches one subscriber and caches its reporting interval.
   *
   * The method rejects null and duplicate subscribers, obtains a shared owner
   * of this integrator, calls `sub->setCharmmContext(m_Context)`, calls
   * `sub->setIntegrator(...)`, appends the subscriber, and appends its current
   * report frequency to the parallel scheduling array.
   *
   * @param[in] sub Non-null shared subscriber owner retained on success.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `sub` is null, is already present,
   * or rejects a second context or integrator attachment.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::NotInitialized` if this integrator is not owned by a
   * `std::shared_ptr` and therefore cannot obtain a shared self-reference.
   * @throws std::bad_alloc If appending either scheduling array cannot
   * allocate.
   *
   * @pre Call @ref setCharmmContext with a non-null context before attaching a
   * subscriber that requires context state.
   * @post On success, the integrator and subscriber retain each other, the
   * subscriber retains the current context, and the current frequency is
   * cached.
   * @warning Attachment is not transactional. Failure after a backlink or first
   * vector append can leave partial state.
   */
  void subscribe(std::shared_ptr<Subscriber> sub);

  /**
   * @brief Attaches subscribers sequentially in vector order.
   *
   * @param[in] sublist Borrowed vector of shared subscriber owners. Each
   * element is passed to the single-subscriber overload and retained on
   * success.
   *
   * @note Exceptions from the single-subscriber overload propagate unchanged,
   * including `ApoCharmmError` and `std::bad_alloc`.
   *
   * @post On success, every element has been attached in input order.
   * @warning The operation is not transactional. Subscribers attached before a
   * failing element remain attached.
   */
  void subscribe(const std::vector<std::shared_ptr<Subscriber>> &sublist);

  /**
   * @brief Removes one subscriber and its cached reporting interval.
   *
   * @param[in] sub Non-null shared subscriber identity to remove. The argument
   * is borrowed for lookup; removing it releases only the integrator's owner.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `sub` is null or is not present.
   *
   * @post On success, the corresponding entries are erased from both parallel
   * arrays.
   * @warning The removed subscriber still retains its context and integrator.
   * Those backlinks are not cleared, and normal resubscription of the same
   * object therefore fails.
   */
  void unsubscribe(std::shared_ptr<Subscriber> sub);

  /**
   * @brief Removes subscribers sequentially in vector order.
   *
   * @param[in] sublist Borrowed vector of subscriber identities. Each element
   * is passed to the single-subscriber overload.
   *
   * @note Exceptions from the single-subscriber overload propagate unchanged.
   *
   * @post On success, every listed subscriber has been removed.
   * @warning The operation is not transactional. Subscribers removed before a
   * failing element remain removed.
   */
  void unsubscribe(const std::vector<std::shared_ptr<Subscriber>> &sublist);

  /**
   * @brief Returns the attached subscriber array.
   *
   * @return Borrowed const alias to the integrator-owned vector. Each element
   * is a shared owner of a subscriber. The alias remains valid until vector
   * reallocation, integrator destruction, or non-const mutation.
   */
  const std::vector<std::shared_ptr<Subscriber>> &getSubscribers(void) const;

  /**
   * @brief Returns mutable access to the attached subscriber array.
   *
   * @return Borrowed mutable alias to the integrator-owned vector.
   * @warning The vector must remain the same length and order as
   * @ref m_ReportFreqList. Direct mutation can violate that invariant and cause
   * out-of-bounds access or callbacks to the wrong subscriber.
   */
  std::vector<std::shared_ptr<Subscriber>> &getSubscribers(void);

  /**
   * @brief Returns the cached subscriber-frequency array.
   *
   * @return Borrowed const alias to positive, dimensionless intervals. Entry
   * `i` schedules subscriber entry `i`. The alias remains valid until vector
   * reallocation, integrator destruction, or non-const mutation.
   */
  const std::vector<int> &getReportFreqList(void) const;

  /**
   * @brief Returns mutable access to cached subscriber frequencies.
   *
   * @return Borrowed mutable alias to the integrator-owned vector.
   * @warning Every value must remain positive and the vector must remain the
   * same length and order as @ref m_Subscribers. Mutation bypasses validation.
   */
  std::vector<int> &getReportFreqList(void);

  void setRemoveCenterOfMassFrequency(const int freq);

  const CudaContainer<double4> &getCoordsDelta(void) const;

  CudaContainer<double4> &getCoordsDelta(void);

  const CudaContainer<double4> &getCoordsDeltaPrevious(void) const;

  CudaContainer<double4> &getCoordsDeltaPrevious(void);
  /** @brief Returns a map of the integrator descriptor. Should be overriden by
   * child classes.*/
  virtual std::map<std::string, std::string> getIntegratorDescriptors(void);

  /** @brief Returns the current step which the integrator is on. */
  int getCurrentPropagatedStep(void) const;

  /** @brief Returns the total number of dynamics steps that have been performed
   * for a "given" trajectory.
   */
  unsigned long long int getTotNumSteps(void) const;

  /** @brief Returns the number of steps most recently requested for dynamics.
   */
  int getNumSteps(void) const;

protected:
  virtual void initializeImpl(void);

  virtual void initializeFromRestartFileImpl(const std::string &rstFileName);

  virtual void propagateOneStepImpl(void);

  void requirePropagationReady(void) const;

protected:
  /**
   * @brief Invokes subscribers scheduled for one local propagation step.
   *
   * The method iterates over the cached frequency array and calls
   * `m_Subscribers[i]->update()` when `istep % m_ReportFreqList[i] == 0`.
   * Calls are synchronous and occur in subscriber-vector order.
   *
   * @param[in] istep Positive step index local to the current propagation call.
   *
   * @pre @ref m_Subscribers and @ref m_ReportFreqList have identical lengths
   * and corresponding order, and every cached frequency is positive.
   * @note Any exception raised by a subscriber propagates immediately. Later
   * subscribers scheduled for the same step are not called.
   */
  void reportIfNeeded(const int istep);

  /** @brief Checks that the kinetic and potential energy are not nans. Throws
   * an error otherwise. Called every min(10^7, min(reportFreqList)) steps.
   * @todo Unittest this  */
  void checkForNanEnergy(void);

  int wrapCurrentPropagatedStep(const unsigned long long int propagatedStep);
  void incrementCurrentPropagatedStep(void);

protected:
  /** @brief Integrator time step, in AKMA units (ps / timfac) */
  double m_TimeStep;

  /** @brief Time unit conversion factor (t (AKMA) = t (psf) / timfac ) */
  double m_Timfac;

  /** @brief If not 0, frequency at which integrator should print debug infos */
  int m_DebugPrintFrequency;

  /** @brief CharmmContext to which Integrator is attached */
  std::shared_ptr<CharmmContext> m_Context;

  /** @todo Pick a better name : this counts the number of steps and is used to
   * compare with nblupdate frequency, rather than counting how many steps have
   * happened since last NBLupdate */
  int m_StepsSinceNeighborListUpdate;

  /** @brief Allows subscribers to have knowledge of which step the integrator
   * has propgated*/
  int m_CurrentPropagatedStep;

  unsigned long long int m_TotNumSteps;

  int m_NumSteps;

  std::shared_ptr<CudaHolonomicConstraint> m_HolonomicConstraint;

  /**
   * @todo  document this
   */
  CudaContainer<double4> m_CoordsRef;
  CudaContainer<double4> m_CoordsDelta;
  CudaContainer<double4> m_CoordsDeltaPrevious;

  std::shared_ptr<cudaStream_t> m_IntegratorStream;
  std::shared_ptr<cudaStream_t> m_IntegratorMemcpyStream;
  bool m_UsingHolonomicConstraints;

  /**
   * @brief Retains subscribers in callback order.
   *
   * Entry `i` must correspond to entry `i` in @ref m_ReportFreqList. Each
   * shared pointer owns the subscriber until successful unsubscription or
   * integrator destruction.
   */
  std::vector<std::shared_ptr<Subscriber>> m_Subscribers;

  /**
   * @brief Stores the frequency snapshot for each attached subscriber.
   *
   * Entry `i` is the positive, dimensionless reporting interval copied from
   * subscriber `i` during attachment. Later changes to the subscriber property
   * do not update this array.
   */
  std::vector<int> m_ReportFreqList;

  int m_NonbondedListUpdateFrequency;

  int m_RemoveCenterOfMassFrequency;

  /** @brief Describe the integrator type. Useful to discriminate methods to
   * use, e.g. for restart subscribers. Could/should be a trait */
  std::string m_IntegratorTypeName;
};

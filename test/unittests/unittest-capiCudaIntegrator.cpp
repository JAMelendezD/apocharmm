// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "CharmmContext.h"
#include "CharmmCrd.h"
#include "CharmmPSF.h"
#include "CharmmParameters.h"
#include "CudaIntegrator.h"
#include "apo_test_helpers.h"
#include "apocharmm_c/CharmmContext.h"
#include "apocharmm_c/CudaIntegrator.h"
#include "apocharmm_c/Error.h"
#include "apocharmm_c/Status.h"
#include "apocharmm_c/detail/CharmmContextHandle.h"
#include "apocharmm_c/detail/CudaIntegratorHandle.h"
#include "catch.hpp"
#include "test_paths.h"

#include <memory>
#include <string>

namespace {

constexpr int RANDOM_SEED = 314159;
constexpr double TEMPERATURE = 300.0;
constexpr double TIME_STEP = 0.002;

class CudaIntegratorProbe final : public CudaIntegrator {
public:
  using CudaIntegrator::CudaIntegrator;

  void initialize(void) override { return; }

  void propagateOneStep(void) override { return; }
};

std::shared_ptr<CharmmContext> MakeInitializedContext(void) {
  const std::string dataPath = getDataPath();

  auto prm =
      std::make_shared<CharmmParameters>(dataPath + "toppar_water_ions.str");
  auto psf = std::make_shared<CharmmPSF>(dataPath + "nacl_pair.psf");
  auto crd = std::make_shared<CharmmCrd>(dataPath + "nacl_pair.cor");

  auto context = std::make_shared<CharmmContext>(psf, prm);
  context->setBoxDimensions({50.0, 50.0, 50.0});
  context->setCoordinates(crd);
  context->useHolonomicConstraints(false);
  context->setRandomSeed(RANDOM_SEED);
  context->assignVelocitiesAtTemperature(TEMPERATURE);

  return context;
}

apo_cuda_integrator MakeIntegratorHandle(void) {
  apo_cuda_integrator integrator;
  integrator.object = std::make_shared<CudaIntegratorProbe>(TIME_STEP);
  return integrator;
}

apo_charmm_context MakeContextHandle(void) {
  apo_charmm_context context;
  context.object = MakeInitializedContext();
  return context;
}

} // namespace

TEST_CASE("CapiCudaIntegratorSetTimeStep") {
  apo_cuda_integrator integrator = MakeIntegratorHandle();
  apo_cuda_integrator emptyIntegrator;

  CHECK(apo_cuda_integrator_set_time_step(&integrator, 0.001) == APO_STATUS_OK);
  CHECK(integrator.object->getTimeStep() == Approx(0.001));
  CHECK(std::string(apo_last_error()).empty() == true);

  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_set_time_step(nullptr, TIME_STEP),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_set_time_step: CudaIntegrator is NULL");
  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_set_time_step(&emptyIntegrator, TIME_STEP),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_set_time_step: CudaIntegrator object is NULL");
}

TEST_CASE("CapiCudaIntegratorContextAttachment") {
  apo_cuda_integrator integrator = MakeIntegratorHandle();
  apo_cuda_integrator emptyIntegrator;
  apo_charmm_context context = MakeContextHandle();
  apo_charmm_context emptyContext;

  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_set_charmm_context(nullptr, &context),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_set_charmm_context: CudaIntegrator is NULL");
  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_set_charmm_context(&emptyIntegrator, &context),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_set_charmm_context: CudaIntegrator object is NULL");
  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_set_charmm_context(&integrator, nullptr),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_set_charmm_context: CharmmContext is NULL");
  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_set_charmm_context(&integrator, &emptyContext),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_set_charmm_context: CharmmContext object is NULL");

  apo_charmm_context incompleteContext;
  incompleteContext.object = std::make_shared<CharmmContext>();
  apo_status status = APO_STATUS_OK;

  CHECK_NOTHROW((status = apo_cuda_integrator_set_charmm_context(
                     &integrator, &incompleteContext)));
  apo_test::CheckNativeError(
      status, APO_STATUS_NOT_INITIALIZED, "NotInitialized",
      "apo_cuda_integrator_set_charmm_context",
      "CharmmContext atom count is not initialized; observed -1",
      "src/CudaIntegrator.cu", "setCharmmContext");
  CHECK(integrator.context == nullptr);
  CHECK(integrator.object->getCharmmContext() == nullptr);

  CHECK(apo_cuda_integrator_set_charmm_context(&integrator, &context) ==
        APO_STATUS_OK);
  CHECK(integrator.context == context.object);
  CHECK(integrator.object->getCharmmContext() == context.object);

  CHECK_NOTHROW(
      (status = apo_cuda_integrator_set_charmm_context(&integrator, &context)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_cuda_integrator_set_charmm_context",
      "A CharmmContext object was already set for this CudaIntegrator.",
      "src/CudaIntegrator.cu", "setCharmmContext");
  CHECK(integrator.context == context.object);
  CHECK(integrator.object->getCharmmContext() == context.object);
}

TEST_CASE("CapiCudaIntegratorPropagation") {
  apo_cuda_integrator integrator = MakeIntegratorHandle();
  apo_cuda_integrator emptyIntegrator;
  apo_status status = APO_STATUS_OK;

  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_propagate(nullptr, 1), APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_propagate: CudaIntegrator is NULL");
  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_propagate(&emptyIntegrator, 1),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_propagate: CudaIntegrator object is NULL");

  CHECK_NOTHROW((status = apo_cuda_integrator_propagate(&integrator, 0)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_cuda_integrator_propagate",
      "Number of propagation steps must be positive; observed 0",
      "src/CudaIntegrator.cu", "propagate");
  CHECK(integrator.object->getTotNumSteps() == 0ULL);

  CHECK_NOTHROW((status = apo_cuda_integrator_propagate(&integrator, -1)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_cuda_integrator_propagate",
      "Number of propagation steps must be positive; observed -1",
      "src/CudaIntegrator.cu", "propagate");
  CHECK(integrator.object->getTotNumSteps() == 0ULL);

  CHECK_NOTHROW((status = apo_cuda_integrator_propagate(&integrator, 1)));
  apo_test::CheckNativeError(status, APO_STATUS_NOT_INITIALIZED,
                             "NotInitialized", "apo_cuda_integrator_propagate",
                             "CharmmContext must be set before propagation",
                             "src/CudaIntegrator.cu", "propagate");

  apo_charmm_context context = MakeContextHandle();
  REQUIRE(apo_cuda_integrator_set_charmm_context(&integrator, &context) ==
          APO_STATUS_OK);

  CHECK(apo_cuda_integrator_propagate(&integrator, 1) == APO_STATUS_OK);
  CHECK(integrator.object->getNumSteps() == 1);
  CHECK(integrator.object->getTotNumSteps() == 1ULL);
  CHECK(integrator.object->getCurrentPropagatedStep() == 1);
  CHECK(std::string(apo_last_error()).empty() == true);
}

TEST_CASE("CapiCudaIntegratorRestartInitialization") {
  apo_cuda_integrator integrator = MakeIntegratorHandle();
  apo_cuda_integrator emptyIntegrator;
  apo_status status = APO_STATUS_OK;

  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_initialize_from_restart_file(nullptr, "restart.rst"),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_initialize_from_restart_file: CudaIntegrator is "
      "NULL");
  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_initialize_from_restart_file(&emptyIntegrator,
                                                       "restart.rst"),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_initialize_from_restart_file: CudaIntegrator "
      "object is NULL");
  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_initialize_from_restart_file(&integrator, nullptr),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_initialize_from_restart_file: Restart file path is "
      "NULL");

  CHECK_NOTHROW((status = apo_cuda_integrator_initialize_from_restart_file(
                     &integrator, "restart.rst")));
  apo_test::CheckNativeError(
      status, APO_STATUS_NOT_IMPLEMENTED, "NotImplemented",
      "apo_cuda_integrator_initialize_from_restart_file",
      "CudaIntegrator::initializeFromRestartFile is not implemented by the "
      "base class",
      "src/CudaIntegrator.cu", "initializeFromRestartFile");
}

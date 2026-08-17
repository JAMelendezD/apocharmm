// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "apo_test_helpers.h"
#include "apocharmm_c/CudaLangevinThermostatIntegrator.h"
#include "apocharmm_c/Error.h"
#include "apocharmm_c/Status.h"
#include "apocharmm_c/detail/CudaLangevinThermostatIntegratorHandle.h"
#include "catch.hpp"

#include <limits>
#include <string>

namespace {

constexpr double TIME_STEP = 0.001;

apo_cuda_langevin_thermostat_integrator *CreateIntegrator(void) {
  apo_cuda_langevin_thermostat_integrator *integrator = nullptr;

  const apo_status status =
      apo_cuda_langevin_thermostat_integrator_create(&integrator, TIME_STEP);

  REQUIRE(status == APO_STATUS_OK);
  REQUIRE(integrator != nullptr);

  return integrator;
}

} // namespace

TEST_CASE("CapiCudaLangevinThermostatConstructionAndDestruction") {
  SECTION("NullOutputPointer") {
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW(status = apo_cuda_langevin_thermostat_integrator_create(
                      nullptr, TIME_STEP));

    apo_test::CheckStatusAndDiagnostic(
        status, APO_STATUS_INVALID_ARGUMENT,
        "apo_cuda_langevin_thermostat_integrator_create: out pointer is NULL");
  }

  SECTION("ZeroTimeStep") {
    apo_cuda_langevin_thermostat_integrator *integrator = nullptr;
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW(status = apo_cuda_langevin_thermostat_integrator_create(
                      &integrator, 0.0));

    CHECK(integrator == nullptr);

    apo_test::CheckNativeError(status, APO_STATUS_INVALID_ARGUMENT,
                               "InvalidArgument",
                               "apo_cuda_langevin_thermostat_integrator_create",
                               "Time step must be positive; observed 0.000000",
                               "src/CudaIntegrator.cu", "setTimeStep");
  }

  SECTION("InfiniteTimeStep") {
    const double infinity = std::numeric_limits<double>::infinity();

    apo_cuda_langevin_thermostat_integrator *integrator = nullptr;
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW(status = apo_cuda_langevin_thermostat_integrator_create(
                      &integrator, infinity));

    CHECK(integrator == nullptr);

    apo_test::CheckNativeError(
        status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
        "apo_cuda_langevin_thermostat_integrator_create",
        "Time step must be finite; observed " + std::to_string(infinity),
        "src/CudaIntegrator.cu", "setTimeStep");
  }

  SECTION("DestroyIsNoThrow") {
    apo_cuda_langevin_thermostat_integrator *integrator = CreateIntegrator();

    CHECK_NOTHROW(apo_cuda_langevin_thermostat_integrator_destroy(integrator));
    CHECK_NOTHROW(apo_cuda_langevin_thermostat_integrator_destroy(nullptr));
  }
}

TEST_CASE("CapiCudaLangevinThermostatValidation") {
  const double infinity = std::numeric_limits<double>::infinity();

  apo_cuda_langevin_thermostat_integrator *integrator = CreateIntegrator();
  apo_status status = APO_STATUS_OK;

  CHECK_NOTHROW(
      status =
          apo_cuda_langevin_thermostat_integrator_set_reference_temperature(
              nullptr, 300.0));

  apo_test::CheckStatusAndDiagnostic(
      status, APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_langevin_thermostat_integrator_set_reference_temperature: "
      "CudaLangevinThermostatIntegrator is NULL");

  CHECK_NOTHROW(
      status =
          apo_cuda_langevin_thermostat_integrator_set_reference_temperature(
              integrator, infinity));

  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_cuda_langevin_thermostat_integrator_set_reference_temperature",
      "Reference temperature must be finite; observed " +
          std::to_string(infinity),
      "src/CudaLangevinThermostatIntegrator.cu", "setReferenceTemperature");

  CHECK_NOTHROW(
      status =
          apo_cuda_langevin_thermostat_integrator_set_reference_temperature(
              integrator, -1.0));

  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_cuda_langevin_thermostat_integrator_set_reference_temperature",
      "Reference temperature must be non-negative; observed -1.000000",
      "src/CudaLangevinThermostatIntegrator.cu", "setReferenceTemperature");

  CHECK_NOTHROW(
      status = apo_cuda_langevin_thermostat_integrator_set_thermostat_friction(
          integrator, infinity));

  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_cuda_langevin_thermostat_integrator_set_thermostat_friction",
      "Thermostat friction must be finite; observed " +
          std::to_string(infinity),
      "src/CudaLangevinThermostatIntegrator.cu", "setThermostatFriction");

  CHECK_NOTHROW(
      status = apo_cuda_langevin_thermostat_integrator_set_thermostat_friction(
          integrator, -1.0));

  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_cuda_langevin_thermostat_integrator_set_thermostat_friction",
      "Thermostat friction must be non-negative; observed -1.000000",
      "src/CudaLangevinThermostatIntegrator.cu", "setThermostatFriction");

  double temperature = 0.0;

  CHECK_NOTHROW(
      status =
          apo_cuda_langevin_thermostat_integrator_get_instantaneous_temperature(
              &temperature, integrator));

  apo_test::CheckNativeError(
      status, APO_STATUS_NOT_INITIALIZED, "NotInitialized",
      "apo_cuda_langevin_thermostat_integrator_get_instantaneous_temperature",
      "CharmmContext must be set before computing instantaneous temperature",
      "src/CudaLangevinThermostatIntegrator.cu", "getInstantaneousTemperature");

  CHECK_NOTHROW(apo_cuda_langevin_thermostat_integrator_destroy(integrator));
}

TEST_CASE("CapiCudaLangevinThermostatGetterInvariants") {
  apo_cuda_langevin_thermostat_integrator *integrator = CreateIntegrator();

  integrator->object->getAverageTemperature().resize(0);

  double temperature = 0.0;
  apo_status status = APO_STATUS_OK;

  CHECK_NOTHROW(
      status = apo_cuda_langevin_thermostat_integrator_get_average_temperature(
          &temperature, integrator));

  apo_test::CheckStatusAndDiagnostic(
      status, APO_STATUS_RUNTIME_ERROR,
      "apo_cuda_langevin_thermostat_integrator_get_average_temperature: "
      "Average temperature does not contain exactly 2 elements");

  CHECK_NOTHROW(apo_cuda_langevin_thermostat_integrator_destroy(integrator));
}

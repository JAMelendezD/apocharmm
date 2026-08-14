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
#include "apocharmm_c/CudaLangevinPistonIntegrator.h"
#include "apocharmm_c/Error.h"
#include "apocharmm_c/Status.h"
#include "apocharmm_c/detail/CudaLangevinPistonIntegratorHandle.h"
#include "catch.hpp"

#include <limits>
#include <string>

namespace {

constexpr double TIME_STEP = 0.001;

apo_cuda_langevin_piston_integrator *CreateIntegrator(void) {
  apo_cuda_langevin_piston_integrator *integrator = nullptr;

  const apo_status status =
      apo_cuda_langevin_piston_integrator_create(&integrator, TIME_STEP);

  REQUIRE(status == APO_STATUS_OK);
  REQUIRE(integrator != nullptr);

  return integrator;
}

} // namespace

TEST_CASE("CapiCudaLangevinPistonConstructionAndDestruction") {
  SECTION("NullOutputPointer") {
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW(status = apo_cuda_langevin_piston_integrator_create(
                      nullptr, TIME_STEP));

    apo_test::CheckStatusAndDiagnostic(
        status, APO_STATUS_INVALID_ARGUMENT,
        "apo_cuda_langevin_piston_integrator_create: out pointer is NULL");
  }

  SECTION("ZeroTimeStep") {
    apo_cuda_langevin_piston_integrator *integrator = nullptr;
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW(
        status = apo_cuda_langevin_piston_integrator_create(&integrator, 0.0));

    CHECK(integrator == nullptr);

    apo_test::CheckNativeError(status, APO_STATUS_INVALID_ARGUMENT,
                               "InvalidArgument",
                               "apo_cuda_langevin_piston_integrator_create",
                               "Time step must be positive; observed 0.000000",
                               "src/CudaIntegrator.cu", "setTimeStep");
  }

  SECTION("InfiniteTimeStep") {
    const double infinity = std::numeric_limits<double>::infinity();

    apo_cuda_langevin_piston_integrator *integrator = nullptr;
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW(status = apo_cuda_langevin_piston_integrator_create(
                      &integrator, infinity));

    CHECK(integrator == nullptr);

    apo_test::CheckNativeError(
        status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
        "apo_cuda_langevin_piston_integrator_create",
        "Time step must be finite; observed " + std::to_string(infinity),
        "src/CudaIntegrator.cu", "setTimeStep");
  }

  SECTION("DestroyIsNoThrow") {
    apo_cuda_langevin_piston_integrator *integrator = CreateIntegrator();

    CHECK_NOTHROW(apo_cuda_langevin_piston_integrator_destroy(integrator));
    CHECK_NOTHROW(apo_cuda_langevin_piston_integrator_destroy(nullptr));
  }
}

TEST_CASE("CapiCudaLangevinPistonValidation") {
  const double infinity = std::numeric_limits<double>::infinity();

  apo_cuda_langevin_piston_integrator *integrator = CreateIntegrator();

  apo_status status = APO_STATUS_OK;

  CHECK_NOTHROW(status =
                    apo_cuda_langevin_piston_integrator_set_reference_pressure(
                        integrator, nullptr, 9));

  apo_test::CheckStatusAndDiagnostic(
      status, APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_langevin_piston_integrator_set_reference_pressure: "
      "pressure_tensor is NULL");

  double pressure[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  CHECK_NOTHROW(status =
                    apo_cuda_langevin_piston_integrator_set_reference_pressure(
                        integrator, pressure, 8));

  apo_test::CheckStatusAndDiagnostic(
      status, APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_langevin_piston_integrator_set_reference_pressure: "
      "pressure_tensor must contain exactly 9 elements");

  pressure[4] = infinity;

  CHECK_NOTHROW(status =
                    apo_cuda_langevin_piston_integrator_set_reference_pressure(
                        integrator, pressure, 9));

  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_cuda_langevin_piston_integrator_set_reference_pressure",
      "Reference pressure tensor values must be finite; non-finite value at "
      "index 4",
      "src/CudaLangevinPistonIntegrator.cu", "setReferencePressure");

  double mass[2] = {1.0, 2.0};

  CHECK_NOTHROW(
      status = apo_cuda_langevin_piston_integrator_set_langevin_piston_mass(
          integrator, mass, 1));

  apo_test::CheckNativeError(
      status, APO_STATUS_NOT_INITIALIZED, "NotInitialized",
      "apo_cuda_langevin_piston_integrator_set_langevin_piston_mass",
      "Crystal type must be set before Langevin piston mass",
      "src/CudaLangevinPistonIntegrator.cu", "setLangevinPistonMass");

  CHECK_NOTHROW(status = apo_cuda_langevin_piston_integrator_set_crystal_type(
                    integrator, APO_CRYSTAL_TYPE_CUBIC));
  REQUIRE(status == APO_STATUS_OK);

  CHECK_NOTHROW(
      status = apo_cuda_langevin_piston_integrator_set_langevin_piston_mass(
          integrator, mass, 2));

  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_cuda_langevin_piston_integrator_set_langevin_piston_mass",
      "Langevin piston mass must contain 1 element(s); observed 2",
      "src/CudaLangevinPistonIntegrator.cu", "setLangevinPistonMass");

  mass[0] = -1.0;

  CHECK_NOTHROW(
      status = apo_cuda_langevin_piston_integrator_set_langevin_piston_mass(
          integrator, mass, 1));

  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_cuda_langevin_piston_integrator_set_langevin_piston_mass",
      "Langevin piston mass values must be non-negative; observed "
      "-1.000000 at index 0",
      "src/CudaLangevinPistonIntegrator.cu", "setLangevinPistonMass");

  CHECK_NOTHROW(
      status = apo_cuda_langevin_piston_integrator_set_langevin_piston_friction(
          integrator, -1.0));

  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_cuda_langevin_piston_integrator_set_langevin_piston_friction",
      "Langevin piston friction must be non-negative; observed -1.000000",
      "src/CudaLangevinPistonIntegrator.cu", "setLangevinPistonFriction");

  double temperature = 0.0;

  CHECK_NOTHROW(
      status =
          apo_cuda_langevin_piston_integrator_get_instantaneous_temperature(
              &temperature, integrator));

  apo_test::CheckNativeError(
      status, APO_STATUS_NOT_INITIALIZED, "NotInitialized",
      "apo_cuda_langevin_piston_integrator_get_instantaneous_temperature",
      "CharmmContext must be set before computing instantaneous temperature",
      "src/CudaLangevinPistonIntegrator.cu", "getInstantaneousTemperature");

  CHECK_NOTHROW(apo_cuda_langevin_piston_integrator_destroy(integrator));
}

TEST_CASE("CapiCudaLangevinPistonGetterInvariants") {
  apo_cuda_langevin_piston_integrator *integrator = CreateIntegrator();

  apo_status status = APO_STATUS_OK;

  SECTION("NoseHooverPistonMass") {
    integrator->object->getNoseHooverPistonMass().resize(0);
    double value = 0.0;

    CHECK_NOTHROW(
        status =
            apo_cuda_langevin_piston_integrator_get_nose_hoover_piston_mass(
                &value, integrator));

    apo_test::CheckStatusAndDiagnostic(
        status, APO_STATUS_RUNTIME_ERROR,
        "apo_cuda_langevin_piston_integrator_get_nose_hoover_piston_mass: "
        "Nose-Hoover piston mass does not contain exactly 1 element");
  }

  SECTION("AverageTemperature") {
    integrator->object->getAverageTemperature().resize(1);
    double value = 0.0;

    CHECK_NOTHROW(
        status = apo_cuda_langevin_piston_integrator_get_average_temperature(
            &value, integrator));

    apo_test::CheckStatusAndDiagnostic(
        status, APO_STATUS_RUNTIME_ERROR,
        "apo_cuda_langevin_piston_integrator_get_average_temperature: "
        "Average temperature does not contain exactly 2 elements");
  }

  SECTION("InstantaneousPressureScalar") {
    integrator->object->getInstantaneousPressureScalar().resize(0);
    double value = 0.0;

    CHECK_NOTHROW(
        status =
            apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_scalar(
                &value, integrator));

    apo_test::CheckStatusAndDiagnostic(
        status, APO_STATUS_RUNTIME_ERROR,
        "apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_scalar:"
        " "
        "Instantaneous pressure scalar does not contain exactly 1 element");
  }

  SECTION("AveragePressureScalar") {
    integrator->object->getAveragePressureScalar().resize(0);
    double value = 0.0;

    CHECK_NOTHROW(
        status =
            apo_cuda_langevin_piston_integrator_get_average_pressure_scalar(
                &value, integrator));

    apo_test::CheckStatusAndDiagnostic(
        status, APO_STATUS_RUNTIME_ERROR,
        "apo_cuda_langevin_piston_integrator_get_average_pressure_scalar: "
        "Average pressure scalar does not contain exactly 1 element");
  }

  CHECK_NOTHROW(apo_cuda_langevin_piston_integrator_destroy(integrator));
}

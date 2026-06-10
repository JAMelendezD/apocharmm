// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "apocharmm_c/CudaLangevinThermostatIntegrator.h"
#include "apocharmm_c/detail/CudaLangevinThermostatIntegratorHandle.h"
#include "apocharmm_c/detail/ErrorInternal.h"
#include "apocharmm_c/detail/SubscriberHandle.h"
#include "apocharmm_c/detail/Validation.h"

#include <memory>

extern "C" apo_status apo_cuda_langevin_thermostat_integrator_create(
    apo_cuda_langevin_thermostat_integrator **out, const double time_step) {
  const char *function_name = "apo_cuda_langevin_thermostat_integrator_create";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::prepare_output_pointer<
                                    apo_cuda_langevin_thermostat_integrator>(
            out, function_name, "out"));

        std::unique_ptr<apo_cuda_langevin_thermostat_integrator> handle(
            new apo_cuda_langevin_thermostat_integrator());
        handle->object =
            std::make_shared<CudaLangevinThermostatIntegrator>(time_step);
        handle->base.object = handle->object;

        *out = handle.release();

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" void apo_cuda_langevin_thermostat_integrator_destroy(
    apo_cuda_langevin_thermostat_integrator *integrator) {
  delete integrator;
  return;
}

extern "C" apo_status
apo_cuda_langevin_thermostat_integrator_set_reference_temperature(
    apo_cuda_langevin_thermostat_integrator *integrator,
    const double temperature) {
  const char *function_name =
      "apo_cuda_langevin_thermostat_integrator_set_reference_temperature";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_handle_object<
                                    apo_cuda_langevin_thermostat_integrator>(
            integrator, function_name, "CudaLangevinThermostatIntegrator"));

        integrator->object->setReferenceTemperature(temperature);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_cuda_langevin_thermostat_integrator_set_thermostat_friction(
    apo_cuda_langevin_thermostat_integrator *integrator,
    const double friction) {
  const char *function_name =
      "apo_cuda_langevin_thermostat_integrator_set_thermostat_friction";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_handle_object<
                                    apo_cuda_langevin_thermostat_integrator>(
            integrator, function_name, "CudaLangevinThermostatIntegrator"));

        integrator->object->setThermostatFriction(friction);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_cuda_langevin_thermostat_integrator_set_thermostat_rng_seed(
    apo_cuda_langevin_thermostat_integrator *integrator, const uint64_t seed) {
  const char *function_name =
      "apo_cuda_langevin_thermostat_integrator_set_thermostat_rng_seed";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_handle_object<
                                    apo_cuda_langevin_thermostat_integrator>(
            integrator, function_name, "CudaLangevinThermostatIntegrator"));

        integrator->object->setThermostatRngSeed(seed);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_cuda_langevin_thermostat_integrator_reset_average_temperature(
    apo_cuda_langevin_thermostat_integrator *integrator) {
  const char *function_name =
      "apo_cuda_langevin_thermostat_integrator_reset_average_temperature";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_handle_object<
                                    apo_cuda_langevin_thermostat_integrator>(
            integrator, function_name, "CudaLangevinThermostatIntegrator"));

        integrator->object->resetAverageTemperature();

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_cuda_langevin_thermostat_integrator_get_reference_temperature(
    double *temperature,
    const apo_cuda_langevin_thermostat_integrator *integrator) {
  const char *function_name =
      "apo_cuda_langevin_thermostat_integrator_get_reference_temperature";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_handle_object<
                                    apo_cuda_langevin_thermostat_integrator>(
            integrator, function_name, "CudaLangevinThermostatIntegrator"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<double>(
            temperature, function_name, "temperature"));

        *temperature = integrator->object->getReferenceTemperature();

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_cuda_langevin_thermostat_integrator_get_thermostat_friction(
    double *friction,
    const apo_cuda_langevin_thermostat_integrator *integrator) {
  const char *function_name =
      "apo_cuda_langevin_thermostat_integrator_get_thermostat_friction";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_handle_object<
                                    apo_cuda_langevin_thermostat_integrator>(
            integrator, function_name, "CudaLangevinThermostatIntegrator"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<double>(
            friction, function_name, "friction"));

        *friction = integrator->object->getThermostatFriction();

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_cuda_langevin_thermostat_integrator_get_thermostat_rng_seed(
    uint64_t *seed, const apo_cuda_langevin_thermostat_integrator *integrator) {
  const char *function_name =
      "apo_cuda_langevin_thermostat_integrator_get_thermostat_rng_seed";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_handle_object<
                                    apo_cuda_langevin_thermostat_integrator>(
            integrator, function_name, "CudaLangevinThermostatIntegrator"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<uint64_t>(
            seed, function_name, "seed"));

        *seed = integrator->object->getThermostatRngSeed();

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_cuda_langevin_thermostat_integrator_get_average_temperature(
    double *temperature,
    const apo_cuda_langevin_thermostat_integrator *integrator) {
  const char *function_name =
      "apo_cuda_langevin_thermostat_integrator_get_average_temperature";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_handle_object<
                                    apo_cuda_langevin_thermostat_integrator>(
            integrator, function_name, "CudaLangevinThermostatIntegrator"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<double>(
            temperature, function_name, "temperature"));

        integrator->object->getAverageTemperature().transferToHost();

        *temperature = integrator->object->getAverageTemperature()[0];

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_cuda_langevin_thermostat_integrator_get_instantaneous_temperature(
    double *temperature,
    const apo_cuda_langevin_thermostat_integrator *integrator) {
  const char *function_name =
      "apo_cuda_langevin_thermostat_integrator_get_instantaneous_temperature";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_handle_object<
                                    apo_cuda_langevin_thermostat_integrator>(
            integrator, function_name, "CudaLangevinThermostatIntegrator"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<double>(
            temperature, function_name, "temperature"));

        *temperature = integrator->object->getInstantaneousTemperature();

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_cuda_langevin_thermostat_integrator_as_cuda_integrator(
    apo_cuda_integrator **out,
    apo_cuda_langevin_thermostat_integrator *integrator) {
  const char *function_name =
      "apo_cuda_langevin_thermostat_integrator_as_cuda_integrator";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::prepare_output_pointer<apo_cuda_integrator>(
                out, function_name, "out"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_handle_object<
                                    apo_cuda_langevin_thermostat_integrator>(
            integrator, function_name, "CudaLangevinThermostatIntegrator"));

        *out = &integrator->base;

        return APO_STATUS_OK;
      },
      function_name);
}

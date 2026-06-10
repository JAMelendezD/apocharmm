// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_CUDA_LANGEVIN_PISTON_INTEGRATOR_H__
#define __APOCHARMM_C_CUDA_LANGEVIN_PISTON_INTEGRATOR_H__

#include "apocharmm_c/CudaIntegrator.h"
#include "apocharmm_c/Enums.h"
#include "apocharmm_c/Export.h"
#include "apocharmm_c/Status.h"

#include <stddef.h>
#include <stdint.h>

#ifndef __cplusplus
#include <stdbool.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct apo_cuda_langevin_piston_integrator
    apo_cuda_langevin_piston_integrator;

APOCHARMM_C_API apo_status apo_cuda_langevin_piston_integrator_create(
    apo_cuda_langevin_piston_integrator **out, const double time_step);

APOCHARMM_C_API void apo_cuda_langevin_piston_integrator_destroy(
    apo_cuda_langevin_piston_integrator *integrator);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_use_nose_hoover_thermostat(
    apo_cuda_langevin_piston_integrator *integrator, const bool flag);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_set_reference_temperature(
    apo_cuda_langevin_piston_integrator *integrator, const double temperature);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_set_nose_hoover_piston_mass(
    apo_cuda_langevin_piston_integrator *integrator, const double mass);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_use_old_temperature(
    apo_cuda_langevin_piston_integrator *integrator, const bool flag);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_set_reference_pressure(
    apo_cuda_langevin_piston_integrator *integrator,
    const double *pressure_tensor, const size_t pressure_tensor_len);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_set_constant_surface_tension(
    apo_cuda_langevin_piston_integrator *integrator, const bool flag);

APOCHARMM_C_API apo_status apo_cuda_langevin_piston_integrator_set_crystal_type(
    apo_cuda_langevin_piston_integrator *integrator,
    const apo_crystal_type crystal_type);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_set_langevin_piston_mass(
    apo_cuda_langevin_piston_integrator *integrator, const double *mass,
    const size_t mass_len);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_set_langevin_piston_friction_seed(
    apo_cuda_langevin_piston_integrator *integrator, const uint64_t seed);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_set_langevin_piston_friction(
    apo_cuda_langevin_piston_integrator *integrator, const double friction);

APOCHARMM_C_API apo_status apo_cuda_langevin_piston_integrator_reset_averages(
    apo_cuda_langevin_piston_integrator *integrator);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_get_reference_temperature(
    double *temperature, const apo_cuda_langevin_piston_integrator *integrator);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_get_nose_hoover_piston_mass(
    double *mass, const apo_cuda_langevin_piston_integrator *integrator);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_get_average_temperature(
    double *temperature, const apo_cuda_langevin_piston_integrator *integrator);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_get_reference_pressure_tensor(
    double *pressure_tensor, const size_t pressure_tensor_len,
    const apo_cuda_langevin_piston_integrator *integrator);

APOCHARMM_C_API apo_status apo_cuda_langevin_piston_integrator_get_crystal_type(
    apo_crystal_type *crystal_type,
    const apo_cuda_langevin_piston_integrator *integrator);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_get_langevin_piston_mass(
    size_t *num_mass, double *mass, const size_t mass_len,
    const apo_cuda_langevin_piston_integrator *integrator);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_tensor(
    double *pressure_tensor, const size_t pressure_tensor_len,
    const apo_cuda_langevin_piston_integrator *integrator);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_scalar(
    double *pressure_scalar,
    const apo_cuda_langevin_piston_integrator *integrator);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_get_average_pressure_tensor(
    double *pressure_tensor, const size_t pressure_tensor_len,
    const apo_cuda_langevin_piston_integrator *integrator);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_get_average_pressure_scalar(
    double *pressure_scalar,
    const apo_cuda_langevin_piston_integrator *integrator);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_get_langevin_piston_friction_seed(
    uint64_t *seed, const apo_cuda_langevin_piston_integrator *integrator);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_get_instantaneous_temperature(
    double *temperature, const apo_cuda_langevin_piston_integrator *integrator);

APOCHARMM_C_API apo_status
apo_cuda_langevin_piston_integrator_as_cuda_integrator(
    apo_cuda_integrator **out, apo_cuda_langevin_piston_integrator *integrator);

#ifdef __cplusplus
}
#endif

#endif

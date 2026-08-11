// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_CHARMM_CONTEXT_H__
#define __APOCHARMM_C_CHARMM_CONTEXT_H__

#include "apocharmm_c/CharmmCrd.h"
#include "apocharmm_c/CharmmParameters.h"
#include "apocharmm_c/CharmmPsf.h"
#include "apocharmm_c/Export.h"
#include "apocharmm_c/ForceManager.h"
#include "apocharmm_c/Status.h"

#include <stddef.h>
#include <stdint.h>

#ifndef __cplusplus
#include <stdbool.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct apo_charmm_context apo_charmm_context;

APOCHARMM_C_API apo_status apo_charmm_context_create(
    apo_charmm_context **out, const apo_force_manager *force_manager);

APOCHARMM_C_API apo_status apo_charmm_context_create_from_psf_parameters(
    apo_charmm_context **out, const apo_charmm_psf *psf,
    const apo_charmm_parameters *parameters);

APOCHARMM_C_API void apo_charmm_context_destroy(apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_set_prm(
    apo_charmm_context *context, apo_charmm_parameters *parameters);

APOCHARMM_C_API apo_status
apo_charmm_context_set_psf(apo_charmm_context *context, apo_charmm_psf *psf);

APOCHARMM_C_API apo_status apo_charmm_context_set_force_manager(
    apo_charmm_context *context, apo_force_manager *force_manager);

APOCHARMM_C_API apo_status apo_charmm_context_set_coordinates_charges(
    apo_charmm_context *context, const double *xyzq, const size_t xyzq_len);

APOCHARMM_C_API apo_status apo_charmm_context_set_coordinates_from_array(
    apo_charmm_context *context, const double *xyz, const size_t xyz_len);

APOCHARMM_C_API apo_status apo_charmm_context_set_coordinates(
    apo_charmm_context *context, const apo_charmm_crd *crd);

APOCHARMM_C_API apo_status
apo_charmm_context_set_charges(apo_charmm_context *context,
                               const double *charges, const size_t charges_len);

APOCHARMM_C_API apo_status apo_charmm_context_set_velocities_inverse_masses(
    apo_charmm_context *context, const double *xyzm, const size_t xyzm_len);

APOCHARMM_C_API apo_status apo_charmm_context_set_velocities(
    apo_charmm_context *context, const double *xyz, const size_t xyz_len);

APOCHARMM_C_API apo_status
apo_charmm_context_set_velocities_from_charmm_velocity_file(
    apo_charmm_context *context, const char *path);

APOCHARMM_C_API apo_status apo_charmm_context_set_masses(
    apo_charmm_context *context, const double *masses, const size_t masses_len);

APOCHARMM_C_API apo_status apo_charmm_context_set_temperature(
    apo_charmm_context *context, const double temperature);

APOCHARMM_C_API apo_status apo_charmm_context_set_periodic_boundary_condition(
    apo_charmm_context *context, const apo_pbc pbc);

APOCHARMM_C_API apo_status apo_charmm_context_set_box_dimensions(
    apo_charmm_context *context, const double *box_dimensions,
    const size_t box_dimensions_len);

APOCHARMM_C_API apo_status apo_charmm_context_set_random_seed(
    apo_charmm_context *context, const uint64_t seed);

APOCHARMM_C_API
apo_status apo_charmm_context_use_holonomic_constraints(
    apo_charmm_context *context, const bool useHolonomicConstraints);

APOCHARMM_C_API apo_status
apo_charmm_context_set_kappa(apo_charmm_context *context, const double kappa);

APOCHARMM_C_API apo_status
apo_charmm_context_set_cutoff(apo_charmm_context *context, const double cutoff);

APOCHARMM_C_API apo_status
apo_charmm_context_set_ctonnb(apo_charmm_context *context, const double ctonnb);

APOCHARMM_C_API apo_status
apo_charmm_context_set_ctofnb(apo_charmm_context *context, const double ctofnb);

APOCHARMM_C_API apo_status apo_charmm_context_set_fft_grid(
    apo_charmm_context *context, const int *grid, const size_t grid_len);

APOCHARMM_C_API apo_status apo_charmm_context_set_pme_spline_order(
    apo_charmm_context *context, const int order);

APOCHARMM_C_API apo_status apo_charmm_context_set_vdw_type(
    apo_charmm_context *context, const int vdw_type);

APOCHARMM_C_API apo_status apo_charmm_context_get_num_atoms(
    int *num_atoms, const apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_get_num_degrees_of_freedom(
    int *ndegf, const apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_get_coordinates_charges(
    double *xyzq, const size_t xyzq_len, const apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_get_velocity_mass(
    double *xyzm, const size_t xyzm_len, const apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_get_periodic_boundary_condition(
    apo_pbc *pbc, const apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_get_box_dimensions(
    double *box_dimensions, const size_t box_dimensions_len,
    const apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_get_random_seed(
    uint64_t *seed, const apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_get_volume(
    double *volume, const apo_charmm_context *context);

APOCHARMM_C_API apo_status
apo_charmm_context_get_kappa(double *kappa, const apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_get_cutoff(
    double *cutoff, const apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_get_ctonnb(
    double *ctonnb, const apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_get_ctofnb(
    double *ctofnb, const apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_get_fft_grid(
    int *grid, const size_t grid_len, const apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_get_pme_spline_order(
    int *order, const apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_get_vdw_type(
    int *vdw_type, const apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_get_force_manager(
    apo_force_manager **out, const apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_assign_velocities_at_temperature(
    apo_charmm_context *context, const double temperature);

APOCHARMM_C_API apo_status apo_charmm_context_compute_temperature(
    double *temperature, apo_charmm_context *context);

APOCHARMM_C_API apo_status apo_charmm_context_calculate_potential_energy(
    apo_charmm_context *context, const bool reset, const bool print);

#ifdef __cplusplus
}
#endif

#endif

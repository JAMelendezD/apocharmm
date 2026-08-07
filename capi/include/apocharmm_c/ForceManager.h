// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_FORCE_MANAGER_H__
#define __APOCHARMM_C_FORCE_MANAGER_H__

#include "apocharmm_c/CharmmParameters.h"
#include "apocharmm_c/CharmmPsf.h"
#include "apocharmm_c/Enums.h"
#include "apocharmm_c/Export.h"
#include "apocharmm_c/Status.h"

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct apo_force_manager apo_force_manager;

APOCHARMM_C_API apo_status
apo_force_manager_create(apo_force_manager **out, const apo_charmm_psf *psf,
                         const apo_charmm_parameters *parameters);

APOCHARMM_C_API void
apo_force_manager_destroy(apo_force_manager *force_manager);

APOCHARMM_C_API apo_status apo_force_manager_set_box_dimensions(
    apo_force_manager *force_manager, const double *box_dimensions,
    const size_t box_dimensions_len);

APOCHARMM_C_API apo_status apo_force_manager_set_kappa(
    apo_force_manager *force_manager, const double kappa);

APOCHARMM_C_API apo_status apo_force_manager_set_cutoff(
    apo_force_manager *force_manager, const double cutoff);

APOCHARMM_C_API apo_status apo_force_manager_set_ctonnb(
    apo_force_manager *force_manager, const double ctonnb);

APOCHARMM_C_API apo_status apo_force_manager_set_ctofnb(
    apo_force_manager *force_manager, const double ctofnb);

APOCHARMM_C_API apo_status apo_force_manager_set_fft_grid(
    apo_force_manager *force_manager, const int *grid, const size_t grid_len);

APOCHARMM_C_API apo_status apo_force_manager_set_pme_spline_order(
    apo_force_manager *force_manager, const int order);

APOCHARMM_C_API apo_status apo_force_manager_set_periodic_boundary_condition(
    apo_force_manager *force_manager, const apo_pbc pbc);

APOCHARMM_C_API apo_status apo_force_manager_set_vdw_type(
    apo_force_manager *force_manager, const int vdw_type);

APOCHARMM_C_API apo_status apo_force_manager_set_print_energy_decomposition(
    apo_force_manager *force_manager, const bool flag);

APOCHARMM_C_API apo_status apo_force_manager_get_num_atoms(
    int *num_atoms, const apo_force_manager *force_manager);

APOCHARMM_C_API apo_status apo_force_manager_is_initialized(
    bool *flag, const apo_force_manager *force_manager);

APOCHARMM_C_API apo_status apo_force_manager_get_box_dimensions(
    double *box_dimensions, const size_t box_dimensions_len,
    const apo_force_manager *force_manager);

APOCHARMM_C_API apo_status apo_force_manager_get_kappa(
    double *kappa, const apo_force_manager *force_manager);

APOCHARMM_C_API apo_status apo_force_manager_get_cutoff(
    double *cutoff, const apo_force_manager *force_manager);

APOCHARMM_C_API apo_status apo_force_manager_get_ctonnb(
    double *ctonnb, const apo_force_manager *force_manager);

APOCHARMM_C_API apo_status apo_force_manager_get_ctofnb(
    double *ctofnb, const apo_force_manager *force_manager);

APOCHARMM_C_API apo_status apo_force_manager_get_fft_grid(
    int *grid, const size_t grid_len, const apo_force_manager *force_manager);

APOCHARMM_C_API apo_status apo_force_manager_get_pme_spline_order(
    int *order, const apo_force_manager *force_manager);

APOCHARMM_C_API apo_status apo_force_manager_get_periodic_boundary_condition(
    apo_pbc *pbc, const apo_force_manager *force_manager);

APOCHARMM_C_API apo_status apo_force_manager_get_vdw_type(
    int *vdw_type, const apo_force_manager *force_manager);

APOCHARMM_C_API apo_status apo_force_manager_is_composite(
    bool *flag, const apo_force_manager *force_manager);

#ifdef __cplusplus
}
#endif

#endif

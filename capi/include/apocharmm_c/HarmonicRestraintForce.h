// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_HARMONIC_RESTRAINT_FORCE_H__
#define __APOCHARMM_C_HARMONIC_RESTRAINT_FORCE_H__

#include "apocharmm_c/AtomSelection.h"
#include "apocharmm_c/Export.h"
#include "apocharmm_c/ForceManager.h"
#include "apocharmm_c/Status.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct apo_harmonic_restraint_force apo_harmonic_restraint_force;

APOCHARMM_C_API apo_status apo_harmonic_restraint_force_create(
    apo_harmonic_restraint_force **out, const int num_atoms);

APOCHARMM_C_API void
apo_harmonic_restraint_force_destroy(apo_harmonic_restraint_force *restraint);

APOCHARMM_C_API apo_status apo_harmonic_restraint_force_set_selection(
    apo_harmonic_restraint_force *restraint,
    const apo_atom_selection *selection);

APOCHARMM_C_API apo_status apo_harmonic_restraint_force_set_force_constant(
    apo_harmonic_restraint_force *restraint, const double force_constant);

APOCHARMM_C_API apo_status apo_harmonic_restraint_force_set_force_constants(
    apo_harmonic_restraint_force *restraint, const double *force_constants,
    const size_t force_constants_len);

APOCHARMM_C_API apo_status
apo_harmonic_restraint_force_set_reference_coordinates(
    apo_harmonic_restraint_force *restraint,
    const double *reference_coordinates,
    const size_t reference_coordinates_len);

APOCHARMM_C_API apo_status apo_harmonic_restraint_force_set_masses(
    apo_harmonic_restraint_force *restraint, const double *masses,
    const size_t masses_len);

APOCHARMM_C_API apo_status apo_harmonic_restraint_force_set_box_dimensions(
    apo_harmonic_restraint_force *restraint, const double *box_dimensions,
    const size_t box_dimensions_len);

APOCHARMM_C_API apo_status apo_force_manager_subscribe_harmonic_restraint_force(
    apo_force_manager *force_manager, apo_harmonic_restraint_force *restraint,
    const char *force_tag);

APOCHARMM_C_API apo_status
apo_force_manager_unsubscribe_harmonic_restraint_force(
    apo_force_manager *force_manager, apo_harmonic_restraint_force *restraint);

#ifdef __cplusplus
}
#endif

#endif

// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_ATOM_SELECTOR_H__
#define __APOCHARMM_C_ATOM_SELECTOR_H__

#include "apocharmm_c/AtomSelection.h"
#include "apocharmm_c/CharmmPsf.h"
#include "apocharmm_c/Export.h"
#include "apocharmm_c/Status.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct apo_atom_selector apo_atom_selector;

APOCHARMM_C_API apo_status apo_atom_selector_create(apo_atom_selector **out,
                                                    const apo_charmm_psf *psf);

APOCHARMM_C_API void apo_atom_selector_destroy(apo_atom_selector *selector);

APOCHARMM_C_API apo_status apo_atom_selector_select(
    apo_atom_selection **out, const apo_atom_selector *selector,
    const char *selection_string);

#ifdef __cplusplus
}
#endif

#endif

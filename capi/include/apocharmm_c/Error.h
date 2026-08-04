// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_ERROR_H__
#define __APOCHARMM_C_ERROR_H__

#include "apocharmm_c/Export.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Return the diagnostic for the most recent failed apoCHARMM C API call on the
 * current thread.
 *
 * The returned pointer is owned by apoCHARMM and must not be freed or modified.
 * It remains valid until a later apoCHARMM C API call on the same thread clears
 * or replaces the thread-local diagnostic. An empty string means that no error
 * is currently recorded.
 */
APOCHARMM_C_API const char *apo_last_error(void);

#ifdef __cplusplus
}
#endif

#endif

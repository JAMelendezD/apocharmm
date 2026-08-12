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
 * @brief Returns the current thread's most recent C ABI diagnostic.
 *
 * Call this function immediately after a status-returning apoCHARMM C function
 * returns a value other than `APO_STATUS_OK`. The function does not clear or
 * otherwise modify the diagnostic. A non-NULL pointer to an empty string means
 * that no diagnostic is currently recorded.
 *
 * The returned pointer is borrowed from apoCHARMM and refers to a
 * null-terminated byte string. A later status-returning apoCHARMM C API call on
 * the same thread clears or replaces the storage, even when that later call
 * succeeds. Successful destruction preserves an existing diagnostic, while a
 * destruction failure may replace it.
 *
 * @return A non-NULL borrowed pointer to storage owned by apoCHARMM. The
 * pointer remains valid until the diagnostic is cleared or replaced on the same
 * thread, or until that thread exits.
 *
 * @note Diagnostics are thread-local; calls on another thread neither expose
 * nor invalidate this thread's diagnostic.
 * @warning Do not free or modify the returned bytes. Copy them before another C
 * ABI call can invalidate the pointer, and do not use the pointer after the
 * originating thread exits.
 * @see apo_status
 */
APOCHARMM_C_API const char *apo_last_error(void);

#ifdef __cplusplus
}
#endif

#endif

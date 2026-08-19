// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

/**
 * @file
 * @brief Declares the C ABI CHARMM parameter-set interface.
 */

#ifndef __APOCHARMM_C_CHARMM_PARAMETERS_H__
#define __APOCHARMM_C_CHARMM_PARAMETERS_H__

#include "apocharmm_c/Export.h"
#include "apocharmm_c/Status.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Represents an opaque owning CHARMM parameter-set handle.
 *
 * A handle created by @ref apo_charmm_parameters_create or
 * @ref apo_charmm_parameters_create_from_files owns shared native parameter
 * state. Release the handle exactly once with
 * @ref apo_charmm_parameters_destroy. Inputs named `parameters` in other C API
 * functions are borrowed for the duration of the call.
 *
 * A ForceManager created from this handle retains independent shared ownership
 * of the native parameter set. Destroying this source handle therefore does not
 * invalidate an already created manager.
 *
 * Status-returning functions clear the calling thread's previous native
 * diagnostic on entry. A successful call leaves that diagnostic empty. A
 * failed call leaves nonempty thread-local text available through
 * @ref apo_last_error until the next guarded C ABI call on the same thread or
 * until thread exit. Handles provide no internal synchronization; callers must
 * not destroy a handle concurrently with another operation that uses it.
 */
typedef struct apo_charmm_parameters apo_charmm_parameters;

/**
 * @brief Creates an owning parameter-set handle from one file.
 *
 * The input path is copied during the call. The native object parses the file
 * on the host before the function returns and retains parsed data rather than
 * the borrowed C string. The file-name extension is not validated.
 *
 * @param[out] out Non-null output slot. It is set to `NULL` before path
 * validation and receives a newly owned handle on success.
 * @param[in] path Borrowed non-null, non-empty, null-terminated path to one
 * CHARMM `.prm` or `.str` file.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `out` is `NULL`, or if `path` is
 * `NULL` or empty.
 * @return `APO_STATUS_RUNTIME_ERROR` if the file cannot be opened or read, a
 * required TOPPAR parameter block is absent, a NONBONDED header continuation
 * is unterminated, a recognized record has an invalid shape or numeric field,
 * native allocation fails, or another uncategorized C++ exception is caught.
 *
 * @post On success, `*out` is non-null and owned by the caller.
 * @post On failure after validating `out`, `*out` remains `NULL`.
 * @post Success clears the calling thread's diagnostic; failure leaves a
 * nonempty diagnostic available through @ref apo_last_error.
 */
APOCHARMM_C_API apo_status
apo_charmm_parameters_create(apo_charmm_parameters **out, const char *path);

/**
 * @brief Creates an owning parameter-set handle from an ordered file list.
 *
 * The function copies all `num_paths` strings and parses the files in array
 * order. File-name extensions are not validated. Later files are merged
 * according to the native section-specific
 * duplicate rules described on @ref charmm_parameters. Neither the pointer
 * array nor its strings are retained after return.
 *
 * @param[out] out Non-null output slot. It is set to `NULL` before file-list
 * validation and receives a newly owned handle on success.
 * @param[in] paths Borrowed array of `num_paths` pointers. When `num_paths` is
 * positive, the array and every element must be non-null; every element must
 * address a non-empty, null-terminated path.
 * @param[in] num_paths Number of elements in `paths`; it must be greater than
 * zero.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `out` is `NULL`, `num_paths` is
 * zero, `paths` is `NULL` for a positive count, or any path element is `NULL`
 * or empty.
 * @return `APO_STATUS_RUNTIME_ERROR` if any file cannot be opened or read, a
 * required TOPPAR parameter block is absent, a NONBONDED header continuation
 * is unterminated, a recognized record has an invalid shape or numeric field,
 * native allocation fails, or another uncategorized C++ exception is caught.
 *
 * @post On success, `*out` is non-null and owned by the caller.
 * @post On failure after validating `out`, `*out` remains `NULL`.
 * @post Success clears the calling thread's diagnostic; failure leaves a
 * nonempty diagnostic available through @ref apo_last_error.
 */
APOCHARMM_C_API apo_status apo_charmm_parameters_create_from_files(
    apo_charmm_parameters **out, const char *const *paths,
    const size_t num_paths);

/**
 * @brief Destroys an owning parameter-set handle.
 *
 * @param[in] parameters Owned handle to destroy, or `NULL`. Passing `NULL` is
 * accepted and has no effect.
 *
 * @post Native destruction is contained within a non-throwing C ABI boundary.
 * @post A normal call preserves the calling thread's previous
 * @ref apo_last_error diagnostic rather than clearing it.
 *
 * @note Do not use the handle after this call. Native parameter state can
 * remain alive when another C handle, such as a ForceManager, retains shared
 * ownership.
 */
APOCHARMM_C_API void
apo_charmm_parameters_destroy(apo_charmm_parameters *parameters);

#ifdef __cplusplus
}
#endif

#endif

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
 * @brief Declares the C ABI CHARMM coordinate-file interface.
 */

#ifndef __APOCHARMM_C_CHARMM_CRD_H__
#define __APOCHARMM_C_CHARMM_CRD_H__

#include "apocharmm_c/Export.h"
#include "apocharmm_c/Status.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Represents an opaque owning handle to parsed CHARMM coordinates.
 *
 * A handle returned by @ref apo_charmm_crd_create owns one shared native
 * CharmmCrd object. Release the handle exactly once with
 * @ref apo_charmm_crd_destroy. Getter inputs named `crd` are borrowed for the
 * duration of the call and are never retained.
 *
 * @ref apo_charmm_context_set_coordinates also borrows the handle only for its
 * call and copies the coordinate values into the context. Destroying this
 * handle afterward does not invalidate coordinates already copied into a
 * context.
 *
 * Every status-returning function clears the calling thread's previous
 * diagnostic on entry. Success leaves an empty diagnostic. Failure leaves
 * thread-local text available through @ref apo_last_error until the next
 * guarded C ABI call on the same thread or thread exit. Copy the diagnostic
 * before making another status-returning call.
 *
 * Handles provide no internal synchronization. Do not destroy a handle
 * concurrently with a call that uses it.
 *
 * @see coordinates
 */
typedef struct apo_charmm_crd apo_charmm_crd;

/**
 * @brief Creates an owning coordinate handle from one CHARMM CRD or COR file.
 *
 * Parsing occurs synchronously and stores both double- and single-precision
 * host coordinate representations. The file-name extension is not validated.
 * The path is copied into native temporary storage and is not retained after
 * construction.
 *
 * @param[out] out Non-NULL output slot. It is set to `NULL` before path
 * validation and receives a newly owned handle on success.
 * @param[in] path Borrowed non-NULL, nonempty, null-terminated path. The
 * pointer and pointed-to storage are not retained. An embedded null byte
 * terminates the path at that position.
 *
 * @retval APO_STATUS_OK Parsing succeeded and `*out` received an owned
 * handle.
 * @retval APO_STATUS_INVALID_ARGUMENT `out` is `NULL`, or `path` is `NULL` or
 * empty.
 * @retval APO_STATUS_RUNTIME_ERROR The file cannot be opened or read; the
 * atom count, record layout, or a coordinate field is invalid; the atom count
 * exceeds `INT_MAX`; host allocation fails; or another uncategorized C++
 * exception is caught.
 *
 * @post On success, `*out` is non-NULL and owned by the caller.
 * @post On failure after `out` is validated, `*out` remains `NULL`.
 * @post The returned status and @ref apo_last_error follow the diagnostic
 * contract documented for @ref apo_charmm_crd.
 */
APOCHARMM_C_API apo_status apo_charmm_crd_create(apo_charmm_crd **out,
                                                 const char *path);

/**
 * @brief Destroys an owning coordinate handle.
 *
 * The private shared owner is released. The native object is destroyed when
 * its final shared owner is released.
 *
 * @param[in] crd Owned handle to destroy, or `NULL`. Passing `NULL` is accepted
 * and has no effect.
 *
 * @post No exception crosses the C ABI boundary.
 * @post A normally returning call preserves the calling thread's previous
 * @ref apo_last_error diagnostic instead of clearing it. A caught destruction
 * exception may replace that diagnostic.
 * @warning Do not use a non-NULL handle after this call.
 */
APOCHARMM_C_API void apo_charmm_crd_destroy(apo_charmm_crd *crd);

/**
 * @brief Returns the stored atom count.
 *
 * @param[out] num_atoms Non-NULL output slot. It is set to zero before handle
 * validation and receives the non-negative dimensionless atom count on
 * success.
 * @param[in] crd Borrowed live coordinate handle containing a native object.
 *
 * @retval APO_STATUS_OK `*num_atoms` received the atom count.
 * @retval APO_STATUS_INVALID_ARGUMENT `num_atoms` is `NULL`, `crd` is `NULL`,
 * or the handle contains no native object.
 * @retval APO_STATUS_RUNTIME_ERROR The native object reports a negative atom
 * count.
 *
 * @post On success, `*num_atoms` contains the non-negative atom count.
 * @post On failure after `num_atoms` is validated, `*num_atoms` remains zero.
 * @post The returned status and @ref apo_last_error follow the diagnostic
 * contract documented for @ref apo_charmm_crd.
 */
APOCHARMM_C_API apo_status
apo_charmm_crd_get_num_atoms(size_t *num_atoms, const apo_charmm_crd *crd);

/**
 * @brief Copies double-precision Cartesian coordinates into a caller buffer.
 *
 * The output contains `N` consecutive `[x, y, z]` records in atom order, where
 * `N` is the coordinate-vector length. Components use angstroms. The function
 * performs a host-to-host copy and does not allocate, transfer CUDA memory, or
 * synchronize a CUDA stream.
 *
 * @param[out] xyz Borrowed contiguous `double` output buffer. It must provide
 * at least `3 * N` elements. It may be `NULL` only when `N` is zero. The caller
 * retains ownership.
 * @param[in] xyz_len Number of `double` elements available in `xyz`, not a byte
 * count. Values greater than `3 * N` are accepted.
 * @param[in] crd Borrowed live coordinate handle containing a native object.
 *
 * @retval APO_STATUS_OK The first `3 * N` elements were written.
 * @retval APO_STATUS_INVALID_ARGUMENT `crd` is `NULL`, the handle contains no
 * native object, `N` is nonzero and `xyz` is `NULL`, or `xyz_len` is less than
 * `3 * N`.
 *
 * @post On success, the first `3 * N` elements contain atom-major Cartesian
 * coordinates. Any trailing elements remain unchanged.
 * @post On failure, the caller buffer remains unchanged and may still contain
 * stale values.
 * @post The returned status and @ref apo_last_error follow the diagnostic
 * contract documented for @ref apo_charmm_crd.
 */
APOCHARMM_C_API apo_status apo_charmm_crd_get_coordinates(
    double *xyz, const size_t xyz_len, const apo_charmm_crd *crd);

#ifdef __cplusplus
}
#endif

#endif

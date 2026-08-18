// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_ATOM_SELECTION_H__
#define __APOCHARMM_C_ATOM_SELECTION_H__

#include "apocharmm_c/Export.h"
#include "apocharmm_c/Status.h"

#include <stddef.h>

#ifndef __cplusplus
#include <stdbool.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Represents an owned, immutable atom-selection result in the C ABI.
 *
 * A successful @ref apo_atom_selector_select call creates one handle. The
 * handle owns shared native selection storage until
 * @ref apo_atom_selection_destroy is called. Accessors borrow the live handle
 * for the duration of one call and never transfer ownership.
 *
 * Atom counts and indices are dimensionless. Indices are zero-based and are
 * returned in ascending order. The current C API exposes no mutating operation,
 * CUDA storage, transfer, stream, or synchronization for a selection.
 *
 * A non-NULL pointer must designate a live handle created by apoCHARMM. Passing
 * a dangling, already-destroyed, or foreign pointer has undefined behavior.
 * Do not overlap destruction with an accessor on another thread.
 *
 * Status-returning accessors clear the calling thread's previous diagnostic at
 * entry. On failure, call @ref apo_last_error immediately on the same thread;
 * its borrowed pointer remains valid only until another diagnostic-changing C
 * ABI call on that thread.
 *
 * @see atom_selection
 */
typedef struct apo_atom_selection apo_atom_selection;

/**
 * @brief Destroys an owned atom-selection handle.
 *
 * @param[in] selection Owned handle to release. `NULL` is accepted and is a
 * no-op. A non-NULL pointer is invalid after this call returns.
 *
 * @post No C++ exception escapes the C ABI boundary.
 * @note A normal destruction preserves the calling thread's existing
 * @ref apo_last_error diagnostic. An internal destruction failure cannot be
 * returned by this void API and may replace that diagnostic.
 * @warning The caller must not destroy the same handle twice or use it after
 * destruction.
 */
APOCHARMM_C_API void apo_atom_selection_destroy(apo_atom_selection *selection);

/**
 * @brief Returns the number of atoms represented by a selection.
 *
 * @param[out] num_atoms Non-NULL output pointer. When this pointer is valid,
 * the function stores zero before validating `selection`, then stores the
 * dimensionless atom count on success.
 * @param[in] selection Borrowed live selection handle. The pointer may not be
 * `NULL` and is not retained.
 * @retval APO_STATUS_OK The atom count was written to `*num_atoms`.
 * @retval APO_STATUS_INVALID_ARGUMENT `num_atoms` is `NULL`, `selection` is
 * `NULL`, or the handle contains no native selection object.
 * @retval APO_STATUS_RUNTIME_ERROR An unexpected standard or nonstandard C++
 * exception crossed the accessor boundary.
 *
 * @post On failure after a valid `num_atoms` pointer is accepted,
 * `*num_atoms == 0`.
 * @note The function clears the previous thread-local diagnostic at entry.
 * Success leaves it empty; failure leaves text available through
 * @ref apo_last_error.
 */
APOCHARMM_C_API apo_status apo_atom_selection_get_num_atoms(
    size_t *num_atoms, const apo_atom_selection *selection);

/**
 * @brief Returns the number of selected atoms.
 *
 * @param[out] num_selected Non-NULL output pointer. When this pointer is valid,
 * the function stores zero before validating `selection`, then stores a
 * dimensionless value in `[0, atom_count]` on success.
 * @param[in] selection Borrowed live selection handle. The pointer may not be
 * `NULL` and is not retained.
 * @retval APO_STATUS_OK The selected-atom count was written to
 * `*num_selected`.
 * @retval APO_STATUS_INVALID_ARGUMENT `num_selected` is `NULL`, `selection` is
 * `NULL`, or the handle contains no native selection object.
 * @retval APO_STATUS_RUNTIME_ERROR An unexpected standard or nonstandard C++
 * exception crossed the accessor boundary.
 *
 * @post On failure after a valid `num_selected` pointer is accepted,
 * `*num_selected == 0`.
 * @note The function clears the previous thread-local diagnostic at entry.
 * Success leaves it empty; failure leaves text available through
 * @ref apo_last_error.
 */
APOCHARMM_C_API apo_status apo_atom_selection_get_num_selected(
    size_t *num_selected, const apo_atom_selection *selection);

/**
 * @brief Copies all selected atom indices into a caller-owned buffer.
 *
 * The required element count is the value returned by
 * @ref apo_atom_selection_get_num_selected. Exactly that many `int` elements
 * are written in strictly ascending, zero-based order. Elements at offsets
 * greater than or equal to the required count are left unchanged.
 *
 * @param[out] atom_indices Buffer receiving zero-based atom indices. A
 * non-NULL buffer must contain at least `len` writable `int` elements. The
 * pointer may be `NULL` only when the selection contains zero atoms.
 * @param[in] len Dimensionless capacity of `atom_indices`, measured in `int`
 * elements. The value may exceed the required count.
 * @param[in] selection Borrowed live selection handle. The pointer may not be
 * `NULL` and is not retained.
 * @retval APO_STATUS_OK Every selected index was copied. For an empty
 * selection, `atom_indices == NULL` is accepted for any `len` value and no
 * element is written.
 * @retval APO_STATUS_INVALID_ARGUMENT `selection` is `NULL`, the handle
 * contains no native object, a nonempty result has a `NULL` buffer, or `len` is
 * smaller than the selected-atom count.
 * @retval APO_STATUS_RUNTIME_ERROR Allocating the temporary native index vector
 * failed, or another unexpected standard or nonstandard C++ exception crossed
 * the boundary.
 *
 * @post On any verified failure path, caller buffer contents remain unchanged.
 * @note The function clears the previous thread-local diagnostic at entry.
 * Success leaves it empty; failure leaves text available through
 * @ref apo_last_error.
 */
APOCHARMM_C_API apo_status apo_atom_selection_get_atom_indices(
    int *atom_indices, const size_t len, const apo_atom_selection *selection);

/**
 * @brief Tests whether one zero-based atom index is selected.
 *
 * @param[out] is_selected Non-NULL output pointer. The function stores `false`
 * before validating the handle or index, then stores the query result on
 * success.
 * @param[in] selection Borrowed live selection handle. The pointer may not be
 * `NULL` and is not retained.
 * @param[in] atom_index Zero-based, dimensionless atom index. The value must be
 * in `[0, atom_count)`.
 * @retval APO_STATUS_OK `*is_selected` contains the selection state.
 * @retval APO_STATUS_INVALID_ARGUMENT `is_selected` is `NULL`, `selection` is
 * `NULL`, the handle contains no native object, or `atom_index` is outside the
 * represented atom range.
 * @retval APO_STATUS_RUNTIME_ERROR An unexpected standard or nonstandard C++
 * exception crossed the query boundary.
 *
 * @post On failure after a valid `is_selected` pointer is accepted,
 * `*is_selected == false`.
 * @note The function clears the previous thread-local diagnostic at entry.
 * Success leaves it empty; failure leaves text available through
 * @ref apo_last_error.
 */
APOCHARMM_C_API apo_status apo_atom_selection_contains(
    bool *is_selected, const apo_atom_selection *selection,
    const int atom_index);

#ifdef __cplusplus
}
#endif

#endif

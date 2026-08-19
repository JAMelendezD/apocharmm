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
 * @brief Declares the C ABI CHARMM PSF interface.
 */

#ifndef __APOCHARMM_C_CHARMM_PSF_H__
#define __APOCHARMM_C_CHARMM_PSF_H__

#include "apocharmm_c/Export.h"
#include "apocharmm_c/Status.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Represents an opaque owning handle to parsed CHARMM PSF state.
 *
 * A handle returned by @ref apo_charmm_psf_create owns one shared native
 * CharmmPSF object. Release the handle exactly once with
 * @ref apo_charmm_psf_destroy. Inputs named `psf` in other C ABI functions are
 * borrowed for the duration of the call and are never retained by those getter
 * functions.
 *
 * ForceManager and CharmmContext handles created from this PSF copy native
 * shared ownership. Destroying the source PSF handle therefore does not
 * invalidate those already-created collaborators.
 *
 * Every status-returning function clears the calling thread's previous
 * diagnostic on entry. Success leaves an empty diagnostic. Failure leaves
 * nonempty thread-local text available through @ref apo_last_error until the
 * next guarded C ABI call on the same thread or thread exit. Copy the
 * diagnostic before making another status-returning call.
 *
 * Handles provide no internal synchronization. Do not destroy a handle
 * concurrently with a call that uses it.
 *
 * @see charmm_psf
 */
typedef struct apo_charmm_psf apo_charmm_psf;

/**
 * @brief Creates an owning PSF handle from one file.
 *
 * The borrowed C string is copied into a native string before return. Parsing
 * occurs synchronously. The file-name extension is not validated. The native
 * parser requires the supported PSF sections described on @ref charmm_psf and
 * constructs derived residue, water, connected-component, and exclusion state.
 *
 * @param[out] out Non-NULL output slot. It is set to `NULL` before path
 * validation and receives a newly owned handle on success.
 * @param[in] path Borrowed non-NULL, non-empty, null-terminated PSF path. The
 * pointer and its storage are not retained.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `out` is `NULL`, or if `path` is
 * `NULL` or empty.
 * @return `APO_STATUS_RUNTIME_ERROR` if the file cannot be opened or read, a
 * required section or record is missing, a numeric field or record shape is
 * invalid, a topology index or count is unsupported, host allocation fails, or
 * another uncategorized C++ exception is caught.
 * @return `APO_STATUS_CUDA_ERROR` if construction of derived CudaContainer
 * state encounters a CUDA allocation, copy, launch, or cleanup failure.
 *
 * @post On success, `*out` is non-NULL and owned by the caller.
 * @post On failure after `out` is validated, `*out` remains `NULL`.
 * @post The returned status and the calling thread's @ref apo_last_error follow
 * the diagnostic contract documented for @ref apo_charmm_psf.
 */
APOCHARMM_C_API apo_status apo_charmm_psf_create(apo_charmm_psf **out,
                                                 const char *path);

/**
 * @brief Destroys an owning PSF handle.
 *
 * @param[in] psf Owned handle to destroy, or `NULL`. Passing `NULL` is accepted
 * and has no effect.
 *
 * @post No exception crosses the C ABI boundary.
 * @post A normally returning call preserves the calling thread's previous
 * @ref apo_last_error diagnostic instead of clearing it.
 *
 * @note Do not use the handle after this call. The underlying native PSF can
 * remain alive while a previously created ForceManager or CharmmContext retains
 * shared ownership.
 */
APOCHARMM_C_API void apo_charmm_psf_destroy(apo_charmm_psf *psf);

/**
 * @brief Returns the atom count.
 *
 * @param[out] num_atoms Non-NULL output slot. It is set to zero before handle
 * validation and receives the dimensionless atom count on success.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `num_atoms` is `NULL`, `psf` is
 * `NULL`, or the handle has no native object.
 * @return `APO_STATUS_RUNTIME_ERROR` if the native atom count is negative or an
 * uncategorized C++ exception is caught.
 *
 * @post On success, `*num_atoms` contains the nonnegative atom count.
 * @post On failure after `num_atoms` is validated, `*num_atoms` remains zero.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 */
APOCHARMM_C_API apo_status
apo_charmm_psf_get_num_atoms(size_t *num_atoms, const apo_charmm_psf *psf);

/**
 * @brief Returns the bond count.
 *
 * @param[out] num_bonds Non-NULL output slot. It is set to zero before handle
 * validation and receives the dimensionless bond count on success.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `num_bonds` is `NULL`, `psf` is
 * `NULL`, or the handle has no native object.
 * @return `APO_STATUS_RUNTIME_ERROR` if the native bond count is negative or an
 * uncategorized C++ exception is caught.
 *
 * @post On success, `*num_bonds` contains the nonnegative bond count.
 * @post On failure after `num_bonds` is validated, `*num_bonds` remains zero.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 */
APOCHARMM_C_API apo_status
apo_charmm_psf_get_num_bonds(size_t *num_bonds, const apo_charmm_psf *psf);

/**
 * @brief Returns the angle count.
 *
 * @param[out] num_angles Non-NULL output slot. It is set to zero before handle
 * validation and receives the dimensionless angle count on success.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `num_angles` is `NULL`, `psf` is
 * `NULL`, or the handle has no native object.
 * @return `APO_STATUS_RUNTIME_ERROR` if the native angle count is negative or
 * an uncategorized C++ exception is caught.
 *
 * @post On success, `*num_angles` contains the nonnegative angle count.
 * @post On failure after `num_angles` is validated, `*num_angles` remains zero.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 */
APOCHARMM_C_API apo_status
apo_charmm_psf_get_num_angles(size_t *num_angles, const apo_charmm_psf *psf);

/**
 * @brief Returns the proper-dihedral count.
 *
 * @param[out] num_dihedrals Non-NULL output slot. It is set to zero before
 * handle validation and receives the dimensionless proper-dihedral count on
 * success.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `num_dihedrals` is `NULL`, `psf` is
 * `NULL`, or the handle has no native object.
 * @return `APO_STATUS_RUNTIME_ERROR` if the native proper-dihedral count is
 * negative or an uncategorized C++ exception is caught.
 *
 * @post On success, `*num_dihedrals` contains the nonnegative count.
 * @post On failure after `num_dihedrals` is validated, `*num_dihedrals` remains
 * zero.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 */
APOCHARMM_C_API apo_status apo_charmm_psf_get_num_dihedrals(
    size_t *num_dihedrals, const apo_charmm_psf *psf);

/**
 * @brief Returns the improper-dihedral count.
 *
 * @param[out] num_impropers Non-NULL output slot. It is set to zero before
 * handle validation and receives the dimensionless improper-dihedral count on
 * success.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `num_impropers` is `NULL`, `psf` is
 * `NULL`, or the handle has no native object.
 * @return `APO_STATUS_RUNTIME_ERROR` if the native improper-dihedral count is
 * negative or an uncategorized C++ exception is caught.
 *
 * @post On success, `*num_impropers` contains the nonnegative count.
 * @post On failure after `num_impropers` is validated, `*num_impropers` remains
 * zero.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 */
APOCHARMM_C_API apo_status apo_charmm_psf_get_num_impropers(
    size_t *num_impropers, const apo_charmm_psf *psf);

/**
 * @brief Returns the CMAP cross-term count.
 *
 * @param[out] num_cross_terms Non-NULL output slot. It is set to zero before
 * handle validation and receives the dimensionless cross-term count on success.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `num_cross_terms` is `NULL`, `psf`
 * is `NULL`, or the handle has no native object.
 * @return `APO_STATUS_RUNTIME_ERROR` if the native cross-term count is negative
 * or an uncategorized C++ exception is caught.
 *
 * @post On success, `*num_cross_terms` contains the nonnegative count.
 * @post On failure after `num_cross_terms` is validated, `*num_cross_terms`
 * remains zero.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 */
APOCHARMM_C_API apo_status apo_charmm_psf_get_num_cross_terms(
    size_t *num_cross_terms, const apo_charmm_psf *psf);

/**
 * @brief Copies fixed-width per-atom segment identifiers.
 *
 * @param[out] segis Borrowed caller-owned output buffer containing at least
 * `8 * num_atoms` `char` elements. Atom `i` occupies
 * `segis[8 * i]` through `segis[8 * i + 7]`. Fields are right-padded with ASCII
 * spaces, are not null-terminated, and are silently truncated to eight bytes.
 * `NULL` is accepted only when the native vector is empty.
 * @param[in] len Number of `char` elements available in `segis`.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `psf` is invalid, or if a nonempty
 * result is requested with `segis == NULL` or `len < 8 * num_atoms`.
 * @return `APO_STATUS_RUNTIME_ERROR` if construction of a validation diagnostic
 * or another uncategorized C++ operation throws.
 *
 * @post On success, the required prefix contains the fixed-width fields; any
 * elements beyond that prefix are unchanged.
 * @post On failure, the caller's buffer contents remain unchanged and can
 * therefore contain stale values.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 */
APOCHARMM_C_API apo_status apo_charmm_psf_get_segment_identifiers(
    char *segis, const size_t len, const apo_charmm_psf *psf);

/**
 * @brief Copies per-atom residue identifiers.
 *
 * @param[out] resis Borrowed caller-owned output buffer containing at least
 * `num_atoms` `int` elements. Values are dimensionless and are written in
 * atom-record order. `NULL` is accepted only when the native vector is empty.
 * @param[in] len Number of `int` elements available in `resis`.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `psf` is invalid, or if a nonempty
 * result is requested with `resis == NULL` or `len < num_atoms`.
 * @return `APO_STATUS_RUNTIME_ERROR` if construction of a validation diagnostic
 * or another uncategorized C++ operation throws.
 *
 * @post On success, the first `num_atoms` elements contain the identifiers; any
 * additional elements are unchanged.
 * @post On failure, the caller's buffer contents remain unchanged and can
 * therefore contain stale values.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 */
APOCHARMM_C_API apo_status apo_charmm_psf_get_residue_identifiers(
    int *resis, const size_t len, const apo_charmm_psf *psf);

/**
 * @brief Copies fixed-width per-atom residue names.
 *
 * @param[out] resns Borrowed caller-owned output buffer containing at least
 * `8 * num_atoms` `char` elements. Atom `i` occupies one eight-byte field.
 * Fields are right-padded with ASCII spaces, are not null-terminated, and are
 * silently truncated to eight bytes. `NULL` is accepted only when the native
 * vector is empty.
 * @param[in] len Number of `char` elements available in `resns`.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `psf` is invalid, or if a nonempty
 * result is requested with `resns == NULL` or `len < 8 * num_atoms`.
 * @return `APO_STATUS_RUNTIME_ERROR` if construction of a validation diagnostic
 * or another uncategorized C++ operation throws.
 *
 * @post On success, the required prefix contains the fixed-width fields; any
 * elements beyond that prefix are unchanged.
 * @post On failure, the caller's buffer contents remain unchanged.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 */
APOCHARMM_C_API apo_status apo_charmm_psf_get_residue_names(
    char *resns, const size_t len, const apo_charmm_psf *psf);

/**
 * @brief Copies fixed-width per-atom atom names.
 *
 * @param[out] anams Borrowed caller-owned output buffer containing at least
 * `8 * num_atoms` `char` elements. Atom `i` occupies one eight-byte field.
 * Fields are right-padded with ASCII spaces, are not null-terminated, and are
 * silently truncated to eight bytes. `NULL` is accepted only when the native
 * vector is empty.
 * @param[in] len Number of `char` elements available in `anams`.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `psf` is invalid, or if a nonempty
 * result is requested with `anams == NULL` or `len < 8 * num_atoms`.
 * @return `APO_STATUS_RUNTIME_ERROR` if construction of a validation diagnostic
 * or another uncategorized C++ operation throws.
 *
 * @post On success, the required prefix contains the fixed-width fields; any
 * elements beyond that prefix are unchanged.
 * @post On failure, the caller's buffer contents remain unchanged.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 */
APOCHARMM_C_API apo_status apo_charmm_psf_get_atom_names(
    char *anams, const size_t len, const apo_charmm_psf *psf);

/**
 * @brief Copies fixed-width per-atom CHARMM atom types.
 *
 * @param[out] atyps Borrowed caller-owned output buffer containing at least
 * `8 * num_atoms` `char` elements. Atom `i` occupies one eight-byte field.
 * Fields are right-padded with ASCII spaces, are not null-terminated, and are
 * silently truncated to eight bytes. `NULL` is accepted only when the native
 * vector is empty.
 * @param[in] len Number of `char` elements available in `atyps`.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `psf` is invalid, or if a nonempty
 * result is requested with `atyps == NULL` or `len < 8 * num_atoms`.
 * @return `APO_STATUS_RUNTIME_ERROR` if construction of a validation diagnostic
 * or another uncategorized C++ operation throws.
 *
 * @post On success, the required prefix contains the fixed-width fields; any
 * elements beyond that prefix are unchanged.
 * @post On failure, the caller's buffer contents remain unchanged.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 */
APOCHARMM_C_API apo_status apo_charmm_psf_get_atom_types(
    char *atyps, const size_t len, const apo_charmm_psf *psf);

/**
 * @brief Copies all per-atom partial charges.
 *
 * @param[out] charges Borrowed caller-owned output buffer containing at least
 * `num_atoms` `double` elements. Values are written in atom-record order and
 * use elementary-charge units. `NULL` is accepted only when the native vector
 * is empty.
 * @param[in] len Number of `double` elements available in `charges`.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `psf` is invalid, or if a nonempty
 * result is requested with `charges == NULL` or `len < num_atoms`.
 * @return `APO_STATUS_RUNTIME_ERROR` if construction of a validation diagnostic
 * or another uncategorized C++ operation throws.
 *
 * @post On success, the first `num_atoms` elements contain copied charges; any
 * additional elements are unchanged.
 * @post On failure, the caller's buffer remains unchanged and can contain stale
 * values.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 */
APOCHARMM_C_API apo_status apo_charmm_psf_get_charges(
    double *charges, const size_t len, const apo_charmm_psf *psf);

/**
 * @brief Copies all per-atom masses.
 *
 * @param[out] masses Borrowed caller-owned output buffer containing at least
 * `num_atoms` `double` elements. Values are written in atom-record order and
 * use atomic mass units. `NULL` is accepted only when the native vector is
 * empty.
 * @param[in] len Number of `double` elements available in `masses`.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `psf` is invalid, or if a nonempty
 * result is requested with `masses == NULL` or `len < num_atoms`.
 * @return `APO_STATUS_RUNTIME_ERROR` if construction of a validation diagnostic
 * or another uncategorized C++ operation throws.
 *
 * @post On success, the first `num_atoms` elements contain copied masses; any
 * additional elements are unchanged.
 * @post On failure, the caller's buffer remains unchanged and can contain stale
 * values.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 */
APOCHARMM_C_API apo_status apo_charmm_psf_get_masses(double *masses,
                                                     const size_t len,
                                                     const apo_charmm_psf *psf);

/**
 * @brief Computes and returns the net partial charge.
 *
 * @param[out] net_charge Non-NULL output slot. It is set to `0.0` before handle
 * validation and receives the sum in elementary-charge units on success.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `net_charge` is `NULL`, `psf` is
 * `NULL`, or the handle has no native object.
 * @return `APO_STATUS_NOT_INITIALIZED` if the native atom count is negative.
 * @return `APO_STATUS_RUNTIME_ERROR` if the native charge-vector length differs
 * from the atom count or another uncategorized C++ exception is caught.
 *
 * @post On success, `*net_charge` contains the computed net charge.
 * @post On failure after `net_charge` is validated, `*net_charge` remains
 * `0.0`.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 */
APOCHARMM_C_API apo_status
apo_charmm_psf_get_net_charge(double *net_charge, const apo_charmm_psf *psf);

/**
 * @brief Computes and returns the total mass.
 *
 * @param[out] total_mass Non-NULL output slot. It is set to `0.0` before handle
 * validation and receives the sum in atomic mass units on success.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `total_mass` is `NULL`, `psf` is
 * `NULL`, or the handle has no native object.
 * @return `APO_STATUS_NOT_INITIALIZED` if the native atom count is negative.
 * @return `APO_STATUS_RUNTIME_ERROR` if the native mass-vector length differs
 * from the atom count or another uncategorized C++ exception is caught.
 *
 * @post On success, `*total_mass` contains the computed total mass.
 * @post On failure after `total_mass` is validated, `*total_mass` remains
 * `0.0`.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 */
APOCHARMM_C_API apo_status
apo_charmm_psf_get_total_mass(double *total_mass, const apo_charmm_psf *psf);

/**
 * @brief Copies the stored PSF path into a caller-owned character buffer.
 *
 * @param[out] file_name Borrowed caller-owned buffer containing at least
 * `stored_path_length` `char` elements. The function writes exactly `len`
 * elements: the path bytes first and ASCII spaces in all remaining positions.
 * No null terminator is written. `NULL` is accepted only when the stored path
 * is empty.
 * @param[in] len Number of `char` elements available in `file_name`. The value
 * must be at least the stored path length.
 * @param[in] psf Borrowed live PSF handle.
 *
 * @return `APO_STATUS_OK` on success.
 * @return `APO_STATUS_INVALID_ARGUMENT` if `psf` is invalid, or if a nonempty
 * path is requested with `file_name == NULL` or `len` smaller than the stored
 * path length.
 * @return `APO_STATUS_RUNTIME_ERROR` if construction of a validation diagnostic
 * or another uncategorized C++ operation throws.
 *
 * @post On success, all `len` output elements have been overwritten.
 * @post On failure, the caller's buffer remains unchanged and can contain stale
 * values.
 * @post Status and @ref apo_last_error follow the @ref apo_charmm_psf
 * diagnostic contract.
 *
 * @warning The current C ABI exposes no companion function that returns the
 * required path-buffer length.
 */
APOCHARMM_C_API apo_status apo_charmm_psf_get_file_name(
    char *file_name, const size_t len, const apo_charmm_psf *psf);

#ifdef __cplusplus
}
#endif

#endif

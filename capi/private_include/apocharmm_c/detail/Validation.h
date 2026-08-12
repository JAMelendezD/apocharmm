// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_DETAIL_VALIDATION_H__
#define __APOCHARMM_C_DETAIL_VALIDATION_H__

#include "apocharmm_c/Status.h"
#include "apocharmm_c/detail/ErrorInternal.h"

#include <cmath>
#include <cstddef>
#include <string>

namespace apocharmm_c {

/**
 * @brief Records an invalid-argument diagnostic for a C ABI entry point.
 *
 * @param[in] function_name Borrowed null-terminated C function label. The
 * pointer may be `nullptr` and is not retained.
 * @param[in] message Human-readable diagnostic copied before return.
 * @return `APO_STATUS_INVALID_ARGUMENT`.
 *
 * @post The current thread's diagnostic contains the selected function label
 * and message.
 */
inline apo_status invalid_argument(const char *function_name,
                                   const std::string &message) {
  return set_last_error(APO_STATUS_INVALID_ARGUMENT, function_name, message);
}

/**
 * @brief Validates and clears an owned-handle output slot before construction.
 *
 * @tparam T Opaque C ABI handle type to be produced by the caller.
 * @param[out] out Borrowed pointer to the caller's handle slot. The pointer
 * must be non-NULL; the slot is set to `nullptr` when validation succeeds.
 * @param[in] function_name Borrowed null-terminated C function label. The
 * pointer may be `nullptr` and is not retained.
 * @param[in] argument_name Borrowed non-NULL, null-terminated argument label
 * used only to format a diagnostic.
 * @return `APO_STATUS_OK` after setting `*out` to `nullptr`, or
 * `APO_STATUS_INVALID_ARGUMENT` when `out` is NULL.
 *
 * @post On success, `*out` is `nullptr`. On failure, no output location exists
 * and no write is performed.
 */
template <typename T>
apo_status prepare_output_pointer(T **out, const char *function_name,
                                  const char *argument_name) {
  if (out == nullptr) {
    return invalid_argument(function_name,
                            std::string(argument_name) + " pointer is NULL");
  }

  *out = nullptr;

  return APO_STATUS_OK;
}

/**
 * @brief Validates that a borrowed C string is non-NULL and nonempty.
 *
 * @param[in] value Borrowed null-terminated string to validate. The pointer may
 * be NULL and is never retained.
 * @param[in] function_name Borrowed null-terminated C function label. The
 * pointer may be `nullptr` and is not retained.
 * @param[in] argument_name Borrowed non-NULL, null-terminated argument label
 * used only to format a diagnostic.
 * @return `APO_STATUS_OK` when `value[0]` is not `\0`, or
 * `APO_STATUS_INVALID_ARGUMENT` when `value` is NULL or empty.
 */
inline apo_status require_c_string(const char *value, const char *function_name,
                                   const char *argument_name) {
  if ((value != nullptr) && (value[0] != '\0'))
    return APO_STATUS_OK;

  return invalid_argument(function_name,
                          std::string(argument_name) + " is NULL or empty");
}

/**
 * @brief Validates that an integer argument is strictly positive.
 *
 * @param[in] value Integer value to validate.
 * @param[in] function_name Borrowed null-terminated C function label. The
 * pointer may be `nullptr` and is not retained.
 * @param[in] argument_name Borrowed non-NULL, null-terminated argument label
 * used only to format a diagnostic.
 * @return `APO_STATUS_OK` when `value` is greater than zero, or
 * `APO_STATUS_INVALID_ARGUMENT` otherwise.
 */
inline apo_status require_positive_int(const int value,
                                       const char *function_name,
                                       const char *argument_name) {
  if (value > 0)
    return APO_STATUS_OK;

  return invalid_argument(function_name,
                          std::string(argument_name) + " must be positive");
}

/**
 * @brief Validates that a borrowed pointer is non-NULL.
 *
 * @tparam T Pointed-to element or handle type.
 * @param[in] pointer Borrowed pointer to validate. The pointer may be NULL and
 * is never retained or dereferenced.
 * @param[in] function_name Borrowed null-terminated C function label. The
 * pointer may be `nullptr` and is not retained.
 * @param[in] argument_name Borrowed non-NULL, null-terminated argument label
 * used only to format a diagnostic.
 * @return `APO_STATUS_OK` when `pointer` is non-NULL, or
 * `APO_STATUS_INVALID_ARGUMENT` otherwise.
 */
template <typename T>
apo_status require_pointer(const T *pointer, const char *function_name,
                           const char *argument_name) {
  if (pointer != nullptr)
    return APO_STATUS_OK;

  return invalid_argument(function_name,
                          std::string(argument_name) + " is NULL");
}

/**
 * @brief Validates an opaque handle and its retained native object.
 *
 * @tparam T C ABI handle struct exposing an `object` member comparable to
 * `nullptr`.
 * @param[in] handle Borrowed handle to validate. The pointer may be NULL and is
 * not retained; its native object is neither copied nor transferred.
 * @param[in] function_name Borrowed null-terminated C function label. The
 * pointer may be `nullptr` and is not retained.
 * @param[in] handle_name Borrowed non-NULL, null-terminated handle label used
 * only to format a diagnostic.
 * @return `APO_STATUS_OK` when both `handle` and `handle->object` are non-NULL,
 * or `APO_STATUS_INVALID_ARGUMENT` otherwise.
 *
 * @note This helper encodes the C ABI invariant that every wrapper handle
 * stores its retained native object in a member named `object`.
 */
template <typename T>
apo_status require_handle_object(const T *handle, const char *function_name,
                                 const char *handle_name) {
  apo_status status = require_pointer<T>(handle, function_name, handle_name);

  if (status != APO_STATUS_OK)
    return status;

  if (handle->object != nullptr)
    return APO_STATUS_OK;

  return invalid_argument(function_name,
                          std::string(handle_name) + " object is NULL");
}

/**
 * @brief Validates the capacity of a caller-owned output buffer.
 *
 * @tparam T Buffer element type.
 * @param[in] buffer Borrowed contiguous candidate output buffer. The pointer
 * may be NULL only when `required_len` is zero. This helper validates but does
 * not write it.
 * @param[in] provided_len Number of `T` elements available in `buffer`.
 * @param[in] required_len Minimum number of `T` elements required by the
 * caller.
 * @param[in] function_name Borrowed null-terminated C function label. The
 * pointer may be `nullptr` and is not retained.
 * @param[in] buffer_name Borrowed non-NULL, null-terminated buffer label used
 * only to format a diagnostic.
 * @return `APO_STATUS_OK` when no elements are required or the non-NULL buffer
 * has sufficient capacity; `APO_STATUS_INVALID_ARGUMENT` otherwise.
 *
 * @post The buffer and its contents are unchanged.
 */
template <typename T>
apo_status
require_output_buffer(const T *buffer, const std::size_t provided_len,
                      const std::size_t required_len, const char *function_name,
                      const char *buffer_name) {
  if (required_len == 0)
    return APO_STATUS_OK;

  if (buffer == nullptr) {
    return invalid_argument(function_name,
                            std::string(buffer_name) + " is NULL");
  }

  if (provided_len < required_len) {
    return invalid_argument(function_name,
                            std::string(buffer_name) + " is too small");
  }

  return APO_STATUS_OK;
}

/**
 * @brief Validates a nonempty flat-array length against an element stride.
 *
 * @param[in] length Number of scalar elements in the flat array.
 * @param[in] stride Number of scalar elements in one logical record. The value
 * must be greater than zero.
 * @param[in] argument_name Borrowed non-NULL, null-terminated array label used
 * only to format a diagnostic.
 * @param[in] function_name Borrowed null-terminated C function label. The
 * pointer may be `nullptr` and is not retained.
 * @return `APO_STATUS_OK` when `length` is nonzero, `stride` is nonzero, and
 * `length` is an exact multiple of `stride`; `APO_STATUS_INVALID_ARGUMENT`
 * otherwise.
 */
inline apo_status require_flat_array_length(const size_t length,
                                            const size_t stride,
                                            const char *argument_name,
                                            const char *function_name) {
  if (length == 0) {
    return apocharmm_c::invalid_argument(
        function_name, std::string(argument_name) + " must not be empty");
  }

  if ((stride == 0) || (length % stride != 0)) {
    return apocharmm_c::invalid_argument(
        function_name, std::string(argument_name) +
                           " length must be a positive multiple of " +
                           std::to_string(stride));
  }

  return APO_STATUS_OK;
}

} // namespace apocharmm_c

/**
 * @brief Returns immediately from a C ABI body when an expression fails.
 *
 * The macro evaluates `expression` exactly once, stores its @ref apo_status,
 * and returns that status from the enclosing function or lambda when it is not
 * `APO_STATUS_OK`.
 *
 * @param[in] expression Expression yielding an @ref apo_status.
 * @post An `APO_STATUS_OK` result continues execution immediately after the
 * macro.
 */
#define APOCHARMM_C_RETURN_IF_ERROR(expression)                                \
  do {                                                                         \
    const apo_status apocharmm_c_status__ = (expression);                      \
    if (apocharmm_c_status__ != APO_STATUS_OK)                                 \
      return apocharmm_c_status__;                                             \
  } while (false);

#endif

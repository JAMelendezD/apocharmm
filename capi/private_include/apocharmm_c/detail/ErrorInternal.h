// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_ERROR_INTERNAL_H__
#define __APOCHARMM_C_ERROR_INTERNAL_H__

#include "apocharmm_c/Export.h"
#include "apocharmm_c/Status.h"

#include "ApoCharmmError.h"

#include <exception>
#include <stdexcept>
#include <string_view>

#if defined(APOCHARMM_C_TESTING)
#define APOCHARMM_C_INTERNAL_TEST_API APOCHARMM_C_API
#else
#define APOCHARMM_C_INTERNAL_TEST_API
#endif

namespace apocharmm_c {

/**
 * @brief Clears the current thread's stored C ABI diagnostic.
 *
 * @post @ref apo_last_error returns a non-NULL pointer to an empty string on
 * the calling thread.
 * @note The operation does not affect another thread's diagnostic.
 */
APOCHARMM_C_INTERNAL_TEST_API void clear_last_error(void) noexcept;

/**
 * @brief Stores a C-string diagnostic and returns a status unchanged.
 *
 * @param[in] status Status to return without modification.
 * @param[in] message Borrowed null-terminated diagnostic. The pointer may be
 * `nullptr`; non-NULL text is copied before return and is not retained.
 * @return The exact value of `status`.
 *
 * @post The current thread has a nonempty diagnostic. A NULL or empty message
 * is stored as `"Unknown error"`; a storage failure uses a static fallback.
 */
apo_status set_last_error(apo_status status, const char *message) noexcept;

/**
 * @brief Stores a length-delimited diagnostic and returns a status unchanged.
 *
 * @param[in] status Status to return without modification.
 * @param[in] message Borrowed diagnostic view. Its data need not be
 * null-terminated and is copied before return.
 * @return The exact value of `status`.
 *
 * @post The current thread has a nonempty diagnostic. An empty message is
 * stored as `"Unknown error"`; a storage failure uses a static fallback.
 */
apo_status set_last_error(apo_status status,
                          const std::string_view message) noexcept;

/**
 * @brief Stores a function-prefixed diagnostic and returns a status unchanged.
 *
 * @param[in] status Status to return without modification.
 * @param[in] function_name Borrowed null-terminated C function label. The
 * pointer may be `nullptr` and is not retained.
 * @param[in] message Borrowed diagnostic view copied before return. An empty
 * view is rendered as `"Unknown error"`.
 * @return The exact value of `status`.
 *
 * @post The current thread has a nonempty diagnostic beginning with the chosen
 * function label and `": "`.
 * @note If dynamic formatting fails, the no-throw fallback buffer is used and
 * the diagnostic may be truncated while remaining null-terminated.
 */
APOCHARMM_C_INTERNAL_TEST_API
apo_status set_last_error(apo_status status, const char *function_name,
                          const std::string_view message) noexcept;

/**
 * @brief Ensures that a failed status has a function-prefixed diagnostic.
 *
 * An existing diagnostic already beginning with `function_name` followed by a
 * colon is preserved. An empty diagnostic receives the standard missing-text
 * fallback; any other diagnostic is prefixed once.
 *
 * @param[in] status Failed status to return without modification.
 * @param[in] function_name Borrowed null-terminated C function label. The
 * pointer may be `nullptr` and is not retained.
 * @return The exact value of `status`.
 *
 * @pre `status` is not `APO_STATUS_OK`.
 * @post The current thread has a nonempty diagnostic associated with the
 * selected function label.
 */
APOCHARMM_C_INTERNAL_TEST_API
apo_status ensure_last_error(apo_status status,
                             const char *function_name) noexcept;

/**
 * @brief Maps a native error category to its public C ABI status.
 *
 * @param[in] code Native category to map.
 * @return The corresponding declared @ref apo_status. Any value outside the
 * declared @ref ApoCharmmErrorCode enumerators maps to
 * `APO_STATUS_RUNTIME_ERROR`.
 */
inline apo_status
status_from_error_code(const ApoCharmmErrorCode code) noexcept {
  switch (code) {
  case ApoCharmmErrorCode::InvalidArgument:
    return APO_STATUS_INVALID_ARGUMENT;
  case ApoCharmmErrorCode::Runtime:
    return APO_STATUS_RUNTIME_ERROR;
  case ApoCharmmErrorCode::Cuda:
    return APO_STATUS_CUDA_ERROR;
  case ApoCharmmErrorCode::NotInitialized:
    return APO_STATUS_NOT_INITIALIZED;
  case ApoCharmmErrorCode::NotImplemented:
    return APO_STATUS_NOT_IMPLEMENTED;
  }
  return APO_STATUS_RUNTIME_ERROR;
}

/**
 * @brief Executes a status-returning C ABI body behind a nonthrowing boundary.
 *
 * The guard clears the calling thread's previous diagnostic before invoking
 * `function` exactly once. An `APO_STATUS_OK` result is returned unchanged. Any
 * other returned status is preserved and receives a diagnostic when necessary.
 *
 * An @ref ApoCharmmError maps through `status_from_error_code()`. A legacy
 * `std::invalid_argument` maps to `APO_STATUS_INVALID_ARGUMENT`. Every other
 * `std::exception` and every nonstandard exception map to
 * `APO_STATUS_RUNTIME_ERROR`. No exception escapes the boundary.
 *
 * @tparam Function Callable type invocable with no arguments and returning an
 * @ref apo_status.
 * @param[in] function Callable borrowed for this invocation and not retained.
 * @param[in] function_name Borrowed null-terminated public C function label.
 * The pointer may be `nullptr` and is not retained.
 * @return The callable's returned status, or the exact translated status for a
 * caught exception.
 *
 * @post The diagnostic is cleared before `function` is invoked. If `function`
 * returns `APO_STATUS_OK` without recording another diagnostic, it remains
 * empty. Failure leaves a nonempty diagnostic prefixed with the selected C
 * function label.
 */
template <typename Function>
apo_status guard(Function &&function, const char *function_name) noexcept {
  clear_last_error();

  try {
    const apo_status status = function();

    if (status == APO_STATUS_OK)
      return status;

    return ensure_last_error(status, function_name);
  } catch (const ApoCharmmError &e) {
    return set_last_error(status_from_error_code(e.getCode()), function_name,
                          e.what());
  } catch (const std::invalid_argument &e) {
    return set_last_error(APO_STATUS_INVALID_ARGUMENT, function_name, e.what());
  } catch (const std::exception &e) {
    return set_last_error(APO_STATUS_RUNTIME_ERROR, function_name, e.what());
  } catch (...) {
    return set_last_error(APO_STATUS_RUNTIME_ERROR, function_name,
                          "Unknown C++ exception");
  }
}

/**
 * @brief Executes a void destroy body behind a nonthrowing C ABI boundary.
 *
 * Unlike @ref guard, this function does not clear the calling thread's
 * diagnostic before invocation. A normally returning destroy body that does
 * not modify the diagnostic therefore preserves a preceding operation's
 * diagnostic. A caught exception replaces the diagnostic using the same
 * exception-to-status category mapping as
 * @ref guard, although the mapped status cannot be returned by the void destroy
 * API.
 *
 * @tparam Function Callable type invocable with no arguments and returning
 * `void`.
 * @param[in] function Destroy callable borrowed for this invocation and not
 * retained.
 * @param[in] function_name Borrowed null-terminated public C destroy-function
 * label. The pointer may be `nullptr` and is not retained.
 *
 * @post No exception escapes. If `function` returns normally without changing
 * the diagnostic, the previous diagnostic is preserved. A caught exception
 * leaves a nonempty diagnostic describing the destruction failure.
 */
template <typename Function>
void guard_destroy(Function &&function, const char *function_name) noexcept {
  try {
    function();
  } catch (const ApoCharmmError &e) {
    (void)set_last_error(status_from_error_code(e.getCode()), function_name,
                         e.what());
  } catch (const std::invalid_argument &e) {
    (void)set_last_error(APO_STATUS_INVALID_ARGUMENT, function_name, e.what());
  } catch (const std::exception &e) {
    (void)set_last_error(APO_STATUS_RUNTIME_ERROR, function_name, e.what());
  } catch (...) {
    (void)set_last_error(APO_STATUS_RUNTIME_ERROR, function_name,
                         "Unknown C++ exception");
  }

  return;
}

} // namespace apocharmm_c

#undef APOCHARMM_C_INTERNAL_TEST_API

#endif

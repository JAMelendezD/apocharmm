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

APOCHARMM_C_INTERNAL_TEST_API void clear_last_error(void) noexcept;

apo_status set_last_error(apo_status status, const char *message) noexcept;

apo_status set_last_error(apo_status status,
                          const std::string_view message) noexcept;

APOCHARMM_C_INTERNAL_TEST_API
apo_status set_last_error(apo_status status, const char *function_name,
                          const std::string_view message) noexcept;

APOCHARMM_C_INTERNAL_TEST_API
apo_status ensure_last_error(apo_status status,
                             const char *function_name) noexcept;

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

template <typename Function>
void guard_destroy(Function &&function, const char *function_name) noexcept {
  // JEG260810: A destroy function can be called implicitly by a
  // language-runtime finalizer when the diagnostic from a preceding
  // status-returning call is being consumed. Successful destruction must not
  // clear that diagnostic.
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

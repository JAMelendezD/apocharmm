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

#include "apocharmm_c/Status.h"

#include "ApoCharmmError.h"

#include <exception>
#include <stdexcept>
#include <string_view>

namespace apocharmm_c {

void clear_last_error(void) noexcept;

apo_status set_last_error(apo_status status, const char *message) noexcept;

apo_status set_last_error(apo_status status,
                          const std::string_view message) noexcept;

apo_status set_last_error(apo_status status, const char *function_name,
                          const std::string_view message) noexcept;

apo_status ensure_last_error(apo_status status,
                             const char *function_name) noexcept;

template <typename Function>
apo_status guard(Function &&function, const char *function_name) noexcept {
  clear_last_error();

  try {
    const apo_status status = function();

    if (status == APO_STATUS_OK)
      return status;

    return ensure_last_error(status, function_name);
  } catch (const ApoCharmmError &e) {
    apo_status status = APO_STATUS_RUNTIME_ERROR;

    switch (e.getCode()) {
    case ApoCharmmErrorCode::InvalidArgument:
      status = APO_STATUS_INVALID_ARGUMENT;
      break;
    case ApoCharmmErrorCode::Runtime:
      status = APO_STATUS_RUNTIME_ERROR;
      break;
    case ApoCharmmErrorCode::Cuda:
      status = APO_STATUS_CUDA_ERROR;
      break;
    case ApoCharmmErrorCode::NotInitialized:
      status = APO_STATUS_NOT_INITIALIZED;
      break;
    case ApoCharmmErrorCode::NotImplemented:
      status = APO_STATUS_NOT_IMPLEMENTED;
      break;
    }

    return set_last_error(status, function_name, e.what());
  } catch (const std::invalid_argument &e) {
    return set_last_error(APO_STATUS_INVALID_ARGUMENT, function_name, e.what());
  } catch (const std::exception &e) {
    return set_last_error(APO_STATUS_RUNTIME_ERROR, function_name, e.what());
  } catch (...) {
    return set_last_error(APO_STATUS_RUNTIME_ERROR, function_name,
                          "Unknown C++ exception");
  }
}

} // namespace apocharmm_c

#endif

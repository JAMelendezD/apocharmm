// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "apocharmm_c/Error.h"
#include "apocharmm_c/detail/ErrorInternal.h"

#include <cstddef>
#include <string>
#include <string_view>
#include <utility>

// Anonymous namespace to prevent usage outside of this file
namespace {

constexpr char FALLBACK_DIAGNOSTIC[] =
    "Failed to store or format the intended error diagnostic";
constexpr std::size_t FALLBACK_DIAGNOSTIC_CAPACITY = 1024;
constexpr std::string_view UNKNOWN_FUNCTION_NAME = "unknown C function";
constexpr std::string_view UNKNOWN_ERROR_MESSAGE = "Unknown error";

/**
 * @brief Owns one thread's current C ABI diagnostic storage.
 *
 * After an error helper returns, `message` aliases the empty string literal,
 * `storage.c_str()`, `fallback`, or the static fallback diagnostic. The fixed
 * buffer provides a null-terminated, allocation-free path when prefix
 * formatting fails.
 */
struct LastErrorState {
  /** @brief Owns a normally formatted diagnostic. */
  std::string storage;
  /** @brief Owns an emergency function-prefixed diagnostic. */
  char fallback[FALLBACK_DIAGNOSTIC_CAPACITY] = {};
  /** @brief Points at the diagnostic exposed by @ref apo_last_error. */
  const char *message = "";
};

/** @brief Stores the independent diagnostic state for the current thread. */
thread_local LastErrorState g_last_error;

/**
 * @brief Returns the function label used to prefix a C ABI diagnostic.
 *
 * @param[in] function_name Borrowed null-terminated label. The pointer may be
 * `nullptr` and is not retained.
 * @return A borrowed view into `function_name`, or static fallback storage when
 * the pointer is NULL or its first character is a newline.
 */
std::string_view get_function_name(const char *function_name) noexcept {
  if ((function_name == nullptr) || (function_name[0] == '\n'))
    return UNKNOWN_FUNCTION_NAME;

  return function_name;
}

} // namespace

extern "C" const char *apo_last_error(void) { return g_last_error.message; }

void apocharmm_c::clear_last_error(void) noexcept {
  g_last_error.storage.clear();
  g_last_error.fallback[0] = '\0';
  g_last_error.message = "";
  return;
}

apo_status apocharmm_c::set_last_error(apo_status status,
                                       const char *message) noexcept {
  return set_last_error(status, (message != nullptr) ? std::string_view(message)
                                                     : std::string_view());
}

apo_status
apocharmm_c::set_last_error(apo_status status,
                            const std::string_view message) noexcept {
  const std::string_view safe_message =
      message.empty() ? UNKNOWN_ERROR_MESSAGE : message;

  try {
    std::string diagnostic(safe_message.data(), safe_message.size());
    g_last_error.storage = std::move(diagnostic);
    g_last_error.message = g_last_error.storage.c_str();
  } catch (...) {
    g_last_error.message = FALLBACK_DIAGNOSTIC;
  }

  return status;
}

apo_status
apocharmm_c::set_last_error(apo_status status, const char *function_name,
                            const std::string_view message) noexcept {
  const std::string_view safe_function_name = get_function_name(function_name);
  const std::string_view safe_message =
      message.empty() ? UNKNOWN_ERROR_MESSAGE : message;

  try {
    std::string diagnostic;
    diagnostic.reserve(safe_function_name.size() + 2 + safe_message.size());
    diagnostic.append(safe_function_name.data(), safe_function_name.size());
    diagnostic.append(": ");
    diagnostic.append(safe_message.data(), safe_message.size());

    g_last_error.storage = std::move(diagnostic);
    g_last_error.message = g_last_error.storage.c_str();
  } catch (...) {
    std::size_t offset = 0;

    const auto append = [&](const std::string_view text) noexcept -> void {
      for (const char c : text) {
        if (offset + 1 >= FALLBACK_DIAGNOSTIC_CAPACITY)
          break;

        g_last_error.fallback[offset++] = c;
      }

      return;
    };

    append(safe_function_name);
    append(": ");
    append(FALLBACK_DIAGNOSTIC);
    g_last_error.fallback[offset] = '\0';
    g_last_error.message = g_last_error.fallback;
  }

  return status;
}

apo_status apocharmm_c::ensure_last_error(apo_status status,
                                          const char *function_name) noexcept {
  const std::string_view safe_function_name = get_function_name(function_name);
  const std::string_view diagnostic = g_last_error.message;

  if ((diagnostic.size() > safe_function_name.size()) &&
      (diagnostic.compare(0, safe_function_name.size(), safe_function_name) ==
       0) &&
      (diagnostic[safe_function_name.size()] == ':')) {
    return status;
  }

  if (diagnostic.empty()) {
    return set_last_error(status, function_name,
                          "C ABI call failed without a diagnostic");
  }

  return set_last_error(status, function_name, diagnostic);
}

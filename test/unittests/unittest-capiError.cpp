// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "ApoCharmmError.h"
#include "apocharmm_c/CharmmParameters.h"
#include "apocharmm_c/Error.h"
#include "apocharmm_c/Status.h"
#include "apocharmm_c/detail/ErrorInternal.h"
#include "catch.hpp"

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>

namespace {

struct ErrorMapping {
  ApoCharmmErrorCode code;
  apo_status status;
  const char *message;
};

constexpr std::array<ErrorMapping, 5> ERROR_MAPPINGS = {
    {{ApoCharmmErrorCode::InvalidArgument, APO_STATUS_INVALID_ARGUMENT,
      "mapped invalid argument"},
     {ApoCharmmErrorCode::Runtime, APO_STATUS_RUNTIME_ERROR, "mapped runtime"},
     {ApoCharmmErrorCode::Cuda, APO_STATUS_CUDA_ERROR, "mapped CUDA"},
     {ApoCharmmErrorCode::NotInitialized, APO_STATUS_NOT_INITIALIZED,
      "mapped not initialized"},
     {ApoCharmmErrorCode::NotImplemented, APO_STATUS_NOT_IMPLEMENTED,
      "mapped not implemented"}}};

constexpr std::array<apo_status, 5> FAILED_STATUSES = {
    {APO_STATUS_INVALID_ARGUMENT, APO_STATUS_RUNTIME_ERROR,
     APO_STATUS_CUDA_ERROR, APO_STATUS_NOT_INITIALIZED,
     APO_STATUS_NOT_IMPLEMENTED}};

std::size_t CountOccurrences(const std::string_view text,
                             const std::string_view needle) noexcept {
  if (needle.empty())
    return 0;

  std::size_t count = 0;
  std::size_t position = 0;

  while ((position = text.find(needle, position)) != std::string_view::npos) {
    count++;
    position += needle.size();
  }

  return count;
}

} // namespace

TEST_CASE("CapiErrorMapsEveryApoCharmmErrorCode") {
  constexpr char FUNCTION_NAME[] = "apo_test_error_mapping";

  for (const ErrorMapping &mapping : ERROR_MAPPINGS) {
    INFO("error code: " << GetApoCharmmErrorCodeName(mapping.code));

    apo_status status = APO_STATUS_OK;
    CHECK_NOTHROW((status = apocharmm_c::guard(
                       [&mapping](void) -> apo_status {
                         throw ApoCharmmError(mapping.code, mapping.message,
                                              "unittest-capiEror.cpp",
                                              "MappedThrowSite", 41);
                       },
                       FUNCTION_NAME)));

    CHECK(status == mapping.status);

    const std::string diagnostic(apo_last_error());
    CHECK(diagnostic.empty() == false);
    CHECK(CountOccurrences(diagnostic, FUNCTION_NAME) == 1);
    CHECK(diagnostic.find(mapping.message) != std::string::npos);
  }
}

TEST_CASE("CapiErrorMapsLegacyAndUnknownExceptions") {
  SECTION("InvalidArgument") {
    constexpr char FUNCTION_NAME[] = "apo_test_legacy_invalid_argument";
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW((status = apocharmm_c::guard(
                       [](void) -> apo_status {
                         throw std::invalid_argument("legacy invalid argument");
                       },
                       FUNCTION_NAME)));

    CHECK(status == APO_STATUS_INVALID_ARGUMENT);

    const std::string diagnostic(apo_last_error());
    CHECK(diagnostic.empty() == false);
    CHECK(diagnostic ==
          "apo_test_legacy_invalid_argument: legacy invalid argument");
    CHECK(CountOccurrences(diagnostic, FUNCTION_NAME) == 1);
  }

  SECTION("RuntimeError") {
    constexpr char FUNCTION_NAME[] = "apo_test_runtime_error";
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW((status = apocharmm_c::guard(
                       [](void) -> apo_status {
                         throw std::runtime_error("generic runtime error");
                       },
                       FUNCTION_NAME)));

    CHECK(status == APO_STATUS_RUNTIME_ERROR);

    const std::string diagnostic(apo_last_error());
    CHECK(diagnostic.empty() == false);
    CHECK(diagnostic == "apo_test_runtime_error: generic runtime error");
    CHECK(CountOccurrences(diagnostic, FUNCTION_NAME) == 1);
  }

  SECTION("UnknownException") {
    constexpr char FUNCTION_NAME[] = "apo_test_unknown_exception";
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW((status = apocharmm_c::guard(
                       [](void) -> apo_status { throw 7; }, FUNCTION_NAME)));

    CHECK(status == APO_STATUS_RUNTIME_ERROR);

    const std::string diagnostic(apo_last_error());
    CHECK(diagnostic.empty() == false);
    CHECK(diagnostic == "apo_test_unknown_exception: Unknown C++ exception");
    CHECK(CountOccurrences(diagnostic, FUNCTION_NAME) == 1);
  }
}

TEST_CASE("CapiErrorGuardCompletesFailedDiagnostic") {
  SECTION("EveryFailedStatusGetsADiagnostic") {
    constexpr char FUNCTION_NAME[] = "apo_test_bare_failure";

    for (const apo_status failedStatus : FAILED_STATUSES) {
      INFO("apo_status: " << static_cast<int>(failedStatus));

      apo_status status = APO_STATUS_OK;
      CHECK_NOTHROW(
          (status = apocharmm_c::guard(
               [failedStatus](void) -> apo_status { return failedStatus; },
               FUNCTION_NAME)));

      CHECK(status == failedStatus);

      const std::string diagnostic(apo_last_error());
      CHECK(diagnostic.empty() == false);
      CHECK(diagnostic ==
            "apo_test_bare_failure: C ABI call failed without a diagnostic");
      CHECK(CountOccurrences(diagnostic, FUNCTION_NAME) == 1);
    }
  }

  SECTION("ExistingDiagnosticIsNotPrefixedTwice") {
    constexpr char FUNCTION_NAME[] = "apo_test_existing_diagnostic";
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW((status = apocharmm_c::guard(
                       [&](void) -> apo_status {
                         return apocharmm_c::set_last_error(
                             APO_STATUS_CUDA_ERROR, FUNCTION_NAME,
                             "existing diagnostic");
                       },
                       FUNCTION_NAME)));

    CHECK(status == APO_STATUS_CUDA_ERROR);

    const std::string diagnostic(apo_last_error());
    CHECK(diagnostic.empty() == false);
    CHECK(diagnostic == "apo_test_existing_diagnostic: existing diagnostic");
    CHECK(CountOccurrences(diagnostic, FUNCTION_NAME) == 1);
  }
}

TEST_CASE("CapiErrorSuccessfulGuardClearsStaleDiagnostic") {
  const apo_status failedStatus = apocharmm_c::guard(
      [](void) -> apo_status { throw std::runtime_error("stale diagnostic"); },
      "apo_test_stale_error");

  REQUIRE(failedStatus == APO_STATUS_RUNTIME_ERROR);
  REQUIRE(std::string(apo_last_error()).empty() == false);

  const apo_status successStatus =
      apocharmm_c::guard([](void) -> apo_status { return APO_STATUS_OK; },
                         "apo_test_successful_call");

  CHECK(successStatus == APO_STATUS_OK);
  CHECK(std::string(apo_last_error()).empty() == true);
}

TEST_CASE("CapiLastErrorPointerFollowsDocumentedLifetime") {
  constexpr char FIRST_FUNCTION_NAME[] = "apo_charmm_parameters_create";
  constexpr char SECOND_FUNCTION_NAME[] =
      "apo_charmm_parameters_create_from_files";

  apo_status status = APO_STATUS_OK;
  CHECK_NOTHROW((status = apo_charmm_parameters_create(nullptr, "unused.prm")));
  REQUIRE(status == APO_STATUS_INVALID_ARGUMENT);

  const char *const firstPointer = apo_last_error();
  REQUIRE(firstPointer != nullptr);

  const std::string firstDiagnostic(firstPointer);
  REQUIRE(firstDiagnostic.empty() == false);
  CHECK(firstDiagnostic == "apo_charmm_parameters_create: out pointer is NULL");
  CHECK(CountOccurrences(firstDiagnostic, FIRST_FUNCTION_NAME) == 1);

  const char *const repeatedPointer = apo_last_error();
  REQUIRE(repeatedPointer != nullptr);
  CHECK(std::string(repeatedPointer) == firstDiagnostic);
  CHECK(std::string(firstPointer) == firstDiagnostic);

  // Do not dereference firstPointer after this call. Its documented lifetime
  // ends when the thread-local diagnostic is cleared or replaced.
  CHECK_NOTHROW(
      (status = apo_charmm_parameters_create_from_files(nullptr, nullptr, 0)));
  REQUIRE(status == APO_STATUS_INVALID_ARGUMENT);

  const char *const replacementPointer = apo_last_error();
  REQUIRE(replacementPointer != nullptr);

  const std::string replacementDiagnostic(replacementPointer);
  CHECK(replacementDiagnostic ==
        "apo_charmm_parameters_create_from_files: out pointer is NULL");
  CHECK(CountOccurrences(replacementDiagnostic, SECOND_FUNCTION_NAME) == 1);
}

TEST_CASE("CapiLastErrorIsThreadLocal") {
  REQUIRE(apocharmm_c::guard([](void) -> apo_status { return APO_STATUS_OK; },
                             "apo_test_main_thread_clear") == APO_STATUS_OK);

  REQUIRE(std::string(apo_last_error()).empty() == true);

  std::atomic<int> readyThreads{0};
  std::atomic<bool> readDiagnostics{false};

  apo_status firstStatus = APO_STATUS_OK;
  apo_status secondStatus = APO_STATUS_OK;
  std::string firstDiagnostic;
  std::string secondDiagnostic;

  const auto runThread = [&](const char *functionName, const char *message,
                             apo_status &status,
                             std::string &diagnostic) -> void {
    status = apocharmm_c::guard(
        [message](void) -> apo_status { throw std::runtime_error(message); },
        functionName);

    readyThreads.fetch_add(1, std::memory_order_release);
    while (!readDiagnostics.load(std::memory_order_acquire))
      std::this_thread::yield();

    diagnostic = apo_last_error();
    return;
  };

  std::thread firstThread(runThread, "apo_test_thread_one",
                          "first thread diagnostic", std::ref(firstStatus),
                          std::ref(firstDiagnostic));
  std::thread secondThread(runThread, "apo_test_thread_two",
                           "second thread diagnostic", std::ref(secondStatus),
                           std::ref(secondDiagnostic));

  while (readyThreads.load(std::memory_order_acquire) != 2)
    std::this_thread::yield();

  readDiagnostics.store(true, std::memory_order_release);

  firstThread.join();
  secondThread.join();

  CHECK(firstStatus == APO_STATUS_RUNTIME_ERROR);
  CHECK(secondStatus == APO_STATUS_RUNTIME_ERROR);
  CHECK(firstDiagnostic == "apo_test_thread_one: first thread diagnostic");
  CHECK(secondDiagnostic == "apo_test_thread_two: second thread diagnostic");
  CHECK(firstDiagnostic != secondDiagnostic);
  CHECK(CountOccurrences(firstDiagnostic, "apo_test_thread_one") == 1);
  CHECK(CountOccurrences(secondDiagnostic, "apo_test_thread_two") == 1);
  CHECK(std::string(apo_last_error()).empty() == true);
}

TEST_CASE("CapiPublicFunctionRejectsInvalidPointerAndClearsOnRecovery") {
  constexpr char PARAMETER_FILE[] = "tmp_capi_error_boundary.prm";
  constexpr char FUNCTION_NAME[] = "apo_charmm_parameters_create";

  {
    std::ofstream parameterFile(PARAMETER_FILE);
    REQUIRE(parameterFile.good());

    parameterFile << "* generated C ABI error boundary unit test\n"
                  << "*\n"
                  << "END\n";
    parameterFile.close();

    REQUIRE(parameterFile.good());
  }

  apo_status failedStatus = APO_STATUS_OK;
  CHECK_NOTHROW(
      (failedStatus = apo_charmm_parameters_create(nullptr, PARAMETER_FILE)));

  CHECK(failedStatus == APO_STATUS_INVALID_ARGUMENT);

  const std::string failedDiagnostic(apo_last_error());
  CHECK(failedDiagnostic.empty() == false);
  CHECK(failedDiagnostic ==
        "apo_charmm_parameters_create: out pointer is NULL");
  CHECK(CountOccurrences(failedDiagnostic, FUNCTION_NAME) == 1);

  apo_charmm_parameters *parameters = nullptr;
  apo_status successStatus = APO_STATUS_RUNTIME_ERROR;
  CHECK_NOTHROW((successStatus = apo_charmm_parameters_create(&parameters,
                                                              PARAMETER_FILE)));

  CHECK(successStatus == APO_STATUS_OK);
  CHECK(parameters != nullptr);
  CHECK(std::string(apo_last_error()).empty() == true);

  apo_charmm_parameters_destroy(parameters);
  CHECK(std::remove(PARAMETER_FILE) == 0);
}

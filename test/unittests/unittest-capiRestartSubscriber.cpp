// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "apo_test_helpers.h"
#include "apocharmm_c/Error.h"
#include "apocharmm_c/RestartSubscriber.h"
#include "apocharmm_c/Status.h"
#include "apocharmm_c/Subscriber.h"
#include "apocharmm_c/detail/RestartSubscriberHandle.h"
#include "catch.hpp"

#include <memory>
#include <string>
#include <string_view>

namespace {

struct RestartSubscriberDeleter {
  void operator()(apo_restart_subscriber *subscriber) const noexcept {
    apo_restart_subscriber_destroy(subscriber);
    return;
  }
};

using RestartSubscriberHandle =
    std::unique_ptr<apo_restart_subscriber, RestartSubscriberDeleter>;

RestartSubscriberHandle CreateRestartSubscriber(const std::string &path) {
  apo_restart_subscriber *rawSubscriber = nullptr;
  REQUIRE(apo_restart_subscriber_create(&rawSubscriber, path.c_str()) ==
          APO_STATUS_OK);
  REQUIRE(rawSubscriber != nullptr);
  return RestartSubscriberHandle(rawSubscriber);
}

RestartSubscriberHandle CreateRestartSubscriber(const std::string &path,
                                                const int reportFrequency) {
  apo_restart_subscriber *rawSubscriber = nullptr;
  REQUIRE(apo_restart_subscriber_create_with_report_frequency(
              &rawSubscriber, path.c_str(), reportFrequency) == APO_STATUS_OK);
  REQUIRE(rawSubscriber != nullptr);
  return RestartSubscriberHandle(rawSubscriber);
}

} // namespace

TEST_CASE("CapiRestartSubscriberConstructionAndBaseConversion") {
  const std::string defaultPath = "tmpCapiRestartSubscriberDefault.rst";
  const std::string frequencyPath = "tmpCapiRestartSubscriberFrequency.rst";
  apo_test::RemoveIfExists(defaultPath);
  apo_test::RemoveIfExists(frequencyPath);

  SECTION("DefaultReportFrequency") {
    RestartSubscriberHandle subscriber = CreateRestartSubscriber(defaultPath);

    apo_subscriber *baseSubscriber = nullptr;
    CHECK(apo_restart_subscriber_as_subscriber(
              &baseSubscriber, subscriber.get()) == APO_STATUS_OK);
    REQUIRE(baseSubscriber != nullptr);
    CHECK(baseSubscriber == &subscriber->base);
    CHECK(baseSubscriber->object == subscriber->object);

    int reportFrequency = 0;
    CHECK(apo_subscriber_get_report_frequency(&reportFrequency,
                                              baseSubscriber) == APO_STATUS_OK);
    CHECK(reportFrequency == 1000);
    CHECK(std::string(apo_last_error()).empty() == true);
  }

  SECTION("ExplicitReportFrequency") {
    constexpr int REPORT_FREQUENCY = 17;
    RestartSubscriberHandle subscriber =
        CreateRestartSubscriber(frequencyPath, REPORT_FREQUENCY);

    apo_subscriber *baseSubscriber = nullptr;
    CHECK(apo_restart_subscriber_as_subscriber(
              &baseSubscriber, subscriber.get()) == APO_STATUS_OK);
    REQUIRE(baseSubscriber != nullptr);

    int reportFrequency = 0;
    CHECK(apo_subscriber_get_report_frequency(&reportFrequency,
                                              baseSubscriber) == APO_STATUS_OK);
    CHECK(reportFrequency == REPORT_FREQUENCY);
    CHECK(std::string(apo_last_error()).empty() == true);
  }

  apo_test::RemoveIfExists(defaultPath);
  apo_test::RemoveIfExists(frequencyPath);
}

TEST_CASE("CapiRestartSubscriberValidatesCreationArguments") {
  constexpr const char *PATH = "tmpCapiRestartSubscriberValidation.rst";
  apo_restart_subscriber sentinel;

  apo_test::CheckStatusAndDiagnostic(
      apo_restart_subscriber_create(nullptr, PATH), APO_STATUS_INVALID_ARGUMENT,
      "apo_restart_subscriber_create: out pointer is NULL");

  apo_restart_subscriber *subscriber = &sentinel;
  apo_test::CheckStatusAndDiagnostic(
      apo_restart_subscriber_create(&subscriber, nullptr),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_restart_subscriber_create: path is NULL");
  CHECK(subscriber == nullptr);

  apo_status status = APO_STATUS_OK;

  subscriber = &sentinel;
  CHECK_NOTHROW((status = apo_restart_subscriber_create(&subscriber, "")));
  apo_test::CheckNativeError(status, APO_STATUS_INVALID_ARGUMENT,
                             "InvalidArgument", "apo_restart_subscriber_create",
                             "Output file name must not be empty",
                             "src/Subscriber.cu", "setFilePath");
  CHECK(subscriber == nullptr);

  apo_test::CheckStatusAndDiagnostic(
      apo_restart_subscriber_create_with_report_frequency(nullptr, PATH, 1),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_restart_subscriber_create_with_report_frequency: out pointer is "
      "NULL");

  subscriber = &sentinel;
  apo_test::CheckStatusAndDiagnostic(
      apo_restart_subscriber_create_with_report_frequency(&subscriber, nullptr,
                                                          1),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_restart_subscriber_create_with_report_frequency: path is NULL");
  CHECK(subscriber == nullptr);

  subscriber = &sentinel;
  CHECK_NOTHROW((status = apo_restart_subscriber_create_with_report_frequency(
                     &subscriber, "", 1)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_restart_subscriber_create_with_report_frequency",
      "Output file name must not be empty", "src/Subscriber.cu", "setFilePath");
  CHECK(subscriber == nullptr);

  subscriber = &sentinel;
  CHECK_NOTHROW((status = apo_restart_subscriber_create_with_report_frequency(
                     &subscriber, PATH, 0)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_restart_subscriber_create_with_report_frequency",
      "Subscriber report frequency must be positive; observed 0",
      "src/Subscriber.cu", "setReportFrequency");
  CHECK(subscriber == nullptr);

  subscriber = &sentinel;
  CHECK_NOTHROW((status = apo_restart_subscriber_create_with_report_frequency(
                     &subscriber, PATH, -1)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_restart_subscriber_create_with_report_frequency",
      "Subscriber report frequency must be positive; observed -1",
      "src/Subscriber.cu", "setReportFrequency");
  CHECK(subscriber == nullptr);
}

TEST_CASE("CapiRestartSubscriberMapsNativeConstructionFailure") {
  constexpr const char *PATH = "missing_capi_restart_subscriber_dir/output.rst";
  apo_restart_subscriber sentinel;
  apo_restart_subscriber *subscriber = &sentinel;
  apo_status status = APO_STATUS_OK;

  CHECK_NOTHROW((status = apo_restart_subscriber_create(&subscriber, PATH)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_restart_subscriber_create",
      "Output directory does not exist: missing_capi_restart_subscriber_dir",
      "src/Subscriber.cu", "checkPath");
  CHECK(subscriber == nullptr);
}

TEST_CASE("CapiRestartSubscriberValidatesBaseConversion") {
  const std::string fileName = "tmpCapiRestartSubscriberConversion.rst";
  apo_test::RemoveIfExists(fileName);

  {
    RestartSubscriberHandle subscriber = CreateRestartSubscriber(fileName);

    apo_test::CheckStatusAndDiagnostic(
        apo_restart_subscriber_as_subscriber(nullptr, subscriber.get()),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_restart_subscriber_as_subscriber: out pointer is NULL");

    apo_subscriber baseSentinel;
    apo_subscriber *baseSubscriber = &baseSentinel;
    apo_test::CheckStatusAndDiagnostic(
        apo_restart_subscriber_as_subscriber(&baseSubscriber, nullptr),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_restart_subscriber_as_subscriber: RestartSubscriber is NULL");
    CHECK(baseSubscriber == nullptr);

    apo_restart_subscriber emptySubscriber;
    baseSubscriber = &baseSentinel;
    apo_test::CheckStatusAndDiagnostic(
        apo_restart_subscriber_as_subscriber(&baseSubscriber, &emptySubscriber),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_restart_subscriber_as_subscriber: RestartSubscriber object is "
        "NULL");
    CHECK(baseSubscriber == nullptr);

    CHECK(apo_restart_subscriber_as_subscriber(
              &baseSubscriber, subscriber.get()) == APO_STATUS_OK);
    CHECK(baseSubscriber == &subscriber->base);
    CHECK(std::string(apo_last_error()).empty() == true);
  }

  apo_test::RemoveIfExists(fileName);
}

TEST_CASE("CapiRestartSubscriberSuccessClearsStaleDiagnostic") {
  const std::string fileName = "tmpCapiRestartSubscriberClearError.rst";
  apo_test::RemoveIfExists(fileName);

  apo_restart_subscriber *rawSubscriber = nullptr;
  REQUIRE(apo_restart_subscriber_create(&rawSubscriber, nullptr) ==
          APO_STATUS_INVALID_ARGUMENT);
  REQUIRE(rawSubscriber == nullptr);
  REQUIRE(std::string(apo_last_error()).empty() == false);

  REQUIRE(apo_restart_subscriber_create(&rawSubscriber, fileName.c_str()) ==
          APO_STATUS_OK);
  RestartSubscriberHandle subscriber(rawSubscriber);
  CHECK(std::string(apo_last_error()).empty() == true);

  subscriber.reset();
  apo_test::RemoveIfExists(fileName);
}

TEST_CASE("CapiRestartSubscriberDestroyIsNoexcept") {
  const std::string fileName = "tmpCapiRestartSubscriberDestroy.rst";
  apo_test::RemoveIfExists(fileName);

  RestartSubscriberHandle subscriber = CreateRestartSubscriber(fileName);
  apo_restart_subscriber *const rawSubscriber = subscriber.release();

  CHECK_NOTHROW(apo_restart_subscriber_destroy(rawSubscriber));
  CHECK(std::string(apo_last_error()).empty() == true);

  CHECK_NOTHROW(apo_restart_subscriber_destroy(nullptr));
  CHECK(std::string(apo_last_error()).empty() == true);

  apo_test::RemoveIfExists(fileName);
}

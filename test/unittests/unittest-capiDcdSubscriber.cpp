// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "apo_test_helpers.h"
#include "apocharmm_c/DcdSubscriber.h"
#include "apocharmm_c/Error.h"
#include "apocharmm_c/Status.h"
#include "apocharmm_c/Subscriber.h"
#include "apocharmm_c/detail/DcdSubscriberHandle.h"
#include "catch.hpp"

#include <memory>
#include <string>

namespace {

struct DcdSubscriberDeleter {
  void operator()(apo_dcd_subscriber *subscriber) const noexcept {
    apo_dcd_subscriber_destroy(subscriber);
    return;
  }
};

using DcdSubscriberHandle =
    std::unique_ptr<apo_dcd_subscriber, DcdSubscriberDeleter>;

DcdSubscriberHandle CreateDcdSubscriber(const std::string &path) {
  apo_dcd_subscriber *rawSubscriber = nullptr;
  REQUIRE(apo_dcd_subscriber_create(&rawSubscriber, path.c_str()) ==
          APO_STATUS_OK);
  REQUIRE(rawSubscriber != nullptr);
  return DcdSubscriberHandle(rawSubscriber);
}

} // namespace

TEST_CASE("CapiDcdSubscriberConstructionAndBaseConversion") {
  const std::string defaultPath = "tmpCapiDcdSubscriberDefault.dcd";
  const std::string frequencyPath = "tmpCapiDcdSubscriberFrequency.dcd";
  apo_test::RemoveIfExists(defaultPath);
  apo_test::RemoveIfExists(frequencyPath);

  SECTION("DefaultReportFrequency") {
    DcdSubscriberHandle subscriber = CreateDcdSubscriber(defaultPath);

    apo_subscriber *baseSubscriber = nullptr;
    CHECK(apo_dcd_subscriber_as_subscriber(&baseSubscriber, subscriber.get()) ==
          APO_STATUS_OK);
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
    apo_dcd_subscriber *rawSubscriber = nullptr;
    REQUIRE(apo_dcd_subscriber_create_with_report_frequency(
                &rawSubscriber, frequencyPath.c_str(), REPORT_FREQUENCY) ==
            APO_STATUS_OK);
    REQUIRE(rawSubscriber != nullptr);
    DcdSubscriberHandle subscriber(rawSubscriber);

    apo_subscriber *baseSubscriber = nullptr;
    CHECK(apo_dcd_subscriber_as_subscriber(&baseSubscriber, subscriber.get()) ==
          APO_STATUS_OK);
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

TEST_CASE("CapiDcdSubscriberValidatesCreationArguments") {
  constexpr const char *PATH = "tmpCapiDcdSubscriberValidation.dcd";
  apo_dcd_subscriber sentinel;

  apo_test::CheckStatusAndDiagnostic(
      apo_dcd_subscriber_create(nullptr, PATH), APO_STATUS_INVALID_ARGUMENT,
      "apo_dcd_subscriber_create: out pointer is NULL");

  apo_dcd_subscriber *subscriber = &sentinel;
  apo_test::CheckStatusAndDiagnostic(
      apo_dcd_subscriber_create(&subscriber, nullptr),
      APO_STATUS_INVALID_ARGUMENT, "apo_dcd_subscriber_create: path is NULL");
  CHECK(subscriber == nullptr);

  apo_status status = APO_STATUS_OK;

  subscriber = &sentinel;
  CHECK_NOTHROW((status = apo_dcd_subscriber_create(&subscriber, "")));
  apo_test::CheckNativeError(status, APO_STATUS_INVALID_ARGUMENT,
                             "InvalidArgument", "apo_dcd_subscriber_create",
                             "Output file name must not be empty",
                             "src/Subscriber.cu", "setFileName");
  CHECK(subscriber == nullptr);

  apo_test::CheckStatusAndDiagnostic(
      apo_dcd_subscriber_create_with_report_frequency(nullptr, PATH, 1),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_dcd_subscriber_create_with_report_frequency: out pointer is NULL");

  subscriber = &sentinel;
  apo_test::CheckStatusAndDiagnostic(
      apo_dcd_subscriber_create_with_report_frequency(&subscriber, nullptr, 1),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_dcd_subscriber_create_with_report_frequency: path is NULL");
  CHECK(subscriber == nullptr);

  subscriber = &sentinel;
  CHECK_NOTHROW((status = apo_dcd_subscriber_create_with_report_frequency(
                     &subscriber, "", 1)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_dcd_subscriber_create_with_report_frequency",
      "Output file name must not be empty", "src/Subscriber.cu", "setFileName");
  CHECK(subscriber == nullptr);

  subscriber = &sentinel;
  CHECK_NOTHROW((status = apo_dcd_subscriber_create_with_report_frequency(
                     &subscriber, PATH, 0)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_dcd_subscriber_create_with_report_frequency",
      "Subscriber report frequency must be positive; observed 0",
      "src/Subscriber.cu", "setReportFrequency");
  CHECK(subscriber == nullptr);

  subscriber = &sentinel;
  CHECK_NOTHROW((status = apo_dcd_subscriber_create_with_report_frequency(
                     &subscriber, PATH, -1)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_dcd_subscriber_create_with_report_frequency",
      "Subscriber report frequency must be positive; observed -1",
      "src/Subscriber.cu", "setReportFrequency");
  CHECK(subscriber == nullptr);
}

TEST_CASE("CapiDcdSubscriberMapsNativeConstructionFailure") {
  constexpr const char *PATH = "missing_capi_dcd_subscriber_dir/output.dcd";
  apo_dcd_subscriber sentinel;
  apo_status status = APO_STATUS_OK;

  SECTION("DefaultReportFrequency") {
    apo_dcd_subscriber *subscriber = &sentinel;

    CHECK_NOTHROW((status = apo_dcd_subscriber_create(&subscriber, PATH)));
    apo_test::CheckNativeError(
        status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
        "apo_dcd_subscriber_create",
        "Output directory does not exist: missing_capi_dcd_subscriber_dir",
        "src/Subscriber.cu", "checkPath");
    CHECK(subscriber == nullptr);
  }

  SECTION("ExplicitReportFrequency") {
    apo_dcd_subscriber *subscriber = &sentinel;

    CHECK_NOTHROW((status = apo_dcd_subscriber_create_with_report_frequency(
                       &subscriber, PATH, 1)));
    apo_test::CheckNativeError(
        status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
        "apo_dcd_subscriber_create_with_report_frequency",
        "Output directory does not exist: missing_capi_dcd_subscriber_dir",
        "src/Subscriber.cu", "checkPath");
    CHECK(subscriber == nullptr);
  }
}

TEST_CASE("CapiDcdSubscriberMapsNativeOpenFailure") {
  apo_dcd_subscriber sentinel;
  apo_status status = APO_STATUS_OK;

  SECTION("DefaultReportFrequency") {
    apo_dcd_subscriber *subscriber = &sentinel;

    CHECK_NOTHROW((status = apo_dcd_subscriber_create(&subscriber, ".")));
    apo_test::CheckNativeError(status, APO_STATUS_RUNTIME_ERROR, "Runtime",
                               "apo_dcd_subscriber_create",
                               "Failed to open DCD file for writing: .",
                               "src/DcdSubscriber.cu", "openFile");
    CHECK(subscriber == nullptr);
  }

  SECTION("ExplicitReportFrequency") {
    apo_dcd_subscriber *subscriber = &sentinel;

    CHECK_NOTHROW((status = apo_dcd_subscriber_create_with_report_frequency(
                       &subscriber, ".", 1)));
    apo_test::CheckNativeError(
        status, APO_STATUS_RUNTIME_ERROR, "Runtime",
        "apo_dcd_subscriber_create_with_report_frequency",
        "Failed to open DCD file for writing: .", "src/DcdSubscriber.cu",
        "openFile");
    CHECK(subscriber == nullptr);
  }
}

TEST_CASE("CapiDcdSubscriberValidatesBaseConversion") {
  const std::string fileName = "tmpCapiDcdSubscriberConversion.dcd";
  apo_test::RemoveIfExists(fileName);

  {
    DcdSubscriberHandle subscriber = CreateDcdSubscriber(fileName);

    apo_test::CheckStatusAndDiagnostic(
        apo_dcd_subscriber_as_subscriber(nullptr, subscriber.get()),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_dcd_subscriber_as_subscriber: out pointer is NULL");

    apo_subscriber baseSentinel;
    apo_subscriber *baseSubscriber = &baseSentinel;
    apo_test::CheckStatusAndDiagnostic(
        apo_dcd_subscriber_as_subscriber(&baseSubscriber, nullptr),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_dcd_subscriber_as_subscriber: DcdSubscriber is NULL");
    CHECK(baseSubscriber == nullptr);

    apo_dcd_subscriber emptySubscriber;
    baseSubscriber = &baseSentinel;
    apo_test::CheckStatusAndDiagnostic(
        apo_dcd_subscriber_as_subscriber(&baseSubscriber, &emptySubscriber),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_dcd_subscriber_as_subscriber: DcdSubscriber object is NULL");
    CHECK(baseSubscriber == nullptr);

    CHECK(apo_dcd_subscriber_as_subscriber(&baseSubscriber, subscriber.get()) ==
          APO_STATUS_OK);
    CHECK(baseSubscriber == &subscriber->base);
    CHECK(std::string(apo_last_error()).empty() == true);
  }

  apo_test::RemoveIfExists(fileName);
}

TEST_CASE("CapiDcdSubscriberSuccessClearsStaleDiagnostic") {
  const std::string fileName = "tmpCapiDcdSubscriberClearError.dcd";
  apo_test::RemoveIfExists(fileName);

  apo_dcd_subscriber *rawSubscriber = nullptr;
  REQUIRE(apo_dcd_subscriber_create(&rawSubscriber, nullptr) ==
          APO_STATUS_INVALID_ARGUMENT);
  REQUIRE(rawSubscriber == nullptr);
  REQUIRE(std::string(apo_last_error()).empty() == false);

  REQUIRE(apo_dcd_subscriber_create(&rawSubscriber, fileName.c_str()) ==
          APO_STATUS_OK);
  DcdSubscriberHandle subscriber(rawSubscriber);
  CHECK(std::string(apo_last_error()).empty() == true);

  subscriber.reset();
  apo_test::RemoveIfExists(fileName);
}

TEST_CASE("CapiDcdSubscriberDestroyIsNoexcept") {
  const std::string fileName = "tmpCapiDcdSubscriberDestroy.dcd";
  apo_test::RemoveIfExists(fileName);

  DcdSubscriberHandle subscriber = CreateDcdSubscriber(fileName);
  apo_dcd_subscriber *const rawSubscriber = subscriber.release();

  CHECK_NOTHROW(apo_dcd_subscriber_destroy(rawSubscriber));
  CHECK(std::string(apo_last_error()).empty() == true);

  CHECK_NOTHROW(apo_dcd_subscriber_destroy(nullptr));
  CHECK(std::string(apo_last_error()).empty() == true);

  apo_test::RemoveIfExists(fileName);
}

// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "CudaIntegrator.h"
#include "Subscriber.h"
#include "apo_test_helpers.h"
#include "apocharmm_c/CudaIntegrator.h"
#include "apocharmm_c/DcdSubscriber.h"
#include "apocharmm_c/Error.h"
#include "apocharmm_c/Status.h"
#include "apocharmm_c/Subscriber.h"
#include "apocharmm_c/detail/CudaIntegratorHandle.h"
#include "apocharmm_c/detail/SubscriberHandle.h"
#include "catch.hpp"

#include <fstream>
#include <memory>
#include <string>
#include <string_view>

namespace {

class TestSubscriber final : public Subscriber {
public:
  void update(void) override { return; }
};

apo_subscriber MakeSubscriberHandle(void) {
  apo_subscriber subscriber;
  subscriber.object = std::make_shared<TestSubscriber>();
  return subscriber;
}

apo_cuda_integrator MakeIntegratorHandle(void) {
  apo_cuda_integrator integrator;
  integrator.object = std::make_shared<CudaIntegrator>(0.002);
  return integrator;
}

} // namespace

TEST_CASE("CapiSubscriberReportFrequency") {
  apo_subscriber subscriber = MakeSubscriberHandle();

  int reportFrequency = 0;
  CHECK(apo_subscriber_get_report_frequency(&reportFrequency, &subscriber) ==
        APO_STATUS_OK);
  CHECK(reportFrequency == 1000);

  CHECK(apo_subscriber_set_report_frequency(&subscriber, 17) == APO_STATUS_OK);
  CHECK(apo_subscriber_get_report_frequency(&reportFrequency, &subscriber) ==
        APO_STATUS_OK);
  CHECK(reportFrequency == 17);
  CHECK(std::string(apo_last_error()).empty() == true);
}

TEST_CASE("CapiSubscriberValidatesPointersAndHandles") {
  apo_subscriber subscriber = MakeSubscriberHandle();
  apo_subscriber emptySubscriber;

  apo_test::CheckStatusAndDiagnostic(
      apo_subscriber_set_report_frequency(nullptr, 1),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_subscriber_set_report_frequency: Subscriber is NULL");
  apo_test::CheckStatusAndDiagnostic(
      apo_subscriber_set_report_frequency(&emptySubscriber, 1),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_subscriber_set_report_frequency: Subscriber object is NULL");

  apo_test::CheckStatusAndDiagnostic(
      apo_subscriber_get_report_frequency(nullptr, &subscriber),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_subscriber_get_report_frequency: report_frequency is NULL");

  int reportFrequency = 17;
  apo_test::CheckStatusAndDiagnostic(
      apo_subscriber_get_report_frequency(&reportFrequency, nullptr),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_subscriber_get_report_frequency: Subscriber is NULL");
  CHECK(reportFrequency == 17);

  apo_test::CheckStatusAndDiagnostic(
      apo_subscriber_get_report_frequency(&reportFrequency, &emptySubscriber),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_subscriber_get_report_frequency: Subscriber object is NULL");
  CHECK(reportFrequency == 17);
}

TEST_CASE("CapiSubscriberMapsNativeReportFrequencyValidation") {
  apo_subscriber subscriber = MakeSubscriberHandle();
  apo_status status = APO_STATUS_OK;

  CHECK_NOTHROW((status = apo_subscriber_set_report_frequency(&subscriber, 0)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_subscriber_set_report_frequency",
      "Subscriber report frequency must be positive; observed 0",
      "src/Subscriber.cu", "setReportFrequency");
  CHECK(subscriber.object->getReportFrequency() == 1000);

  CHECK_NOTHROW(
      (status = apo_subscriber_set_report_frequency(&subscriber, -1)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_subscriber_set_report_frequency",
      "Subscriber report frequency must be positive; observed -1",
      "src/Subscriber.cu", "setReportFrequency");
  CHECK(subscriber.object->getReportFrequency() == 1000);
}

TEST_CASE("CapiDcdSubscriberMapsNativeConstructionValidation") {
  constexpr const char *PATH = "tmpCapiDcdSubscriberValidation.dcd";
  apo_test::RemoveIfExists(PATH);

  apo_dcd_subscriber *subscriber = nullptr;
  apo_status status = APO_STATUS_OK;

  apo_test::CheckStatusAndDiagnostic(
      apo_dcd_subscriber_create(&subscriber, nullptr),
      APO_STATUS_INVALID_ARGUMENT, "apo_dcd_subscriber_create: path is NULL");
  CHECK(subscriber == nullptr);

  CHECK_NOTHROW((status = apo_dcd_subscriber_create(&subscriber, "")));
  apo_test::CheckNativeError(status, APO_STATUS_INVALID_ARGUMENT,
                             "InvalidArgument", "apo_dcd_subscriber_create",
                             "Output file name must not be empty",
                             "src/Subscriber.cu", "setFileName");
  CHECK(subscriber == nullptr);

  CHECK_NOTHROW((status = apo_dcd_subscriber_create_with_report_frequency(
                     &subscriber, PATH, 0)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_dcd_subscriber_create_with_report_frequency",
      "Subscriber report frequency must be positive; observed 0",
      "src/Subscriber.cu", "setReportFrequency");
  CHECK(subscriber == nullptr);

  CHECK_NOTHROW((status = apo_dcd_subscriber_create_with_report_frequency(
                     &subscriber, PATH, -1)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_dcd_subscriber_create_with_report_frequency",
      "Subscriber report frequency must be positive; observed -1",
      "src/Subscriber.cu", "setReportFrequency");
  CHECK(subscriber == nullptr);

  std::ifstream input(PATH);
  CHECK(input.good() == false);

  CHECK_NOTHROW(
      (status = apo_dcd_subscriber_create(
           &subscriber, "missing_capi_dcd_subscriber_dir/output.dcd")));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_dcd_subscriber_create",
      "Output directory does not exist: missing_capi_dcd_subscriber_dir",
      "src/Subscriber.cu", "checkPath");
  CHECK(subscriber == nullptr);

  CHECK_NOTHROW((status = apo_dcd_subscriber_create(&subscriber, ".")));
  apo_test::CheckNativeError(status, APO_STATUS_RUNTIME_ERROR, "Runtime",
                             "apo_dcd_subscriber_create",
                             "Failed to open DCD file for writing: .",
                             "src/DcdSubscriber.cu", "openFile");
  CHECK(subscriber == nullptr);

  apo_test::RemoveIfExists(PATH);
}

TEST_CASE("CapiSubscriberAttachmentValidatesHandles") {
  apo_cuda_integrator integrator = MakeIntegratorHandle();
  apo_cuda_integrator emptyIntegrator;
  apo_subscriber subscriber = MakeSubscriberHandle();
  apo_subscriber emptySubscriber;

  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_subscribe(nullptr, &subscriber),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_subscribe: CudaIntegrator is NULL");
  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_subscribe(&emptyIntegrator, &subscriber),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_subscribe: CudaIntegrator object is NULL");
  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_subscribe(&integrator, nullptr),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_subscribe: Subscriber is NULL");
  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_subscribe(&integrator, &emptySubscriber),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_subscribe: Subscriber object is NULL");

  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_unsubscribe(nullptr, &subscriber),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_unsubscribe: CudaIntegrator is NULL");
  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_unsubscribe(&emptyIntegrator, &subscriber),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_unsubscribe: CudaIntegrator object is NULL");
  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_unsubscribe(&integrator, nullptr),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_unsubscribe: Subscriber is NULL");
  apo_test::CheckStatusAndDiagnostic(
      apo_cuda_integrator_unsubscribe(&integrator, &emptySubscriber),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_cuda_integrator_unsubscribe: Subscriber object is NULL");
}

TEST_CASE("CapiSubscriberMapsNativeAttachmentValidation") {
  apo_cuda_integrator integrator = MakeIntegratorHandle();
  apo_subscriber subscriber = MakeSubscriberHandle();
  subscriber.object->setFileName("capiSubscriber.out");

  apo_status status = APO_STATUS_OK;

  CHECK_NOTHROW(
      (status = apo_cuda_integrator_unsubscribe(&integrator, &subscriber)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_cuda_integrator_unsubscribe",
      "Subscriber not found (file \"capiSubscriber.out\")",
      "src/CudaIntegrator.cu", "unsubscribe");

  CHECK(apo_cuda_integrator_subscribe(&integrator, &subscriber) ==
        APO_STATUS_OK);
  REQUIRE(integrator.object->getSubscribers().size() == 1);
  REQUIRE(integrator.object->getReportFreqList().size() == 1);

  CHECK_NOTHROW(
      (status = apo_cuda_integrator_subscribe(&integrator, &subscriber)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_cuda_integrator_subscribe",
      "Subscriber is already subscribed to this CudaIntegrator",
      "src/CudaIntegrator.cu", "subscribe");
  CHECK(integrator.object->getSubscribers().size() == 1);
  CHECK(integrator.object->getReportFreqList().size() == 1);

  CHECK(apo_cuda_integrator_unsubscribe(&integrator, &subscriber) ==
        APO_STATUS_OK);
  CHECK(integrator.object->getSubscribers().empty() == true);
  CHECK(integrator.object->getReportFreqList().empty() == true);

  CHECK_NOTHROW(
      (status = apo_cuda_integrator_unsubscribe(&integrator, &subscriber)));
  apo_test::CheckNativeError(
      status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
      "apo_cuda_integrator_unsubscribe",
      "Subscriber not found (file \"capiSubscriber.out\")",
      "src/CudaIntegrator.cu", "unsubscribe");
}

TEST_CASE("CapiSubscriberSuccessClearsStaleDiagnostic") {
  apo_subscriber subscriber = MakeSubscriberHandle();

  REQUIRE(apo_subscriber_set_report_frequency(nullptr, 1) ==
          APO_STATUS_INVALID_ARGUMENT);
  REQUIRE(std::string(apo_last_error()).empty() == false);

  CHECK(apo_subscriber_set_report_frequency(&subscriber, 7) == APO_STATUS_OK);
  CHECK(std::string(apo_last_error()).empty() == true);
}

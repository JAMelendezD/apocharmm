// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "ApoCharmmError.h"
#include "CudaIntegrator.h"
#include "Subscriber.h"
#include "apo_test_helpers.h"
#include "catch.hpp"

namespace {

constexpr double TIME_STEP = 0.002;

class CountingSubscriber : public Subscriber {
public:
  explicit CountingSubscriber(const int reportFrequency) : Subscriber() {
    this->setReportFrequency(reportFrequency);
  }

  void update(void) override {
    updateCount++;
    return;
  }

  int updateCount = 0;
};

} // namespace

TEST_CASE("CudaIntegratorBaseConstructionAndTimeStep") {
  SECTION("DefaultConstructor") {
    CudaIntegrator integrator;

    CHECK(integrator.getTimeStep() == Approx(0.0));
    CHECK(integrator.getCharmmContext() == nullptr);
    CHECK(integrator.getSubscribers().empty() == true);
    CHECK(integrator.getReportFreqList().empty() == true);
    CHECK(integrator.getCurrentPropagatedStep() == 0);
    CHECK(integrator.getTotNumSteps() == 0);
    CHECK(integrator.getNumSteps() == 0);
    CHECK(integrator.getCoordsDelta().size() == 0);
    CHECK(integrator.getCoordsDeltaPrevious().size() == 0);
  }

  SECTION("TimeStepConstructor") {
    CudaIntegrator integrator(TIME_STEP);

    CHECK(integrator.getTimeStep() == Approx(TIME_STEP));
  }

  SECTION("TimeStepAndDebugConstructor") {
    CudaIntegrator integrator(TIME_STEP, 7);

    CHECK(integrator.getTimeStep() == Approx(TIME_STEP));
  }

  SECTION("SetTimeStep") {
    CudaIntegrator integrator;

    integrator.setTimeStep(TIME_STEP);
    CHECK(integrator.getTimeStep() == Approx(TIME_STEP));

    integrator.setTimeStep(0.001);
    CHECK(integrator.getTimeStep() == Approx(0.001));
  }

  SECTION("SetFrequenciesDoNotThrow") {
    CudaIntegrator integrator;

    CHECK_NOTHROW(integrator.setDebugPrintFrequency(10));
    CHECK_NOTHROW(integrator.setNonbondedListUpdateFrequency(5));
    CHECK_NOTHROW(integrator.setRemoveCenterOfMassFrequency(100));
  }
}

TEST_CASE("CudaIntegratorBaseInvalidOperations") {
  SECTION("PropagateWithoutCharmmContextThrows") {
    CudaIntegrator integrator(TIME_STEP);

    CHECK_THROWS_AS(integrator.propagate(1), std::invalid_argument);
  }

  SECTION("PropagateZeroStepsWithoutCharmmContextStillThrows") {
    CudaIntegrator integrator(TIME_STEP);

    CHECK_THROWS_AS(integrator.propagate(0), std::invalid_argument);
  }

  SECTION("InitializeFromEmptyRestartFileThrows") {
    CudaIntegrator integrator(TIME_STEP);

    CHECK_THROWS_AS(integrator.initializeFromRestartFile(""),
                    std::invalid_argument);
  }
}

TEST_CASE("CudaIntegratorBaseSubscriberManagement") {
  auto integrator = std::make_shared<CudaIntegrator>(TIME_STEP);

  auto sub2 = std::make_shared<CountingSubscriber>(2);
  auto sub3 = std::make_shared<CountingSubscriber>(3);

  SECTION("SubscribeNullSubscriberThrows") {
    apo_test::CheckApoCharmmError([&]() { integrator->subscribe(nullptr); },
                                  ApoCharmmErrorCode::InvalidArgument,
                                  "Subscriber must not be null");

    CHECK(integrator->getSubscribers().empty() == true);
    CHECK(integrator->getReportFreqList().empty() == true);
  }

  SECTION("SubscribeRequiresSharedOwnership") {
    CudaIntegrator stackIntegrator(TIME_STEP);

    apo_test::CheckApoCharmmError(
        [&]() { stackIntegrator.subscribe(sub2); },
        ApoCharmmErrorCode::NotInitialized,
        "CudaIntegrator must be owned by std::shared_ptr before subscribing");

    CHECK(stackIntegrator.getSubscribers().empty() == true);
    CHECK(stackIntegrator.getReportFreqList().empty() == true);
  }

  SECTION("SubscribeSingleSubscriber") {
    integrator->subscribe(sub2);

    REQUIRE(integrator->getSubscribers().size() == 1);
    REQUIRE(integrator->getReportFreqList().size() == 1);

    CHECK(integrator->getSubscribers()[0] == sub2);
    CHECK(integrator->getReportFreqList()[0] == 2);
  }

  SECTION("SubscribeDuplicateSubscriberThrowsWithoutDuplicatingLists") {
    integrator->subscribe(sub2);

    apo_test::CheckApoCharmmError(
        [&]() { integrator->subscribe(sub2); },
        ApoCharmmErrorCode::InvalidArgument,
        "Subscriber is already subscribed to this CudaIntegrator");

    REQUIRE(integrator->getSubscribers().size() == 1);
    REQUIRE(integrator->getReportFreqList().size() == 1);
    CHECK(integrator->getSubscribers()[0] == sub2);
    CHECK(integrator->getReportFreqList()[0] == 2);
  }

  SECTION("SubscribeSubscriberList") {
    integrator->subscribe(std::vector<std::shared_ptr<Subscriber>>{sub2, sub3});

    REQUIRE(integrator->getSubscribers().size() == 2);
    REQUIRE(integrator->getReportFreqList().size() == 2);

    CHECK(integrator->getSubscribers()[0] == sub2);
    CHECK(integrator->getSubscribers()[1] == sub3);
    CHECK(integrator->getReportFreqList() == std::vector<int>{2, 3});
  }

  SECTION("UnsubscribeSingleSubscriber") {
    integrator->subscribe(std::vector<std::shared_ptr<Subscriber>>{sub2, sub3});

    integrator->unsubscribe(sub2);

    REQUIRE(integrator->getSubscribers().size() == 1);
    REQUIRE(integrator->getReportFreqList().size() == 1);

    CHECK(integrator->getSubscribers()[0] == sub3);
    CHECK(integrator->getReportFreqList()[0] == 3);
  }

  SECTION("UnsubscribeSubscriberList") {
    integrator->subscribe(std::vector<std::shared_ptr<Subscriber>>{sub2, sub3});

    integrator->unsubscribe(
        std::vector<std::shared_ptr<Subscriber>>{sub2, sub3});

    CHECK(integrator->getSubscribers().empty() == true);
    CHECK(integrator->getReportFreqList().empty() == true);
  }

  SECTION("UnsubscribeNullSubscriberThrows") {
    apo_test::CheckApoCharmmError([&]() { integrator->unsubscribe(nullptr); },
                                  ApoCharmmErrorCode::InvalidArgument,
                                  "Subscriber must not be null");

    CHECK(integrator->getSubscribers().empty() == true);
    CHECK(integrator->getReportFreqList().empty() == true);
  }

  SECTION("UnsubscribeMissingSubscriberThrows") {
    apo_test::CheckApoCharmmError([&]() { integrator->unsubscribe(sub2); },
                                  ApoCharmmErrorCode::InvalidArgument,
                                  "Subscriber not found (file \"\")");
  }
}

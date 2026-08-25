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
#include "Subscriber.h"
#include "apo_test_helpers.h"
#include "catch.hpp"

#include <filesystem>
#include <fstream>
#include <string>

namespace {

class TestSubscriber final : public Subscriber {
public:
  using Subscriber::Subscriber;

  void update(void) override { return; }
};

} // namespace

TEST_CASE("SubscriberConstructionAndReportFrequency") {
  SECTION("DefaultConstruction") {
    TestSubscriber subscriber;

    CHECK(subscriber.getFilePath().empty() == true);
    CHECK(subscriber.getReportFrequency() == 1000);
  }

  SECTION("SetPositiveReportFrequency") {
    TestSubscriber subscriber;

    subscriber.setReportFrequency(17);
    CHECK(subscriber.getReportFrequency() == 17);
  }

  SECTION("RejectsZeroReportFrequency") {
    TestSubscriber subscriber;

    apo_test::CheckApoCharmmError(
        [&]() { subscriber.setReportFrequency(0); },
        ApoCharmmErrorCode::InvalidArgument,
        "Subscriber report frequency must be positive; observed 0");
    CHECK(subscriber.getReportFrequency() == 1000);
  }

  SECTION("RejectsNegativeReportFrequency") {
    TestSubscriber subscriber;

    apo_test::CheckApoCharmmError(
        [&]() { subscriber.setReportFrequency(-7); },
        ApoCharmmErrorCode::InvalidArgument,
        "Subscriber report frequency must be positive; observed -7");
    CHECK(subscriber.getReportFrequency() == 1000);
  }

  SECTION("ConstructorValidatesReportFrequencyBeforeOpeningFile") {
    const std::filesystem::path filePath = "tmpSubscriberInvalidFrequency.txt";
    apo_test::RemoveIfExists(filePath);

    apo_test::CheckApoCharmmError(
        [&]() { (void)TestSubscriber(filePath, 0); },
        ApoCharmmErrorCode::InvalidArgument,
        "Subscriber report frequency must be positive; observed 0");

    std::ifstream input(filePath);
    CHECK(input.good() == false);

    apo_test::RemoveIfExists(filePath);
  }
}

TEST_CASE("SubscriberValidatesOutputFile") {
  SECTION("RejectsEmptyFileName") {
    TestSubscriber subscriber;

    apo_test::CheckApoCharmmError(
        [&]() { subscriber.setFilePath(std::filesystem::path{}); },
        ApoCharmmErrorCode::InvalidArgument,
        "Output file name must not be empty");
  }

  SECTION("RejectsMissingOutputDirectory") {
    const std::filesystem::path filePath =
        std::filesystem::path{"missing_subscriber_dir"} / "output.txt";

    apo_test::CheckApoCharmmError(
        [&]() { (void)TestSubscriber(filePath); },
        ApoCharmmErrorCode::InvalidArgument,
        "Output directory does not exist: missing_subscriber_dir");
  }

  SECTION("OpenRequiresFileName") {
    TestSubscriber subscriber;

    apo_test::CheckApoCharmmError([&]() { subscriber.openFile(); },
                                  ApoCharmmErrorCode::NotInitialized,
                                  "Subscriber output file name is not set");
  }

  SECTION("ReportsOutputOpenFailure") {
    apo_test::CheckApoCharmmError(
        []() { (void)TestSubscriber("."); }, ApoCharmmErrorCode::Runtime,
        "Failed to open subscriber output file for writing: .");
  }
}

TEST_CASE("SubscriberCommentSections") {
  SECTION("RejectsEmptyComment") {
    TestSubscriber subscriber;

    apo_test::CheckApoCharmmError([&]() { subscriber.addCommentSection(""); },
                                  ApoCharmmErrorCode::InvalidArgument,
                                  "Comment section must not be empty");
  }

  SECTION("RequiresOpenOutputFile") {
    TestSubscriber subscriber;

    apo_test::CheckApoCharmmError(
        [&]() { subscriber.addCommentSection("comment"); },
        ApoCharmmErrorCode::NotInitialized,
        "Subscriber output file is not open");
  }

  SECTION("WritesNewlineTerminatedComments") {
    const std::filesystem::path filePath = "tmpSubscriberComments.txt";
    apo_test::RemoveIfExists(filePath);

    {
      TestSubscriber subscriber(filePath);
      subscriber.addCommentSection("first");
      subscriber.addCommentSection("second\n");
    }

    std::ifstream input(filePath);
    REQUIRE(input.good());

    std::string firstLine;
    std::string secondLine;
    REQUIRE(static_cast<bool>(std::getline(input, firstLine)));
    REQUIRE(static_cast<bool>(std::getline(input, secondLine)));
    CHECK(firstLine == "first");
    CHECK(secondLine == "second");

    apo_test::RemoveIfExists(filePath);
  }
}

#if defined(__linux__)
TEST_CASE("SubscriberReportsCommentWriteFailure") {
  TestSubscriber subscriber("/dev/full");

  apo_test::CheckApoCharmmError(
      [&]() { subscriber.addCommentSection("comment"); },
      ApoCharmmErrorCode::Runtime,
      "Failed to write comment section to subscriber output file: /dev/full");
}
#endif

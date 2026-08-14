// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "apo_test_helpers.h"
#include "catch.hpp"
#include "str_utils.h"

#include <limits>
#include <sstream>
#include <string>
#include <string_view>

TEST_CASE("StrUtilsTryGetLineHandlesFinalAndEmptyLines") {
  const std::string fileData = "first\n\nthird";
  std::size_t position = 0;
  std::size_t lineNumber = 0;
  std::string line = "stale";

  CHECK(apo::try_get_line(line, position, lineNumber, fileData));
  CHECK(line == "first");
  CHECK(lineNumber == 1);

  CHECK(apo::try_get_line(line, position, lineNumber, fileData));
  CHECK(line.empty());
  CHECK(lineNumber == 2);

  CHECK(apo::try_get_line(line, position, lineNumber, fileData));
  CHECK(line == "third");
  CHECK(lineNumber == 3);
  CHECK(position == fileData.size());

  CHECK_FALSE(apo::try_get_line(line, position, lineNumber, fileData));
  CHECK(line.empty());
  CHECK(position == fileData.size());
  CHECK(lineNumber == 3);
}

TEST_CASE("StrUtilsGetLineReportsUnexpectedEndOfFile") {
  std::size_t position = 0;
  std::size_t lineNumber = 0;
  std::string line;

  apo_test::CheckApoCharmmError(
      [&line, &position, &lineNumber](void) {
        apo::get_line(line, position, lineNumber, "", "ATOM records",
                      "PSF \"test.psf\"");
      },
      ApoCharmmErrorCode::Runtime,
      "Unexpected end of file while reading ATOM records in PSF \"test.psf\"");
}

TEST_CASE("StrUtilsExtractsFixedWidthFieldsSafely") {
  const std::string line = "0123456789";
  std::string_view field = "stale";

  CHECK(apo::try_get_fixed_width_field(field, line, 2, 4));
  CHECK(field == "2345");

  CHECK(apo::try_get_fixed_width_field(field, line, line.size(), 0));
  CHECK(field.empty());

  CHECK_FALSE(apo::try_get_fixed_width_field(field, line, 8, 3));
  CHECK(field.empty());

  CHECK_FALSE(apo::try_get_fixed_width_field(
      field, line, std::numeric_limits<std::size_t>::max(), 1));
  CHECK(field.empty());
}

TEST_CASE("StrUtilsParsesIntegersStrictly") {
  int value = 99;

  CHECK(apo::try_parse_int(value, " -17 "));
  CHECK(value == -17);

  value = 99;
  CHECK_FALSE(apo::try_parse_int(value, ""));
  CHECK(value == 99);

  CHECK_FALSE(apo::try_parse_int(value, "12BAD"));
  CHECK(value == 99);

  CHECK_FALSE(apo::try_parse_int(value, "999999999999999999999"));
  CHECK(value == 99);

  apo_test::CheckApoCharmmError(
      [](void) {
        static_cast<void>(apo::parse_int(
            "BAD", "atom index", "BOND section of PSF \"test.psf\" at line 7"));
      },
      ApoCharmmErrorCode::Runtime,
      "Invalid atom index value \"BAD\" in BOND section of PSF \"test.psf\" "
      "at line 7");
}

TEST_CASE("StrUtilsParsesUnsignedIntegersStrictly") {
  unsigned long long int value = 99;

  CHECK(apo::try_parse_ull(value, " 17 "));
  CHECK(value == 17);

  value = 99;
  CHECK_FALSE(apo::try_parse_ull(value, "-1"));
  CHECK(value == 99);

  CHECK_FALSE(apo::try_parse_ull(value, "12BAD"));
  CHECK(value == 99);

  const std::string outOfRange =
      std::to_string(std::numeric_limits<unsigned long long int>::max()) + "0";
  CHECK_FALSE(apo::try_parse_ull(value, outOfRange));
  CHECK(value == 99);

  apo_test::CheckApoCharmmError(
      [](void) {
        static_cast<void>(apo::parse_ull(
            "-1", "count", "ATOM section of PSF \"test.psf\" at line 4"));
      },
      ApoCharmmErrorCode::Runtime,
      "Invalid count value \"-1\" in ATOM section of PSF \"test.psf\" at "
      "line 4");
}

TEST_CASE("StrUtilsParsesFloatingPointValuesStrictly") {
  double value = 99.0;

  CHECK(apo::try_parse_finite_double(value, "1.25"));
  CHECK(value == Approx(1.25));

  CHECK(apo::try_parse_finite_double(value, "-2.5E+01"));
  CHECK(value == Approx(-25.0));

  CHECK(apo::try_parse_finite_double(value, "1.25D+02"));
  CHECK(value == Approx(125.0));

  CHECK(apo::try_parse_finite_double(value, "2.5d-01"));
  CHECK(value == Approx(0.25));

  value = 99.0;
  CHECK_FALSE(apo::try_parse_finite_double(value, "1.0BAD"));
  CHECK(value == Approx(99.0));

  CHECK_FALSE(apo::try_parse_finite_double(value, "1E9999"));
  CHECK(value == Approx(99.0));

  CHECK_FALSE(apo::try_parse_finite_double(value, "NAN"));
  CHECK(value == Approx(99.0));

  CHECK_FALSE(apo::try_parse_finite_double(value, "INF"));
  CHECK(value == Approx(99.0));

  double nonfiniteValue = 0.0;
  CHECK(apo::try_parse_double(nonfiniteValue, "NAN"));
  CHECK(std::isnan(nonfiniteValue));

  apo_test::CheckApoCharmmError(
      [](void) {
        static_cast<void>(apo::parse_double(
            "BAD", "charge", "ATOM section of PSF \"test.psf\" at line 5"));
      },
      ApoCharmmErrorCode::Runtime,
      "Invalid charge value \"BAD\" in ATOM section of PSF \"test.psf\" at "
      "line 5");
}

TEST_CASE("StrUtilsReadsStreamsAndNormalizesLineEndings") {
  std::istringstream input("first\r\nsecond\nthird");
  std::size_t lineNumber = 0;
  std::string line = "stale";

  CHECK(apo::try_get_line(line, lineNumber, input));
  CHECK(line == "first");
  CHECK(lineNumber == 1);

  apo::get_line(line, lineNumber, input, "second record", "test stream");
  CHECK(line == "second");
  CHECK(lineNumber == 2);

  apo::get_line(line, lineNumber, input, "third record", "test stream");
  CHECK(line == "third");
  CHECK(lineNumber == 3);

  CHECK_FALSE(apo::try_get_line(line, lineNumber, input));
  CHECK(line.empty());
  CHECK(lineNumber == 3);
}

TEST_CASE("StrUtilsFindsRequiredStreamLines") {
  std::istringstream input("header\n"
                           "ignored\n"
                           " !TARGET\r\n"
                           "payload\n");

  std::size_t lineNumber = 0;
  std::string line;

  apo::find_required_line(input, lineNumber, " !TARGET", "target section",
                          "test stream");
  CHECK(lineNumber == 3);

  apo::get_line(line, lineNumber, input, "payload record", "test stream");
  CHECK(line == "payload");
  CHECK(lineNumber == 4);
}

TEST_CASE("StrUtilsReportsStreamReadFailures") {
  SECTION("UnexpectedEndOfFile") {
    std::istringstream input("");
    std::size_t lineNumber = 0;
    std::string line;

    apo_test::CheckApoCharmmError(
        [&]() {
          apo::get_line(line, lineNumber, input, "ATOM record", "test stream");
        },
        ApoCharmmErrorCode::Runtime,
        "Unexpected end of file while reading ATOM record in test stream");
  }

  SECTION("MissingRequiredLine") {
    std::istringstream input("header\npayload\n");
    std::size_t lineNumber = 0;

    apo_test::CheckApoCharmmError(
        [&]() {
          apo::find_required_line(input, lineNumber, " !TARGET",
                                  "target section", "test stream");
        },
        ApoCharmmErrorCode::Runtime,
        "Could not find target section in test stream");
  }
}

TEST_CASE("StrUtilsParsesFixedWidthFields") {
  CHECK(apo::get_fixed_width_field("0123456789", 2, 4, "digits",
                                   "test record") == "2345");

  CHECK(apo::parse_fixed_width_int("  -17  ", 0, 7, "count", "test record") ==
        -17);

  CHECK(apo::parse_fixed_width_ull("  42  ", 0, 6, "count", "test record") ==
        42ULL);

  CHECK(apo::parse_fixed_width_double(" 1.25D+02 ", 0, 10, "value",
                                      "test record") == Approx(125.0));

  apo_test::CheckApoCharmmError(
      [](void) {
        static_cast<void>(
            apo::get_fixed_width_field("123", 2, 2, "count", "test record"));
      },
      ApoCharmmErrorCode::Runtime,
      "Field \"count\" is truncated in test record");

  apo_test::CheckApoCharmmError(
      [](void) {
        static_cast<void>(
            apo::parse_fixed_width_int(" BAD ", 0, 5, "count", "test record"));
      },
      ApoCharmmErrorCode::Runtime,
      "Invalid count value \"BAD\" in test record");

  apo_test::CheckApoCharmmError(
      [](void) {
        static_cast<void>(
            apo::parse_fixed_width_ull(" -1 ", 0, 4, "count", "test record"));
      },
      ApoCharmmErrorCode::Runtime, "Invalid count value \"-1\" in test record");

  apo_test::CheckApoCharmmError(
      [](void) {
        static_cast<void>(apo::parse_fixed_width_double(
            " INF ", 0, 5, "floating-point", "test record"));
      },
      ApoCharmmErrorCode::Runtime,
      "Invalid floating-point value \"INF\" in test record");
}

TEST_CASE("StrUtilsSplitRejectsEmptyDelimiter") {
  apo_test::CheckApoCharmmError(
      [](void) { static_cast<void>(apo::split("ABC", "")); },
      ApoCharmmErrorCode::InvalidArgument, "Split delimiter must not be empty");
}

TEST_CASE("StrUtilsSplitPreservesCurrentWhitespaceBehavior") {
  CHECK(apo::split("  A   B C  ") == std::vector<std::string>{"A", "B", "C"});

  CHECK(apo::split("A,,B", ",") == std::vector<std::string>{"A", "", "B"});

  CHECK(apo::split("A,   B", ",") == std::vector<std::string>{"A", "B"});
}

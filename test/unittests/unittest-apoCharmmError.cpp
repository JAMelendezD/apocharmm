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
#include "catch.hpp"

#include <exception>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace {

constexpr std::string_view EXPECTED_SKULL =
    R"ASCII(                            /---------\
                           /           \
                          /             \
                         /               \
                         !         XXXX  !
                         !  ----   XXXX  !
                         !          XXX  !
                         !       X       !
                          --\   XXX   /--
                           ! !  XXX  ! !
                           ! !       ! !
                           ! I I I I I !
                           !  I I I I  !
                            \         /
                             --     --
                               \---/
                        XXX             XXX
                       XXXX             XXXX
                       XXXXX           XXXXX
                          XXX         XXX
                            XXX     XXX
                               XXXXX
                              XXX XXX
                            XXX     XXX
                          XXX         XXX
                       XXXXX           XXXXX
                       XXXX             XXXX
                        XXX             XXX)ASCII";

static_assert(
    noexcept(GetApoCharmmErrorCodeName(ApoCharmmErrorCode::InvalidArgument)));
static_assert(noexcept(std::declval<const ApoCharmmError &>().getCode()));
static_assert(noexcept(std::declval<const ApoCharmmError &>().getMessage()));
static_assert(noexcept(std::declval<const ApoCharmmError &>().getSourceFile()));
static_assert(
    noexcept(std::declval<const ApoCharmmError &>().getSourceFunction()));
static_assert(noexcept(std::declval<const ApoCharmmError &>().getSourceLine()));
static_assert(noexcept(std::declval<const ApoCharmmError &>().what()));
static_assert(std::is_copy_constructible_v<ApoCharmmError>);
static_assert(std::is_move_constructible_v<ApoCharmmError>);

[[noreturn]] void ThrowApoCharmmErrorForTest(const ApoCharmmErrorCode code,
                                             const std::string_view message,
                                             int &expectedSourceLine) {
  expectedSourceLine = __LINE__ + 1;
  APOCHARMM_THROW(code, message);
}

} // namespace

TEST_CASE("ApoCharmmErrorCodeValuesAndNames") {
  CHECK(static_cast<int>(ApoCharmmErrorCode::InvalidArgument) == 0);
  CHECK(static_cast<int>(ApoCharmmErrorCode::Runtime) == 1);
  CHECK(static_cast<int>(ApoCharmmErrorCode::Cuda) == 2);
  CHECK(static_cast<int>(ApoCharmmErrorCode::NotInitialized) == 3);
  CHECK(static_cast<int>(ApoCharmmErrorCode::NotImplemented) == 4);

  CHECK(GetApoCharmmErrorCodeName(ApoCharmmErrorCode::InvalidArgument) ==
        "InvalidArgument");
  CHECK(GetApoCharmmErrorCodeName(ApoCharmmErrorCode::Runtime) == "Runtime");
  CHECK(GetApoCharmmErrorCodeName(ApoCharmmErrorCode::Cuda) == "Cuda");
  CHECK(GetApoCharmmErrorCodeName(ApoCharmmErrorCode::NotInitialized) ==
        "NotInitialized");
  CHECK(GetApoCharmmErrorCodeName(ApoCharmmErrorCode::NotImplemented) ==
        "NotImplemented");
  CHECK(GetApoCharmmErrorCodeName(static_cast<ApoCharmmErrorCode>(-1)) ==
        "Unknown");
}

TEST_CASE("ApoCharmmErrorAccessorsAndFormatting") {
  const ApoCharmmError error(ApoCharmmErrorCode::Runtime, "integration failed",
                             "error_test.cpp", "FormatTest", 37);

  CHECK(error.getCode() == ApoCharmmErrorCode::Runtime);
  CHECK(error.getMessage() == "integration failed");
  CHECK(error.getSourceFile() == "error_test.cpp");
  CHECK(error.getSourceFunction() == "FormatTest");
  CHECK(error.getSourceLine() == 37);

  const std::string expected =
      "\n\n" + std::string(EXPECTED_SKULL) +
      "\n\napoCHARMM error [Runtime]: integration failed\n  source: "
      "error_test.cpp:37\n  function: FormatTest";
  CHECK(std::string(error.what()) == expected);
}

TEST_CASE("ApoCharmmErrorMessageEdgeCases") {
  SECTION("EmptyHumanMessage") {
    const ApoCharmmError error(ApoCharmmErrorCode::InvalidArgument, "",
                               "empty.cpp", "EmptyMessage", 11);
    const std::string diagnostic(error.what());

    CHECK(error.getMessage().empty() == true);
    CHECK(diagnostic.find("apoCHARMM error [InvalidArgument]: \n  source: "
                          "empty.cpp:11\n  function: EmptyMessage") !=
          std::string::npos);
  }

  SECTION("MultilineHumanMessage") {
    constexpr std::string_view message = "first line\nsecond line";
    const ApoCharmmError error(ApoCharmmErrorCode::Cuda, message,
                               "multiline.cpp", "MultilineMessage", 19);
    const std::string diagnostic(error.what());

    CHECK(error.getMessage() == message);
    CHECK(diagnostic.find(
              "apoCHARMM error [Cuda]: first line\nsecond line\n  source: "
              "multiline.cpp:19\n  function: MultilineMessage") !=
          std::string_view::npos);
  }

  SECTION("TrailingNewlinesRemainWellFormed") {
    constexpr std::string_view message = "trailing line\n\n";
    const ApoCharmmError error(ApoCharmmErrorCode::NotInitialized, message,
                               "trailing.cpp", "TrailingNewlines", 23);
    const std::string diagnostic(error.what());

    CHECK(error.getMessage() == message);
    CHECK(diagnostic.find(
              "apoCHARMM error [NotInitialized]: trailing line\n\n\n  source: "
              "trailing.cpp:23\n  function: TrailingNewlines") !=
          std::string_view::npos);
  }
}

TEST_CASE("ApoCharmmErrorCopyandMoveConstruction") {
  SECTION("CopyConstruction") {
    const ApoCharmmError original(ApoCharmmErrorCode::Cuda, "copy message",
                                  "copy.cpp", "CopyConstruction", 29);
    const ApoCharmmError copy(original);

    CHECK(copy.getCode() == original.getCode());
    CHECK(copy.getMessage() == original.getMessage());
    CHECK(copy.getSourceFile() == original.getSourceFile());
    CHECK(copy.getSourceFunction() == original.getSourceFunction());
    CHECK(copy.getSourceLine() == original.getSourceLine());
    CHECK(copy.what() == original.what());
  }

  SECTION("MoveConstruction") {
    ApoCharmmError original(ApoCharmmErrorCode::NotImplemented, "move message",
                            "move.cpp", "MoveConstruction", 31);
    const std::string expectedWhat(original.what());

    const ApoCharmmError moved(std::move(original));

    CHECK(moved.getCode() == ApoCharmmErrorCode::NotImplemented);
    CHECK(moved.getMessage() == "move message");
    CHECK(moved.getSourceFile() == "move.cpp");
    CHECK(moved.getSourceFunction() == "MoveConstruction");
    CHECK(moved.getSourceLine() == 31);
    CHECK(moved.what() == expectedWhat);
  }
}

TEST_CASE("ApoCharmmThrowCapturesSourceAndInheritance") {
  SECTION("CapturesSourceMetadata") {
    int expectedSourceLine = 0;
    bool caught = false;

    try {
      ThrowApoCharmmErrorForTest(ApoCharmmErrorCode::NotImplemented,
                                 "macro failure", expectedSourceLine);
    } catch (const ApoCharmmError &error) {
      caught = true;

      const std::string sourceFile(error.getSourceFile());
      const std::string diagnostic(error.what());
      const std::string expectedSource =
          "  source: " + sourceFile + ':' + std::to_string(expectedSourceLine);

      CHECK(error.getCode() == ApoCharmmErrorCode::NotImplemented);
      CHECK(error.getMessage() == "macro failure");
      REQUIRE(sourceFile.empty() == false);
      CHECK(sourceFile.find("unittest-apoCharmmError.cpp") !=
            std::string::npos);
      CHECK(sourceFile.front() != '/');
      CHECK(sourceFile.front() != '\\');
      CHECK(sourceFile.find(":/") == std::string::npos);
      CHECK(sourceFile.find(":\\") == std::string::npos);
      CHECK(error.getSourceFunction() == "ThrowApoCharmmErrorForTest");
      CHECK(error.getSourceLine() == expectedSourceLine);
      CHECK(
          diagnostic.find("apoCHARMM error [NotImplemented]: macro failure") !=
          std::string::npos);
      CHECK(diagnostic.find(expectedSource) != std::string::npos);
      CHECK(diagnostic.find("  function: ThrowApoCharmmErrorForTest") !=
            std::string::npos);
    }

    CHECK(caught == true);
  }

  SECTION("CatchableAsApoCharmmError") {
    int expectedSourceLine = 0;
    CHECK_THROWS_AS(ThrowApoCharmmErrorForTest(ApoCharmmErrorCode::Runtime,
                                               "derived", expectedSourceLine),
                    ApoCharmmError);
  }

  SECTION("CatchableAsRuntimeError") {
    int expectedSourceLine = 0;
    CHECK_THROWS_AS(ThrowApoCharmmErrorForTest(ApoCharmmErrorCode::Runtime,
                                               "runtime", expectedSourceLine),
                    std::runtime_error);
  }

  SECTION("CatchableAsException") {
    int expectedSourceLine = 0;
    CHECK_THROWS_AS(ThrowApoCharmmErrorForTest(ApoCharmmErrorCode::Runtime,
                                               "exception", expectedSourceLine),
                    std::exception);
  }
}

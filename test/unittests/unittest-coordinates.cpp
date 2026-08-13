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
#include "Coordinates.h"
#include "apo_test_helpers.h"
#include "catch.hpp"

#include <type_traits>

static_assert(std::is_nothrow_destructible_v<Coordinates>);

namespace {

constexpr double DOUBLE_TOLERANCE = 1.0e-12;
constexpr double FLOAT_TOLERANCE = 1.0e-6;

} // namespace

TEST_CASE("CoordinatesDefaultConstructor") {
  Coordinates coordinates;

  CHECK(coordinates.getNumAtoms() == -1);
  CHECK(coordinates.getCoordinatesDP().empty() == true);
  CHECK(coordinates.getCoordinatesSP().empty() == true);
}

TEST_CASE("CoordinatesFromDouble3Vector") {
  const std::vector<double3> input = {
      {1.25, -2.5, 3.75}, {0.0, 4.5, -6.25}, {-7.125, 8.0, 9.5}};
  const std::vector<double3> expectedDP = {
      {1.25, -2.5, 3.75}, {0.0, 4.5, -6.25}, {-7.125, 8.0, 9.5}};
  const std::vector<float3> expectedSP = {
      {1.25f, -2.5f, 3.75f}, {0.0f, 4.5f, -6.25f}, {-7.125f, 8.0f, 9.5f}};

  Coordinates coordinates(input);

  apo_test::CheckVectorsClose<double3>(coordinates.getCoordinatesDP(),
                                       expectedDP, DOUBLE_TOLERANCE);
  apo_test::CheckVectorsClose<float3>(coordinates.getCoordinatesSP(),
                                      expectedSP, FLOAT_TOLERANCE);
}

TEST_CASE("CoordinatesFromFloat3Vector") {
  const std::vector<float3> input = {
      {1.25f, -2.5f, 3.75f}, {0.0f, 4.5f, -6.25f}, {-7.125f, 8.0f, 9.5f}};
  const std::vector<double3> expectedDP = {
      {1.25, -2.5, 3.75}, {0.0, 4.5, -6.25}, {-7.125, 8.0, 9.5}};
  const std::vector<float3> expectedSP = {
      {1.25f, -2.5f, 3.75f}, {0.0f, 4.5f, -6.25f}, {-7.125f, 8.0f, 9.5f}};

  Coordinates coordinates(input);

  apo_test::CheckVectorsClose<double3>(coordinates.getCoordinatesDP(),
                                       expectedDP, DOUBLE_TOLERANCE);
  apo_test::CheckVectorsClose<float3>(coordinates.getCoordinatesSP(),
                                      expectedSP, FLOAT_TOLERANCE);
}

TEST_CASE("CoordinatesFromDoubleNestedVector") {
  const std::vector<std::vector<double>> input = {
      {1.25, -2.5, 3.75}, {0.0, 4.5, -6.25}, {-7.125, 8.0, 9.5}};
  const std::vector<double3> expectedDP = {
      {1.25, -2.5, 3.75}, {0.0, 4.5, -6.25}, {-7.125, 8.0, 9.5}};
  const std::vector<float3> expectedSP = {
      {1.25f, -2.5f, 3.75f}, {0.0f, 4.5f, -6.25f}, {-7.125f, 8.0f, 9.5f}};

  Coordinates coordinates(input);

  apo_test::CheckVectorsClose<double3>(coordinates.getCoordinatesDP(),
                                       expectedDP, DOUBLE_TOLERANCE);
  apo_test::CheckVectorsClose<float3>(coordinates.getCoordinatesSP(),
                                      expectedSP, FLOAT_TOLERANCE);
}

TEST_CASE("CoordinatesFromFloatNestedVector") {
  const std::vector<std::vector<float>> input = {
      {1.25f, -2.5f, 3.75f}, {0.0f, 4.5f, -6.25f}, {-7.125f, 8.0f, 9.5f}};
  const std::vector<double3> expectedDP = {
      {1.25, -2.5, 3.75}, {0.0, 4.5, -6.25}, {-7.125, 8.0, 9.5}};
  const std::vector<float3> expectedSP = {
      {1.25f, -2.5f, 3.75f}, {0.0f, 4.5f, -6.25f}, {-7.125f, 8.0f, 9.5f}};

  Coordinates coordinates(input);

  apo_test::CheckVectorsClose<double3>(coordinates.getCoordinatesDP(),
                                       expectedDP, DOUBLE_TOLERANCE);
  apo_test::CheckVectorsClose<float3>(coordinates.getCoordinatesSP(),
                                      expectedSP, FLOAT_TOLERANCE);
}

TEST_CASE("CoordinatesEmptyInputVectors") {
  SECTION("EmptyDouble3Vector") {
    const std::vector<double3> input;
    Coordinates coordinates(input);

    CHECK(coordinates.getNumAtoms() == 0);
    CHECK(coordinates.getCoordinatesDP().empty() == true);
    CHECK(coordinates.getCoordinatesSP().empty() == true);
  }

  SECTION("EmptyFloat3Vector") {
    const std::vector<float3> input;
    Coordinates coordinates(input);

    CHECK(coordinates.getNumAtoms() == 0);
    CHECK(coordinates.getCoordinatesDP().empty() == true);
    CHECK(coordinates.getCoordinatesSP().empty() == true);
  }

  SECTION("EmptyDoubleNestedVector") {
    const std::vector<std::vector<double>> input;
    Coordinates coordinates(input);

    CHECK(coordinates.getNumAtoms() == 0);
    CHECK(coordinates.getCoordinatesDP().empty() == true);
    CHECK(coordinates.getCoordinatesSP().empty() == true);
  }

  SECTION("EmptyFloatNestedVector") {
    const std::vector<std::vector<float>> input;
    Coordinates coordinates(input);

    CHECK(coordinates.getNumAtoms() == 0);
    CHECK(coordinates.getCoordinatesDP().empty() == true);
    CHECK(coordinates.getCoordinatesSP().empty() == true);
  }
}

TEST_CASE("CoordinatesInputValidationUsesApoCharmmError") {
  SECTION("DoubleNestedVectorEntryWidth") {
    const std::vector<std::vector<double>> input = {{1.0, 2.0, 3.0},
                                                    {4.0, 5.0}};

    apo_test::CheckApoCharmmError(
        [&input](void) {
          Coordinates coordinates(input);
          static_cast<void>(coordinates);
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Coordinate entry 1 must contain exactly 3 values; observed 2");
  }

  SECTION("FloatNestedVectorEntryWidth") {
    const std::vector<std::vector<float>> input = {{1.0f, 2.0f, 3.0f},
                                                   {4.0f, 5.0f, 6.0f, 7.0f}};

    apo_test::CheckApoCharmmError(
        [&input](void) {
          Coordinates coordinates(input);
          static_cast<void>(coordinates);
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Coordinate entry 1 must contain exactly 3 values; observed 4");
  }

  SECTION("NegativeAtomCount") {
    Coordinates coordinates;
    coordinates.setNumAtoms(2);

    apo_test::CheckApoCharmmError(
        [&coordinates](void) { coordinates.setNumAtoms(-1); },
        ApoCharmmErrorCode::InvalidArgument,
        "Atom count must be non-negative; observed -1");

    CHECK(coordinates.getNumAtoms() == 2);
    CHECK(coordinates.getCoordinatesDP().size() == 2);
    CHECK(coordinates.getCoordinatesSP().size() == 2);
  }
}

TEST_CASE("CoordinatesSetNumAtoms") {
  Coordinates coordinates;

  coordinates.setNumAtoms(3);

  CHECK(coordinates.getNumAtoms() == 3);
  CHECK(coordinates.getCoordinatesDP().size() == 3);
  CHECK(coordinates.getCoordinatesSP().size() == 3);

  coordinates.setNumAtoms(1);

  CHECK(coordinates.getNumAtoms() == 1);
  CHECK(coordinates.getCoordinatesDP().size() == 1);
  CHECK(coordinates.getCoordinatesSP().size() == 1);

  coordinates.setNumAtoms(0);

  CHECK(coordinates.getNumAtoms() == 0);
  CHECK(coordinates.getCoordinatesDP().empty() == true);
  CHECK(coordinates.getCoordinatesSP().empty() == true);
}

TEST_CASE("CoordinatesMutableGetters") {
  Coordinates coordinates;
  coordinates.setNumAtoms(2);

  std::vector<double3> &coordinatesDP = coordinates.getCoordinatesDP();
  std::vector<float3> &coordinatesSP = coordinates.getCoordinatesSP();

  coordinatesDP[0] = {1.0, 2.0, 3.0};
  coordinatesDP[1] = {-1.0, -2.0, -3.0};
  coordinatesSP[0] = {5.0f, 6.0f, 7.0f};
  coordinatesSP[1] = {-5.0f, -6.0f, -7.0f};

  const Coordinates &constCoordinates = coordinates;
  const std::vector<double3> &constCoordinatesDP =
      constCoordinates.getCoordinatesDP();
  const std::vector<float3> &constCoordinatesSP =
      constCoordinates.getCoordinatesSP();

  REQUIRE(constCoordinatesDP.size() == 2);
  REQUIRE(constCoordinatesSP.size() == 2);

  CHECK(constCoordinatesDP[0].x == Approx(1.0));
  CHECK(constCoordinatesDP[0].y == Approx(2.0));
  CHECK(constCoordinatesDP[0].z == Approx(3.0));
  CHECK(constCoordinatesDP[1].x == Approx(-1.0));
  CHECK(constCoordinatesDP[1].y == Approx(-2.0));
  CHECK(constCoordinatesDP[1].z == Approx(-3.0));

  CHECK(constCoordinatesSP[0].x == Approx(5.0f));
  CHECK(constCoordinatesSP[0].y == Approx(6.0f));
  CHECK(constCoordinatesSP[0].z == Approx(7.0f));
  CHECK(constCoordinatesSP[1].x == Approx(-5.0f));
  CHECK(constCoordinatesSP[1].y == Approx(-6.0f));
  CHECK(constCoordinatesSP[1].z == Approx(-7.0f));
}

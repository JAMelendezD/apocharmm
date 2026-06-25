// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "Coordinates.h"
#include "apo_test_helpers.h"
#include "catch.hpp"

namespace {

constexpr double DOUBLE_TOLERANCE = 1.0e-12;
constexpr double FLOAT_TOLERANCE = 1.0e-6;

} // namespace

TEST_CASE("CoordinatesDefaultConstructor") {
  Coordinates coordinates;

  CHECK(coordinates.getNumAtoms() == -1);
  CHECK(coordinates.getCoordinatesD().empty() == true);
  CHECK(coordinates.getCoordinatesF().empty() == true);
}

TEST_CASE("CoordinatesFromDouble3Vector") {
  const std::vector<double3> input = {
      {1.25, -2.5, 3.75}, {0.0, 4.5, -6.25}, {-7.125, 8.0, 9.5}};
  const std::vector<double4> expectedD = {
      {1.25, -2.5, 3.75, 0.0}, {0.0, 4.5, -6.25, 0.0}, {-7.125, 8.0, 9.5, 0.0}};
  const std::vector<float4> expectedF = {{1.25f, -2.5f, 3.75f, 0.0f},
                                         {0.0f, 4.5f, -6.25f, 0.0f},
                                         {-7.125f, 8.0f, 9.5f, 0.0f}};

  Coordinates coordinates(input);

  apo_test::CheckVectorsClose<double4>(coordinates.getCoordinatesD(),
                                             expectedD, DOUBLE_TOLERANCE);
  apo_test::CheckVectorsClose<float4>(coordinates.getCoordinatesF(),
                                            expectedF, FLOAT_TOLERANCE);
}

TEST_CASE("CoordinatesFromFloat3Vector") {
  const std::vector<float3> input = {
      {1.25f, -2.5f, 3.75f}, {0.0f, 4.5f, -6.25f}, {-7.125f, 8.0f, 9.5f}};
  const std::vector<double4> expectedD = {
      {1.25, -2.5, 3.75, 0.0}, {0.0, 4.5, -6.25, 0.0}, {-7.125, 8.0, 9.5, 0.0}};
  const std::vector<float4> expectedF = {{1.25f, -2.5f, 3.75f, 0.0f},
                                         {0.0f, 4.5f, -6.25f, 0.0f},
                                         {-7.125f, 8.0f, 9.5f, 0.0f}};

  Coordinates coordinates(input);

  apo_test::CheckVectorsClose<double4>(coordinates.getCoordinatesD(),
                                             expectedD, DOUBLE_TOLERANCE);
  apo_test::CheckVectorsClose<float4>(coordinates.getCoordinatesF(),
                                            expectedF, FLOAT_TOLERANCE);
}

TEST_CASE("CoordinatesFromDoubleNestedVector") {
  const std::vector<std::vector<double>> input = {
      {1.25, -2.5, 3.75}, {0.0, 4.5, -6.25}, {-7.125, 8.0, 9.5}};
  const std::vector<double4> expectedD = {
      {1.25, -2.5, 3.75, 0.0}, {0.0, 4.5, -6.25, 0.0}, {-7.125, 8.0, 9.5, 0.0}};
  const std::vector<float4> expectedF = {{1.25f, -2.5f, 3.75f, 0.0f},
                                         {0.0f, 4.5f, -6.25f, 0.0f},
                                         {-7.125f, 8.0f, 9.5f, 0.0f}};

  Coordinates coordinates(input);

  apo_test::CheckVectorsClose<double4>(coordinates.getCoordinatesD(),
                                             expectedD, DOUBLE_TOLERANCE);
  apo_test::CheckVectorsClose<float4>(coordinates.getCoordinatesF(),
                                            expectedF, FLOAT_TOLERANCE);
}

TEST_CASE("CoordinatesFromFloatNestedVector") {
  const std::vector<std::vector<float>> input = {
      {1.25f, -2.5f, 3.75f}, {0.0f, 4.5f, -6.25f}, {-7.125f, 8.0f, 9.5f}};
  const std::vector<double4> expectedD = {
      {1.25, -2.5, 3.75, 0.0}, {0.0, 4.5, -6.25, 0.0}, {-7.125, 8.0, 9.5, 0.0}};
  const std::vector<float4> expectedF = {{1.25f, -2.5f, 3.75f, 0.0f},
                                         {0.0f, 4.5f, -6.25f, 0.0f},
                                         {-7.125f, 8.0f, 9.5f, 0.0f}};

  Coordinates coordinates(input);

  apo_test::CheckVectorsClose<double4>(coordinates.getCoordinatesD(),
                                             expectedD, DOUBLE_TOLERANCE);
  apo_test::CheckVectorsClose<float4>(coordinates.getCoordinatesF(),
                                            expectedF, FLOAT_TOLERANCE);
}

TEST_CASE("CoordinatesEmptyInputVectors") {
  SECTION("EmptyDouble3Vector") {
    const std::vector<double3> input;
    Coordinates coordinates(input);

    CHECK(coordinates.getNumAtoms() == 0);
    CHECK(coordinates.getCoordinatesD().empty() == true);
    CHECK(coordinates.getCoordinatesF().empty() == true);
  }

  SECTION("EmptyFloat3Vector") {
    const std::vector<float3> input;
    Coordinates coordinates(input);

    CHECK(coordinates.getNumAtoms() == 0);
    CHECK(coordinates.getCoordinatesD().empty() == true);
    CHECK(coordinates.getCoordinatesF().empty() == true);
  }

  SECTION("EmptyDoubleNestedVector") {
    const std::vector<std::vector<double>> input;
    Coordinates coordinates(input);

    CHECK(coordinates.getNumAtoms() == 0);
    CHECK(coordinates.getCoordinatesD().empty() == true);
    CHECK(coordinates.getCoordinatesF().empty() == true);
  }

  SECTION("EmptyFloatNestedVector") {
    const std::vector<std::vector<float>> input;
    Coordinates coordinates(input);

    CHECK(coordinates.getNumAtoms() == 0);
    CHECK(coordinates.getCoordinatesD().empty() == true);
    CHECK(coordinates.getCoordinatesF().empty() == true);
  }
}

TEST_CASE("CoordinatesSetNumAtoms") {
  Coordinates coordinates;

  coordinates.setNumAtoms(3);

  CHECK(coordinates.getNumAtoms() == 3);
  CHECK(coordinates.getCoordinatesD().size() == 3);
  CHECK(coordinates.getCoordinatesF().size() == 3);

  coordinates.setNumAtoms(1);

  CHECK(coordinates.getNumAtoms() == 1);
  CHECK(coordinates.getCoordinatesD().size() == 1);
  CHECK(coordinates.getCoordinatesF().size() == 1);

  coordinates.setNumAtoms(0);

  CHECK(coordinates.getNumAtoms() == 0);
  CHECK(coordinates.getCoordinatesD().empty() == true);
  CHECK(coordinates.getCoordinatesF().empty() == true);
}

TEST_CASE("CoordinatesMutableGetters") {
  Coordinates coordinates;
  coordinates.setNumAtoms(2);

  std::vector<double4> &coordinatesD = coordinates.getCoordinatesD();
  std::vector<float4> &coordinatesF = coordinates.getCoordinatesF();

  coordinatesD[0] = {1.0, 2.0, 3.0, 4.0};
  coordinatesD[1] = {-1.0, -2.0, -3.0, -4.0};
  coordinatesF[0] = {5.0f, 6.0f, 7.0f, 8.0f};
  coordinatesF[1] = {-5.0f, -6.0f, -7.0f, -8.0f};

  const Coordinates &constCoordinates = coordinates;
  const std::vector<double4> &constCoordinatesD =
      constCoordinates.getCoordinatesD();
  const std::vector<float4> &constCoordinatesF =
      constCoordinates.getCoordinatesF();

  REQUIRE(constCoordinatesD.size() == 2);
  REQUIRE(constCoordinatesF.size() == 2);

  CHECK(constCoordinatesD[0].x == Approx(1.0));
  CHECK(constCoordinatesD[0].y == Approx(2.0));
  CHECK(constCoordinatesD[0].z == Approx(3.0));
  CHECK(constCoordinatesD[0].w == Approx(4.0));
  CHECK(constCoordinatesD[1].x == Approx(-1.0));
  CHECK(constCoordinatesD[1].y == Approx(-2.0));
  CHECK(constCoordinatesD[1].z == Approx(-3.0));
  CHECK(constCoordinatesD[1].w == Approx(-4.0));

  CHECK(constCoordinatesF[0].x == Approx(5.0f));
  CHECK(constCoordinatesF[0].y == Approx(6.0f));
  CHECK(constCoordinatesF[0].z == Approx(7.0f));
  CHECK(constCoordinatesF[0].w == Approx(8.0f));
  CHECK(constCoordinatesF[1].x == Approx(-5.0f));
  CHECK(constCoordinatesF[1].y == Approx(-6.0f));
  CHECK(constCoordinatesF[1].z == Approx(-7.0f));
  CHECK(constCoordinatesF[1].w == Approx(-8.0f));
}

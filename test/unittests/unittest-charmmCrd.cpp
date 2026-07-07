// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "CharmmCrd.h"
#include "apo_test_helpers.h"
#include "catch.hpp"
#include "test_paths.h"

#include <fstream>
#include <iomanip>

namespace {

constexpr double DOUBLE_TOLERANCE = 1.0e-12;
constexpr double FLOAT_TOLERANCE = 1.0e-6;

void WriteStandardCharmmCrd(const std::string &fileName,
                            const std::vector<double3> &coords) {
  std::ofstream file(fileName);
  REQUIRE(file.good());

  file << "* generated standard CharmmCrd unit test\n";
  file << "*\n";
  file << std::setw(10) << coords.size() << "\n";

  for (std::size_t i = 0; i < coords.size(); i++) {
    const int atomNumber = static_cast<int>(i + 1);
    const int residueNumber = (i < 2) ? 1 : 2;
    const std::string residueName = (residueNumber == 1) ? "RS1" : "RS2";
    const std::string atomName = "A" + std::to_string(atomNumber);

    file << std::setw(5) << atomNumber << std::setw(5) << residueNumber << " "
         << std::left << std::setw(4) << residueName << " " << std::setw(4)
         << atomName << std::right << std::fixed << std::setprecision(5)
         << std::setw(10) << coords[i].x << std::setw(10) << coords[i].y
         << std::setw(10) << coords[i].z << " SEG1" << std::setw(5)
         << residueNumber << std::setw(10) << 0.0 << "\n";
  }

  file.close();
  REQUIRE(file.good());

  return;
}

void WriteExtendedCharmmCrd(const std::string &fileName,
                            const std::vector<double3> &coords) {
  std::ofstream file(fileName);
  REQUIRE(file.good());

  file << "* generated extended CharmmCrd unit test\n";
  file << "*\n";
  file << std::setw(10) << coords.size() << " EXT\n";

  for (std::size_t i = 0; i < coords.size(); i++) {
    const int atomNumber = static_cast<int>(i + 1);
    const int residueNumber = (i < 2) ? 1 : 2;
    const std::string residueName = (residueNumber == 1) ? "RES1" : "RES2";
    const std::string atomName = "ATOM" + std::to_string(atomNumber);

    file << std::setw(10) << atomNumber << std::setw(10) << residueNumber << " "
         << std::left << std::setw(8) << residueName << " " << std::setw(8)
         << atomName << "  " << std::right << std::fixed
         << std::setprecision(10) << std::setw(20) << coords[i].x
         << std::setw(20) << coords[i].y << std::setw(20) << coords[i].z
         << " SEG1" << std::setw(10) << residueNumber << std::setw(20) << 0.0
         << "\n";
  }

  file.close();
  REQUIRE(file.good());

  return;
}

void WriteMalformedCharmmCrd(const std::string &fileName) {
  std::ofstream file(fileName);
  REQUIRE(file.good());

  file << "* generated malformed CharmmCrd unit test\n";
  file << "*\n";
  file << std::setw(10) << 1 << "\n";
  file << "short\n";

  file.close();
  REQUIRE(file.good());

  return;
}

} // namespace

TEST_CASE("CharmmCrdParsesStandardFormat") {
  const std::string fileName = "tmp_charmm_crd_standard.cor";
  const std::vector<double3> expectedDP = {
      {1.25, -2.5, 3.75}, {0.0, 4.5, -6.25}, {-7.125, 8.0, 9.5}};
  const std::vector<float3> expectedSP = {
      {1.25f, -2.5f, 3.75f}, {0.0f, 4.5f, -6.25f}, {-7.125f, 8.0f, 9.5f}};

  WriteStandardCharmmCrd(fileName, expectedDP);

  CharmmCrd crd(fileName);

  apo_test::CheckVectorsClose<double3>(crd.getCoordinatesDP(), expectedDP,
                                       DOUBLE_TOLERANCE);
  apo_test::CheckVectorsClose<float3>(crd.getCoordinatesSP(), expectedSP,
                                      FLOAT_TOLERANCE);

  apo_test::RemoveIfExists(fileName);
}

TEST_CASE("CharmmCrdParsesExtendedFormat") {
  const std::string fileName = "tmp_charmm_crd_extended.cor";
  const std::vector<double3> expectedDP = {{-10.123456789, 20.25, -30.5},
                                           {1000.0, -2000.5, 3000.75},
                                           {0.125, 0.25, 0.5}};
  const std::vector<float3> expectedSP = {{-10.123456789f, 20.25f, -30.5f},
                                          {1000.0f, -2000.5f, 3000.75f},
                                          {0.125f, 0.25f, 0.5f}};

  WriteExtendedCharmmCrd(fileName, expectedDP);

  CharmmCrd crd(fileName);

  apo_test::CheckVectorsClose<double3>(crd.getCoordinatesDP(), expectedDP,
                                       DOUBLE_TOLERANCE);
  apo_test::CheckVectorsClose<float3>(crd.getCoordinatesSP(), expectedSP,
                                      FLOAT_TOLERANCE);

  apo_test::RemoveIfExists(fileName);
}

TEST_CASE("CharmmCrdParsesRepositoryNaclPair") {
  const std::string fileName = getDataPath() + "nacl_pair.cor";
  const std::vector<double3> expectedDP = {{0.0, 0.0, 0.0}, {2.82, 2.82, 2.82}};
  const std::vector<float3> expectedSP = {{0.0f, 0.0f, 0.0f},
                                          {2.82f, 2.82f, 2.82f}};

  CharmmCrd crd(fileName);

  apo_test::CheckVectorsClose<double3>(crd.getCoordinatesDP(), expectedDP,
                                       DOUBLE_TOLERANCE);
  apo_test::CheckVectorsClose<float3>(crd.getCoordinatesSP(), expectedSP,
                                      FLOAT_TOLERANCE);
}

TEST_CASE("CharmmCrdRejectsMissingFile") {
  const std::string fileName = "tmp_charmm_crd_missing.cor";

  apo_test::RemoveIfExists(fileName);

  CHECK_THROWS_AS(CharmmCrd(fileName), std::runtime_error);
}

TEST_CASE("CharmmCrdRejectsMalformedFile") {
  const std::string fileName = "tmp_charmm_crd_malformed.cor";

  WriteMalformedCharmmCrd(fileName);

  CHECK_THROWS(CharmmCrd(fileName));

  apo_test::RemoveIfExists(fileName);
}

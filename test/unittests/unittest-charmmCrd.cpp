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

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

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

std::string MakeStandardCharmmCrdRecord(const std::string &x,
                                        const std::string &y,
                                        const std::string &z) {
  std::ostringstream record;
  record << std::setw(5) << 1 << std::setw(5) << 1 << " " << std::left
         << std::setw(4) << "RS1" << " " << std::setw(4) << "A1" << std::right
         << std::setw(10) << x << std::setw(10) << y << std::setw(10) << z
         << "\n";
  return record.str();
}

std::string MakeExtendedCharmmCrdRecord(const std::string &x,
                                        const std::string &y,
                                        const std::string &z) {
  std::ostringstream record;
  record << std::setw(10) << 1 << std::setw(10) << 1 << " " << std::left
         << std::setw(8) << "RES1" << " " << std::setw(8) << "ATOM1" << "  "
         << std::right << std::setw(20) << x << std::setw(20) << y
         << std::setw(20) << z << "\n";
  return record.str();
}

std::string MakeCharmmCrdText(const std::string &atomCountLine,
                              const std::string &coordinateRecords) {
  return "* generated invalid CharmmCrd unit test\n*\n" + atomCountLine + "\n" +
         coordinateRecords;
}

void CheckCharmmCrdFileError(const std::string &fileName,
                             const std::string &contents,
                             const std::string_view expectedMessage) {
  apo_test::WriteTextFile(fileName, contents);

  apo_test::CheckApoCharmmError(
      [&fileName](void) {
        CharmmCrd crd(fileName);
        static_cast<void>(crd);
      },
      ApoCharmmErrorCode::Runtime, expectedMessage);

  apo_test::RemoveIfExists(fileName);
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
  const std::vector<double3> expectedDP = {{0.0, 0.0, 0.0}, {2.82, 2.82, 2.82}};
  const std::vector<float3> expectedSP = {{0.0f, 0.0f, 0.0f},
                                          {2.82f, 2.82f, 2.82f}};

  CharmmCrd crd(apo_test::GetDataDir() / "nacl_pair.cor");

  apo_test::CheckVectorsClose<double3>(crd.getCoordinatesDP(), expectedDP,
                                       DOUBLE_TOLERANCE);
  apo_test::CheckVectorsClose<float3>(crd.getCoordinatesSP(), expectedSP,
                                      FLOAT_TOLERANCE);
}

TEST_CASE("CharmmCrdValidationUsesApoCharmmError") {
  SECTION("EmptyPath") {
    apo_test::CheckApoCharmmError(
        [](void) {
          CharmmCrd crd("");
          static_cast<void>(crd);
        },
        ApoCharmmErrorCode::InvalidArgument,
        "CHARMM coordinate file path must not be empty");
  }

  SECTION("MissingFile") {
    const std::string fileName = "tmp_charmm_crd_missing.cor";
    apo_test::RemoveIfExists(fileName);

    apo_test::CheckApoCharmmError(
        [&fileName](void) {
          CharmmCrd crd(fileName);
          static_cast<void>(crd);
        },
        ApoCharmmErrorCode::Runtime,
        "Failed to open file \"" + fileName + "\"");
  }

  SECTION("MissingAtomCountAtEndOfFile") {
    const std::string fileName = "tmp_charmm_crd_missing_atom_count.cor";

    CheckCharmmCrdFileError(
        fileName, "* generated title-only CharmmCrd unit test\n",
        "Unexpected end of file while reading atom count in CHARMM coordinate "
        "file \"" +
            fileName + "\"");
  }

  SECTION("BlankAtomCount") {
    const std::string fileName = "tmp_charmm_crd_blank_atom_count.cor";

    CheckCharmmCrdFileError(fileName, MakeCharmmCrdText("", ""),
                            "Missing atom count in CHARMM coordinate file \"" +
                                fileName + "\" at line 3");
  }

  SECTION("InvalidAtomCount") {
    const std::string fileName = "tmp_charmm_crd_invalid_atom_count.cor";

    CheckCharmmCrdFileError(
        fileName, MakeCharmmCrdText("BAD", ""),
        "Invalid atom count value \"BAD\" in CHARMM coordinate file \"" +
            fileName + "\" at line 3");
  }

  SECTION("AtomCountWithTrailingCharacters") {
    const std::string fileName = "tmp_charmm_crd_trailing_atom_count.cor";

    CheckCharmmCrdFileError(
        fileName, MakeCharmmCrdText("1BAD", ""),
        "Invalid atom count value \"1BAD\" in CHARMM coordinate file \"" +
            fileName + "\" at line 3");
  }

  SECTION("NegativeAtomCount") {
    const std::string fileName = "tmp_charmm_crd_negative_atom_count.cor";

    CheckCharmmCrdFileError(
        fileName, MakeCharmmCrdText("-1", ""),
        "Invalid atom count value \"-1\" in CHARMM coordinate file \"" +
            fileName + "\" at line 3");
  }

  SECTION("OutOfRangeAtomCount") {
    const std::string fileName = "tmp_charmm_crd_out_of_range_atom_count.cor";
    const std::string atomCount =
        std::to_string(std::numeric_limits<unsigned long long int>::max()) +
        "0";

    CheckCharmmCrdFileError(fileName, MakeCharmmCrdText(atomCount, ""),
                            "Invalid atom count value \"" + atomCount +
                                "\" in CHARMM coordinate file \"" + fileName +
                                "\" at line 3");
  }

  SECTION("AtomCountExceedsSupportedRange") {
    const std::string fileName = "tmp_charmm_crd_large_atom_count.cor";
    const unsigned long long int unsupportedCount =
        static_cast<unsigned long long int>(std::numeric_limits<int>::max()) +
        1ULL;

    CheckCharmmCrdFileError(
        fileName, MakeCharmmCrdText(std::to_string(unsupportedCount), ""),
        "Atom count exceeds supported range in CHARMM coordinate file \"" +
            fileName + "\" at line 3");
  }

  SECTION("MissingCoordinateRecord") {
    const std::string fileName = "tmp_charmm_crd_missing_record.cor";

    CheckCharmmCrdFileError(
        fileName, MakeCharmmCrdText("1", ""),
        "Unexpected end of file while reading coordinate record 1 in CHARMM "
        "coordinate file \"" +
            fileName + "\"");
  }

  SECTION("TruncatedCoordinateRecord") {
    const std::string fileName = "tmp_charmm_crd_truncated_record.cor";

    CheckCharmmCrdFileError(
        fileName, MakeCharmmCrdText("1", "short\n"),
        "Coordinate record 1 is truncated in CHARMM coordinate file \"" +
            fileName + "\" at line 4");
  }

  SECTION("TruncatedExtendedCoordinateRecord") {
    const std::string fileName = "tmp_charmm_crd_truncated_extended_record.cor";

    CheckCharmmCrdFileError(
        fileName, MakeCharmmCrdText("1 EXT", "short\n"),
        "Coordinate record 1 is truncated in CHARMM coordinate file \"" +
            fileName + "\" at line 4");
  }

  SECTION("InvalidXCoordinate") {
    const std::string fileName = "tmp_charmm_crd_invalid_x.cor";

    CheckCharmmCrdFileError(
        fileName,
        MakeCharmmCrdText("1",
                          MakeStandardCharmmCrdRecord("BAD", "2.0", "3.0")),
        "Invalid X coordinate value \"BAD\" in coordinate record 1 of CHARMM "
        "coordinate file \"" +
            fileName + "\" at line 4");
  }

  SECTION("InvalidYCoordinate") {
    const std::string fileName = "tmp_charmm_crd_invalid_y.cor";

    CheckCharmmCrdFileError(
        fileName,
        MakeCharmmCrdText("1",
                          MakeStandardCharmmCrdRecord("1.0", "BAD", "3.0")),
        "Invalid Y coordinate value \"BAD\" in coordinate record 1 of CHARMM "
        "coordinate file \"" +
            fileName + "\" at line 4");
  }

  SECTION("InvalidZCoordinate") {
    const std::string fileName = "tmp_charmm_crd_invalid_z.cor";

    CheckCharmmCrdFileError(
        fileName,
        MakeCharmmCrdText("1",
                          MakeStandardCharmmCrdRecord("1.0", "2.0", "BAD")),
        "Invalid Z coordinate value \"BAD\" in coordinate record 1 of CHARMM "
        "coordinate file \"" +
            fileName + "\" at line 4");
  }

  SECTION("InvalidExtendedXCoordinate") {
    const std::string fileName = "tmp_charmm_crd_invalid_extended_x.cor";

    CheckCharmmCrdFileError(
        fileName,
        MakeCharmmCrdText("1 EXT",
                          MakeExtendedCharmmCrdRecord("BAD", "2.0", "3.0")),
        "Invalid X coordinate value \"BAD\" in coordinate record 1 of "
        "CHARMM coordinate file \"" +
            fileName + "\" at line 4");
  }

  SECTION("InvalidExtendedYCoordinate") {
    const std::string fileName = "tmp_charmm_crd_invalid_extended_y.cor";

    CheckCharmmCrdFileError(
        fileName,
        MakeCharmmCrdText("1 EXT",
                          MakeExtendedCharmmCrdRecord("1.0", "BAD", "3.0")),
        "Invalid Y coordinate value \"BAD\" in coordinate record 1 of "
        "CHARMM coordinate file \"" +
            fileName + "\" at line 4");
  }

  SECTION("InvalidExtendedZCoordinate") {
    const std::string fileName = "tmp_charmm_crd_invalid_extended_z.cor";

    CheckCharmmCrdFileError(
        fileName,
        MakeCharmmCrdText("1 EXT",
                          MakeExtendedCharmmCrdRecord("1.0", "2.0", "BAD")),
        "Invalid Z coordinate value \"BAD\" in coordinate record 1 of "
        "CHARMM coordinate file \"" +
            fileName + "\" at line 4");
  }

  SECTION("CoordinateWithTrailingCharacters") {
    const std::string fileName = "tmp_charmm_crd_trailing_coordinate.cor";

    CheckCharmmCrdFileError(
        fileName,
        MakeCharmmCrdText("1",
                          MakeStandardCharmmCrdRecord("1.0BAD", "2.0", "3.0")),
        "Invalid X coordinate value \"1.0BAD\" in coordinate record 1 of "
        "CHARMM coordinate file \"" +
            fileName + "\" at line 4");
  }

  SECTION("OutOfRangeCoordinate") {
    const std::string fileName = "tmp_charmm_crd_out_of_range_coordinate.cor";

    CheckCharmmCrdFileError(
        fileName,
        MakeCharmmCrdText("1",
                          MakeStandardCharmmCrdRecord("1E9999", "2.0", "3.0")),
        "Invalid X coordinate value \"1E9999\" in coordinate record 1 of "
        "CHARMM coordinate file \"" +
            fileName + "\" at line 4");
  }

  SECTION("NonfiniteCoordinate") {
    const std::string fileName = "tmp_charmm_crd_nonfinite_coordinate.cor";

    CheckCharmmCrdFileError(
        fileName,
        MakeCharmmCrdText("1",
                          MakeStandardCharmmCrdRecord("NAN", "2.0", "3.0")),
        "Invalid X coordinate value \"NAN\" in coordinate record 1 of CHARMM "
        "coordinate file \"" +
            fileName + "\" at line 4");
  }
}

// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: Samarjeet Prasad, James E. Gonzales II
//
// ENDLICENSE

#include "CharmmCrd.h"

#include "ApoCharmmError.h"
#include "str_utils.h"

#include <array>
#include <cstddef>
#include <filesystem>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

CharmmCrd::CharmmCrd(const std::filesystem::path &filePath) : Coordinates() {
  APOCHARMM_REQUIRE(!filePath.empty(), ApoCharmmErrorCode::InvalidArgument,
                    "CHARMM coordinate file path must not be empty");

  this->readCharmmCrdFile(filePath);
}

void CharmmCrd::readCharmmCrdFile(const std::filesystem::path &filePath) {
  std::string fileData = "";
  apo::read_file_into_string(fileData, filePath);

  const std::string sourceName =
      "CHARMM coordinate file \"" + filePath.string() + "\"";

  std::size_t pos = 0;
  std::size_t lineNumber = 0;
  std::string line = "";
  std::vector<std::string> tokens;

  // Parse TITLE
  do {
    apo::get_line(line, pos, lineNumber, fileData, "atom count", sourceName);
  } while (!line.empty() && (line.front() == '*'));

  tokens = apo::split(line);

  APOCHARMM_REQUIRE(!tokens.empty() && !tokens.front().empty(),
                    ApoCharmmErrorCode::Runtime,
                    "Missing atom count in " + sourceName + " at line " +
                        std::to_string(lineNumber));

  const std::string atomCountContext =
      sourceName + " at line " + std::to_string(lineNumber);
  const unsigned long long int natom =
      apo::parse_ull(tokens.front(), "atom count", atomCountContext);

  APOCHARMM_REQUIRE(natom <= static_cast<unsigned long long int>(
                                 std::numeric_limits<int>::max()),
                    ApoCharmmErrorCode::Runtime,
                    "Atom count exceeds supported range in " + sourceName +
                        " at line " + std::to_string(lineNumber));

  const bool isExt = ((tokens.size() >= 2) && (tokens[1] == "EXT"));
  const std::array<std::size_t, 3> fieldOffsets =
      isExt ? std::array<std::size_t, 3>{40, 60, 80}
            : std::array<std::size_t, 3>{20, 30, 40};
  const std::size_t fieldWidth = isExt ? 20 : 10;
  const std::array<std::string_view, 3> axisNames = {"X", "Y", "Z"};

  this->setNumAtoms(static_cast<int>(natom));

  for (unsigned long long int i = 0; i < natom; i++) {
    const std::size_t recordNumber = static_cast<std::size_t>(i + 1);

    apo::get_line(line, pos, lineNumber, fileData,
                  "coordinate record " + std::to_string(recordNumber),
                  sourceName);

    const std::string recordContext =
        "coordinate record " + std::to_string(recordNumber) + " of " +
        sourceName + " at line " + std::to_string(lineNumber);
    const std::string truncatedMessage =
        "Coordinate record " + std::to_string(recordNumber) +
        " is truncated in " + sourceName + " at line " +
        std::to_string(lineNumber);

    std::array<double, 3> xyz = {0.0, 0.0, 0.0};

    for (std::size_t axis = 0; axis < axisNames.size(); axis++) {
      std::string_view field;

      APOCHARMM_REQUIRE(apo::try_get_fixed_width_field(
                            field, line, fieldOffsets[axis], fieldWidth),
                        ApoCharmmErrorCode::Runtime, truncatedMessage);

      xyz[axis] = apo::parse_double(
          field, std::string(axisNames[axis]) + " coordinate", recordContext);
    }

    m_CoordinatesDP[i] = make_double3(xyz[0], xyz[1], xyz[2]);
    m_CoordinatesSP[i] =
        make_float3(static_cast<float>(xyz[0]), static_cast<float>(xyz[1]),
                    static_cast<float>(xyz[2]));
  }

  return;
}

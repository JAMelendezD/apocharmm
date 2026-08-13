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

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

void ReadCharmmCrdLine(std::string &line, std::size_t &position,
                       std::size_t &lineNumber, const std::string_view fileData,
                       const std::string_view recordName,
                       const std::string &fileName) {
  APOCHARMM_REQUIRE(position < fileData.size(), ApoCharmmErrorCode::Runtime,
                    "Unexpected end of file while reading " +
                        std::string(recordName) +
                        " in CHARMM coordinate file \"" + fileName + "\"");

  const std::size_t lineEnd = fileData.find('\n', position);
  if (lineEnd == std::string::npos) {
    line.assign(fileData.substr(position));
    position = fileData.size();
  } else {
    line.assign(fileData.substr(position, lineEnd - position));
    position = lineEnd + 1;
  }

  lineNumber++;

  return;
}

double ParseCharmmCrdCoordinate(
    const std::string &line, const std::size_t offset, const std::size_t width,
    const std::string_view axisName, const std::size_t recordNumber,
    const std::string &fileName, const std::size_t lineNumber) {
  APOCHARMM_REQUIRE(line.size() >= offset + width, ApoCharmmErrorCode::Runtime,
                    "Coordinate record " + std::to_string(recordNumber) +
                        " is truncated in CHARMM coordinate file \"" +
                        fileName + "\" at line " + std::to_string(lineNumber));

  const std::string token = apo::trim(line.substr(offset, width));
  std::size_t parsedCharacters = 0;
  double value = 0.0;
  bool conversionSucceeded = false;

  try {
    value = std::stod(token, &parsedCharacters);
    conversionSucceeded = true;
  } catch (const std::invalid_argument &) {
  } catch (const std::out_of_range &) {
  }

  const bool isValid = conversionSucceeded &&
                       (parsedCharacters == token.size()) &&
                       std::isfinite(value);

  APOCHARMM_REQUIRE(
      isValid, ApoCharmmErrorCode::Runtime,
      "Invalid " + std::string(axisName) + " coordinate value \"" + token +
          "\" in coordinate record " + std::to_string(recordNumber) +
          " of CHARMM coordinate file \"" + fileName + "\" at line " +
          std::to_string(lineNumber));

  return value;
}

} // namespace

CharmmCrd::CharmmCrd(const std::string &fileName) : Coordinates() {
  APOCHARMM_REQUIRE(!fileName.empty(), ApoCharmmErrorCode::InvalidArgument,
                    "CHARMM coordinate file path must not be empty");

  this->readCharmmCrdFile(fileName);
}

void CharmmCrd::readCharmmCrdFile(const std::string &fileName) {
  std::string fileData = "";
  apo::read_file_into_string(fileData, fileName);

  std::size_t pos = 0;
  std::size_t lineNumber = 0;
  std::string line = "";
  std::vector<std::string> tokens;

  // Parse TITLE
  do {
    ReadCharmmCrdLine(line, pos, lineNumber, fileData, "atom count", fileName);
  } while (!line.empty() && (line.front() == '*'));

  tokens = apo::split(line);
  APOCHARMM_REQUIRE(!tokens.empty() && !tokens.front().empty(),
                    ApoCharmmErrorCode::Runtime,
                    "Missing atom count in CHARMM coordinate file \"" +
                        fileName + "\" at line " + std::to_string(lineNumber));

  // Read coordinate data
  const std::string &atomCountToken = tokens.front();
  std::size_t parsedCharacters = 0;
  unsigned long long int natom = 0;
  bool conversionSucceeded = false;

  try {
    natom = std::stoull(atomCountToken, &parsedCharacters);
    conversionSucceeded = true;
  } catch (const std::invalid_argument &) {
  } catch (const std::out_of_range &) {
  }

  const bool atomCountIsValid = conversionSucceeded &&
                                (parsedCharacters == atomCountToken.size()) &&
                                (atomCountToken.front() != '-');

  APOCHARMM_REQUIRE(atomCountIsValid, ApoCharmmErrorCode::Runtime,
                    "Invalid atom count value \"" + atomCountToken +
                        "\" in CHARMM coordinate file \"" + fileName +
                        "\" at line " + std::to_string(lineNumber));

  APOCHARMM_REQUIRE(
      natom <=
          static_cast<unsigned long long int>(std::numeric_limits<int>::max()),
      ApoCharmmErrorCode::Runtime,
      "Atom count exceeds supported range in CHARMM coordinate file \"" +
          fileName + "\" at line " + std::to_string(lineNumber));

  const bool isExt = ((tokens.size() >= 2) && (tokens[1] == "EXT"));
  this->setNumAtoms(static_cast<int>(natom));
  for (unsigned long long int i = 0; i < natom; i++) {
    ReadCharmmCrdLine(line, pos, lineNumber, fileData,
                      "coordinate record " + std::to_string(i + 1), fileName);
    double x = -9999.9999, y = -9999.9999, z = -9999.9999;
    if (isExt) {
      x = ParseCharmmCrdCoordinate(line, 40, 20, "X", i + 1, fileName,
                                   lineNumber);
      y = ParseCharmmCrdCoordinate(line, 60, 20, "Y", i + 1, fileName,
                                   lineNumber);
      z = ParseCharmmCrdCoordinate(line, 80, 20, "Z", i + 1, fileName,
                                   lineNumber);
    } else {
      x = ParseCharmmCrdCoordinate(line, 20, 10, "X", i + 1, fileName,
                                   lineNumber);
      y = ParseCharmmCrdCoordinate(line, 30, 10, "Y", i + 1, fileName,
                                   lineNumber);
      z = ParseCharmmCrdCoordinate(line, 40, 10, "Z", i + 1, fileName,
                                   lineNumber);
    }
    m_CoordinatesDP[i] = make_double3(x, y, z);
    m_CoordinatesSP[i] = make_float3(
        static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));
  }

  return;
}

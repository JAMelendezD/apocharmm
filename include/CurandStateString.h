// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#pragma once

#include "ApoCharmmError.h"
#include "str_utils.h"

#include <curand_kernel.h>

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace apo {

// struct curand_state_string {
//   unsigned long long int pos = 0;
//   std::vector<curandStatePhilox4_32_10_t> states;
// };

inline char to_hex_digit(const unsigned char value) {
  return static_cast<char>((value < 10) ? ('0' + value) : ('A' + value - 10));
}

inline unsigned char from_hex_digit(const char value) {
  if ((value >= '0') && (value <= '9'))
    return static_cast<unsigned char>(value - '0');
  if ((value >= 'A') && (value <= 'F'))
    return static_cast<unsigned char>(value - 'A' + 10);
  if ((value >= 'a') && (value <= 'f'))
    return static_cast<unsigned char>(value - 'a' + 10);

  APOCHARMM_THROW(ApoCharmmErrorCode::InvalidArgument,
                  "Invalid hex digit in curand state string");
}

inline std::string
curand_states_to_string(const unsigned long long int pos,
                        const std::vector<curandStatePhilox4_32_10_t> &states) {
  std::string encoded = "APO_PHILOX_V1:";
  encoded += std::to_string(pos);
  encoded += ":";
  encoded += std::to_string(states.size());
  encoded += ":";
  encoded += std::to_string(sizeof(curandStatePhilox4_32_10_t));
  encoded += ":";

  const unsigned char *bytes =
      reinterpret_cast<const unsigned char *>(states.data());
  const std::size_t numBytes =
      states.size() * sizeof(curandStatePhilox4_32_10_t);

  encoded.reserve(encoded.size() + 2 * numBytes);

  for (std::size_t i = 0; i < numBytes; i++) {
    encoded.push_back(to_hex_digit(static_cast<unsigned char>(bytes[i] >> 4)));
    encoded.push_back(
        to_hex_digit(static_cast<unsigned char>(bytes[i] & 0x0F)));
  }

  return encoded;
}

inline void
curand_states_from_string(unsigned long long int &pos,
                          std::vector<curandStatePhilox4_32_10_t> &states,
                          const std::string &value) {
  constexpr const char *PREFIX = "APO_PHILOX_V1:";
  const std::string trimmed = apo::trim(value);
  unsigned long long int parsedPosition = 0;

  APOCHARMM_REQUIRE(trimmed.rfind(PREFIX, 0) == 0,
                    ApoCharmmErrorCode::InvalidArgument,
                    "Curand state string has invalid prefix");

  std::size_t begin = std::string(PREFIX).size();
  std::size_t end = trimmed.find(':', begin);
  APOCHARMM_REQUIRE(end != std::string::npos,
                    ApoCharmmErrorCode::InvalidArgument,
                    "Curand state string is missing sequence position");

  APOCHARMM_REQUIRE(
      apo::try_parse_ull(parsedPosition, trimmed.substr(begin, end - begin)),
      ApoCharmmErrorCode::InvalidArgument,
      "Curand state string has invalid sequence position");

  begin = end + 1;
  end = trimmed.find(':', begin);
  APOCHARMM_REQUIRE(end != std::string::npos,
                    ApoCharmmErrorCode::InvalidArgument,
                    "Curand state string is missing state count");

  unsigned long long int numStatesValue = 0;
  APOCHARMM_REQUIRE(
      apo::try_parse_ull(numStatesValue, trimmed.substr(begin, end - begin)),
      ApoCharmmErrorCode::InvalidArgument,
      "Curand state string has invalid state count");

  if constexpr (sizeof(std::size_t) < sizeof(unsigned long long int)) {
    APOCHARMM_REQUIRE(numStatesValue <= std::numeric_limits<std::size_t>::max(),
                      ApoCharmmErrorCode::InvalidArgument,
                      "Curand state string state count is out of range");
  }

  const std::size_t numStates = static_cast<std::size_t>(numStatesValue);

  begin = end + 1;
  end = trimmed.find(':', begin);
  APOCHARMM_REQUIRE(end != std::string::npos,
                    ApoCharmmErrorCode::InvalidArgument,
                    "Curand state string is missing state size");

  unsigned long long int stateSizeValue = 0;
  APOCHARMM_REQUIRE(
      apo::try_parse_ull(stateSizeValue, trimmed.substr(begin, end - begin)),
      ApoCharmmErrorCode::InvalidArgument,
      "Curand state string has invalid state size");

  if constexpr (sizeof(std::size_t) < sizeof(unsigned long long int)) {
    APOCHARMM_REQUIRE(stateSizeValue <= std::numeric_limits<std::size_t>::max(),
                      ApoCharmmErrorCode::InvalidArgument,
                      "Curand state string state size is out of range");
  }

  const std::size_t stateSize = static_cast<std::size_t>(stateSizeValue);

  APOCHARMM_REQUIRE(
      stateSize == sizeof(curandStatePhilox4_32_10_t),
      ApoCharmmErrorCode::Runtime,
      "Curand state string was written with a different curand state size");

  APOCHARMM_REQUIRE(numStates <=
                        std::numeric_limits<std::size_t>::max() / stateSize,
                    ApoCharmmErrorCode::InvalidArgument,
                    "Curand state string byte count is out of range");

  const std::size_t numBytes = numStates * stateSize;

  APOCHARMM_REQUIRE(numBytes <= std::numeric_limits<std::size_t>::max() / 2,
                    ApoCharmmErrorCode::InvalidArgument,
                    "Curand state string hex length is out of range");

  const std::string hex = trimmed.substr(end + 1);
  const std::size_t expectedHexLength = 2 * numBytes;
  APOCHARMM_REQUIRE(hex.size() == expectedHexLength,
                    ApoCharmmErrorCode::InvalidArgument,
                    "Curand state string has invalid hex length");

  std::vector<curandStatePhilox4_32_10_t> parsedStates(numStates);

  unsigned char *bytes = reinterpret_cast<unsigned char *>(parsedStates.data());
  for (std::size_t i = 0; i < numBytes; i++) {
    const unsigned char hi = from_hex_digit(hex[2 * i + 0]);
    const unsigned char lo = from_hex_digit(hex[2 * i + 1]);
    bytes[i] = static_cast<unsigned char>((hi << 4) | lo);
  }

  pos = parsedPosition;
  states.swap(parsedStates);

  return;
}

} // namespace apo

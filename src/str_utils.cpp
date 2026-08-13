// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "str_utils.h"

#include "ApoCharmmError.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace {

template <typename T>
bool try_parse_number(T &value, const std::string &token) {
  static_assert(std::is_same_v<T, int> ||
                    std::is_same_v<T, unsigned long long int> ||
                    std::is_same_v<T, double>,
                "Unsupported numeric value type");

  std::size_t parsed_characters = 0;
  T parsed_value{};

  try {
    if constexpr (std::is_same_v<T, int>)
      parsed_value = std::stoi(token, &parsed_characters);
    else if constexpr (std::is_same_v<T, unsigned long long int>)
      parsed_value = std::stoull(token, &parsed_characters);
    else
      parsed_value = std::stod(token, &parsed_characters);
  } catch (const std::invalid_argument &) {
    return false;
  } catch (const std::out_of_range &) {
    return false;
  }

  if (parsed_characters != token.size())
    return false;

  value = parsed_value;

  return true;
}

std::string make_invalid_value_message(const std::string_view token,
                                       const std::string_view field_name,
                                       const std::string_view context) {
  return "Invalid " + std::string(field_name) + " value \"" +
         std::string(token) + "\" in " + std::string(context);
}

} // namespace

void apo::ltrim_ip(std::string &str) {
  str.erase(str.begin(),
            std::find_if(str.begin(), str.end(),
                         [](unsigned char c) { return !std::isspace(c); }));
  return;
}

void apo::rtrim_ip(std::string &str) {
  str.erase(std::find_if(str.rbegin(), str.rend(),
                         [](unsigned char c) { return !std::isspace(c); })
                .base(),
            str.end());
  return;
}

void apo::trim_ip(std::string &str) {
  apo::rtrim_ip(str);
  apo::ltrim_ip(str);
  return;
}

std::string apo::ltrim(const std::string_view str) {
  std::string s{str};
  apo::ltrim_ip(s);
  return s;
}

std::string apo::rtrim(const std::string_view str) {
  std::string s{str};
  apo::rtrim_ip(s);
  return s;
}

std::string apo::trim(const std::string_view str) {
  std::string s{str};
  apo::trim_ip(s);
  return s;
}

void apo::to_lower_ip(std::string &str) {
  std::transform(
      str.begin(), str.end(), str.begin(),
      [](const unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return;
}

void apo::to_upper_ip(std::string &str) {
  std::transform(
      str.begin(), str.end(), str.begin(),
      [](const unsigned char c) { return static_cast<char>(std::toupper(c)); });
  return;
}

std::string apo::to_lower(const std::string_view str) {
  std::string s{str};
  apo::to_lower_ip(s);
  return s;
}

std::string apo::to_upper(const std::string_view str) {
  std::string s{str};
  apo::to_upper_ip(s);
  return s;
}

std::vector<std::string> apo::split(const std::string_view str,
                                    const std::string_view delimiter) {
  APOCHARMM_REQUIRE(!delimiter.empty(), ApoCharmmErrorCode::InvalidArgument,
                    "Split delimiter must not be empty");

  const std::string trimmed = apo::trim(str);
  std::string_view remaining(trimmed);
  std::vector<std::string> tokens;

  std::size_t delimiter_position = 0;

  while ((delimiter_position = remaining.find(delimiter)) !=
         std::string_view::npos) {
    tokens.emplace_back(remaining.substr(0, delimiter_position));
    remaining.remove_prefix(delimiter_position + delimiter.size());

    while (!remaining.empty() &&
           std::isspace(static_cast<unsigned char>(remaining.front())) != 0) {
      remaining.remove_prefix(1);
    }
  }

  tokens.emplace_back(remaining);

  return tokens;
}

bool apo::try_get_line(std::string &line, std::size_t &pos,
                       const std::string_view file_data) {
  if (pos >= file_data.size()) {
    line.clear();
    pos = file_data.size();
    return false;
  }

  const std::size_t line_end = file_data.find('\n', pos);

  if (line_end == std::string_view::npos) {
    line.assign(file_data.substr(pos));
    pos = file_data.size();
  } else {
    line.assign(file_data.substr(pos, line_end - pos));
    pos = line_end + 1;
  }

  return true;
}

bool apo::try_get_line(std::string &line, std::size_t &pos,
                       std::size_t &line_number,
                       const std::string_view file_data) {
  if (!apo::try_get_line(line, pos, file_data))
    return false;

  line_number++;

  return true;
}

void apo::get_line(std::string &line, std::size_t &pos,
                   const std::string_view file_data) {
  static_cast<void>(apo::try_get_line(line, pos, file_data));
  return;
}

void apo::get_line(std::string &line, std::size_t &pos,
                   std::size_t &line_number, const std::string_view file_data,
                   const std::string_view record_name,
                   const std::string_view source_name) {
  APOCHARMM_REQUIRE(apo::try_get_line(line, pos, line_number, file_data),
                    ApoCharmmErrorCode::Runtime,
                    "Unexpected end of file while reading " +
                        std::string(record_name) + " in " +
                        std::string(source_name));
  return;
}

bool apo::try_get_fixed_width_field(std::string_view &field,
                                    const std::string_view line,
                                    const std::size_t offset,
                                    const std::size_t width) noexcept {
  field = std::string_view{};

  if ((offset > line.size()) || (width > line.size() - offset))
    return false;

  field = line.substr(offset, width);

  return true;
}

void apo::read_file_into_string(std::string &file_data,
                                const std::string &file_name) {
  std::ifstream ifs(file_name, std::ios::in | std::ios::binary);

  APOCHARMM_REQUIRE(ifs.is_open(), ApoCharmmErrorCode::Runtime,
                    "Failed to open file \"" + file_name + "\"");

  ifs.seekg(0, std::ios::end);
  const std::streampos end_position = ifs.tellg();

  APOCHARMM_REQUIRE(end_position != std::streampos(-1),
                    ApoCharmmErrorCode::Runtime,
                    "Failed to determine size of file \"" + file_name + "\"");

  // Store the size of the file
  const std::streamoff file_size = end_position - std::streampos(0);

  APOCHARMM_REQUIRE(file_size >= 0, ApoCharmmErrorCode::Runtime,
                    "File size is invalid for \"" + file_name + "\"");

  const std::uintmax_t unsigned_file_size =
      static_cast<std::uintmax_t>(file_size);

  APOCHARMM_REQUIRE(
      (unsigned_file_size <=
       static_cast<std::uintmax_t>(std::numeric_limits<std::size_t>::max())) &&
          (unsigned_file_size <=
           static_cast<std::uintmax_t>(
               std::numeric_limits<std::streamsize>::max())),
      ApoCharmmErrorCode::Runtime,
      "File is too large to read into memory: \"" + file_name + "\"");

  // Initialize a std::string with length fileSize, filled with null characters
  std::string contents(static_cast<std::size_t>(unsigned_file_size), '\0');

  // Go back to the beginning of the file and store the contents in the string
  ifs.seekg(0, std::ios::beg);

  APOCHARMM_REQUIRE(ifs.good(), ApoCharmmErrorCode::Runtime,
                    "Failed to seek in file \"" + file_name + "\"");

  if (!contents.empty()) {
    const std::streamsize expected_size =
        static_cast<std::streamsize>(contents.size());

    ifs.read(contents.data(), expected_size);

    APOCHARMM_REQUIRE(!ifs.bad() && (ifs.gcount() == expected_size),
                      ApoCharmmErrorCode::Runtime,
                      "Failed while reading file \"" + file_name + "\"");
  }

  file_data.swap(contents);

  return;
}

bool apo::try_parse_int(int &value, const std::string_view token) {
  const std::string normalized_token = apo::trim(token);
  int parsed_value = 0;

  if (!try_parse_number<int>(parsed_value, normalized_token))
    return false;

  value = parsed_value;

  return true;
}

bool apo::try_parse_ull(unsigned long long int &value,
                        const std::string_view token) {
  const std::string normalized_token = apo::trim(token);

  if (normalized_token.empty() || (normalized_token.front() == '-'))
    return false;

  unsigned long long int parsed_value = 0;

  if (!try_parse_number<unsigned long long int>(parsed_value, normalized_token))
    return false;

  value = parsed_value;

  return true;
}

bool apo::try_parse_double(double &value, const std::string_view token) {
  std::string normalized_token = apo::trim(token);

  std::replace(normalized_token.begin(), normalized_token.end(), 'D', 'E');
  std::replace(normalized_token.begin(), normalized_token.end(), 'd', 'E');

  double parsed_value = 0.0;

  if (!try_parse_number<double>(parsed_value, normalized_token))
    return false;

  value = parsed_value;

  return true;
}

bool apo::try_parse_finite_double(double &value, const std::string_view token) {
  double parsed_value = 0.0;

  if (!apo::try_parse_double(parsed_value, token) ||
      !std::isfinite(parsed_value)) {
    return false;
  }

  value = parsed_value;

  return true;
}

int apo::parse_int(const std::string_view token,
                   const std::string_view field_name,
                   const std::string_view context) {
  int value = 0;

  APOCHARMM_REQUIRE(
      apo::try_parse_int(value, token), ApoCharmmErrorCode::Runtime,
      make_invalid_value_message(apo::trim(token), field_name, context));

  return value;
}

unsigned long long int apo::parse_ull(const std::string_view token,
                                      const std::string_view field_name,
                                      const std::string_view context) {
  unsigned long long int value = 0;

  APOCHARMM_REQUIRE(
      apo::try_parse_ull(value, token), ApoCharmmErrorCode::Runtime,
      make_invalid_value_message(apo::trim(token), field_name, context));

  return value;
}

double apo::parse_double(const std::string_view token,
                         const std::string_view field_name,
                         const std::string_view context) {
  double value = 0.0;

  APOCHARMM_REQUIRE(
      apo::try_parse_finite_double(value, token), ApoCharmmErrorCode::Runtime,
      make_invalid_value_message(apo::trim(token), field_name, context));

  return value;
}

bool apo::contains_wildcard(const std::string_view str) {
  return str.find_first_of("*%#+") != std::string_view::npos;
}

std::string apo::cDoubleToFortSciStr(const double val, const int prec) {
  APOCHARMM_REQUIRE(prec >= 0, ApoCharmmErrorCode::InvalidArgument,
                    "Precision must be non-negative; observed " +
                        std::to_string(prec));

  std::ostringstream oss;
  oss.imbue(std::locale::classic());
  oss << std::scientific << std::setprecision(prec) << val;
  std::string str = oss.str();
  std::replace(str.begin(), str.end(), 'e', 'D');
  return str;
}

double apo::fortSciStrToCDouble(const std::string_view str) {
  return apo::parse_double(str, "floating-point", "Fortran scientific value");
}

std::string apo::get_rst_field(const std::string &line,
                               const std::size_t offset,
                               const std::size_t width,
                               const std::string_view field_name,
                               const std::string &rst_name) {
  std::string_view field;

  APOCHARMM_REQUIRE(apo::try_get_fixed_width_field(field, line, offset, width),
                    ApoCharmmErrorCode::Runtime,
                    "Restart field \"" + std::string(field_name) +
                        "\" is truncated in file \"" + rst_name + "\"");

  return std::string(field);
}

double apo::parse_rst_double(const std::string &line, const std::size_t offset,
                             const std::size_t width,
                             const std::string_view field_name,
                             const std::string &rst_name) {
  const std::string field =
      get_rst_field(line, offset, width, field_name, rst_name);

  double value = 0.0;

  APOCHARMM_REQUIRE(apo::try_parse_double(value, field),
                    ApoCharmmErrorCode::Runtime,
                    "Restart field \"" + std::string(field_name) +
                        "\" is not a valid floating-point value in file \"" +
                        rst_name + "\"");

  APOCHARMM_REQUIRE(std::isfinite(value), ApoCharmmErrorCode::Runtime,
                    "Restart field \"" + std::string(field_name) +
                        "\" must be finite in file \"" + rst_name +
                        "\"; observed " + std::to_string(value));

  return value;
}

int apo::parse_rst_int(const std::string &line, const std::size_t offset,
                       const std::size_t width,
                       const std::string_view field_name,
                       const std::string &rst_name) {
  const std::string field =
      get_rst_field(line, offset, width, field_name, rst_name);

  int value = 0;

  APOCHARMM_REQUIRE(
      apo::try_parse_int(value, field), ApoCharmmErrorCode::Runtime,
      "Restart field \"" + std::string(field_name) +
          "\" is not a valid integer value in file \"" + rst_name + "\"");

  return value;
}

unsigned long long int apo::parse_rst_ull(const std::string &line,
                                          const std::size_t offset,
                                          const std::size_t width,
                                          const std::string_view field_name,
                                          const std::string &rst_name) {
  const std::string field =
      get_rst_field(line, offset, width, field_name, rst_name);

  unsigned long long int value = 0;

  APOCHARMM_REQUIRE(
      apo::try_parse_ull(value, field), ApoCharmmErrorCode::Runtime,
      "Restart field \"" + std::string(field_name) +
          "\" is not a valid integer value in file \"" + rst_name + "\"");

  return value;
}

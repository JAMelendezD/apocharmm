// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#pragma once

#include <cstddef>
#include <filesystem>
#include <iosfwd>
#include <string>
#include <string_view>
#include <vector>

namespace apo {

void ltrim_ip(std::string &str);
void rtrim_ip(std::string &str);
void trim_ip(std::string &str);

std::string ltrim(const std::string_view str);
std::string rtrim(const std::string_view str);
std::string trim(const std::string_view str);

void to_lower_ip(std::string &str);
void to_upper_ip(std::string &str);

std::string to_lower(const std::string_view str);
std::string to_upper(const std::string_view str);

std::vector<std::string> split(const std::string_view str,
                               const std::string_view delimiter = " ");

bool try_get_line(std::string &line, std::size_t &pos,
                  const std::string_view file_data);

bool try_get_line(std::string &line, std::size_t &pos, std::size_t &line_number,
                  const std::string_view file_data);

bool try_get_line(std::string &line, std::size_t &line_number,
                  std::istream &input);

void get_line(std::string &line, std::size_t &pos,
              const std::string_view file_data);

void get_line(std::string &line, std::size_t &pos, std::size_t &line_number,
              const std::string_view file_data,
              const std::string_view record_name,
              const std::string_view source_name);

void get_line(std::string &line, std::size_t &line_number, std::istream &input,
              const std::string_view record_name,
              const std::string_view source_name);

void find_required_line(std::istream &input, std::size_t &line_number,
                        const std::string_view target,
                        const std::string_view record_name,
                        const std::string_view source_name);

bool try_get_fixed_width_field(std::string_view &field,
                               const std::string_view line,
                               const std::size_t offset,
                               const std::size_t width) noexcept;

void read_file_into_string(std::string &file_data,
                           const std::filesystem::path &file_path);

bool try_parse_int(int &value, const std::string_view token);

bool try_parse_ull(unsigned long long int &value, const std::string_view token);

bool try_parse_double(double &value, const std::string_view token);

bool try_parse_finite_double(double &value, const std::string_view token);

int parse_int(const std::string_view token, const std::string_view field_name,
              const std::string_view context);

unsigned long long int parse_ull(const std::string_view token,
                                 const std::string_view field_name,
                                 const std::string_view context);

double parse_double(const std::string_view token,
                    const std::string_view field_name,
                    const std::string_view context);

bool contains_wildcard(const std::string_view str);

std::string cDoubleToFortSciStr(const double val, const int prec);

double fortSciStrToCDouble(const std::string_view str);

std::string get_fixed_width_field(const std::string_view line,
                                  const std::size_t offset,
                                  const std::size_t width,
                                  const std::string_view field_name,
                                  const std::string_view context);

double parse_fixed_width_double(const std::string_view line,
                                const std::size_t offset,
                                const std::size_t width,
                                const std::string_view field_name,
                                const std::string_view context);

int parse_fixed_width_int(const std::string_view line, const std::size_t offset,
                          const std::size_t width,
                          const std::string_view field_name,
                          const std::string_view context);

unsigned long long int parse_fixed_width_ull(const std::string_view line,
                                             const std::size_t offset,
                                             const std::size_t width,
                                             const std::string_view field_name,
                                             const std::string_view context);

} // namespace apo

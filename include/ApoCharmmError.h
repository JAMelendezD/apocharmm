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

#include <stdexcept>
#include <string>
#include <string_view>

enum class ApoCharmmErrorCode {
  InvalidArgument,
  Runtime,
  Cuda,
  NotInitialized,
  NotImplemented
};

std::string_view
GetApoCharmmErrorCodeName(const ApoCharmmErrorCode code) noexcept;

class ApoCharmmError : public std::runtime_error {
public:
  ApoCharmmError(const ApoCharmmErrorCode code, const std::string_view message,
                 const std::string_view sourceFile,
                 const std::string_view sourceFunction, const int sourceLine);

public:
  ApoCharmmErrorCode getCode(void) const noexcept;
  std::string_view getMessage(void) const noexcept;
  std::string_view getSourceFile(void) const noexcept;
  std::string_view getSourceFunction(void) const noexcept;
  int getSourceLine(void) const noexcept;

private:
  ApoCharmmErrorCode m_Code;
  std::string m_Message;
  std::string m_SourceFile;
  std::string m_SourceFunction;
  int m_SourceLine;
};

#define APOCHARMM_THROW(errorCode, message)                                    \
  throw ::ApoCharmmError((errorCode), (message), __FILE__, __func__, __LINE__)

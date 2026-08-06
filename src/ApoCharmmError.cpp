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

#include <cstddef>
#include <sstream>
#include <utility>

namespace {

constexpr std::string_view APOCHARMM_SKULL =
    R"ASCII(                            /---------\
                           /           \
                          /             \
                         /               \
                         !         XXXX  !
                         !  ----   XXXX  !
                         !          XXX  !
                         !       X       !
                          --\   XXX   /--
                           ! !  XXX  ! !
                           ! !       ! !
                           ! I I I I I !
                           !  I I I I  !
                            \         /
                             --     --
                               \---/
                        XXX             XXX
                       XXXX             XXXX
                       XXXXX           XXXXX
                          XXX         XXX
                            XXX     XXX
                               XXXXX
                              XXX XXX
                            XXX     XXX
                          XXX         XXX
                       XXXXX           XXXXX
                       XXXX             XXXX
                        XXX             XXX)ASCII";

std::string_view
GetProjectRelativeSourceFile(const std::string_view sourceFile) noexcept {
  constexpr std::string_view SOURCE_ROOT{APOCHARMM_SOURCE_ROOT};

  if (sourceFile.compare(0, SOURCE_ROOT.size(), SOURCE_ROOT) == 0)
    return sourceFile.substr(SOURCE_ROOT.size());

  const std::size_t separator = sourceFile.find_last_of("/\\");
  if (separator != std::string_view::npos)
    return sourceFile.substr(separator + 1);

  return sourceFile;
}

} // namespace

std::string_view
GetApoCharmmErrorCodeName(const ApoCharmmErrorCode code) noexcept {
  switch (code) {
  case ApoCharmmErrorCode::InvalidArgument:
    return "InvalidArgument";
  case ApoCharmmErrorCode::Runtime:
    return "Runtime";
  case ApoCharmmErrorCode::Cuda:
    return "Cuda";
  case ApoCharmmErrorCode::NotInitialized:
    return "NotInitialized";
  case ApoCharmmErrorCode::NotImplemented:
    return "NotImplemented";
  }
  return "Unknown";
}

ApoCharmmError::ApoCharmmError(const ApoCharmmErrorCode code,
                               const std::string_view message,
                               const std::string_view sourceFile,
                               const std::string_view sourceFunction,
                               const int sourceLine)
    : std::runtime_error([code, message, sourceFile, sourceFunction,
                          sourceLine](void) -> std::string {
        const std::string_view projectRelativeSourceFile =
            GetProjectRelativeSourceFile(sourceFile);

        std::ostringstream diagnostic;
        diagnostic << '\n'
                   << '\n'
                   << APOCHARMM_SKULL << '\n'
                   << '\n'
                   << "apoCHARMM error [" << GetApoCharmmErrorCodeName(code)
                   << "]: " << message << '\n'
                   << "  source: " << projectRelativeSourceFile << ':'
                   << sourceLine << '\n'
                   << "  function: " << sourceFunction;

        return diagnostic.str();
      }()),
      m_Code(code), m_Message(message),
      m_SourceFile(GetProjectRelativeSourceFile(sourceFile)),
      m_SourceFunction(sourceFunction), m_SourceLine(sourceLine) {}

ApoCharmmErrorCode ApoCharmmError::getCode(void) const noexcept {
  return m_Code;
}

std::string_view ApoCharmmError::getMessage(void) const noexcept {
  return m_Message;
}

std::string_view ApoCharmmError::getSourceFile(void) const noexcept {
  return m_SourceFile;
}

std::string_view ApoCharmmError::getSourceFunction(void) const noexcept {
  return m_SourceFunction;
}

int ApoCharmmError::getSourceLine(void) const noexcept { return m_SourceLine; }

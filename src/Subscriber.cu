// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author:  Samarjeet Prasad, James E. Gonzales II
//
// ENDLICENSE

#include "Subscriber.h"

#include "ApoCharmmError.h"
#include "CudaIntegrator.h"

#include <fstream>
#include <iostream>
#include <sys/stat.h>

Subscriber::Subscriber(void)
    : m_ReportFrequency(1000), m_FileName(""), m_FileStream(),
      m_CharmmContext(nullptr), m_Integrator(nullptr) {}

Subscriber::Subscriber(const std::string &fileName) : Subscriber() {
  this->setFileName(fileName);
  this->openFile();
}

Subscriber::Subscriber(const std::string &fileName, const int reportFrequency)
    : Subscriber() {
  this->setReportFrequency(reportFrequency);
  this->setFileName(fileName);
  this->openFile();
}

void Subscriber::setReportFrequency(const int reportFrequency) {
  APOCHARMM_REQUIRE(reportFrequency > 0, ApoCharmmErrorCode::InvalidArgument,
                    "Subscriber report frequency must be positive; observed " +
                        std::to_string(reportFrequency));

  m_ReportFrequency = reportFrequency;

  return;
}

void Subscriber::setFileName(const std::string &fileName) {
  APOCHARMM_REQUIRE(!fileName.empty(), ApoCharmmErrorCode::InvalidArgument,
                    "Output file name must not be empty");

  this->checkPath(fileName);

  m_FileName = fileName;

  return;
}

void Subscriber::setCharmmContext(std::shared_ptr<CharmmContext> ctx) {
  APOCHARMM_REQUIRE(m_CharmmContext == nullptr,
                    ApoCharmmErrorCode::InvalidArgument,
                    "Subscriber already has a CharmmContext");

  m_CharmmContext = ctx;

  return;
}

void Subscriber::setIntegrator(std::shared_ptr<CudaIntegrator> integrator) {
  APOCHARMM_REQUIRE(m_Integrator == nullptr,
                    ApoCharmmErrorCode::InvalidArgument,
                    "Subscriber already has an Integrator");

  m_Integrator = integrator;

  return;
}

int Subscriber::getReportFrequency(void) const { return m_ReportFrequency; }

const std::string &Subscriber::getFileName(void) const { return m_FileName; }

std::string &Subscriber::getFileName(void) { return m_FileName; }

void Subscriber::checkPath(const std::string &fileName) {
  std::string totalFilePath = fileName;

  // If no "/" character in pathandfile, nothing to check
  if (totalFilePath.find('/') == std::string::npos)
    return;

  struct stat sb;
  std::size_t botDirPos = totalFilePath.find_last_of('/');
  std::string dirName = totalFilePath.substr(0, botDirPos);
  if (dirName != "") {
    // Check that the directory actually exists
    APOCHARMM_REQUIRE(stat(dirName.c_str(), &sb) == 0,
                      ApoCharmmErrorCode::InvalidArgument,
                      "Output directory does not exist: " + dirName);
  }

  return;
}

void Subscriber::openFile(void) {
  APOCHARMM_REQUIRE(!m_FileName.empty(), ApoCharmmErrorCode::NotInitialized,
                    "Subscriber output file name is not set");

  this->checkPath(m_FileName);

  if (m_FileStream.is_open())
    m_FileStream.close();

  m_FileStream.clear();
  m_FileStream.open(m_FileName, std::ios::out);

  APOCHARMM_REQUIRE(m_FileStream.is_open() && m_FileStream.good(),
                    ApoCharmmErrorCode::Runtime,
                    "Failed to open subscriber output file for writing: " +
                        m_FileName);

  return;
}

void Subscriber::addCommentSection(const std::string &commentLines) {
  APOCHARMM_REQUIRE(!commentLines.empty(), ApoCharmmErrorCode::InvalidArgument,
                    "Comment section must not be empty");

  APOCHARMM_REQUIRE(m_FileStream.is_open(), ApoCharmmErrorCode::NotInitialized,
                    "Subscriber output file is not open");

  // Make sure that last char is a line break
  std::string commentSection = commentLines;
  if (commentSection.back() != '\n')
    commentSection += "\n";

  m_FileStream << commentSection;
  m_FileStream.flush();

  APOCHARMM_REQUIRE(
      m_FileStream.good(), ApoCharmmErrorCode::Runtime,
      "Failed to write comment section to subscriber output file: " +
          m_FileName);

  return;
}

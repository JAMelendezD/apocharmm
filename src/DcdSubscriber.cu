// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: Samarjeet Prasad, James E. Gonzales II
//
// ENDLICENSE

#include "DcdSubscriber.h"

#include "ApoCharmmError.h"
#include "CharmmContext.h"
#include "CudaIntegrator.h"

#include <array>
#include <cstddef>
#include <fstream>
#include <limits>
#include <vector>

DcdSubscriber::DcdSubscriber(const std::string &fileName) : Subscriber() {
  m_NumFramesWritten = 0;
  m_IsHeaderWritten = false;
  this->setFileName(fileName);
  this->openFile();
}

DcdSubscriber::DcdSubscriber(const std::string &fileName,
                             const int reportFrequency)
    : Subscriber() {
  m_NumFramesWritten = 0;
  m_IsHeaderWritten = false;
  this->setReportFrequency(reportFrequency);
  this->setFileName(fileName);
  this->openFile();
}

void DcdSubscriber::update(void) {
  APOCHARMM_REQUIRE(!m_FileName.empty(), ApoCharmmErrorCode::NotInitialized,
                    "DcdSubscriber requires an output file before update");

  APOCHARMM_REQUIRE(m_CharmmContext != nullptr,
                    ApoCharmmErrorCode::NotInitialized,
                    "DcdSubscriber requires a CharmmContext before update");

  APOCHARMM_REQUIRE(m_Integrator != nullptr, ApoCharmmErrorCode::NotInitialized,
                    "DcdSubscriber requires an integrator before update");

  APOCHARMM_REQUIRE(m_FileStream.is_open(), ApoCharmmErrorCode::NotInitialized,
                    "DcdSubscriber output file is not open for writing: " +
                        m_FileName);

  const std::vector<double> boxDimensions = m_CharmmContext->getBoxDimensions();

  APOCHARMM_REQUIRE(
      (boxDimensions.size() == 3) && (boxDimensions[0] > 0.0) &&
          (boxDimensions[1] > 0.0) && (boxDimensions[2] > 0.0),
      ApoCharmmErrorCode::NotInitialized,
      "DcdSubscriber requires three positive box dimensions before update");

  const int numAtoms = m_CharmmContext->getNumAtoms();

  APOCHARMM_REQUIRE(numAtoms > 0, ApoCharmmErrorCode::NotInitialized,
                    "DcdSubscriber requires at least one atom before update");

  const unsigned long long int nextNumFramesWritten =
      static_cast<unsigned long long int>(m_NumFramesWritten) + 1ULL;
  const unsigned long long int nextNumStepsWritten =
      nextNumFramesWritten *
      static_cast<unsigned long long int>(m_ReportFrequency);

  APOCHARMM_REQUIRE(
      (nextNumFramesWritten <=
       static_cast<unsigned long long int>(std::numeric_limits<int>::max())) &&
          (nextNumStepsWritten <= static_cast<unsigned long long int>(
                                      std::numeric_limits<int>::max())),
      ApoCharmmErrorCode::Runtime,
      "DcdSubscriber frame metadata exceeds DCD integer range");

  // Write header if needed
  if (!m_IsHeaderWritten)
    this->writeHeader();

  // JEG250610: Add conditional statement for this
  this->writeXtalData();

  // Write coordinate data
  this->writeCoordData();

  // Update header to account for new data that was written
  const int numFramesWritten = static_cast<int>(nextNumFramesWritten);
  const int numStepsWritten = static_cast<int>(nextNumStepsWritten);

  m_FileStream.seekp(8);
  m_FileStream.write(reinterpret_cast<const char *>(&numFramesWritten),
                     sizeof(int));

  m_FileStream.seekp(20);
  m_FileStream.write(reinterpret_cast<const char *>(&numStepsWritten),
                     sizeof(int));

  m_FileStream.seekp(0, std::ofstream::end);
  m_FileStream.flush();

  APOCHARMM_REQUIRE(m_FileStream.good(), ApoCharmmErrorCode::Runtime,
                    "Failed to write DCD frame: " + m_FileName);

  m_NumFramesWritten = numFramesWritten;

  return;
}

void DcdSubscriber::openFile(void) {
  APOCHARMM_REQUIRE(!m_FileName.empty(), ApoCharmmErrorCode::NotInitialized,
                    "DcdSubscriber output file name is not set");

  this->checkPath(m_FileName);

  if (m_FileStream.is_open())
    m_FileStream.close();

  m_FileStream.clear();
  m_FileStream.open(m_FileName, std::ios::out | std::ios::binary);

  APOCHARMM_REQUIRE(m_FileStream.is_open() && m_FileStream.good(),
                    ApoCharmmErrorCode::Runtime,
                    "Failed to open DCD file for writing: " + m_FileName);

  m_NumFramesWritten = 0;
  m_IsHeaderWritten = false;

  return;
}

void DcdSubscriber::writeHeader(void) {
  std::array<char, 1024> buffer;
  int *iPtr = reinterpret_cast<int *>(buffer.data());
  float *fPtr = reinterpret_cast<float *>(buffer.data());
  int blockSizeBytes = 0;

  int start = 0; // JEG250610: Get this from integrator
  int nstep = 0; // This is updated as it is written
  int nsavv = 0;
  int ndegf = m_CharmmContext->getNumDegreesOfFreedom();
  int namnf = 0;
  float delta = m_Integrator->getTimeStep();
  bool hasXtalData = true;
  bool has4dData = false;
  bool hasCheqData = false;
  bool hasContData = false;
  int version = 49;

  // Fill input control block
  blockSizeBytes = 84;
  buffer.fill(0);
  iPtr[0] = 1146244931;        // 'C','O','R','D'
  iPtr[1] = 0;                 // Number of frames written (NFILE)
  iPtr[2] = start;             // Starting step (START)
  iPtr[3] = m_ReportFrequency; // Reporting frequency (NSAVC)
  iPtr[4] = nstep;  // Report frequency * number of frames written (NSTEP)
  iPtr[5] = nsavv;  // Velocity saving frequency (NSAVV)
  iPtr[6] = 0;      // UNUSED OPTION
  iPtr[7] = 0;      // UNUSED OPTION
  iPtr[8] = ndegf;  // Number of degrees of freedom (NDEGF)
  iPtr[9] = namnf;  // Number of fixed atoms (NAMNF)
  fPtr[10] = delta; // Time step (DELTA)
  iPtr[11] = (hasXtalData) ? 1 : 0;
  iPtr[12] = (has4dData) ? 1 : 0;
  iPtr[13] = (hasCheqData) ? 1 : 0;
  iPtr[14] = (hasContData) ? 1 : 0;
  iPtr[15] = 0;       // UNUSED OPTION
  iPtr[16] = 0;       // UNUSED OPTION
  iPtr[17] = 0;       // UNUSED OPTION
  iPtr[18] = 0;       // UNUSED OPTION
  iPtr[19] = 0;       // UNUSED OPTION
  iPtr[20] = version; // CHARMM version, should always be >= 22

  // Write start size of input control block
  m_FileStream.write(reinterpret_cast<const char *>(&blockSizeBytes),
                     sizeof(int));

  // Write input control block
  m_FileStream.write(buffer.data(), blockSizeBytes);

  // Write end size of input control block
  m_FileStream.write(reinterpret_cast<const char *>(&blockSizeBytes),
                     sizeof(int));

  // Fill title block
  blockSizeBytes = 164;
  buffer.fill(' ');
  iPtr[0] = 2;
  std::string title = "* GENERATED BY APOCHARMM";
  for (std::size_t i = 0; i < title.size(); i++)
    buffer[4 + i] = title[i];
  buffer[84] = '*';

  // Write start size of title block
  m_FileStream.write(reinterpret_cast<const char *>(&blockSizeBytes),
                     sizeof(int));

  // Write title block
  m_FileStream.write(buffer.data(), blockSizeBytes);

  // Write end size of title block
  m_FileStream.write(reinterpret_cast<const char *>(&blockSizeBytes),
                     sizeof(int));

  // Fill atom block
  blockSizeBytes = 4;
  buffer.fill(0);
  iPtr[0] = m_CharmmContext->getNumAtoms();

  // Write start size of atom block
  m_FileStream.write(reinterpret_cast<const char *>(&blockSizeBytes),
                     sizeof(int));

  // Write atom block
  m_FileStream.write(buffer.data(), blockSizeBytes);

  // Write end size of atom block
  m_FileStream.write(reinterpret_cast<const char *>(&blockSizeBytes),
                     sizeof(int));

  m_FileStream.flush();

  APOCHARMM_REQUIRE(m_FileStream.good(), ApoCharmmErrorCode::Runtime,
                    "Failed to write DCD header: " + m_FileName);

  m_IsHeaderWritten = true;

  return;
}

void DcdSubscriber::writeXtalData(void) {
  std::array<char, 1024> buffer;
  double *dPtr = reinterpret_cast<double *>(buffer.data());
  int blockSizeBytes = 6 * static_cast<int>(sizeof(double));
  std::vector<double> boxDims = m_CharmmContext->getBoxDimensions();

  // Fille XTAL block
  buffer.fill(0);
  dPtr[0] = boxDims[0];
  dPtr[2] = boxDims[1];
  dPtr[5] = boxDims[2];

  // Write start size of XTAL block
  m_FileStream.write(reinterpret_cast<const char *>(&blockSizeBytes),
                     sizeof(int));

  // Write XTAL block
  m_FileStream.write(buffer.data(), blockSizeBytes);

  // Write end size of XTAL block
  m_FileStream.write(reinterpret_cast<const char *>(&blockSizeBytes),
                     sizeof(int));

  return;
}

void DcdSubscriber::writeCoordData(void) {
  m_CharmmContext->getCoordinatesChargesSP().transferToHost();
  const float4 *xyzq =
      m_CharmmContext->getCoordinatesChargesSP().getHostArray().data();
  const int numAtoms = m_CharmmContext->getNumAtoms();
  const int blockSizeBytes = numAtoms * static_cast<int>(sizeof(float));

  // Write X coordinate block
  m_FileStream.write(reinterpret_cast<const char *>(&blockSizeBytes),
                     sizeof(int));
  for (int i = 0; i < numAtoms; i++)
    m_FileStream.write(reinterpret_cast<const char *>(&xyzq[i].x),
                       sizeof(float));
  m_FileStream.write(reinterpret_cast<const char *>(&blockSizeBytes),
                     sizeof(int));

  // Write Y coordinate block
  m_FileStream.write(reinterpret_cast<const char *>(&blockSizeBytes),
                     sizeof(int));
  for (int i = 0; i < numAtoms; i++)
    m_FileStream.write(reinterpret_cast<const char *>(&xyzq[i].y),
                       sizeof(float));
  m_FileStream.write(reinterpret_cast<const char *>(&blockSizeBytes),
                     sizeof(int));

  // Write Z coordinate block
  m_FileStream.write(reinterpret_cast<const char *>(&blockSizeBytes),
                     sizeof(int));
  for (int i = 0; i < numAtoms; i++)
    m_FileStream.write(reinterpret_cast<const char *>(&xyzq[i].z),
                       sizeof(float));
  m_FileStream.write(reinterpret_cast<const char *>(&blockSizeBytes),
                     sizeof(int));

  return;
}

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
#include "CharmmContext.h"
#include "CharmmCrd.h"
#include "CharmmPSF.h"
#include "CharmmParameters.h"
#include "CudaIntegrator.h"
#include "DcdSubscriber.h"
#include "ForceManager.h"
#include "apo_test_helpers.h"
#include "catch.hpp"
#include "test_paths.h"

#include <array>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double TIME_STEP = 0.001;
constexpr int REPORT_FREQUENCY = 2;
constexpr double FLOAT_TOLERANCE = 1.0e-6;

const std::vector<double> BOX_DIMENSIONS = {50.0, 50.0, 50.0};

struct DcdFrame {
  std::array<double, 6> xtal = {};
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;
};

struct DcdFile {
  int cord = 0;
  int numFrames = 0;
  int start = 0;
  int reportFrequency = 0;
  int numSteps = 0;
  int numDegreesOfFreedom = 0;
  int hasXtalData = 0;
  int version = 0;
  int numAtoms = 0;
  float timeStep = 0.0f;
  std::vector<DcdFrame> frames;
};

template <typename T> T ReadValue(std::ifstream &input) {
  T value{};
  input.read(reinterpret_cast<char *>(&value), sizeof(T));
  REQUIRE(input.good());

  return value;
}

template <typename T>
T ReadBufferValue(const std::vector<char> &buffer, const std::size_t offset) {
  REQUIRE(offset + sizeof(T) <= buffer.size());

  T value{};
  std::memcpy(static_cast<void *>(&value),
              static_cast<const void *>(buffer.data() + offset), sizeof(T));

  return value;
}

std::vector<char> ReadBytes(std::ifstream &input, const std::size_t count) {
  std::vector<char> buffer(count);
  input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
  REQUIRE(input.good());

  return buffer;
}

std::vector<float> ReadFloatBlock(std::ifstream &input, const int numAtoms) {
  const int blockSizeBytes = ReadValue<int>(input);
  REQUIRE(blockSizeBytes == numAtoms * static_cast<int>(sizeof(float)));

  std::vector<float> values(static_cast<std::size_t>(numAtoms));
  input.read(reinterpret_cast<char *>(values.data()), blockSizeBytes);
  REQUIRE(input.good());

  const int endBlockSizeBytes = ReadValue<int>(input);
  REQUIRE(endBlockSizeBytes == blockSizeBytes);

  return values;
}

DcdFile ReadDcdFile(const std::string &fileName) {
  std::ifstream input(fileName, std::ios::binary);
  REQUIRE(input.good());

  DcdFile dcd;

  const int controlBlockSize = ReadValue<int>(input);
  REQUIRE(controlBlockSize == 84);

  const std::vector<char> controlBlock =
      ReadBytes(input, static_cast<std::size_t>(controlBlockSize));

  const int controlBlockEndSize = ReadValue<int>(input);
  REQUIRE(controlBlockEndSize == controlBlockSize);

  dcd.cord = ReadBufferValue<int>(controlBlock, 0);
  dcd.numFrames = ReadBufferValue<int>(controlBlock, 4);
  dcd.start = ReadBufferValue<int>(controlBlock, 8);
  dcd.reportFrequency = ReadBufferValue<int>(controlBlock, 12);
  dcd.numSteps = ReadBufferValue<int>(controlBlock, 16);
  dcd.numDegreesOfFreedom = ReadBufferValue<int>(controlBlock, 32);
  dcd.timeStep = ReadBufferValue<float>(controlBlock, 40);
  dcd.hasXtalData = ReadBufferValue<int>(controlBlock, 44);
  dcd.version = ReadBufferValue<int>(controlBlock, 80);

  const int titleBlockSize = ReadValue<int>(input);
  REQUIRE(titleBlockSize == 164);
  ReadBytes(input, static_cast<std::size_t>(titleBlockSize));
  const int titleEndBlocksize = ReadValue<int>(input);
  REQUIRE(titleEndBlocksize == titleBlockSize);

  const int atomBlockSize = ReadValue<int>(input);
  REQUIRE(atomBlockSize == 4);
  dcd.numAtoms = ReadValue<int>(input);
  const int atomEndBlockSize = ReadValue<int>(input);
  REQUIRE(atomEndBlockSize == atomBlockSize);

  while (input.peek() != std::ifstream::traits_type::eof()) {
    DcdFrame frame;

    const int xtalBlockSize = ReadValue<int>(input);
    REQUIRE(xtalBlockSize == 6 * static_cast<int>(sizeof(double)));
    for (std::size_t i = 0; i < frame.xtal.size(); i++)
      frame.xtal[i] = ReadValue<double>(input);
    const int xtalEndBlockSize = ReadValue<int>(input);
    REQUIRE(xtalEndBlockSize == xtalBlockSize);

    frame.x = ReadFloatBlock(input, dcd.numAtoms);
    frame.y = ReadFloatBlock(input, dcd.numAtoms);
    frame.z = ReadFloatBlock(input, dcd.numAtoms);

    dcd.frames.push_back(frame);
  }

  return dcd;
}

std::shared_ptr<CharmmContext> CreateContext(void) {
  auto prm = std::make_shared<CharmmParameters>(getDataPath() +
                                                "toppar_water_ions.str");
  auto psf = std::make_shared<CharmmPSF>(getDataPath() + "nacl_pair.psf");
  auto crd = std::make_shared<CharmmCrd>(getDataPath() + "nacl_pair.cor");

  auto ctx = std::make_shared<CharmmContext>(psf, prm);
  ctx->setBoxDimensions(BOX_DIMENSIONS);
  ctx->setCoordinates(crd);
  ctx->useHolonomicConstraints(false);

  return ctx;
}

void CheckFrameMatchesContext(const DcdFrame &frame,
                              std::shared_ptr<CharmmContext> ctx) {
  ctx->getCoordinatesChargesSP().transferToHost();

  const std::vector<float4> &xyzq =
      ctx->getCoordinatesChargesSP().getHostArray();

  REQUIRE(frame.x.size() == xyzq.size());
  REQUIRE(frame.y.size() == xyzq.size());
  REQUIRE(frame.z.size() == xyzq.size());

  for (std::size_t i = 0; i < xyzq.size(); i++) {
    CHECK(frame.x[i] == Approx(xyzq[i].x).margin(FLOAT_TOLERANCE));
    CHECK(frame.y[i] == Approx(xyzq[i].y).margin(FLOAT_TOLERANCE));
    CHECK(frame.z[i] == Approx(xyzq[i].z).margin(FLOAT_TOLERANCE));
  }

  return;
}

} // namespace

TEST_CASE("DcdSubscriberConstructorsAndReportFrequency") {
  SECTION("DefaultReportFrequencyConstructor") {
    const std::string fileName = "tmp_dcd_subscriber_default.dcd";
    apo_test::RemoveIfExists(fileName);

    DcdSubscriber subscriber(fileName);

    CHECK(subscriber.getFileName() == fileName);
    CHECK(subscriber.getReportFrequency() == 1000);

    subscriber.setReportFrequency(REPORT_FREQUENCY);
    CHECK(subscriber.getReportFrequency() == REPORT_FREQUENCY);

    std::ifstream file(fileName, std::ios::binary);
    CHECK(file.good());

    apo_test::RemoveIfExists(fileName);
  }

  SECTION("ExplicitReportFrequencyConstructor") {
    const std::string fileName = "tmp_dcd_subscriber_frequency.dcd";
    apo_test::RemoveIfExists(fileName);

    DcdSubscriber subscriber(fileName, REPORT_FREQUENCY);

    CHECK(subscriber.getFileName() == fileName);
    CHECK(subscriber.getReportFrequency() == REPORT_FREQUENCY);

    apo_test::RemoveIfExists(fileName);
  }

  SECTION("RejectsMissingOutputDirectory") {
    apo_test::CheckApoCharmmError(
        []() { (void)DcdSubscriber("missing_dcd_subscriber_dir/output.dcd"); },
        ApoCharmmErrorCode::InvalidArgument,
        "Output directory does not exist: missing_dcd_subscriber_dir");
  }
}

TEST_CASE("DcdSubscriberUpdateWritesOneFrame") {
  const std::string fileName = "tmp_dcd_subscriber_one_frame.dcd";
  apo_test::RemoveIfExists(fileName);

  auto ctx = CreateContext();
  auto integrator = std::make_shared<CudaIntegrator>(TIME_STEP);

  DcdSubscriber subscriber(fileName, REPORT_FREQUENCY);
  subscriber.setCharmmContext(ctx);
  subscriber.setIntegrator(integrator);

  apo_test::CheckApoCharmmError([&]() { subscriber.setCharmmContext(ctx); },
                                ApoCharmmErrorCode::InvalidArgument,
                                "Subscriber already has a CharmmContext");
  apo_test::CheckApoCharmmError([&]() { subscriber.setIntegrator(integrator); },
                                ApoCharmmErrorCode::InvalidArgument,
                                "Subscriber already has an Integrator");

  subscriber.update();

  const DcdFile dcd = ReadDcdFile(fileName);

  CHECK(dcd.cord == 1146244931);
  CHECK(dcd.numFrames == 1);
  CHECK(dcd.start == 0);
  CHECK(dcd.reportFrequency == REPORT_FREQUENCY);
  CHECK(dcd.numSteps == REPORT_FREQUENCY);
  CHECK(dcd.numDegreesOfFreedom == ctx->getNumDegreesOfFreedom());
  CHECK(dcd.timeStep ==
        Approx(static_cast<float>(TIME_STEP)).margin(FLOAT_TOLERANCE));
  CHECK(dcd.hasXtalData == 1);
  CHECK(dcd.version == 49);
  CHECK(dcd.numAtoms == ctx->getNumAtoms());

  REQUIRE(dcd.frames.size() == 1);
  CHECK(dcd.frames[0].xtal[0] == Approx(BOX_DIMENSIONS[0]));
  CHECK(dcd.frames[0].xtal[2] == Approx(BOX_DIMENSIONS[1]));
  CHECK(dcd.frames[0].xtal[5] == Approx(BOX_DIMENSIONS[2]));
  CheckFrameMatchesContext(dcd.frames[0], ctx);

  apo_test::RemoveIfExists(fileName);
}

TEST_CASE("DcdSubscriberUpdatePatchesHeaderForMultipleFrames") {
  const std::string fileName = "tmp_dcd_subscriber_two_frames.dcd";
  apo_test::RemoveIfExists(fileName);

  auto ctx = CreateContext();
  auto integrator = std::make_shared<CudaIntegrator>(TIME_STEP);

  DcdSubscriber subscriber(fileName, REPORT_FREQUENCY);
  subscriber.setCharmmContext(ctx);
  subscriber.setIntegrator(integrator);

  subscriber.update();
  subscriber.update();

  const DcdFile dcd = ReadDcdFile(fileName);

  CHECK(dcd.numFrames == 2);
  CHECK(dcd.reportFrequency == REPORT_FREQUENCY);
  CHECK(dcd.numSteps == 2 * REPORT_FREQUENCY);
  CHECK(dcd.numAtoms == ctx->getNumAtoms());

  REQUIRE(dcd.frames.size() == 2);
  CheckFrameMatchesContext(dcd.frames[0], ctx);
  CheckFrameMatchesContext(dcd.frames[1], ctx);

  apo_test::RemoveIfExists(fileName);
}

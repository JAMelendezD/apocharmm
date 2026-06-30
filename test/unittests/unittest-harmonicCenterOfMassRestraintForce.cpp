// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "AtomSelection.h"
#include "CharmmContext.h"
#include "CharmmCrd.h"
#include "CharmmPSF.h"
#include "CharmmParameters.h"
#include "CudaContainer.h"
#include "CudaEnergyVirial.h"
#include "CudaLangevinThermostatIntegrator.h"
#include "ForceManager.h"
#include "HarmonicCenterOfMassRestraintForce.h"
#include "apo_test_helpers.h"
#include "catch.hpp"
#include "cuda_utils.h"
#include "test_paths.h"

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <vector_types.h>

namespace {

using HmcmForce = HarmonicCenterOfMassRestraintForce<long long int, float>;

constexpr int NUM_ATOMS = 3;
constexpr int RANDOM_SEED = 314159;
constexpr double TEMPERATURE = 300.0;
constexpr double TIME_STEP = 0.001;
constexpr double LANGEVIN_THERMOSTAT_FRICTION = 0.0;
constexpr double FORCE_TOLERANCE = 0.0;
constexpr double ENERGY_TOLERANCE = 0.0;
constexpr double VIRIAL_TOLERANCE = 0.0;
constexpr double INV_FORCE_SCALE_TEST = 1.0 / static_cast<double>(1LL << 40);

const std::vector<double> BOX_DIMENSIONS = {50.0, 50.0, 50.0};
const std::vector<double> SMALL_BOX_DIMENSIONS = {20.0, 20.0, 20.0};

AtomSelection MakeSelection(const int numAtoms,
                            const std::vector<int> &selectedAtoms) {
  AtomSelection selection(numAtoms, AtomSelection::InitialValue::NONE);

  for (const int atomIndex : selectedAtoms)
    selection.set(atomIndex);

  return selection;
}

std::vector<double3>
CopyFixedPointForceToHost(const std::shared_ptr<Force<long long int>> &force) {
  std::vector<long long int> fx(force->size());
  std::vector<long long int> fy(force->size());
  std::vector<long long int> fz(force->size());

  force->getXYZ(fx.data(), fy.data(), fz.data());

  std::vector<double3> hostForces(force->size());
  for (std::size_t i = 0; i < hostForces.size(); i++) {
    hostForces[i].x = static_cast<double>(fx[i]) * INV_FORCE_SCALE_TEST;
    hostForces[i].y = static_cast<double>(fy[i]) * INV_FORCE_SCALE_TEST;
    hostForces[i].z = static_cast<double>(fz[i]) * INV_FORCE_SCALE_TEST;
  }

  return hostForces;
}

double GetHmcmEnergy(HmcmForce &restraint) {
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));
  restraint.getEnergyVirial()->copyToHost(*restraint.getStream());
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));
  return restraint.getEnergyVirial()->getEnergy("hmcm");
}

std::vector<double> GetHmcmVirial(HmcmForce &restraint,
                                  CudaContainer<float4> &xyzq,
                                  const int numAtoms,
                                  const std::vector<double> &boxDimensions) {
  std::shared_ptr<cudaStream_t> stream = restraint.getStream();

  restraint.getForce()->convert<double>(*stream);
  cudaCheck(cudaStreamSynchronize(*stream));

  restraint.getEnergyVirial()->calcVirial(
      numAtoms, xyzq.getDeviceArray().data(), boxDimensions[0],
      boxDimensions[1], boxDimensions[2], restraint.getForce()->stride(),
      reinterpret_cast<double *>(restraint.getForce()->xyz()), *stream);
  cudaCheck(cudaStreamSynchronize(*stream));

  restraint.getEnergyVirial()->copyToHost(*stream);
  cudaCheck(cudaStreamSynchronize(*stream));

  double virial[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  restraint.getEnergyVirial()->getVirial(virial);

  return std::vector<double>(virial, virial + 9);
}

std::vector<double>
ExpectedVirialFromCenterGradient(const std::vector<double> &center,
                                 const std::vector<double> &gradient) {
  REQUIRE(center.size() == 3);
  REQUIRE(gradient.size() == 3);

  return {-center[0] * gradient[0], -center[0] * gradient[1],
          -center[0] * gradient[2], -center[1] * gradient[0],
          -center[1] * gradient[1], -center[1] * gradient[2],
          -center[2] * gradient[0], -center[2] * gradient[1],
          -center[2] * gradient[2]};
}

} // namespace

TEST_CASE("HarmonicCenterOfMassRestraintForceConstructorAndMetadata") {
  HmcmForce restraint(NUM_ATOMS);

  CHECK(decltype(restraint)::contributesVirial == true);
  CHECK(restraint.getStream() != nullptr);
  CHECK(restraint.getForce() != nullptr);
  CHECK(restraint.getForce()->size() == NUM_ATOMS);
  CHECK(restraint.getForce()->stride() >= NUM_ATOMS);
  CHECK(restraint.getEnergyVirial() != nullptr);
}

TEST_CASE("HarmonicCenterOfMassRestraintForceValidatesInput") {
  HmcmForce restraint(NUM_ATOMS);

  SECTION("ForceConstantValidation") {
    CHECK_THROWS_AS(restraint.setForceConstant(-1.0), std::invalid_argument);
    CHECK_THROWS_AS(
        restraint.setForceConstant(std::numeric_limits<double>::infinity()),
        std::invalid_argument);
    CHECK_NOTHROW(restraint.setForceConstant(0.0));
    CHECK_NOTHROW(restraint.setForceConstant(1.0));
  }

  SECTION("ReferenceDistanceValidation") {
    CHECK_THROWS_AS(restraint.setReferenceDistance(-1.0),
                    std::invalid_argument);
    CHECK_THROWS_AS(
        restraint.setReferenceDistance(std::numeric_limits<double>::infinity()),
        std::invalid_argument);
    CHECK_NOTHROW(restraint.setReferenceDistance(0.0));
    CHECK_NOTHROW(restraint.setReferenceDistance(1.0));
  }

  SECTION("ReferencePositionValidation") {
    CHECK_THROWS_AS(restraint.setReferencePosition({0.0, 0.0}, {1, 1, 1}),
                    std::invalid_argument);
    CHECK_THROWS_AS(restraint.setReferencePosition({0.0, 0.0, 0.0}, {1, 1}),
                    std::invalid_argument);
    CHECK_THROWS_AS(restraint.setReferencePosition({0.0, 0.0, 0.0}, {0, 0, 0}),
                    std::invalid_argument);
    CHECK_THROWS_AS(restraint.setReferencePosition({0.0, 0.0, 0.0}, {1, 2, 1}),
                    std::invalid_argument);
    CHECK_THROWS_AS(
        restraint.setReferencePosition(
            {0.0, std::numeric_limits<double>::infinity(), 0.0}, {1, 1, 1}),
        std::invalid_argument);
    CHECK_NOTHROW(restraint.setReferencePosition({0.0, 0.0, 0.0}));
    CHECK_NOTHROW(restraint.setReferencePosition({0.0, 0.0, 0.0}, {1, 0, 1}));
  }

  SECTION("BoxDimensionValidation") {
    CHECK_THROWS_AS(restraint.setBoxDimensions({10.0, 10.0}),
                    std::invalid_argument);
    CHECK_THROWS_AS(restraint.setBoxDimensions({10.0, -10.0, 10.0}),
                    std::invalid_argument);
    CHECK_THROWS_AS(restraint.setBoxDimensions(
                        {10.0, std::numeric_limits<double>::infinity(), 10.0}),
                    std::invalid_argument);
    CHECK_NOTHROW(restraint.setBoxDimensions({10.0, 10.0, 10.0}));
  }

  SECTION("SelectionValidation") {
    AtomSelection emptySelection(NUM_ATOMS, AtomSelection::InitialValue::NONE);
    CHECK_THROWS_AS(restraint.setSelection(emptySelection),
                    std::invalid_argument);

    AtomSelection wrongSizeSelection(NUM_ATOMS + 1,
                                     AtomSelection::InitialValue::ALL);
    CHECK_THROWS_AS(restraint.setSelection(wrongSizeSelection),
                    std::invalid_argument);

    CHECK_NOTHROW(restraint.setSelection(MakeSelection(NUM_ATOMS, {0})));
  }

  SECTION("MassValidation") {
    CHECK_THROWS_AS(restraint.setMasses({1.0}), std::invalid_argument);
    CHECK_THROWS_AS(restraint.setMasses({0.0, 0.0, 0.0}),
                    std::invalid_argument);
    CHECK_THROWS_AS(restraint.setMasses(
                        {1.0, std::numeric_limits<double>::infinity(), 1.0}),
                    std::invalid_argument);
    CHECK_THROWS_AS(restraint.setMasses({1.0, 1.0, -1.0}),
                    std::invalid_argument);
    CHECK_NOTHROW(restraint.setMasses({1.0, 2.0, 3.0}));
  }

  SECTION("InitializeValidation") {
    CHECK_NOTHROW(restraint.initialize(NUM_ATOMS, {10.0, 10.0, 10.0}));
    CHECK_THROWS_AS(restraint.initialize(4, {10.0, 10.0, 10.0}),
                    std::runtime_error);
  }
}

TEST_CASE("HarmonicCenterOfMassRestraintFroceRequiresBoxForNonzeroForce") {
  HmcmForce restraint(NUM_ATOMS);
  restraint.setSelection(MakeSelection(NUM_ATOMS, {0, 1}));
  restraint.setForceConstant(1.0);
  restraint.setReferencePosition({0.0, 0.0, 0.0});

  CudaContainer<float4> xyzq = std::vector<float4>{
      make_float4(1.0f, 0.0f, 0.0f, 0.0f), make_float4(2.0f, 0.0f, 0.0f, 0.0f),
      make_float4(3.0f, 0.0f, 0.0f, 0.0f)};

  CHECK_THROWS_AS(
      restraint.calcForce(xyzq.getDeviceArray().data(), true, false),
      std::runtime_error);
}

TEST_CASE("HarmonicCenterOfMassRestraintForceNoopsForZeroForceConstant") {
  HmcmForce restraint(NUM_ATOMS);
  restraint.setBoxDimensions(SMALL_BOX_DIMENSIONS);
  restraint.setSelection(MakeSelection(NUM_ATOMS, {0, 1}));
  restraint.setReferencePosition({0.0, 0.0, 0.0});
  restraint.setForceConstant(0.0);

  CudaContainer<float4> xyzq = std::vector<float4>{
      make_float4(1.0f, 0.0f, 0.0f, 0.0f), make_float4(3.0f, 0.0f, 0.0f, 0.0f),
      make_float4(5.0f, 0.0f, 0.0f, 0.0f)};

  restraint.clear();
  restraint.calcForce(xyzq.getDeviceArray().data(), true, true);
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));

  const std::vector<double3> expectedForces = {make_double3(0.0, 0.0, 0.0),
                                               make_double3(0.0, 0.0, 0.0),
                                               make_double3(0.0, 0.0, 0.0)};

  apo_test::CheckVectorsClose<double3>(
      "HMCM zero-force-constant forces",
      CopyFixedPointForceToHost(restraint.getForce()), expectedForces,
      FORCE_TOLERANCE);
  CHECK(GetHmcmEnergy(restraint) == Approx(0.0).margin(ENERGY_TOLERANCE));
}

TEST_CASE("HarmonicCenterOfMassRestraintForceUnweightedPointRestraint") {
  HmcmForce restraint(NUM_ATOMS);
  restraint.setBoxDimensions(SMALL_BOX_DIMENSIONS);
  restraint.setSelection(MakeSelection(NUM_ATOMS, {0, 1}));
  restraint.setForceConstant(2.5);
  restraint.setReferencePosition({1.0, 1.0, 1.0});

  CudaContainer<float4> xyzq = std::vector<float4>{
      make_float4(1.0f, 2.0f, 3.0f, 0.0f), make_float4(3.0f, 4.0f, 5.0f, 0.0f),
      make_float4(9.0f, 9.0f, 9.0f, 0.0f)};

  restraint.clear();
  restraint.calcForce(xyzq.getDeviceArray().data(), true, true);
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));

  const std::vector<double3> expectedForces = {make_double3(2.5, 5.0, 7.5),
                                               make_double3(2.5, 5.0, 7.5),
                                               make_double3(0.0, 0.0, 0.0)};

  apo_test::CheckVectorsClose<double3>(
      "HMCM unweighted point restraint forces",
      CopyFixedPointForceToHost(restraint.getForce()), expectedForces,
      FORCE_TOLERANCE);
  CHECK(GetHmcmEnergy(restraint) == Approx(35.0).margin(ENERGY_TOLERANCE));

  const std::vector<double> expectedVirial =
      ExpectedVirialFromCenterGradient({2.0, 3.0, 4.0}, {5.0, 10.0, 15.0});

  apo_test::CheckVectorsClose<double>(
      "HMCM unweighted point restraint virial",
      GetHmcmVirial(restraint, xyzq, NUM_ATOMS, SMALL_BOX_DIMENSIONS),
      expectedVirial, VIRIAL_TOLERANCE);
}

TEST_CASE("HarmonicCenterOfMassRestraintForceCoordinateMask") {
  HmcmForce restraint(NUM_ATOMS);
  restraint.setBoxDimensions(SMALL_BOX_DIMENSIONS);
  restraint.setSelection(MakeSelection(NUM_ATOMS, {0, 1}));
  restraint.setForceConstant(2.5);
  restraint.setReferencePosition({1.0, 0.0, 0.0}, {1, 0, 0});

  CudaContainer<float4> xyzq = std::vector<float4>{
      make_float4(1.0f, 2.0f, 3.0f, 0.0f), make_float4(3.0f, 4.0f, 5.0f, 0.0f),
      make_float4(9.0f, 9.0f, 9.0f, 0.0f)};

  restraint.clear();
  restraint.calcForce(xyzq.getDeviceArray().data(), true, true);
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));

  const std::vector<double3> expectedForces = {make_double3(2.5, 0.0, 0.0),
                                               make_double3(2.5, 0.0, 0.0),
                                               make_double3(0.0, 0.0, 0.0)};

  apo_test::CheckVectorsClose<double3>(
      "HMCM coordinate-mask forces",
      CopyFixedPointForceToHost(restraint.getForce()), expectedForces,
      FORCE_TOLERANCE);
  CHECK(GetHmcmEnergy(restraint) == Approx(2.5).margin(ENERGY_TOLERANCE));

  const std::vector<double> expectedVirial =
      ExpectedVirialFromCenterGradient({2.0, 3.0, 4.0}, {5.0, 0.0, 0.0});

  apo_test::CheckVectorsClose<double>(
      "HMCM coordinate-mask virial",
      GetHmcmVirial(restraint, xyzq, NUM_ATOMS, SMALL_BOX_DIMENSIONS),
      expectedVirial, VIRIAL_TOLERANCE);
}

TEST_CASE("HarmonicCenterOfMassRestraintForceMassWeighting") {
  HmcmForce restraint(NUM_ATOMS);
  restraint.setBoxDimensions(SMALL_BOX_DIMENSIONS);
  restraint.setSelection(MakeSelection(NUM_ATOMS, {0, 1}));
  restraint.setMasses({1.0, 3.0, 100.0});
  restraint.setForceConstant(2.0);
  restraint.setReferencePosition({1.0, 0.0, 0.0});

  CudaContainer<float4> xyzq = std::vector<float4>{
      make_float4(0.0f, 0.0f, 0.0f, 0.0f), make_float4(4.0f, 0.0f, 0.0f, 0.0f),
      make_float4(9.0f, 9.0f, 9.0f, 0.0f)};

  restraint.clear();
  restraint.calcForce(xyzq.getDeviceArray().data(), true, true);
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));

  const std::vector<double3> expectedForces = {make_double3(2.0, 0.0, 0.0),
                                               make_double3(6.0, 0.0, 0.0),
                                               make_double3(0.0, 0.0, 0.0)};

  apo_test::CheckVectorsClose<double3>(
      "HMCM mass-weighted forces",
      CopyFixedPointForceToHost(restraint.getForce()), expectedForces,
      FORCE_TOLERANCE);
  CHECK(GetHmcmEnergy(restraint) == Approx(8.0).margin(ENERGY_TOLERANCE));

  const std::vector<double> expectedVirial =
      ExpectedVirialFromCenterGradient({3.0, 0.0, 0.0}, {8.0, 0.0, 0.0});

  apo_test::CheckVectorsClose<double>(
      "HMCM mass-weighted virial",
      GetHmcmVirial(restraint, xyzq, NUM_ATOMS, SMALL_BOX_DIMENSIONS),
      expectedVirial, VIRIAL_TOLERANCE);
}

TEST_CASE("HarmonicCenterOfMassRestraintForceCanDisableMassWeighting") {
  HmcmForce restraint(NUM_ATOMS);
  restraint.setBoxDimensions(SMALL_BOX_DIMENSIONS);
  restraint.setSelection(MakeSelection(NUM_ATOMS, {0, 1}));
  restraint.setMasses({1.0, 3.0, 100.0});
  restraint.setMassWeighting(false);
  restraint.setForceConstant(2.0);
  restraint.setReferencePosition({1.0, 0.0, 0.0});

  CudaContainer<float4> xyzq = std::vector<float4>{
      make_float4(0.0f, 0.0f, 0.0f, 0.0f), make_float4(4.0f, 0.0f, 0.0f, 0.0f),
      make_float4(9.0f, 9.0f, 9.0f, 0.0f)};

  restraint.clear();
  restraint.calcForce(xyzq.getDeviceArray().data(), true, true);
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));

  const std::vector<double3> expectedForces = {make_double3(2.0, 0.0, 0.0),
                                               make_double3(2.0, 0.0, 0.0),
                                               make_double3(0.0, 0.0, 0.0)};

  apo_test::CheckVectorsClose<double3>(
      "HMCM unweighted-after-masses forces",
      CopyFixedPointForceToHost(restraint.getForce()), expectedForces,
      FORCE_TOLERANCE);
  CHECK(GetHmcmEnergy(restraint) == Approx(2.0).margin(ENERGY_TOLERANCE));

  const std::vector<double> expectedVirial =
      ExpectedVirialFromCenterGradient({2.0, 0.0, 0.0}, {4.0, 0.0, 0.0});

  apo_test::CheckVectorsClose<double>(
      "HMCM unweighted-after-masses virial",
      GetHmcmVirial(restraint, xyzq, NUM_ATOMS, SMALL_BOX_DIMENSIONS),
      expectedVirial, VIRIAL_TOLERANCE);
}

TEST_CASE("HarmonicCenterOfMassRestraintForceReferenceDistance") {
  HmcmForce restraint(NUM_ATOMS);
  restraint.setBoxDimensions(SMALL_BOX_DIMENSIONS);
  restraint.setSelection(MakeSelection(NUM_ATOMS, {0, 1}));
  restraint.setForceConstant(3.0);
  restraint.setReferencePosition({0.0, 0.0, 0.0});
  restraint.setReferenceDistance(1.0);

  CudaContainer<float4> xyzq = std::vector<float4>{
      make_float4(3.0f, 0.0f, 0.0f, 0.0f), make_float4(5.0f, 0.0f, 0.0f, 0.0f),
      make_float4(9.0f, 9.0f, 9.0f, 0.0f)};

  restraint.clear();
  restraint.calcForce(xyzq.getDeviceArray().data(), true, true);
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));

  const std::vector<double3> expectedForces = {make_double3(9.0, 0.0, 0.0),
                                               make_double3(9.0, 0.0, 0.0),
                                               make_double3(0.0, 0.0, 0.0)};

  apo_test::CheckVectorsClose<double3>(
      "HMCM reference-distance forces",
      CopyFixedPointForceToHost(restraint.getForce()), expectedForces,
      FORCE_TOLERANCE);
  CHECK(GetHmcmEnergy(restraint) == Approx(27.0).margin(ENERGY_TOLERANCE));

  const std::vector<double> expectedVirial =
      ExpectedVirialFromCenterGradient({4.0, 0.0, 0.0}, {18.0, 0.0, 0.0});

  apo_test::CheckVectorsClose<double>(
      "HMCM reference-distance virial",
      GetHmcmVirial(restraint, xyzq, NUM_ATOMS, SMALL_BOX_DIMENSIONS),
      expectedVirial, VIRIAL_TOLERANCE);
}

TEST_CASE("HarmonicCenterOfMassRestraintForcePeriodicWrappingAndVirial") {
  HmcmForce restraint(NUM_ATOMS);
  restraint.setBoxDimensions({10.0, 10.0, 10.0});
  restraint.setSelection(MakeSelection(NUM_ATOMS, {0, 1}));
  restraint.setForceConstant(2.0);
  restraint.setReferencePosition({1.1, 0.0, 0.0});

  CudaContainer<float4> xyzq = std::vector<float4>{
      make_float4(9.8f, 0.0f, 0.0f, 0.0f), make_float4(0.4f, 0.0f, 0.0f, 0.0f),
      make_float4(5.0f, 5.0f, 5.0f, 0.0f)};

  restraint.clear();
  restraint.calcForce(xyzq.getDeviceArray().data(), true, true);
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));

  const std::vector<double3> expectedForces = {make_double3(-2.0, 0.0, 0.0),
                                               make_double3(-2.0, 0.0, 0.0),
                                               make_double3(0.0, 0.0, 0.0)};

  apo_test::CheckVectorsClose<double3>(
      "HMCM periodic-wrapping forces",
      CopyFixedPointForceToHost(restraint.getForce()), expectedForces,
      FORCE_TOLERANCE);
  CHECK(GetHmcmEnergy(restraint) == Approx(2.0).margin(ENERGY_TOLERANCE));

  const std::vector<double> expectedVirial =
      ExpectedVirialFromCenterGradient({0.1, 0.0, 0.0}, {-4.0, 0.0, 0.0});

  apo_test::CheckVectorsClose<double>(
      "HMCM periodic-wrapping virial",
      GetHmcmVirial(restraint, xyzq, NUM_ATOMS, {10.0, 10.0, 10.0}),
      expectedVirial, VIRIAL_TOLERANCE);
}

TEST_CASE("HarmonicCenterOfMassRestraintForceClearResetsForcesAndEnergy") {
  HmcmForce restraint(NUM_ATOMS);
  restraint.setBoxDimensions(SMALL_BOX_DIMENSIONS);
  restraint.setSelection(MakeSelection(NUM_ATOMS, {0, 1}));
  restraint.setReferencePosition({0.0, 0.0, 0.0});
  restraint.setForceConstant(1.0);

  CudaContainer<float4> xyzq = std::vector<float4>{
      make_float4(1.0f, 0.0f, 0.0f, 0.0f), make_float4(3.0f, 0.0f, 0.0f, 0.0f),
      make_float4(9.0f, 9.0f, 9.0f, 0.0f)};

  restraint.clear();
  restraint.calcForce(xyzq.getDeviceArray().data(), true, true);
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));

  CHECK(GetHmcmEnergy(restraint) > 0.0);

  restraint.clear();
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));

  const std::vector<double3> expectedForces = {make_double3(0.0, 0.0, 0.0),
                                               make_double3(0.0, 0.0, 0.0),
                                               make_double3(0.0, 0.0, 0.0)};

  apo_test::CheckVectorsClose<double3>(
      "HMCM cleared forces", CopyFixedPointForceToHost(restraint.getForce()),
      expectedForces, FORCE_TOLERANCE);
  CHECK(GetHmcmEnergy(restraint) == Approx(0.0).margin(ENERGY_TOLERANCE));

  const std::vector<double> expectedVirial =
      ExpectedVirialFromCenterGradient({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});

  apo_test::CheckVectorsClose<double>(
      "HMCM cleared virial",
      GetHmcmVirial(restraint, xyzq, NUM_ATOMS, SMALL_BOX_DIMENSIONS),
      expectedVirial, VIRIAL_TOLERANCE);
}

TEST_CASE("HarmonicCenterOfMassRestraintForceCanBeSubscribedToForceManager") {
  const std::string dataPath = getDataPath();

  auto prm =
      std::make_shared<CharmmParameters>(dataPath + "toppar_water_ions.str");
  auto psf = std::make_shared<CharmmPSF>(dataPath + "nacl_pair.psf");
  auto crd = std::make_shared<CharmmCrd>(dataPath + "nacl_pair.cor");

  auto ctx = std::make_shared<CharmmContext>(psf, prm);
  ctx->setBoxDimensions(BOX_DIMENSIONS);
  ctx->setCoordinates(crd);
  ctx->useHolonomicConstraints(false);
  ctx->setRandomSeedForVelocities(RANDOM_SEED);
  ctx->assignVelocitiesAtTemperature(TEMPERATURE);

  auto restraint = std::make_shared<HmcmForce>(ctx->getNumAtoms());
  restraint->setSelection(MakeSelection(ctx->getNumAtoms(), {0, 1}));
  restraint->setReferencePosition({1.41, 1.41, 1.41});
  restraint->setForceConstant(1.0);

  auto fm = ctx->getForceManager();
  CHECK_NOTHROW(fm->subscribe(restraint, "hmcm", restraint->getStream(),
                              restraint->getForce(),
                              restraint->getEnergyVirial()));

  auto integrator =
      std::make_shared<CudaLangevinThermostatIntegrator>(TIME_STEP);
  integrator->setReferenceTemperature(TEMPERATURE);
  integrator->setThermostatFriction(LANGEVIN_THERMOSTAT_FRICTION);
  integrator->setThermostatRngSeed(RANDOM_SEED);
  integrator->setCharmmContext(ctx);

  CHECK_NOTHROW(integrator->propagate(1));

  apo_test::CheckFiniteTemperature(ctx->computeTemperature());

  CHECK_NOTHROW(fm->unsubscribe(restraint));
}

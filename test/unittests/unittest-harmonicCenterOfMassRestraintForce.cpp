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

#include <cmath>
#include <filesystem>
#include <limits>
#include <memory>
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
  apo_test::CheckApoCharmmError([]() -> void { (void)HmcmForce(0); },
                                ApoCharmmErrorCode::InvalidArgument,
                                "Atom count must be positive; observed 0");

  apo_test::CheckApoCharmmError([]() -> void { (void)HmcmForce(-1); },
                                ApoCharmmErrorCode::InvalidArgument,
                                "Atom count must be positive; observed -1");

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
  const double infinity = std::numeric_limits<double>::infinity();

  SECTION("ForceConstantValidation") {
    apo_test::CheckApoCharmmError(
        [&restraint]() -> void { restraint.setForceConstant(-1.0); },
        ApoCharmmErrorCode::InvalidArgument,
        "Force constant must be non-negative; observed " +
            std::to_string(-1.0));

    apo_test::CheckApoCharmmError(
        [&restraint, infinity]() -> void {
          restraint.setForceConstant(infinity);
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Force constant must be finite; observed " + std::to_string(infinity));

    CHECK_NOTHROW(restraint.setForceConstant(0.0));
    CHECK_NOTHROW(restraint.setForceConstant(1.0));
  }

  SECTION("ReferenceDistanceValidation") {
    apo_test::CheckApoCharmmError(
        [&restraint]() -> void { restraint.setReferenceDistance(-1.0); },
        ApoCharmmErrorCode::InvalidArgument,
        "Reference distance must be non-negative; observed " +
            std::to_string(-1.0));

    apo_test::CheckApoCharmmError(
        [&restraint, infinity]() -> void {
          restraint.setReferenceDistance(infinity);
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Reference distance must be finite; observed " +
            std::to_string(infinity));

    CHECK_NOTHROW(restraint.setReferenceDistance(0.0));
    CHECK_NOTHROW(restraint.setReferenceDistance(1.0));
  }

  SECTION("ReferencePositionValidation") {
    apo_test::CheckApoCharmmError(
        [&restraint]() -> void {
          restraint.setReferencePosition({0.0, 0.0}, {1, 1, 1});
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Reference-position array size mismatch; expected 3, observed 2");

    apo_test::CheckApoCharmmError(
        [&restraint]() -> void {
          restraint.setReferencePosition({0.0, 0.0, 0.0}, {1, 1});
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Reference-mask array size mismatch; expected 3, observed 2");

    apo_test::CheckApoCharmmError(
        [&restraint]() -> void {
          restraint.setReferencePosition({0.0, 0.0, 0.0}, {0, 0, 0});
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Reference mask must activate at least one coordinate");

    apo_test::CheckApoCharmmError(
        [&restraint]() -> void {
          restraint.setReferencePosition({0.0, 0.0, 0.0}, {1, 2, 1});
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Reference mask at index 1 must be 0 or 1; observed 2");

    apo_test::CheckApoCharmmError(
        [&restraint, infinity]() -> void {
          restraint.setReferencePosition({0.0, infinity, 0.0}, {1, 1, 1});
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Reference position at index 1 must be finite; observed " +
            std::to_string(infinity));

    CHECK_NOTHROW(restraint.setReferencePosition({0.0, 0.0, 0.0}));
    CHECK_NOTHROW(restraint.setReferencePosition({0.0, 0.0, 0.0}, {1, 0, 1}));
  }

  SECTION("BoxDimensionValidation") {
    apo_test::CheckApoCharmmError(
        [&restraint]() -> void { restraint.setBoxDimensions({10.0, 10.0}); },
        ApoCharmmErrorCode::InvalidArgument,
        "Box-dimension array size mismatch; expected 3, observed 2");

    apo_test::CheckApoCharmmError(
        [&restraint, infinity]() -> void {
          restraint.setBoxDimensions({10.0, infinity, 10.0});
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Box dimension at index 1 must be finite; observed " +
            std::to_string(infinity));

    apo_test::CheckApoCharmmError(
        [&restraint]() -> void {
          restraint.setBoxDimensions({10.0, 0.0, 10.0});
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Box dimension at index 1 must be positive; observed " +
            std::to_string(0.0));

    apo_test::CheckApoCharmmError(
        [&restraint]() -> void {
          restraint.setBoxDimensions({10.0, -1.0, 10.0});
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Box dimension at index 1 must be positive; observed " +
            std::to_string(-1.0));

    CHECK_NOTHROW(restraint.setBoxDimensions({10.0, 10.0, 10.0}));
  }

  SECTION("SelectionValidation") {
    AtomSelection emptySelection(NUM_ATOMS, AtomSelection::InitialValue::NONE);
    apo_test::CheckApoCharmmError(
        [&restraint, &emptySelection]() -> void {
          restraint.setSelection(emptySelection);
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Selection must contain at least one atom; observed 0");

    AtomSelection wrongSizeSelection(NUM_ATOMS + 1,
                                     AtomSelection::InitialValue::ALL);
    apo_test::CheckApoCharmmError(
        [&restraint, &wrongSizeSelection]() -> void {
          restraint.setSelection(wrongSizeSelection);
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Selection atom count mismatch; expected 3, observed 4");

    CHECK_NOTHROW(restraint.setSelection(MakeSelection(NUM_ATOMS, {0})));
  }

  SECTION("MassValidation") {
    apo_test::CheckApoCharmmError(
        [&restraint]() -> void { restraint.setMasses({1.0}); },
        ApoCharmmErrorCode::InvalidArgument,
        "Mass array size mismatch; expected 3, observed 1");

    apo_test::CheckApoCharmmError(
        [&restraint, infinity]() -> void {
          restraint.setMasses({1.0, infinity, 1.0});
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Mass at index 1 must be finite; observed " + std::to_string(infinity));

    apo_test::CheckApoCharmmError(
        [&restraint]() -> void { restraint.setMasses({1.0, 1.0, -1.0}); },
        ApoCharmmErrorCode::InvalidArgument,
        "Mass at index 2 must be non-negative; observed " +
            std::to_string(-1.0));

    apo_test::CheckApoCharmmError(
        [&restraint]() -> void { restraint.setMasses({0.0, 0.0, 0.0}); },
        ApoCharmmErrorCode::InvalidArgument,
        "Selected atoms must have positive total weight; observed " +
            std::to_string(0.0));

    CHECK_NOTHROW(restraint.setMassWeighting(false));

    apo_test::CheckApoCharmmError(
        [&restraint]() -> void { restraint.setMassWeighting(true); },
        ApoCharmmErrorCode::InvalidArgument,
        "Selected atoms must have positive total weight; observed " +
            std::to_string(0.0));

    CHECK_NOTHROW(restraint.setMasses({1.0, 2.0, 3.0}));
  }

  SECTION("InitializeValidation") {
    CHECK_NOTHROW(restraint.initialize(NUM_ATOMS, {10.0, 10.0, 10.0}));

    apo_test::CheckApoCharmmError(
        [&restraint]() -> void {
          restraint.initialize(NUM_ATOMS + 1, {10.0, 10.0, 10.0});
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Initialization atom count mismatch; expected 3, observed 4");
  }
}

TEST_CASE("HarmonicCenterOfMassRestraintForceRequiresBoxForNonzeroForce") {
  HmcmForce restraint(NUM_ATOMS);
  restraint.setSelection(MakeSelection(NUM_ATOMS, {0, 1}));
  restraint.setForceConstant(1.0);
  restraint.setReferencePosition({0.0, 0.0, 0.0});

  CudaContainer<float4> xyzq = std::vector<float4>{
      make_float4(1.0f, 0.0f, 0.0f, 0.0f), make_float4(2.0f, 0.0f, 0.0f, 0.0f),
      make_float4(3.0f, 0.0f, 0.0f, 0.0f)};

  apo_test::CheckApoCharmmError(
      [&restraint, &xyzq]() -> void {
        restraint.calcForce(xyzq.getDeviceArray().data(), true, false);
      },
      ApoCharmmErrorCode::NotInitialized,
      "Box dimensions must be set before force evaluation");
}

TEST_CASE("HarmonicCenterOfMassRestraintForceChecksEveryKernelLaunchVariant") {
  HmcmForce restraint(NUM_ATOMS);
  restraint.setBoxDimensions(SMALL_BOX_DIMENSIONS);
  restraint.setSelection(MakeSelection(NUM_ATOMS, {0, 1}));
  restraint.setForceConstant(1.0);
  restraint.setReferencePosition({0.0, 0.0, 0.0});

  CudaContainer<float4> xyzq = std::vector<float4>{
      make_float4(1.0f, 0.0f, 0.0f, 0.0f), make_float4(2.0f, 0.0f, 0.0f, 0.0f),
      make_float4(3.0f, 0.0f, 0.0f, 0.0f)};

  const auto checkLaunches = [&restraint,
                              &xyzq](const bool calcEnergy,
                                     const bool calcVirial) -> void {
    restraint.clear();
    restraint.calcForce(xyzq.getDeviceArray().data(), calcEnergy, calcVirial);
    cudaCheck(cudaStreamSynchronize(*restraint.getStream()));
    return;
  };

  CHECK_NOTHROW(checkLaunches(true, true));
  CHECK_NOTHROW(checkLaunches(true, false));
  CHECK_NOTHROW(checkLaunches(false, true));
  CHECK_NOTHROW(checkLaunches(false, false));
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
  auto prm = std::make_shared<CharmmParameters>(apo_test::GetTopparDir() /
                                                "toppar_water_ions.str");
  auto psf =
      std::make_shared<CharmmPSF>(apo_test::GetDataDir() / "nacl_pair.psf");
  auto crd =
      std::make_shared<CharmmCrd>(apo_test::GetDataDir() / "nacl_pair.cor");

  auto ctx = std::make_shared<CharmmContext>(psf, prm);
  ctx->setBoxDimensions(BOX_DIMENSIONS);
  ctx->setCoordinates(crd);
  ctx->useHolonomicConstraints(false);
  ctx->setRandomSeed(RANDOM_SEED);
  ctx->assignVelocitiesAtTemperature(TEMPERATURE);

  auto restraint = std::make_shared<HmcmForce>(ctx->getNumAtoms());
  restraint->setSelection(MakeSelection(ctx->getNumAtoms(), {0, 1}));
  restraint->setReferencePosition({1.41, 1.41, 1.41});
  restraint->setForceConstant(1.0);

  auto fm = ctx->getForceManager();

  apo_test::CheckApoCharmmError(
      [&fm, &restraint]() -> void {
        fm->subscribe(restraint, "", restraint->getStream(),
                      restraint->getForce(), restraint->getEnergyVirial());
      },
      ApoCharmmErrorCode::InvalidArgument, "Force tag must not be empty");

  CHECK_NOTHROW(fm->subscribe(restraint, "hmcm", restraint->getStream(),
                              restraint->getForce(),
                              restraint->getEnergyVirial()));

  apo_test::CheckApoCharmmError(
      [&fm, &restraint]() -> void {
        fm->subscribe(restraint, "hmcm", restraint->getStream(),
                      restraint->getForce(), restraint->getEnergyVirial());
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Force is already subscribed to this ForceManager");

  auto integrator =
      std::make_shared<CudaLangevinThermostatIntegrator>(TIME_STEP);
  integrator->setReferenceTemperature(TEMPERATURE);
  integrator->setThermostatFriction(LANGEVIN_THERMOSTAT_FRICTION);
  integrator->setThermostatRngSeed(RANDOM_SEED);
  integrator->setCharmmContext(ctx);

  CHECK_NOTHROW(integrator->propagate(1));

  apo_test::CheckFiniteTemperature(ctx->computeTemperature());

  CHECK_NOTHROW(fm->unsubscribe(restraint));
  apo_test::CheckApoCharmmError(
      [&fm, &restraint]() -> void { fm->unsubscribe(restraint); },
      ApoCharmmErrorCode::InvalidArgument,
      "Force is not subscribed to this ForceManager");
}

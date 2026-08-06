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
#include "CudaLangevinThermostatIntegrator.h"
#include "ForceManager.h"
#include "HarmonicRestraintForce.h"
#include "apo_test_helpers.h"
#include "catch.hpp"
#include "test_paths.h"

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <vector_types.h>

namespace {

constexpr int NUM_ATOMS = 3;
constexpr int RANDOM_SEED = 314159;
constexpr double TEMPERATURE = 300.0;
constexpr double TIME_STEP = 0.001;
constexpr double LANGEVIN_THERMOSTAT_FRICTION = 0.0;
constexpr double FORCE_TOLERANCE = 0.0;
constexpr double ENERGY_TOLERANCE = 0.0;
constexpr double INV_FORCE_SCALE_TEST = 1.0 / static_cast<double>(1LL << 40);

const std::vector<double> BOX_DIMENSIONS = {50.0, 50.0, 50.0};

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

double
GetHarmonicEnergy(HarmonicRestraintForce<long long int, float> &restraint) {
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));
  restraint.getEnergyVirial()->copyToHost(*restraint.getStream());
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));
  return restraint.getEnergyVirial()->getEnergy("harm");
}

} // namespace

TEST_CASE("HarmonicRestraintForceConstructorAndDefaults") {
  HarmonicRestraintForce<long long int, float> restraint(NUM_ATOMS);

  CHECK(decltype(restraint)::contributesVirial == false);
  CHECK(restraint.getStream() != nullptr);
  CHECK(restraint.getForce() != nullptr);
  CHECK(restraint.getForce()->size() == NUM_ATOMS);
  CHECK(restraint.getForce()->stride() >= NUM_ATOMS);
  CHECK(restraint.getEnergyVirial() != nullptr);

  CudaContainer<float4> xyzq = std::vector<float4>{
      make_float4(1.0f, 2.0f, 3.0f, 0.0f), make_float4(4.0f, 5.0f, 6.0f, 0.0f),
      make_float4(7.0f, 8.0f, 9.0f, 0.0f)};

  restraint.clear();
  restraint.calcForce(xyzq.getDeviceArray().data(), true, false);
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));

  const std::vector<double3> expectedForces = {make_double3(0.0, 0.0, 0.0),
                                               make_double3(0.0, 0.0, 0.0),
                                               make_double3(0.0, 0.0, 0.0)};

  apo_test::CheckVectorsClose<double3>(
      CopyFixedPointForceToHost(restraint.getForce()), expectedForces,
      FORCE_TOLERANCE);
  CHECK(GetHarmonicEnergy(restraint) == Approx(0.0).margin(ENERGY_TOLERANCE));
}

TEST_CASE("HarmonicRestraintForceConstructionAndInitializationValidation") {
  apo_test::CheckApoCharmmError(
      []() -> void { (void)HarmonicRestraintForce<long long int, float>(0); },
      ApoCharmmErrorCode::InvalidArgument,
      "Atom count must be positive; observed 0");

  apo_test::CheckApoCharmmError(
      []() -> void { (void)HarmonicRestraintForce<long long int, float>(-1); },
      ApoCharmmErrorCode::InvalidArgument,
      "Atom count must be positive; observed -1");

  HarmonicRestraintForce<long long int, float> restraint(NUM_ATOMS);

  CudaContainer<float4> xyzq = std::vector<float4>{
      make_float4(0.0f, 0.0f, 0.0f, 0.0f), make_float4(0.0f, 0.0f, 0.0f, 0.0f),
      make_float4(0.0f, 0.0f, 0.0f, 0.0f)};

  CHECK_NOTHROW(restraint.setBoxDimensions(BOX_DIMENSIONS));
  CHECK_NOTHROW(
      restraint.calcForce(xyzq.getDeviceArray().data(), false, false));
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));

  apo_test::CheckApoCharmmError(
      [&restraint]() -> void {
        restraint.initialize(NUM_ATOMS + 1, BOX_DIMENSIONS);
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Initialization atom count mismatch; expected 3, observed 4");

  CHECK_NOTHROW(restraint.initialize(NUM_ATOMS, BOX_DIMENSIONS));
}

TEST_CASE("HarmonicRestraintForceSelectionAndForceConstantValidation") {
  HarmonicRestraintForce<long long int, float> restraint(NUM_ATOMS);
  const double infinity = std::numeric_limits<double>::infinity();

  AtomSelection wrongSelection(NUM_ATOMS + 1);
  apo_test::CheckApoCharmmError(
      [&restraint, &wrongSelection]() -> void {
        restraint.setSelection(wrongSelection);
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Selection atom count mismatch; expected 3, observed 4");

  apo_test::CheckApoCharmmError(
      [&restraint, infinity]() -> void {
        restraint.setForceConstant(infinity);
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Force constant must be finite; observed " + std::to_string(infinity));

  apo_test::CheckApoCharmmError(
      [&restraint]() -> void { restraint.setForceConstant(-1.0); },
      ApoCharmmErrorCode::InvalidArgument,
      "Force constant must be non-negative; observed " + std::to_string(-1.0));

  CHECK_NOTHROW(restraint.setForceConstant(0.0));
  CHECK_NOTHROW(restraint.setForceConstant(1.0));

  apo_test::CheckApoCharmmError(
      [&restraint]() -> void { restraint.setForceConstants({1.0, 2.0}); },
      ApoCharmmErrorCode::InvalidArgument,
      "Force-constant array size mismatch; expected 3, observed 2");

  apo_test::CheckApoCharmmError(
      [&restraint]() -> void {
        restraint.setForceConstants({1.0, 2.0, 3.0, 4.0});
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Force-constant array size mismatch; expected 3, observed 4");

  apo_test::CheckApoCharmmError(
      [&restraint, infinity]() -> void {
        restraint.setForceConstants({1.0, infinity, 3.0});
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Force constant at index 1 must be finite; observed " +
          std::to_string(infinity));

  apo_test::CheckApoCharmmError(
      [&restraint]() -> void { restraint.setForceConstants({1.0, 2.0, -1.0}); },
      ApoCharmmErrorCode::InvalidArgument,
      "Force constant at index 2 must be non-negative; observed " +
          std::to_string(-1.0));

  CHECK_NOTHROW(restraint.setForceConstants({0.0, 1.0, 2.0}));
}

TEST_CASE("HarmonicRestraintForceReferenceCoordinateAndMassValidation") {
  HarmonicRestraintForce<long long int, float> restraint(NUM_ATOMS);
  const double infinity = std::numeric_limits<double>::infinity();

  apo_test::CheckApoCharmmError(
      [&restraint]() -> void {
        restraint.setReferenceCoordinates(std::vector<double3>{
            make_double3(0.0, 0.0, 0.0), make_double3(0.0, 0.0, 0.0)});
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Reference-coordinate array size mismatch; expected 3, observed 2");

  apo_test::CheckApoCharmmError(
      [&restraint]() -> void {
        restraint.setReferenceCoordinates(
            std::vector<std::vector<double>>{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}});
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Reference-coordinate array size mismatch; expected 3, observed 2");

  apo_test::CheckApoCharmmError(
      [&restraint]() -> void {
        restraint.setReferenceCoordinates(std::vector<std::vector<double>>{
            {0.0, 0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0, 0.0}});
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Reference coordinate at atom index 1 has invalid size; expected 3, "
      "observed 2");

  apo_test::CheckApoCharmmError(
      [&restraint, infinity]() -> void {
        restraint.setReferenceCoordinates(std::vector<double3>{
            make_double3(0.0, 0.0, 0.0), make_double3(0.0, infinity, 0.0),
            make_double3(0.0, 0.0, 0.0)});
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Reference coordinate at atom index 1, Y component must be finite; "
      "observed " +
          std::to_string(infinity));

  apo_test::CheckApoCharmmError(
      [&restraint, infinity]() -> void {
        restraint.setReferenceCoordinates(std::vector<std::vector<double>>{
            {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, infinity}});
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Reference coordinate at atom index 2, Z component must be finite; "
      "observed " +
          std::to_string(infinity));

  CHECK_NOTHROW(
      restraint.setReferenceCoordinates(std::vector<std::vector<double>>{
          {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}}));

  apo_test::CheckApoCharmmError(
      [&restraint]() -> void { restraint.setMasses({1.0, 2.0}); },
      ApoCharmmErrorCode::InvalidArgument,
      "Mass array size mismatch; expected 3, observed 2");

  apo_test::CheckApoCharmmError(
      [&restraint, infinity]() -> void {
        restraint.setMasses({1.0, infinity, 3.0});
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Mass at index 1 must be finite; observed " + std::to_string(infinity));

  apo_test::CheckApoCharmmError(
      [&restraint]() -> void { restraint.setMasses({1.0, -1.0, 3.0}); },
      ApoCharmmErrorCode::InvalidArgument,
      "Mass at index 1 must be non-negative; observed " + std::to_string(-1.0));

  CHECK_NOTHROW(restraint.setMasses({0.0, 1.0, 2.0}));
}

TEST_CASE("HarmonicRestraintForceBoxDimensionValidation") {
  HarmonicRestraintForce<long long int, float> restraint(NUM_ATOMS);
  const double infinity = std::numeric_limits<double>::infinity();

  apo_test::CheckApoCharmmError(
      [&restraint]() -> void { restraint.setBoxDimensions({50.0, 50.0}); },
      ApoCharmmErrorCode::InvalidArgument,
      "Box-dimension array size mismatch; expected 3, observed 2");

  apo_test::CheckApoCharmmError(
      [&restraint]() -> void {
        restraint.setBoxDimensions({50.0, 50.0, 50.0, 50.0});
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Box-dimension array size mismatch; expected 3, observed 4");

  apo_test::CheckApoCharmmError(
      [&restraint, infinity]() -> void {
        restraint.setBoxDimensions({50.0, infinity, 50.0});
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Box dimension at index 1 must be finite; observed " +
          std::to_string(infinity));

  apo_test::CheckApoCharmmError(
      [&restraint]() -> void { restraint.setBoxDimensions({50.0, 0.0, 50.0}); },
      ApoCharmmErrorCode::InvalidArgument,
      "Box dimension at index 1 must be positive; observed " +
          std::to_string(0.0));

  apo_test::CheckApoCharmmError(
      [&restraint]() -> void {
        restraint.setBoxDimensions({50.0, -1.0, 50.0});
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Box dimension at index 1 must be positive; observed " +
          std::to_string(-1.0));

  CHECK_NOTHROW(restraint.setBoxDimensions(BOX_DIMENSIONS));
}

TEST_CASE("HarmonicRestraintForceForceManagerSubscriptionValidation") {
  ForceManager forceManager;
  auto restraint =
      std::make_shared<HarmonicRestraintForce<long long int, float>>(NUM_ATOMS);

  apo_test::CheckApoCharmmError(
      [&forceManager, &restraint]() -> void {
        forceManager.subscribe(restraint, "", restraint->getStream(),
                               restraint->getForce(),
                               restraint->getEnergyVirial());
      },
      ApoCharmmErrorCode::InvalidArgument, "Force tag must not be empty");

  CHECK_NOTHROW(forceManager.subscribe(
      restraint, "harm", restraint->getStream(), restraint->getForce(),
      restraint->getEnergyVirial()));

  apo_test::CheckApoCharmmError(
      [&forceManager, &restraint]() -> void {
        forceManager.subscribe(restraint, "harm", restraint->getStream(),
                               restraint->getForce(),
                               restraint->getEnergyVirial());
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Force is already subscribed to this ForceManager");

  CHECK_NOTHROW(forceManager.unsubscribe(restraint));

  apo_test::CheckApoCharmmError(
      [&forceManager, &restraint]() -> void {
        forceManager.unsubscribe(restraint);
      },
      ApoCharmmErrorCode::InvalidArgument,
      "Force is not subscribed to this ForceManager");
}

TEST_CASE("HarmonicRestraintForceCalculatesForceAndEnergy") {
  const AtomSelection selection = MakeSelection(NUM_ATOMS, {0, 1});

  HarmonicRestraintForce<long long int, float> restraint(NUM_ATOMS);
  restraint.setSelection(selection);
  restraint.setForceConstants({2.0, 0.5, 99.0});
  restraint.setMasses({1.5, 2.0, 5.0});
  restraint.setReferenceCoordinates(std::vector<double3>{
      make_double3(0.0, 0.0, 0.0), make_double3(1.0, 1.0, 1.0),
      make_double3(10.0, 10.0, 10.0)});
  restraint.initialize(NUM_ATOMS, BOX_DIMENSIONS);

  CudaContainer<float4> xyzq =
      std::vector<float4>{make_float4(1.0f, -2.0f, 0.5f, 0.0f),
                          make_float4(2.0f, 3.0f, -1.0f, 0.0f),
                          make_float4(9.0f, 8.0f, 7.0f, 0.0f)};

  restraint.clear();
  restraint.calcForce(xyzq.getDeviceArray().data(), true, false);
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));

  const std::vector<double3> expectedForces = {make_double3(6.0, -12.0, 3.0),
                                               make_double3(2.0, 4.0, -4.0),
                                               make_double3(0.0, 0.0, 0.0)};

  apo_test::CheckVectorsClose<double3>(
      CopyFixedPointForceToHost(restraint.getForce()), expectedForces,
      FORCE_TOLERANCE);

  constexpr double expectedEnergy = 24.75;
  CHECK(GetHarmonicEnergy(restraint) ==
        Approx(expectedEnergy).margin(ENERGY_TOLERANCE));
}

TEST_CASE("HarmonicRestraintForceSelectionMakesScalarForceConstant") {
  HarmonicRestraintForce<long long int, float> restraint(NUM_ATOMS);

  restraint.setReferenceCoordinates(std::vector<double3>{
      make_double3(0.0, 0.0, 0.0), make_double3(0.0, 0.0, 0.0),
      make_double3(0.0, 0.0, 0.0)});
  restraint.setSelection(MakeSelection(NUM_ATOMS, {1}));
  restraint.setForceConstant(3.0);

  CudaContainer<float4> xyzq = std::vector<float4>{
      make_float4(1.0f, 0.0f, 0.0f, 0.0f), make_float4(0.0f, 2.0f, 0.0f, 0.0f),
      make_float4(0.0f, 0.0f, 3.0f, 0.0f)};

  restraint.clear();
  restraint.calcForce(xyzq.getDeviceArray().data(), true, false);
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));

  const std::vector<double3> expectedForces = {make_double3(0.0, 0.0, 0.0),
                                               make_double3(0.0, 12.0, 0.0),
                                               make_double3(0.0, 0.0, 0.0)};

  apo_test::CheckVectorsClose<double3>(
      CopyFixedPointForceToHost(restraint.getForce()), expectedForces,
      FORCE_TOLERANCE);

  constexpr double expectedEnergy = 12.0;
  CHECK(GetHarmonicEnergy(restraint) ==
        Approx(expectedEnergy).margin(ENERGY_TOLERANCE));
}

TEST_CASE("HarmonicRestraintForceClearResetsForcesAndEnergy") {
  HarmonicRestraintForce<long long int, float> restraint(NUM_ATOMS);

  restraint.setForceConstant(1.0);
  restraint.setReferenceCoordinates(std::vector<double3>{
      make_double3(0.0, 0.0, 0.0), make_double3(0.0, 0.0, 0.0),
      make_double3(0.0, 0.0, 0.0)});

  CudaContainer<float4> xyzq = std::vector<float4>{
      make_float4(1.0f, 0.0f, 0.0f, 0.0f), make_float4(0.0f, 1.0f, 0.0f, 0.0f),
      make_float4(0.0f, 0.0f, 1.0f, 0.0f)};

  restraint.clear();
  restraint.calcForce(xyzq.getDeviceArray().data(), true, false);
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));

  restraint.clear();
  cudaCheck(cudaStreamSynchronize(*restraint.getStream()));

  const std::vector<double3> expectedForces = {make_double3(0.0, 0.0, 0.0),
                                               make_double3(0.0, 0.0, 0.0),
                                               make_double3(0.0, 0.0, 0.0)};

  apo_test::CheckVectorsClose<double3>(
      CopyFixedPointForceToHost(restraint.getForce()), expectedForces,
      FORCE_TOLERANCE);

  constexpr double expectedEnergy = 0.0;
  CHECK(GetHarmonicEnergy(restraint) ==
        Approx(expectedEnergy).margin(ENERGY_TOLERANCE));
}

TEST_CASE("HarmonicRestraintForceCanBeSubscribedToForceManager") {
  const std::string dataPath = getDataPath();

  auto prm =
      std::make_shared<CharmmParameters>(dataPath + "toppar_water_ions.str");
  auto psf = std::make_shared<CharmmPSF>(dataPath + "nacl_pair.psf");
  auto crd = std::make_shared<CharmmCrd>(dataPath + "nacl_pair.cor");

  auto ctx = std::make_shared<CharmmContext>(psf, prm);
  ctx->setBoxDimensions(BOX_DIMENSIONS);
  ctx->setCoordinates(crd);
  ctx->useHolonomicConstraints(false);
  ctx->setRandomSeed(RANDOM_SEED);
  ctx->assignVelocitiesAtTemperature(TEMPERATURE);

  auto restraint =
      std::make_shared<HarmonicRestraintForce<long long int, float>>(
          ctx->getNumAtoms());
  restraint->setReferenceCoordinates(crd->getCoordinatesDP());
  restraint->setMasses(psf->getMasses());
  restraint->setSelection(MakeSelection(ctx->getNumAtoms(), {0, 1}));
  restraint->setForceConstant(10.0);

  auto fm = ctx->getForceManager();
  CHECK_NOTHROW(fm->subscribe(restraint, "harm", restraint->getStream(),
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

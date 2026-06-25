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
#include "CudaLangevinThermostatIntegrator.h"
#include "ForceManager.h"
#include "HarmonicRestraintForce.h"
#include "apo_test_helpers.h"
#include "catch.hpp"
#include "test_paths.h"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
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

TEST_CASE("HarmonicRestraintForceSettersAndValidation") {
  HarmonicRestraintForce<long long int, float> restraint(NUM_ATOMS);

  SECTION("SelectionSizeMustMatch") {
    AtomSelection wrongSelection(NUM_ATOMS + 1);

    CHECK_THROWS_AS(restraint.setSelection(wrongSelection),
                    std::invalid_argument);
  }

  SECTION("ForceConstantValidation") {
    CHECK_THROWS_AS(restraint.setForceConstant(-1.0), std::invalid_argument);
    CHECK_THROWS_AS(
        restraint.setForceConstant(std::numeric_limits<double>::infinity()),
        std::invalid_argument);
    CHECK_NOTHROW(restraint.setForceConstant(0.0));
    CHECK_NOTHROW(restraint.setForceConstant(1.0));
  }

  SECTION("ForceConstantsSizeMustMatch") {
    CHECK_THROWS_AS(restraint.setForceConstants({1.0, 2.0}),
                    std::invalid_argument);
    CHECK_NOTHROW(restraint.setForceConstants({1.0, 2.0, 3.0}));
    CHECK_THROWS_AS(restraint.setForceConstants({1.0, 2.0, 3.0, 4.0}),
                    std::invalid_argument);
  }

  SECTION("ReferenceCoordinateSizeMustMatch") {
    CHECK_THROWS_AS(restraint.setReferenceCoordinates(
                        std::vector<double4>{make_double4(0.0, 0.0, 0.0, 1.0),
                                             make_double4(0.0, 0.0, 0.0, 1.0)}),
                    std::invalid_argument);
    CHECK_NOTHROW(restraint.setReferenceCoordinates(std::vector<double4>{
        make_double4(0.0, 0.0, 0.0, 1.0), make_double4(0.0, 0.0, 0.0, 1.0),
        make_double4(0.0, 0.0, 0.0, 1.0)}));
    CHECK_THROWS_AS(restraint.setReferenceCoordinates(
                        std::vector<double4>{make_double4(0.0, 0.0, 0.0, 1.0),
                                             make_double4(0.0, 0.0, 0.0, 1.0),
                                             make_double4(0.0, 0.0, 0.0, 1.0),
                                             make_double4(0.0, 0.0, 0.0, 1.0)}),
                    std::invalid_argument);
  }

  SECTION("MassSizeMustMatch") {
    CHECK_THROWS_AS(restraint.setMasses({1.0, 2.0}), std::invalid_argument);
    CHECK_NOTHROW(restraint.setMasses({1.0, 2.0, 3.0}));
    CHECK_THROWS_AS(restraint.setMasses({1.0, 2.0, 3.0, 4.0}),
                    std::invalid_argument);
  }

  SECTION("InitializeValidatesAtomCountAndBoxDimensions") {
    CHECK_NOTHROW(restraint.initialize(NUM_ATOMS, BOX_DIMENSIONS));
    CHECK_THROWS_AS(restraint.initialize(NUM_ATOMS + 1, BOX_DIMENSIONS),
                    std::runtime_error);
    CHECK_THROWS_AS(restraint.setBoxDimensions({1.0, 2.0}),
                    std::invalid_argument);
    CHECK_THROWS_AS(restraint.setBoxDimensions({50.0, -1.0, 50.0}),
                    std::invalid_argument);
    CHECK_THROWS_AS(restraint.setBoxDimensions(
                        {50.0, std::numeric_limits<double>::infinity(), 50.0}),
                    std::invalid_argument);
  }
}

TEST_CASE("HarmonicRestraintForceCalculatesForceAndEnergy") {
  const AtomSelection selection = MakeSelection(NUM_ATOMS, {0, 1});

  HarmonicRestraintForce<long long int, float> restraint(NUM_ATOMS);
  restraint.setSelection(selection);
  restraint.setForceConstants({2.0, 0.5, 99.0});
  restraint.setMasses({1.5, 2.0, 5.0});
  restraint.setReferenceCoordinates(std::vector<double4>{
      make_double4(0.0, 0.0, 0.0, 1.0), make_double4(1.0, 1.0, 1.0, 2.0),
      make_double4(10.0, 10.0, 10.0, 3.0)});
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

  restraint.setReferenceCoordinates(std::vector<double4>{
      make_double4(0.0, 0.0, 0.0, 1.0), make_double4(0.0, 0.0, 0.0, 2.0),
      make_double4(0.0, 0.0, 0.0, 3.0)});
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
  restraint.setReferenceCoordinates(std::vector<double4>{
      make_double4(0.0, 0.0, 0.0, 0.0), make_double4(0.0, 0.0, 0.0, 0.0),
      make_double4(0.0, 0.0, 0.0, 0.0)});

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

  auto fm = std::make_shared<ForceManager>(psf, prm);
  fm->setBoxDimensions(BOX_DIMENSIONS);

  auto ctx = std::make_shared<CharmmContext>(fm);
  ctx->setCoordinates(crd);
  ctx->useHolonomicConstraints(false);
  ctx->setRandomSeedForVelocities(RANDOM_SEED);
  ctx->assignVelocitiesAtTemperature(TEMPERATURE);

  auto restraint =
      std::make_shared<HarmonicRestraintForce<long long int, float>>(
          ctx->getNumAtoms());
  restraint->setReferenceCoordinates(crd->getCoordinatesD());
  restraint->setMasses(psf->getMasses());
  restraint->setSelection(MakeSelection(ctx->getNumAtoms(), {0, 1}));
  restraint->setForceConstant(10.0);

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

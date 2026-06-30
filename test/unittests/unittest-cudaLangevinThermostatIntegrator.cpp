// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "CharmmContext.h"
#include "CharmmCrd.h"
#include "CharmmPSF.h"
#include "CharmmParameters.h"
#include "CudaLangevinThermostatIntegrator.h"
#include "ForceManager.h"
#include "apo_test_helpers.h"
#include "catch.hpp"
#include "test_paths.h"

#include <cstdint>

namespace {

constexpr int RANDOM_SEED = 314159;
constexpr int NUM_STEPS = 5;
constexpr double TIME_STEP = 0.001;
constexpr double TEMPERATURE = 300.0;
constexpr double REFERENCE_TEMPERATURE = 315;
constexpr double THERMOSTAT_FRICTION = 1.0;
constexpr std::uint64_t THERMOSTAT_SEED = 987654321ULL;
constexpr unsigned long long int RNG_SEQUENCE_POS = 8ULL;
constexpr double DETERMINISTIC_TOLERANCE = 0.0;
constexpr double TOLERANCE = 1.0e-10;
constexpr double FINITE_TEMPERATURE_UPPER_BOUND = 1.e8;

const std::vector<double> BOX_DIMENSIONS = {50.0, 50.0, 50.0};

void ConfigureIntegrator(
    std::shared_ptr<CudaLangevinThermostatIntegrator> integrator) {
  integrator->setReferenceTemperature(REFERENCE_TEMPERATURE);
  integrator->setThermostatFriction(THERMOSTAT_FRICTION);
  integrator->setThermostatRngSeed(THERMOSTAT_SEED);
  integrator->setRngSequencePos(0ULL);
  return;
}

void CheckIntegratorStateMatches(
    std::shared_ptr<CudaLangevinThermostatIntegrator> left,
    std::shared_ptr<CudaLangevinThermostatIntegrator> right,
    const double tolerance) {
  CHECK(left->getReferenceTemperature() ==
        Approx(right->getReferenceTemperature()).margin(tolerance));
  CHECK(left->getThermostatFriction() ==
        Approx(right->getThermostatFriction()).margin(tolerance));
  CHECK(left->getThermostatRngSeed() == right->getThermostatRngSeed());
  CHECK(left->getRngSequencePos() == right->getRngSequencePos());
  CHECK(left->getAverageWindowSize() == right->getAverageWindowSize());

  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(left->getKineticEnergy()),
      apo_test::CopyToHost<double>(right->getKineticEnergy()), tolerance);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(left->getAverageTemperature()),
      apo_test::CopyToHost<double>(right->getAverageTemperature()), tolerance);

  return;
}

} // namespace

TEST_CASE("CudaLangevinThermostatIntegratorConstructorDefaults") {
  CudaLangevinThermostatIntegrator integrator(TIME_STEP);

  CHECK(integrator.getTimeStep() == Approx(TIME_STEP));
  CHECK(integrator.getReferenceTemperature() == Approx(300.0));
  CHECK(integrator.getThermostatFriction() == Approx(0.0));
  CHECK(integrator.getRngSequencePos() == 0ULL);
  CHECK(integrator.getAverageWindowSize() == 0);
  CHECK(integrator.getKineticEnergy().size() == 2);
  CHECK(integrator.getAverageTemperature().size() == 2);
}

TEST_CASE("CudaLangevinThermostatIntegratorSettersAndReset") {
  CudaLangevinThermostatIntegrator integrator(TIME_STEP);

  integrator.setReferenceTemperature(REFERENCE_TEMPERATURE);
  integrator.setThermostatFriction(THERMOSTAT_FRICTION);
  integrator.setThermostatRngSeed(THERMOSTAT_SEED);
  integrator.setRngSequencePos(RNG_SEQUENCE_POS);

  CHECK(integrator.getReferenceTemperature() == Approx(REFERENCE_TEMPERATURE));
  CHECK(integrator.getThermostatFriction() == Approx(THERMOSTAT_FRICTION));
  CHECK(integrator.getThermostatRngSeed() == THERMOSTAT_SEED);
  CHECK(integrator.getRngSequencePos() == RNG_SEQUENCE_POS);

  integrator.getAverageTemperature().set(42.0);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(integrator.getAverageTemperature()),
      std::vector<double>{42.0, 42.0}, 0.0);

  integrator.resetAverageTemperature();

  CHECK(integrator.getAverageWindowSize() == 0);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(integrator.getAverageTemperature()),
      std::vector<double>{0.0, 0.0}, 0.0);
}

TEST_CASE("CudaLangevinThermostatIntegratorContextInitialization") {
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

  auto integrator =
      std::make_shared<CudaLangevinThermostatIntegrator>(TIME_STEP);

  ConfigureIntegrator(integrator);

  integrator->setCharmmContext(ctx);

  CHECK(integrator->getCharmmContext() == ctx);
  CHECK(integrator->getNumberOfAtoms() == 2);
  CHECK(integrator->getCoordsDelta().size() == 2);
  CHECK(integrator->getCoordsDeltaPrevious().size() == 2);
  CHECK(integrator->getRngSequencePos() == 0ULL);
  apo_test::CheckFiniteTemperature(ctx->computeTemperature());

  CHECK_THROWS_AS(integrator->setCharmmContext(ctx), std::invalid_argument);
}

TEST_CASE(
    "CudaLangevinThermostatIntegratorFrictionReinitializesBeforePropagation") {
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

  auto integrator =
      std::make_shared<CudaLangevinThermostatIntegrator>(TIME_STEP);
  integrator->setReferenceTemperature(REFERENCE_TEMPERATURE);
  integrator->setThermostatRngSeed(THERMOSTAT_SEED);
  integrator->setCharmmContext(ctx);

  CHECK_NOTHROW(integrator->setThermostatFriction(THERMOSTAT_FRICTION));
  CHECK(integrator->getThermostatFriction() == Approx(THERMOSTAT_FRICTION));
  CHECK(integrator->getTotNumSteps() == 0ULL);
}

TEST_CASE("CudaLangevinThermostatIntegratorShortPropagation") {
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

  auto integrator =
      std::make_shared<CudaLangevinThermostatIntegrator>(TIME_STEP);
  ConfigureIntegrator(integrator);
  integrator->setCharmmContext(ctx);
  integrator->resetAverageTemperature();

  const std::vector<double4> initialCoordinates =
      apo_test::CopyToHost<double4>(ctx->getCoordinatesCharges());

  integrator->propagate(NUM_STEPS);

  const std::vector<double4> finalCoordinates =
      apo_test::CopyToHost<double4>(ctx->getCoordinatesCharges());
  const std::vector<double> averageTemperature =
      apo_test::CopyToHost<double>(integrator->getAverageTemperature());
  const std::vector<double> kineticEnergy =
      apo_test::CopyToHost<double>(integrator->getKineticEnergy());

  CHECK(integrator->getCurrentPropagatedStep() == NUM_STEPS);
  CHECK(integrator->getTotNumSteps() ==
        static_cast<unsigned long long int>(NUM_STEPS));
  CHECK(integrator->getNumSteps() == NUM_STEPS);
  CHECK(integrator->getRngSequencePos() == 4ULL * NUM_STEPS);
  CHECK(integrator->getAverageWindowSize() == NUM_STEPS);

  apo_test::CheckFiniteTemperature(ctx->computeTemperature());
  apo_test::CheckFiniteTemperature(integrator->getInstantaneousTemperature());

  REQUIRE(averageTemperature.size() == 2);
  apo_test::CheckFiniteTemperature(averageTemperature[0]);
  apo_test::CheckFiniteTemperature(averageTemperature[1]);

  REQUIRE(kineticEnergy.size() == 2);
  apo_test::CheckFiniteNonnegativeScalar(kineticEnergy[0]);
  apo_test::CheckFiniteNonnegativeScalar(kineticEnergy[1]);

  REQUIRE(finalCoordinates.size() == initialCoordinates.size());
  CHECK(((finalCoordinates[0].x != initialCoordinates[0].x) ||
         (finalCoordinates[0].y != initialCoordinates[0].y) ||
         (finalCoordinates[0].z != initialCoordinates[0].z)));
}

TEST_CASE("CudaLangevinThermostatIntegratorDeterministicTrajectory") {
  const std::string dataPath = getDataPath();

  auto prm1 =
      std::make_shared<CharmmParameters>(dataPath + "toppar_water_ions.str");
  auto prm2 =
      std::make_shared<CharmmParameters>(dataPath + "toppar_water_ions.str");

  auto psf1 = std::make_shared<CharmmPSF>(dataPath + "nacl_pair.psf");
  auto psf2 = std::make_shared<CharmmPSF>(dataPath + "nacl_pair.psf");

  auto crd1 = std::make_shared<CharmmCrd>(dataPath + "nacl_pair.cor");
  auto crd2 = std::make_shared<CharmmCrd>(dataPath + "nacl_pair.cor");

  auto ctx1 = std::make_shared<CharmmContext>(psf1, prm1);
  ctx1->setBoxDimensions(BOX_DIMENSIONS);
  ctx1->setCoordinates(crd1);
  ctx1->useHolonomicConstraints(false);
  ctx1->setRandomSeedForVelocities(RANDOM_SEED);
  ctx1->assignVelocitiesAtTemperature(TEMPERATURE);

  auto ctx2 = std::make_shared<CharmmContext>(psf2, prm2);
  ctx2->setBoxDimensions(BOX_DIMENSIONS);
  ctx2->setCoordinates(crd2);
  ctx2->useHolonomicConstraints(false);
  ctx2->setRandomSeedForVelocities(RANDOM_SEED);
  ctx2->assignVelocitiesAtTemperature(TEMPERATURE);

  apo_test::CheckVectorsClose<double4>(
      apo_test::CopyToHost<double4>(ctx1->getCoordinatesCharges()),
      apo_test::CopyToHost<double4>(ctx2->getCoordinatesCharges()),
      DETERMINISTIC_TOLERANCE);
  apo_test::CheckVectorsClose<double4>(
      apo_test::CopyToHost<double4>(ctx1->getVelocityMass()),
      apo_test::CopyToHost<double4>(ctx2->getVelocityMass()),
      DETERMINISTIC_TOLERANCE);

  auto integrator1 =
      std::make_shared<CudaLangevinThermostatIntegrator>(TIME_STEP);
  ConfigureIntegrator(integrator1);
  integrator1->setCharmmContext(ctx1);

  auto integrator2 =
      std::make_shared<CudaLangevinThermostatIntegrator>(TIME_STEP);
  ConfigureIntegrator(integrator2);
  integrator2->setCharmmContext(ctx2);

  CheckIntegratorStateMatches(integrator1, integrator2,
                              DETERMINISTIC_TOLERANCE);

  const std::vector<double4> initialCoordinates1 =
      apo_test::CopyToHost<double4>(ctx1->getCoordinatesCharges());
  const std::vector<double4> initialCoordinates2 =
      apo_test::CopyToHost<double4>(ctx2->getCoordinatesCharges());

  integrator1->propagate(NUM_STEPS);
  integrator2->propagate(NUM_STEPS);

  const std::vector<double4> finalCoordinates1 =
      apo_test::CopyToHost<double4>(ctx1->getCoordinatesCharges());
  const std::vector<double4> finalCoordinates2 =
      apo_test::CopyToHost<double4>(ctx2->getCoordinatesCharges());

  apo_test::CheckVectorsClose<double4>(finalCoordinates1, finalCoordinates2,
                                       DETERMINISTIC_TOLERANCE);
  apo_test::CheckVectorsClose<double4>(
      apo_test::CopyToHost<double4>(ctx1->getVelocityMass()),
      apo_test::CopyToHost<double4>(ctx2->getVelocityMass()),
      DETERMINISTIC_TOLERANCE);
  CheckIntegratorStateMatches(integrator1, integrator2,
                              DETERMINISTIC_TOLERANCE);

  CHECK(integrator1->getCurrentPropagatedStep() == NUM_STEPS);
  CHECK(integrator2->getCurrentPropagatedStep() == NUM_STEPS);
  CHECK(integrator1->getRngSequencePos() == 4ULL * NUM_STEPS);
  CHECK(integrator2->getRngSequencePos() == 4ULL * NUM_STEPS);

  CHECK(((finalCoordinates1[0].x != initialCoordinates1[0].x) ||
         (finalCoordinates1[0].y != initialCoordinates1[0].y) ||
         (finalCoordinates1[0].z != initialCoordinates1[0].z)));
  CHECK(((finalCoordinates2[0].x != initialCoordinates2[0].x) ||
         (finalCoordinates2[0].y != initialCoordinates2[0].y) ||
         (finalCoordinates2[0].z != initialCoordinates2[0].z)));
}

TEST_CASE("CudaLangevinThermostatIntegratorRestartValidation") {
  auto integrator =
      std::make_shared<CudaLangevinThermostatIntegrator>(TIME_STEP);

  CHECK_THROWS_AS(integrator->initializeFromRestartFile("missing.rst"),
                  std::runtime_error);

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

  integrator->setCharmmContext(ctx);

  CHECK_THROWS_AS(integrator->initializeFromRestartFile("missing.rst"),
                  std::runtime_error);
}

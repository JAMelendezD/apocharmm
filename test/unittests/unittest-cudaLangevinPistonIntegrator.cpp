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
#include "CudaLangevinPistonIntegrator.h"
#include "ForceManager.h"
#include "apo_test_helpers.h"
#include "catch.hpp"
#include "test_paths.h"

#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

namespace {

constexpr int RANDOM_SEED = 314159;
constexpr int NUM_STEPS = 3;
constexpr double TIME_STEP = 0.001;
constexpr double TEMPERATURE = 300.0;
constexpr double REFERENCE_TEMPERATURE = 315.0;
constexpr double NOSE_HOOVER_PISTON_MASS = 500.0;
constexpr double LANGEVIN_PISTON_MASS = 5.0;
constexpr double LANGEVIN_PISTON_FRICTION = 5.0;
constexpr std::uint64_t LANGEVIN_PISTON_SEED = 987654321ULL;
constexpr int MAX_PREDICTOR_CORRECTOR_ITERATIONS = 4;
constexpr double DETERMINISTIC_TOLERANCE = 0.0;
constexpr double FLOAT_TOLERANCE = 1.0e-10;

const std::vector<double> BOX_DIMENSIONS = {50.0, 50.0, 50.0};
const std::vector<double> REFERENCE_PRESSURE_TENSOR = {1.0, 0.0, 0.0, 0.0, 1.0,
                                                       0.0, 0.0, 0.0, 1.0};

void ConfigureIntegrator(
    std::shared_ptr<CudaLangevinPistonIntegrator> integrator) {
  integrator->useNoseHooverThermostat(true);
  integrator->setReferenceTemperature(REFERENCE_TEMPERATURE);
  integrator->setNoseHooverPistonMass(NOSE_HOOVER_PISTON_MASS);
  integrator->setMaxPredictorCorrectorIterations(
      MAX_PREDICTOR_CORRECTOR_ITERATIONS);
  integrator->useOldTemperature(false);

  integrator->setReferencePressure(REFERENCE_PRESSURE_TENSOR);
  integrator->setConstantSurfaceTension(false);
  integrator->setCrystalType(CRYSTAL::CUBIC);
  integrator->setLangevinPistonMass({LANGEVIN_PISTON_MASS});
  integrator->setLangevinPistonFrictionSeed(LANGEVIN_PISTON_SEED);
  integrator->setRngSequencePos(0ULL);
  integrator->setLangevinPistonFriction(LANGEVIN_PISTON_FRICTION);

  return;
}

void CheckIntegratorStateMatches(
    std::shared_ptr<CudaLangevinPistonIntegrator> left,
    std::shared_ptr<CudaLangevinPistonIntegrator> right,
    const double tolerance) {
  CHECK(left->getReferenceTemperature() ==
        Approx(right->getReferenceTemperature()).margin(tolerance));
  CHECK(left->usingOldTemperature() == right->usingOldTemperature());
  CHECK(left->getMaxPredictorCorrectorIterations() ==
        right->getMaxPredictorCorrectorIterations());
  CHECK(left->getCrystalType() == right->getCrystalType());
  CHECK(left->getLangevinPistonFrictionSeed() ==
        right->getLangevinPistonFrictionSeed());
  CHECK(left->getRngSequencePos() == right->getRngSequencePos());

  CHECK(apo_test::CopyScalarToHost<double>(left->getNoseHooverPistonMass()) ==
        Approx(apo_test::CopyScalarToHost<double>(
                   right->getNoseHooverPistonMass()))
            .margin(tolerance));
  CHECK(
      apo_test::CopyScalarToHost<double>(left->getNoseHooverPistonVelocity()) ==
      Approx(apo_test::CopyScalarToHost<double>(
                 right->getNoseHooverPistonVelocity()))
          .margin(tolerance));
  CHECK(apo_test::CopyScalarToHost<double>(
            left->getNoseHooverPistonVelocityPrevious()) ==
        Approx(apo_test::CopyScalarToHost<double>(
                   right->getNoseHooverPistonVelocityPrevious()))
            .margin(tolerance));
  CHECK(apo_test::CopyScalarToHost<double>(left->getNoseHooverPistonForce()) ==
        Approx(apo_test::CopyScalarToHost<double>(
                   right->getNoseHooverPistonForce()))
            .margin(tolerance));
  CHECK(apo_test::CopyScalarToHost<double>(
            left->getNoseHooverPistonForcePrevious()) ==
        Approx(apo_test::CopyScalarToHost<double>(
                   right->getNoseHooverPistonForcePrevious()))
            .margin(tolerance));

  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(left->getReferencePressureTensor()),
      apo_test::CopyToHost<double>(right->getReferencePressureTensor()),
      tolerance);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(left->getLangevinPistonMass()),
      apo_test::CopyToHost<double>(right->getLangevinPistonMass()), tolerance);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(left->getLangevinPistonOnStepPosition()),
      apo_test::CopyToHost<double>(right->getLangevinPistonOnStepPosition()),
      tolerance);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(left->getLangevinPistonHalfStepPosition()),
      apo_test::CopyToHost<double>(right->getLangevinPistonHalfStepPosition()),
      tolerance);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(left->getLangevinPistonOnStepVelocity()),
      apo_test::CopyToHost<double>(right->getLangevinPistonOnStepVelocity()),
      tolerance);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(left->getLangevinPistonHalfStepVelocity()),
      apo_test::CopyToHost<double>(right->getLangevinPistonHalfStepVelocity()),
      tolerance);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(left->getLangevinPistonDeltaPosition()),
      apo_test::CopyToHost<double>(right->getLangevinPistonDeltaPosition()),
      tolerance);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(
          left->getLangevinPistonDeltaPositionPrevious()),
      apo_test::CopyToHost<double>(
          right->getLangevinPistonDeltaPositionPrevious()),
      tolerance);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(
          left->getLangevinPistonDeltaPositionPredicted()),
      apo_test::CopyToHost<double>(
          right->getLangevinPistonDeltaPositionPredicted()),
      tolerance);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(left->getLangevinPistonDeltaPressure()),
      apo_test::CopyToHost<double>(right->getLangevinPistonDeltaPressure()),
      tolerance);

  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(left->getKineticEnergy()),
      apo_test::CopyToHost<double>(right->getKineticEnergy()), tolerance);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(left->getAverageTemperature()),
      apo_test::CopyToHost<double>(right->getAverageTemperature()), tolerance);

  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(left->getInstantaneousPressureTensor()),
      apo_test::CopyToHost<double>(right->getInstantaneousPressureTensor()),
      tolerance);
  CHECK(apo_test::CopyScalarToHost<double>(
            left->getInstantaneousPressureScalar()) ==
        Approx(apo_test::CopyScalarToHost<double>(
                   right->getInstantaneousPressureScalar()))
            .margin(tolerance));
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(left->getAveragePressureTensor()),
      apo_test::CopyToHost<double>(right->getAveragePressureTensor()),
      tolerance);
  CHECK(apo_test::CopyScalarToHost<double>(left->getAveragePressureScalar()) ==
        Approx(apo_test::CopyScalarToHost<double>(
                   right->getAveragePressureScalar()))
            .margin(tolerance));

  return;
}

} // namespace

TEST_CASE("CudaLangevinPistonIntegratorConstructorDefaults") {
  auto integrator = std::make_shared<CudaLangevinPistonIntegrator>(TIME_STEP);

  CHECK(integrator->getTimeStep() == Approx(TIME_STEP));
  CHECK(integrator->getReferenceTemperature() == Approx(300.0));
  CHECK(integrator->usingOldTemperature() == false);
  CHECK(integrator->getMaxPredictorCorrectorIterations() == 3);
  CHECK(integrator->getCrystalType() == CRYSTAL::NONE);
  CHECK(integrator->getRngSequencePos() == 0ULL);

  CHECK(apo_test::CopyScalarToHost<double>(
            integrator->getNoseHooverPistonMass()) == Approx(-9999.9999));
  CHECK(apo_test::CopyScalarToHost<double>(
            integrator->getNoseHooverPistonVelocity()) == Approx(0.0));
  CHECK(apo_test::CopyScalarToHost<double>(
            integrator->getNoseHooverPistonVelocityPrevious()) == Approx(0.0));
  CHECK(apo_test::CopyScalarToHost<double>(
            integrator->getNoseHooverPistonForce()) == Approx(0.0));
  CHECK(apo_test::CopyScalarToHost<double>(
            integrator->getNoseHooverPistonForcePrevious()) == Approx(0.0));

  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost(integrator->getReferencePressureTensor()),
      REFERENCE_PRESSURE_TENSOR, 0.0);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost(integrator->getKineticEnergy()),
      std::vector<double>{0.0, 0.0}, 0.0);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost(integrator->getAverageTemperature()),
      std::vector<double>{0.0, 0.0}, 0.0);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost(integrator->getAveragePressureTensor()),
      std::vector<double>(9, 0.0), 0.0);
  CHECK(apo_test::CopyScalarToHost<double>(
            integrator->getAveragePressureScalar()) == Approx(0.0).margin(0.0));

  CHECK(integrator->getLangevinPistonMass().size() == 0);
}

TEST_CASE("CudaLangevinPistonIntegratorSettersAndCrystalType") {
  auto integrator = std::make_shared<CudaLangevinPistonIntegrator>(TIME_STEP);

  SECTION("RejectsPressureAndPistonSetterBeforeCrystalType") {
    CHECK_THROWS_AS(integrator->setLangevinPistonMass({1.0}),
                    std::runtime_error);
    CHECK_THROWS_AS(integrator->setLangevinPistonFriction(1.0),
                    std::runtime_error);
    CHECK_THROWS_AS(integrator->setReferencePressure({1.0, 0.0}),
                    std::invalid_argument);
    CHECK_THROWS_AS(integrator->setCrystalType(CRYSTAL::NONE),
                    std::invalid_argument);
  }

  SECTION("CubicCrystal") {
    integrator->setReferenceTemperature(REFERENCE_TEMPERATURE);
    integrator->setNoseHooverPistonMass(NOSE_HOOVER_PISTON_MASS);
    integrator->setNoseHooverPistonVelocity(0.125);
    integrator->setNoseHooverPistonVelocityPrevious(-0.25);
    integrator->setNoseHooverPistonForce(0.5);
    integrator->setNoseHooverPistonForcePrevious(-0.75);
    integrator->setMaxPredictorCorrectorIterations(
        MAX_PREDICTOR_CORRECTOR_ITERATIONS);
    integrator->useOldTemperature(true);
    integrator->setReferencePressure(REFERENCE_PRESSURE_TENSOR);
    integrator->setConstantSurfaceTension(false);
    integrator->setCrystalType(CRYSTAL::CUBIC);
    integrator->setLangevinPistonMass({LANGEVIN_PISTON_MASS});
    integrator->setLangevinPistonFrictionSeed(LANGEVIN_PISTON_SEED);
    integrator->setRngSequencePos(4ULL);
    integrator->setLangevinPistonFriction(LANGEVIN_PISTON_FRICTION);

    CHECK(integrator->getReferenceTemperature() ==
          Approx(REFERENCE_TEMPERATURE));
    CHECK(integrator->usingOldTemperature() == true);
    CHECK(integrator->getMaxPredictorCorrectorIterations() ==
          MAX_PREDICTOR_CORRECTOR_ITERATIONS);
    CHECK(integrator->getCrystalType() == CRYSTAL::CUBIC);
    CHECK(integrator->getLangevinPistonFrictionSeed() == LANGEVIN_PISTON_SEED);
    CHECK(integrator->getRngSequencePos() == 4ULL);

    CHECK(apo_test::CopyScalarToHost<double>(
              integrator->getNoseHooverPistonMass()) ==
          Approx(NOSE_HOOVER_PISTON_MASS));
    apo_test::CheckVectorsClose<double>(
        apo_test::CopyToHost<double>(integrator->getReferencePressureTensor()),
        REFERENCE_PRESSURE_TENSOR, 0.0);
    apo_test::CheckVectorsClose<double>(
        apo_test::CopyToHost<double>(integrator->getLangevinPistonMass()),
        std::vector<double>{LANGEVIN_PISTON_MASS}, 0.0);

    CHECK_THROWS_AS(integrator->setLangevinPistonMass({1.0, 2.0}),
                    std::invalid_argument);
  }

  SECTION("TetragonalCrystal") {
    integrator->setCrystalType(CRYSTAL::TETRAGONAL);
    integrator->setLangevinPistonMass({1.0, 2.0});

    CHECK(integrator->getCrystalType() == CRYSTAL::TETRAGONAL);
    apo_test::CheckVectorsClose<double>(
        apo_test::CopyToHost<double>(integrator->getLangevinPistonMass()),
        std::vector<double>{1.0, 2.0}, 0.0);

    CHECK(integrator->getLangevinPistonOnStepPosition().size() == 2);
    CHECK(integrator->getLangevinPistonDeltaPressure().size() == 2);
  }

  SECTION("OrthorhombicCrystal") {
    integrator->setCrystalType(CRYSTAL::ORTHORHOMBIC);
    integrator->setLangevinPistonMass({1.0, 2.0, 3.0});

    CHECK(integrator->getCrystalType() == CRYSTAL::ORTHORHOMBIC);
    apo_test::CheckVectorsClose<double>(
        apo_test::CopyToHost<double>(integrator->getLangevinPistonMass()),
        std::vector<double>{1.0, 2.0, 3.0}, 0.0);

    CHECK(integrator->getLangevinPistonOnStepPosition().size() == 3);
    CHECK(integrator->getLangevinPistonDeltaPressure().size() == 3);
  }
}

TEST_CASE("CudaLangevinPistonIntegratorResetAverages") {
  auto integrator = std::make_shared<CudaLangevinPistonIntegrator>(TIME_STEP);

  integrator->getAverageTemperature().set(42.0);
  integrator->getAveragePressureTensor().set(2.0);
  integrator->getAveragePressureScalar().set(3.0);

  integrator->resetAverages();

  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(integrator->getAverageTemperature()),
      std::vector<double>{0.0, 0.0}, 0.0);
  apo_test::CheckVectorsClose<double>(
      apo_test::CopyToHost<double>(integrator->getAveragePressureTensor()),
      std::vector<double>(9, 0.0), 0.0);
  CHECK(apo_test::CopyScalarToHost<double>(
            integrator->getAveragePressureScalar()) == Approx(0.0).margin(0.0));
}

TEST_CASE("CudaLangevinPistonIntegratorRequireCrystalTypeBeforeContext") {
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

  auto integrator = std::make_shared<CudaLangevinPistonIntegrator>(TIME_STEP);

  CHECK_THROWS_AS(integrator->setCharmmContext(ctx), std::runtime_error);
}

TEST_CASE("CudaLangevinPistonIntegratorContextInitialization") {
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

  auto integrator = std::make_shared<CudaLangevinPistonIntegrator>(TIME_STEP);
  ConfigureIntegrator(integrator);

  integrator->setCharmmContext(ctx);

  CHECK(integrator->getCharmmContext() == ctx);
  CHECK(integrator->getCoordsDelta().size() == 2);
  CHECK(integrator->getCoordsDeltaPrevious().size() == 2);
  CHECK(integrator->getLangevinPistonMass().size() == 1);
  CHECK(integrator->getRngSequencePos() == 0ULL);

  apo_test::CheckFiniteTemperature(ctx->computeTemperature());

  apo_test::CheckApoCharmmError(
      [&]() { integrator->setCharmmContext(ctx); },
      ApoCharmmErrorCode::InvalidArgument,
      "A CharmmContext object was already set for this CudaIntegrator.");
}

TEST_CASE("CudaLangevinPistonIntegratorAutoPistonMass") {
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

  auto integrator = std::make_shared<CudaLangevinPistonIntegrator>(TIME_STEP);
  integrator->setCrystalType(CRYSTAL::CUBIC);
  integrator->setLangevinPistonFrictionSeed(LANGEVIN_PISTON_SEED);
  integrator->setRngSequencePos(0ULL);
  integrator->setLangevinPistonFriction(0.0);
  integrator->setCharmmContext(ctx);

  constexpr double expectedMass = 0.02 * (22.9898 + 35.45);

  CHECK(
      apo_test::CopyScalarToHost<double>(integrator->getLangevinPistonMass()) ==
      Approx(expectedMass).margin(1.0e-4));
}

TEST_CASE("CudaLangevinPistonIntegratorShortPropagation") {
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

  auto integrator = std::make_shared<CudaLangevinPistonIntegrator>(TIME_STEP);
  ConfigureIntegrator(integrator);
  integrator->setCharmmContext(ctx);
  integrator->resetAverages();

  const std::vector<double4> initialCoordinates =
      apo_test::CopyToHost<double4>(ctx->getCoordinatesChargesDP());
  const std::vector<double> initialBoxDimensions = ctx->getBoxDimensions();

  integrator->propagate(NUM_STEPS);

  const std::vector<double4> finalCoordinates =
      apo_test::CopyToHost<double4>(ctx->getCoordinatesChargesDP());
  const std::vector<double> finalBoxDimensions = ctx->getBoxDimensions();

  CHECK(integrator->getCurrentPropagatedStep() == NUM_STEPS);
  CHECK(integrator->getTotNumSteps() ==
        static_cast<unsigned long long int>(NUM_STEPS));
  CHECK(integrator->getNumSteps() == NUM_STEPS);
  CHECK(
      integrator->getRngSequencePos() ==
      static_cast<unsigned long long int>(MAX_PREDICTOR_CORRECTOR_ITERATIONS) *
          NUM_STEPS);

  apo_test::CheckFiniteTemperature(ctx->computeTemperature());
  apo_test::CheckFiniteTemperature(integrator->getInstantaneousTemperature());

  apo_test::CheckFiniteNonnegativeVector(
      apo_test::CopyToHost<double>(integrator->getKineticEnergy()));
  apo_test::CheckFiniteNonnegativeVector(
      apo_test::CopyToHost<double>(integrator->getAverageTemperature()));
  apo_test::CheckFiniteVector(apo_test::CopyToHost<double>(
      integrator->getInstantaneousPressureTensor()));
  apo_test::CheckFiniteScalar(apo_test::CopyScalarToHost<double>(
      integrator->getInstantaneousPressureScalar()));
  apo_test::CheckFiniteVector(
      apo_test::CopyToHost<double>(integrator->getAveragePressureTensor()));
  apo_test::CheckFiniteScalar(apo_test::CopyScalarToHost<double>(
      integrator->getAveragePressureScalar()));

  REQUIRE(finalCoordinates.size() == initialCoordinates.size());
  CHECK(((finalCoordinates[0].x != initialCoordinates[0].x) ||
         (finalCoordinates[0].y != initialCoordinates[0].y) ||
         (finalCoordinates[0].z != initialCoordinates[0].z)));

  REQUIRE(finalBoxDimensions.size() == initialBoxDimensions.size());
  CHECK(std::isfinite(finalBoxDimensions[0]));
  CHECK(std::isfinite(finalBoxDimensions[1]));
  CHECK(std::isfinite(finalBoxDimensions[2]));
}

TEST_CASE("CudaLangevinPistonIntegratorDeterministicTrajectory") {
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
  ctx1->setRandomSeed(RANDOM_SEED);
  ctx1->assignVelocitiesAtTemperature(TEMPERATURE);

  auto ctx2 = std::make_shared<CharmmContext>(psf2, prm2);
  ctx2->setBoxDimensions(BOX_DIMENSIONS);
  ctx2->setCoordinates(crd2);
  ctx2->useHolonomicConstraints(false);
  ctx2->setRandomSeed(RANDOM_SEED);
  ctx2->assignVelocitiesAtTemperature(TEMPERATURE);

  apo_test::CheckVectorsClose<double4>(
      apo_test::CopyToHost<double4>(ctx1->getCoordinatesChargesDP()),
      apo_test::CopyToHost<double4>(ctx2->getCoordinatesChargesDP()),
      DETERMINISTIC_TOLERANCE);
  apo_test::CheckVectorsClose<double4>(
      apo_test::CopyToHost<double4>(ctx1->getVelocitiesInverseMasses()),
      apo_test::CopyToHost<double4>(ctx2->getVelocitiesInverseMasses()),
      DETERMINISTIC_TOLERANCE);

  auto integrator1 = std::make_shared<CudaLangevinPistonIntegrator>(TIME_STEP);
  ConfigureIntegrator(integrator1);
  integrator1->setCharmmContext(ctx1);

  auto integrator2 = std::make_shared<CudaLangevinPistonIntegrator>(TIME_STEP);
  ConfigureIntegrator(integrator2);
  integrator2->setCharmmContext(ctx2);

  CheckIntegratorStateMatches(integrator1, integrator2,
                              DETERMINISTIC_TOLERANCE);

  const std::vector<double4> initialCoordinates1 =
      apo_test::CopyToHost<double4>(ctx1->getCoordinatesChargesDP());
  const std::vector<double4> initialCoordinates2 =
      apo_test::CopyToHost<double4>(ctx2->getCoordinatesChargesDP());

  integrator1->propagate(NUM_STEPS);
  integrator2->propagate(NUM_STEPS);

  const std::vector<double4> finalCoordinates1 =
      apo_test::CopyToHost<double4>(ctx1->getCoordinatesChargesDP());
  const std::vector<double4> finalCoordinates2 =
      apo_test::CopyToHost<double4>(ctx2->getCoordinatesChargesDP());

  apo_test::CheckVectorsClose<double4>(finalCoordinates1, finalCoordinates2,
                                       DETERMINISTIC_TOLERANCE);
  apo_test::CheckVectorsClose<double4>(
      apo_test::CopyToHost<double4>(ctx1->getVelocitiesInverseMasses()),
      apo_test::CopyToHost<double4>(ctx2->getVelocitiesInverseMasses()),
      DETERMINISTIC_TOLERANCE);
  CheckIntegratorStateMatches(integrator1, integrator2,
                              DETERMINISTIC_TOLERANCE);

  CHECK(integrator1->getCurrentPropagatedStep() == NUM_STEPS);
  CHECK(integrator2->getCurrentPropagatedStep() == NUM_STEPS);
  CHECK(
      integrator1->getRngSequencePos() ==
      static_cast<unsigned long long int>(MAX_PREDICTOR_CORRECTOR_ITERATIONS) *
          NUM_STEPS);
  CHECK(
      integrator2->getRngSequencePos() ==
      static_cast<unsigned long long int>(MAX_PREDICTOR_CORRECTOR_ITERATIONS) *
          NUM_STEPS);

  CHECK(((finalCoordinates1[0].x != initialCoordinates1[0].x) ||
         (finalCoordinates1[0].y != initialCoordinates1[0].y) ||
         (finalCoordinates1[0].z != initialCoordinates1[0].z)));
  CHECK(((finalCoordinates2[0].x != initialCoordinates2[0].x) ||
         (finalCoordinates2[0].y != initialCoordinates2[0].y) ||
         (finalCoordinates2[0].z != initialCoordinates2[0].z)));
}

TEST_CASE("CudaLangevinPistonIntegratorRestartValidation") {
  auto integrator = std::make_shared<CudaLangevinPistonIntegrator>(TIME_STEP);

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
  ctx->setRandomSeed(RANDOM_SEED);
  ctx->assignVelocitiesAtTemperature(TEMPERATURE);

  ConfigureIntegrator(integrator);
  integrator->setCharmmContext(ctx);

  CHECK_THROWS_AS(integrator->initializeFromRestartFile("missing.rst"),
                  std::runtime_error);
}

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
#include "CudaIntegrator.h"
#include "CudaLangevinPistonIntegrator.h"
#include "CudaLangevinThermostatIntegrator.h"
#include "CudaNoseHooverIntegrator.h"
#include "ForceManager.h"
#include "RestartSubscriber.h"
#include "apo_test_helpers.h"
#include "catch.hpp"
#include "test_paths.h"

#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace {

constexpr int NUM_STEPS = 3;
constexpr int REPORT_FREQUENCY = NUM_STEPS;
constexpr int RANDOM_SEED = 314159;
constexpr double TIME_STEP = 0.001;
constexpr double REFERENCE_TEMPERATURE = 300.0;
constexpr double NOSE_HOOVER_PISTON_MASS = 300.0;
constexpr double LANGEVIN_THERMOSTAT_FRICTION = 1.0;
constexpr double LANGEVIN_PISTON_MASS = 5.0;
constexpr double LANGEVIN_PISTON_FRICTION = 5.0;
constexpr std::uint64_t LANGEVIN_RNG_SEED = 987654321ULL;
constexpr double TOLERANCE = 0.0;

const std::vector<double> BOX_DIMENSIONS = {50.0, 50.0, 50.0};
const std::vector<double> REFERENCE_PRESSURE_TENSOR = {1.0, 0.0, 0.0, 0.0, 1.0,
                                                       0.0, 0.0, 0.0, 1.0};

std::shared_ptr<CharmmContext> CreateContext(const bool assignVelocities) {
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

  if (assignVelocities == true) {
    ctx->setRandomSeedForVelocities(RANDOM_SEED);
    ctx->assignVelocitiesAtTemperature(REFERENCE_TEMPERATURE);
  }

  return ctx;
}

std::shared_ptr<CudaNoseHooverIntegrator>
CreateNoseHooverIntegrator(std::shared_ptr<CharmmContext> ctx) {
  auto integrator = std::make_shared<CudaNoseHooverIntegrator>(TIME_STEP);
  integrator->setReferenceTemperature(REFERENCE_TEMPERATURE);
  integrator->setNoseHooverPistonMass(NOSE_HOOVER_PISTON_MASS);
  integrator->useOldTemperature(false);
  integrator->setCharmmContext(ctx);
  return integrator;
}

std::shared_ptr<CudaLangevinThermostatIntegrator>
CreateLangevinThermostatIntegrator(std::shared_ptr<CharmmContext> ctx) {
  auto integrator =
      std::make_shared<CudaLangevinThermostatIntegrator>(TIME_STEP);
  integrator->setReferenceTemperature(REFERENCE_TEMPERATURE);
  integrator->setThermostatFriction(LANGEVIN_THERMOSTAT_FRICTION);
  integrator->setThermostatRngSeed(LANGEVIN_RNG_SEED);
  integrator->setRngSequencePos(0ULL);
  integrator->setCharmmContext(ctx);
  return integrator;
}

std::shared_ptr<CudaLangevinPistonIntegrator>
CreateLangevinPistonIntegrator(std::shared_ptr<CharmmContext> ctx) {
  auto integrator = std::make_shared<CudaLangevinPistonIntegrator>(TIME_STEP);
  integrator->useNoseHooverThermostat(true);
  integrator->setReferenceTemperature(REFERENCE_TEMPERATURE);
  integrator->setNoseHooverPistonMass(NOSE_HOOVER_PISTON_MASS);
  integrator->useOldTemperature(false);
  integrator->setReferencePressure(REFERENCE_PRESSURE_TENSOR);
  integrator->setConstantSurfaceTension(false);
  integrator->setCrystalType(CRYSTAL::CUBIC);
  integrator->setLangevinPistonMass({LANGEVIN_PISTON_MASS});
  integrator->setLangevinPistonFrictionSeed(LANGEVIN_RNG_SEED);
  integrator->setRngSequencePos(0ULL);
  integrator->setLangevinPistonFriction(LANGEVIN_PISTON_FRICTION); // *
  integrator->setCharmmContext(ctx);
  return integrator;
}

void WriteRestartFile(const std::string &fileName,
                      const std::shared_ptr<CharmmContext> &ctx,
                      const std::shared_ptr<CudaIntegrator> &integrator) {
  RestartSubscriber rst(fileName, REPORT_FREQUENCY);
  rst.setCharmmContext(ctx);
  rst.setIntegrator(integrator);
  rst.update();
  return;
}

void CheckRestartFileContainsRequiredSections(const std::string &fileName) {
  std::ifstream fin(fileName);
  REQUIRE(fin.good());

  const std::string contents((std::istreambuf_iterator<char>(fin)),
                             std::istreambuf_iterator<char>());

  CHECK(contents.find("REST") != std::string::npos);
  CHECK(contents.find("APO") != std::string::npos);
  CHECK(contents.find("!CRYSTAL PARAMETERS") != std::string::npos);
  CHECK(
      contents.find("!NATOM,NPRIV,NSTEP,NSAVC,NSAVV,JHSTRT,NDEGF,SEED,NSAVL") !=
      std::string::npos);
  CHECK(contents.find("!ENERGIES and STATISTICS") != std::string::npos);
  CHECK(contents.find("!XOLD, YOLD, ZOLD") != std::string::npos);
  CHECK(contents.find("!VX, VY, VZ") != std::string::npos);
  CHECK(contents.find("!X, Y, Z") != std::string::npos);

  return;
}

void CheckCharmmContextStateMatches(
    const std::shared_ptr<CharmmContext> &ctx1,
    const std::shared_ptr<CharmmContext> &ctx2) {
  apo_test::CheckVectorsClose<double4>(
      "Coordinates/charges",
      apo_test::CopyToHost<double4>(ctx1->getCoordinatesCharges()),
      apo_test::CopyToHost<double4>(ctx2->getCoordinatesCharges()), TOLERANCE);
  apo_test::CheckVectorsClose<double4>(
      "Velocities/masses",
      apo_test::CopyToHost<double4>(ctx1->getVelocityMass()),
      apo_test::CopyToHost<double4>(ctx2->getVelocityMass()), TOLERANCE);
  apo_test::CheckVectorsClose<double>("Box dimensions",
                                      ctx1->getBoxDimensions(),
                                      ctx2->getBoxDimensions(), TOLERANCE);
  return;
}

void CheckBaseIntegratorStateMatches(CudaIntegrator &integrator1,
                                     CudaIntegrator &integrator2) {
  CHECK(integrator1.getTotNumSteps() == integrator2.getTotNumSteps());
  CHECK(integrator1.getTotNumSteps() ==
        static_cast<unsigned long long int>(NUM_STEPS + 1));
  apo_test::CheckVectorsClose<double4>(
      "Coordinates delta",
      apo_test::CopyToHost<double4>(integrator1.getCoordsDelta()),
      apo_test::CopyToHost<double4>(integrator2.getCoordsDelta()), TOLERANCE);
  apo_test::CheckVectorsClose<double4>(
      "Coordinate delta previous",
      apo_test::CopyToHost<double4>(integrator1.getCoordsDeltaPrevious()),
      apo_test::CopyToHost<double4>(integrator2.getCoordsDeltaPrevious()),
      TOLERANCE);
  return;
}

void CheckNoseHooverStateMatches(CudaNoseHooverIntegrator &integrator1,
                                 CudaNoseHooverIntegrator &integrator2) {
  CheckBaseIntegratorStateMatches(integrator1, integrator2);

  CHECK(integrator1.getReferenceTemperature() ==
        Approx(integrator2.getReferenceTemperature()).margin(TOLERANCE));
  CHECK(integrator1.usingOldTemperature() == integrator2.usingOldTemperature());

  apo_test::CheckVectorsClose<double>(
      "Nose-Hoover piston velocity",
      apo_test::CopyToHost<double>(integrator1.getNoseHooverPistonVelocity()),
      apo_test::CopyToHost<double>(integrator2.getNoseHooverPistonVelocity()),
      TOLERANCE);
  apo_test::CheckVectorsClose<double>(
      "Nose-Hoover piston velocity previous",
      apo_test::CopyToHost<double>(
          integrator1.getNoseHooverPistonVelocityPrevious()),
      apo_test::CopyToHost<double>(
          integrator2.getNoseHooverPistonVelocityPrevious()),
      TOLERANCE);
  apo_test::CheckVectorsClose<double>(
      "Nose-Hoover piston force",
      apo_test::CopyToHost<double>(integrator1.getNoseHooverPistonForce()),
      apo_test::CopyToHost<double>(integrator2.getNoseHooverPistonForce()),
      TOLERANCE);
  apo_test::CheckVectorsClose<double>(
      "Nose-Hoover piston force previous",
      apo_test::CopyToHost<double>(
          integrator1.getNoseHooverPistonForcePrevious()),
      apo_test::CopyToHost<double>(
          integrator2.getNoseHooverPistonForcePrevious()),
      TOLERANCE);
  apo_test::CheckVectorsClose<double>(
      "Kinetic energy",
      apo_test::CopyToHost<double>(integrator1.getKineticEnergy()),
      apo_test::CopyToHost<double>(integrator2.getKineticEnergy()), TOLERANCE);

  return;
}

void CheckLangevinThermostatStateMatches(
    CudaLangevinThermostatIntegrator &integrator1,
    CudaLangevinThermostatIntegrator &integrator2) {
  CheckBaseIntegratorStateMatches(integrator1, integrator2);

  CHECK(integrator1.getReferenceTemperature() ==
        Approx(integrator2.getReferenceTemperature()).margin(TOLERANCE));
  CHECK(integrator1.getThermostatFriction() ==
        Approx(integrator2.getThermostatFriction()).margin(TOLERANCE));
  CHECK(integrator1.getThermostatRngSeed() ==
        integrator2.getThermostatRngSeed());
  CHECK(integrator1.getRngSequencePos() == integrator2.getRngSequencePos());

  apo_test::CheckVectorsClose<double>(
      "Kinetic energy",
      apo_test::CopyToHost<double>(integrator1.getKineticEnergy()),
      apo_test::CopyToHost<double>(integrator2.getKineticEnergy()), TOLERANCE);

  return;
}

void CheckLangevinPistonStateMatches(
    CudaLangevinPistonIntegrator &integrator1,
    CudaLangevinPistonIntegrator &integrator2) {
  CheckBaseIntegratorStateMatches(integrator1, integrator2);

  CHECK(integrator1.getReferenceTemperature() ==
        Approx(integrator2.getReferenceTemperature()).margin(TOLERANCE));
  CHECK(integrator1.usingOldTemperature() == integrator2.usingOldTemperature());
  CHECK(integrator1.getCrystalType() == integrator2.getCrystalType());
  CHECK(integrator1.getLangevinPistonFrictionSeed() ==
        integrator2.getLangevinPistonFrictionSeed());
  CHECK(integrator1.getRngSequencePos() == integrator2.getRngSequencePos());

  apo_test::CheckVectorsClose<double>(
      "Nose-Hoover piston velocity",
      apo_test::CopyToHost<double>(integrator1.getNoseHooverPistonVelocity()),
      apo_test::CopyToHost<double>(integrator2.getNoseHooverPistonVelocity()),
      TOLERANCE);
  apo_test::CheckVectorsClose<double>(
      "Nose-Hoover piston velocity previous",
      apo_test::CopyToHost<double>(
          integrator1.getNoseHooverPistonVelocityPrevious()),
      apo_test::CopyToHost<double>(
          integrator2.getNoseHooverPistonVelocityPrevious()),
      TOLERANCE);
  apo_test::CheckVectorsClose<double>(
      "Nose-Hoover piston force",
      apo_test::CopyToHost<double>(integrator1.getNoseHooverPistonForce()),
      apo_test::CopyToHost<double>(integrator2.getNoseHooverPistonForce()),
      TOLERANCE);
  apo_test::CheckVectorsClose<double>(
      "Nose-Hoover piston force previous",
      apo_test::CopyToHost<double>(
          integrator1.getNoseHooverPistonForcePrevious()),
      apo_test::CopyToHost<double>(
          integrator2.getNoseHooverPistonForcePrevious()),
      TOLERANCE);
  apo_test::CheckVectorsClose<double>(
      "Kinetic energy",
      apo_test::CopyToHost<double>(integrator1.getKineticEnergy()),
      apo_test::CopyToHost<double>(integrator2.getKineticEnergy()), TOLERANCE);

  apo_test::CheckVectorsClose<double>(
      "Langevin piston on step position",
      apo_test::CopyToHost<double>(
          integrator1.getLangevinPistonOnStepPosition()),
      apo_test::CopyToHost<double>(
          integrator2.getLangevinPistonOnStepPosition()),
      TOLERANCE);
  apo_test::CheckVectorsClose<double>(
      "Langevin piston half step position",
      apo_test::CopyToHost<double>(
          integrator1.getLangevinPistonHalfStepPosition()),
      apo_test::CopyToHost<double>(
          integrator2.getLangevinPistonHalfStepPosition()),
      TOLERANCE);
  apo_test::CheckVectorsClose<double>(
      "Langevin piston on step velocity",
      apo_test::CopyToHost<double>(
          integrator1.getLangevinPistonOnStepVelocity()),
      apo_test::CopyToHost<double>(
          integrator2.getLangevinPistonOnStepVelocity()),
      TOLERANCE);
  apo_test::CheckVectorsClose<double>(
      "Langevin piston half step velocity",
      apo_test::CopyToHost<double>(
          integrator1.getLangevinPistonHalfStepVelocity()),
      apo_test::CopyToHost<double>(
          integrator2.getLangevinPistonHalfStepVelocity()),
      TOLERANCE);
  apo_test::CheckVectorsClose<double>(
      "Langevin piston delta position",
      apo_test::CopyToHost<double>(
          integrator1.getLangevinPistonDeltaPosition()),
      apo_test::CopyToHost<double>(
          integrator2.getLangevinPistonDeltaPosition()),
      TOLERANCE);
  apo_test::CheckVectorsClose<double>(
      "Langevin piston delta position previous",
      apo_test::CopyToHost<double>(
          integrator1.getLangevinPistonDeltaPositionPrevious()),
      apo_test::CopyToHost<double>(
          integrator2.getLangevinPistonDeltaPositionPrevious()),
      TOLERANCE);
  apo_test::CheckVectorsClose<double>(
      "Instantaneous pressure tensor",
      apo_test::CopyToHost<double>(
          integrator1.getInstantaneousPressureTensor()),
      apo_test::CopyToHost<double>(
          integrator2.getInstantaneousPressureTensor()),
      TOLERANCE);
  CHECK(apo_test::CopyScalarToHost<double>(
            integrator1.getInstantaneousPressureScalar()) ==
        Approx(apo_test::CopyScalarToHost<double>(
                   integrator2.getInstantaneousPressureScalar()))
            .margin(TOLERANCE));

  return;
}

} // namespace

TEST_CASE("RestartSubscriberConstructionAndReportFrequency") {
  SECTION("DefaultConstruction") {
    RestartSubscriber rst;

    CHECK(rst.getFileName().empty() == true);
    CHECK(rst.getReportFrequency() == 1000);
  }

  SECTION("FileConstructor") {
    const std::string fileName = "tmpRestartSubscriberConstructor.rst";
    apo_test::RemoveIfExists(fileName);

    RestartSubscriber rst(fileName);

    CHECK(rst.getFileName() == fileName);
    CHECK(rst.getReportFrequency() == 1000);

    std::ifstream fin(fileName);
    CHECK(fin.good());

    apo_test::RemoveIfExists(fileName);
  }

  SECTION("FileAndReportFrequencyConstructor") {
    const std::string fileName = "tmpRestartSubscriberReportFrequency.rst";
    apo_test::RemoveIfExists(fileName);

    RestartSubscriber rst(fileName, REPORT_FREQUENCY);

    CHECK(rst.getFileName() == fileName);
    CHECK(rst.getReportFrequency() == REPORT_FREQUENCY);

    rst.setReportFrequency(1);
    CHECK(rst.getReportFrequency() == 1);

    apo_test::RemoveIfExists(fileName);
  }

  SECTION("RejectsMissingOutputDirectory") {
    CHECK_THROWS_AS(
        RestartSubscriber("missing_restart_subscriber_dir/output.rst"),
        std::invalid_argument);
  }
}

TEST_CASE("RestartSubscriberRejectsUnsupportedIntegrator") {
  const std::string fileName = "tmpUnsupportedRestartSubscriber.rst";
  apo_test::RemoveIfExists(fileName);

  auto ctx = CreateContext(true);
  auto integrator = std::make_shared<CudaIntegrator>(TIME_STEP);

  RestartSubscriber rst(fileName, REPORT_FREQUENCY);
  rst.setCharmmContext(ctx);
  rst.setIntegrator(integrator);

  CHECK_THROWS_AS(rst.update(), std::runtime_error);

  apo_test::RemoveIfExists(fileName);
}

TEST_CASE("RestartSubscriberRejectsDuplicateContextOrIntegrator") {
  const std::string fileName = "tmpDuplicateRestartSubscriber.rst";
  apo_test::RemoveIfExists(fileName);

  auto ctx = CreateContext(true);
  auto integrator = CreateNoseHooverIntegrator(ctx);

  RestartSubscriber rst(fileName, REPORT_FREQUENCY);

  rst.setCharmmContext(ctx);
  CHECK_THROWS_AS(rst.setCharmmContext(ctx), std::invalid_argument);

  rst.setIntegrator(integrator);
  CHECK_THROWS_AS(rst.setIntegrator(integrator), std::invalid_argument);

  apo_test::RemoveIfExists(fileName);
}

TEST_CASE("RestartSubscriberNoseHooverTrajectoryEquivalence") {
  const std::string fileName = "tmpRestartNoseHoover.rst";
  apo_test::RemoveIfExists(fileName);

  auto ctx1 = CreateContext(true);
  auto ctx2 = CreateContext(true);
  auto ctx3 = CreateContext(false);

  auto integrator1 = CreateNoseHooverIntegrator(ctx1);
  auto integrator2 = CreateNoseHooverIntegrator(ctx2);
  auto integrator3 = CreateNoseHooverIntegrator(ctx3);

  integrator1->propagate(NUM_STEPS + 1);

  integrator2->propagate(NUM_STEPS);
  WriteRestartFile(fileName, ctx2, integrator2);
  CheckRestartFileContainsRequiredSections(fileName);
  integrator2->propagate(1);

  integrator3->initializeFromRestartFile(fileName);
  integrator3->propagate(1);

  CheckCharmmContextStateMatches(ctx1, ctx2);
  CheckCharmmContextStateMatches(ctx1, ctx3);
  CheckCharmmContextStateMatches(ctx2, ctx3);

  CheckNoseHooverStateMatches(*integrator1, *integrator2);
  CheckNoseHooverStateMatches(*integrator1, *integrator3);
  CheckNoseHooverStateMatches(*integrator2, *integrator3);

  apo_test::RemoveIfExists(fileName);
}

TEST_CASE("RestartSubscriberLangevinThermostatTrajectoryEquivalence") {
  const std::string fileName = "tmpRestartLangevinThermostat.rst";
  apo_test::RemoveIfExists(fileName);

  auto ctx1 = CreateContext(true);
  auto ctx2 = CreateContext(true);
  auto ctx3 = CreateContext(false);

  auto integrator1 = CreateLangevinThermostatIntegrator(ctx1);
  auto integrator2 = CreateLangevinThermostatIntegrator(ctx2);
  auto integrator3 = CreateLangevinThermostatIntegrator(ctx3);

  integrator1->propagate(NUM_STEPS + 1);

  integrator2->propagate(NUM_STEPS);
  WriteRestartFile(fileName, ctx2, integrator2);
  CheckRestartFileContainsRequiredSections(fileName);
  integrator2->propagate(1);

  integrator3->initializeFromRestartFile(fileName);
  integrator3->propagate(1);

  CheckCharmmContextStateMatches(ctx1, ctx2);
  CheckCharmmContextStateMatches(ctx1, ctx3);
  CheckCharmmContextStateMatches(ctx2, ctx3);

  CheckLangevinThermostatStateMatches(*integrator1, *integrator2);
  CheckLangevinThermostatStateMatches(*integrator1, *integrator3);
  CheckLangevinThermostatStateMatches(*integrator2, *integrator3);

  apo_test::RemoveIfExists(fileName);
}

TEST_CASE("RestartSubscriberLangevinPistonTrajectoryEquivalence") {
  const std::string fileName = "tmpRestartLangevinPiston.rst";
  apo_test::RemoveIfExists(fileName);

  auto ctx1 = CreateContext(true);
  auto ctx2 = CreateContext(true);
  auto ctx3 = CreateContext(false);

  auto integrator1 = CreateLangevinPistonIntegrator(ctx1);
  auto integrator2 = CreateLangevinPistonIntegrator(ctx2);
  auto integrator3 = CreateLangevinPistonIntegrator(ctx3);

  integrator1->propagate(NUM_STEPS + 1);

  integrator2->propagate(NUM_STEPS);
  WriteRestartFile(fileName, ctx2, integrator2);
  CheckRestartFileContainsRequiredSections(fileName);
  integrator2->propagate(1);

  integrator3->initializeFromRestartFile(fileName);
  integrator3->propagate(1);

  CheckCharmmContextStateMatches(ctx1, ctx2);
  CheckCharmmContextStateMatches(ctx1, ctx3);
  CheckCharmmContextStateMatches(ctx2, ctx3);

  CheckLangevinPistonStateMatches(*integrator1, *integrator2);
  CheckLangevinPistonStateMatches(*integrator1, *integrator3);
  CheckLangevinPistonStateMatches(*integrator2, *integrator3);

  apo_test::RemoveIfExists(fileName);
}

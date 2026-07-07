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
#include "ForceManager.h"
#include "PBC.h"
#include "apo_test_helpers.h"
#include "catch.hpp"
#include "test_paths.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double TOLERANCE = 0.0;

const std::vector<double> BOX_DIMENSIONS = {40.0, 41.0, 42.0};

} // namespace

TEST_CASE("CharmmContextForceManagerConstructorMirrorsBackendState") {
  auto prm = std::make_shared<CharmmParameters>(getDataPath() +
                                                "toppar_water_ions.str");
  auto psf = std::make_shared<CharmmPSF>(getDataPath() + "nacl_pair.psf");

  auto fm = std::make_shared<ForceManager>(psf, prm);
  fm->setBoxDimensions(BOX_DIMENSIONS);

  std::cout << "auto ctx = std::make_shared<CharmmContext>(fm);" << std::endl;
  auto ctx = std::make_shared<CharmmContext>(fm);

  std::cout << "CHECK(ctx->getForceManager() == fm);" << std::endl;
  CHECK(ctx->getForceManager() == fm);
  CHECK(ctx->getPsf() == psf);
  CHECK(ctx->getPrm() == prm);

  apo_test::CheckVectorsClose<double>("CharmmContext box dimensions",
                                      ctx->getBoxDimensions(), BOX_DIMENSIONS,
                                      TOLERANCE);

  apo_test::CheckVectorsClose<double>("ForceManager box dimensions",
                                      fm->getBoxDimensions(), BOX_DIMENSIONS,
                                      TOLERANCE);

  CHECK(ctx->getPeriodicBoundaryCondition() ==
        fm->getPeriodicBoundaryCondition());

  CHECK(fm->isInitialized() == true);
}

TEST_CASE("CharmmContextStagedStateConfiguresForceManagerOnAttach") {
  auto prm = std::make_shared<CharmmParameters>(getDataPath() +
                                                "toppar_water_ions.str");
  auto psf = std::make_shared<CharmmPSF>(getDataPath() + "nacl_pair.psf");

  auto ctx = std::make_shared<CharmmContext>();
  ctx->setPsf(psf);
  ctx->setPrm(prm);
  ctx->setBoxDimensions(BOX_DIMENSIONS);
  ctx->setPeriodicBoundaryCondition(PBC::P21);

  auto fm = std::make_shared<ForceManager>();
  ctx->setForceManager(fm);

  CHECK(ctx->getForceManager() == fm);

  CHECK(ctx->getPsf() == psf);
  CHECK(ctx->getPrm() == prm);
  CHECK(fm->getPsf() == psf);
  CHECK(fm->getPrm() == prm);

  apo_test::CheckVectorsClose<double>("CharmmContext staged box dimensions",
                                      ctx->getBoxDimensions(), BOX_DIMENSIONS,
                                      TOLERANCE);

  apo_test::CheckVectorsClose<double>("ForceManager staged box dimensions",
                                      fm->getBoxDimensions(), BOX_DIMENSIONS,
                                      TOLERANCE);

  CHECK(ctx->getPeriodicBoundaryCondition() == PBC::P21);
  CHECK(fm->getPeriodicBoundaryCondition() == PBC::P21);
}

TEST_CASE("CharmmContextStagedStateCanLoadCoordinatesAfterBackendInitialize") {
  auto prm = std::make_shared<CharmmParameters>(getDataPath() +
                                                "toppar_water_ions.str");
  auto psf = std::make_shared<CharmmPSF>(getDataPath() + "nacl_pair.psf");
  auto crd = std::make_shared<CharmmCrd>(getDataPath() + "nacl_pair.cor");

  auto ctx = std::make_shared<CharmmContext>();
  ctx->setPsf(psf);
  ctx->setPrm(prm);
  ctx->setBoxDimensions(BOX_DIMENSIONS);
  ctx->setPeriodicBoundaryCondition(PBC::P1);

  auto fm = std::make_shared<ForceManager>();
  ctx->setForceManager(fm);

  CHECK(fm->isInitialized() == true);

  ctx->setCoordinates(crd);

  CHECK(ctx->getNumAtoms() == psf->getNumAtoms());

  auto &coordinatesCharges = ctx->getCoordinatesChargesDP();
  auto &velocitiesInverseMasses = ctx->getVelocitiesInverseMasses();

  CHECK(coordinatesCharges.size() ==
        static_cast<std::size_t>(psf->getNumAtoms()));
  CHECK(velocitiesInverseMasses.size() ==
        static_cast<std::size_t>(psf->getNumAtoms()));

  coordinatesCharges.transferFromDevice();
  velocitiesInverseMasses.transferFromDevice();

  const std::vector<double> charges = psf->getCharges();
  const std::vector<double> masses = psf->getMasses();

  REQUIRE(charges.empty() == false);
  REQUIRE(masses.empty() == false);

  CHECK(coordinatesCharges[0].w == Approx(charges[0]).margin(TOLERANCE));
  CHECK(velocitiesInverseMasses[0].w ==
        Approx(1.0 / masses[0]).margin(TOLERANCE));
}

TEST_CASE("CharmmContextRejectsInvalidStagedBoxDimensions") {
  CharmmContext ctx;

  CHECK_THROWS_AS(ctx.setBoxDimensions({40.0, 40.0}), std::invalid_argument);
  CHECK_THROWS_AS(ctx.setBoxDimensions({40.0, 40.0, 40.0, 40.0}),
                  std::invalid_argument);
  CHECK_THROWS_AS(ctx.setBoxDimensions({40.0, 0.0, 40.0}),
                  std::invalid_argument);
  CHECK_THROWS_AS(ctx.setBoxDimensions({40.0, -1.0, 40.0}),
                  std::invalid_argument);
}

TEST_CASE("ForceManagerContextBackPointerDoesNotOwnCharmmContext") {
  auto prm = std::make_shared<CharmmParameters>(getDataPath() +
                                                "toppar_water_ions.str");
  auto psf = std::make_shared<CharmmPSF>(getDataPath() + "nacl_pair.psf");
  auto crd = std::make_shared<CharmmCrd>(getDataPath() + "nacl_pair.cor");

  auto fm = std::make_shared<ForceManager>(psf, prm);
  fm->setBoxDimensions(BOX_DIMENSIONS);

  {
    auto ctx = std::make_shared<CharmmContext>(fm);
    ctx->setCoordinates(crd);

    CHECK(fm->hasCharmmContext() == true);
    CHECK(fm->getContext() == ctx);
  }

  CHECK(fm->hasCharmmContext() == false);
  CHECK(fm->getContext() == nullptr);
}

TEST_CASE("CharmmContextCanAttachForceManagerBeforeStagingState") {
  auto prm = std::make_shared<CharmmParameters>(getDataPath() +
                                                "toppar_water_ions.str");
  auto psf = std::make_shared<CharmmPSF>(getDataPath() + "nacl_pair.psf");

  auto ctx = std::make_shared<CharmmContext>();
  auto fm = std::make_shared<ForceManager>();

  ctx->setForceManager(fm);

  CHECK(fm->isInitialized() == false);

  ctx->setPsf(psf);
  CHECK(fm->isInitialized() == false);

  ctx->setPrm(prm);
  CHECK(fm->isInitialized() == false);

  ctx->setBoxDimensions(BOX_DIMENSIONS);

  CHECK(fm->isInitialized() == true);
  CHECK(ctx->getPsf() == psf);
  CHECK(ctx->getPrm() == prm);
  CHECK(fm->getPsf() == psf);
  CHECK(fm->getPrm() == prm);

  apo_test::CheckVectorsClose<double>("CharmmContext box dimensions",
                                      ctx->getBoxDimensions(), BOX_DIMENSIONS,
                                      TOLERANCE);
  apo_test::CheckVectorsClose<double>("ForceManager box dimensions",
                                      fm->getBoxDimensions(), BOX_DIMENSIONS,
                                      TOLERANCE);
}

TEST_CASE("CharmmContextConstructsFromPsfAndParameters") {
  auto prm = std::make_shared<CharmmParameters>(getDataPath() +
                                                "toppar_water_ions.str");
  auto psf = std::make_shared<CharmmPSF>(getDataPath() + "nacl_pair.psf");
  auto crd = std::make_shared<CharmmCrd>(getDataPath() + "nacl_pair.cor");

  auto ctx = std::make_shared<CharmmContext>(psf, prm);

  REQUIRE(ctx->getForceManager() != nullptr);
  CHECK(ctx->getPsf() == psf);
  CHECK(ctx->getPrm() == prm);
  CHECK(ctx->getForceManager()->getPsf() == psf);
  CHECK(ctx->getForceManager()->getPrm() == prm);

  CHECK(ctx->getForceManager()->isInitialized() == false);

  ctx->setBoxDimensions(BOX_DIMENSIONS);

  CHECK(ctx->getForceManager()->isInitialized() == true);

  apo_test::CheckVectorsClose<double>("CharmmContext box dimensions",
                                      ctx->getBoxDimensions(), BOX_DIMENSIONS,
                                      TOLERANCE);

  apo_test::CheckVectorsClose<double>(
      "ForceManager box dimensions", ctx->getForceManager()->getBoxDimensions(),
      BOX_DIMENSIONS, TOLERANCE);

  ctx->setCoordinates(crd);

  CHECK(ctx->getNumAtoms() == psf->getNumAtoms());
  CHECK(ctx->getCoordinatesChargesDP().size() ==
        static_cast<std::size_t>(psf->getNumAtoms()));
  CHECK(ctx->getVelocitiesInverseMasses().size() ==
        static_cast<std::size_t>(psf->getNumAtoms()));
}

TEST_CASE("CharmmContextForwardsForceManagerConfiguration") {
  auto prm = std::make_shared<CharmmParameters>(getDataPath() +
                                                "toppar_water_ions.str");
  auto psf = std::make_shared<CharmmPSF>(getDataPath() + "nacl_pair.psf");

  auto ctx = std::make_shared<CharmmContext>(psf, prm);
  auto fm = ctx->getForceManager();

  ctx->setKappa(0.45f);
  ctx->setCutoff(9.0f);
  ctx->setCtonnb(7.5f);
  ctx->setCtofnb(8.5f);
  ctx->setFFTGrid(32, 34, 36);
  ctx->setPmeSplineOrder(6);
  ctx->setVdwType(VDW_DBEXP);

  CHECK(ctx->getKappa() == Approx(0.45f));
  CHECK(ctx->getCutoff() == Approx(9.0f));
  CHECK(ctx->getCtonnb() == Approx(7.5f));
  CHECK(ctx->getCtofnb() == Approx(8.5f));
  apo_test::CheckVectorsEqual<int>("context FFT grid", ctx->getFFTGrid(),
                                   std::vector<int>{32, 34, 36});
  CHECK(ctx->getPmeSplineOrder() == 6);
  CHECK(ctx->getVdwType() == VDW_DBEXP);

  CHECK(fm->getKappa() == Approx(0.45f));
  CHECK(fm->getCutoff() == Approx(9.0f));
  CHECK(fm->getCtonnb() == Approx(7.5f));
  CHECK(fm->getCtofnb() == Approx(8.5f));
  apo_test::CheckVectorsEqual<int>("context FFT grid", fm->getFFTGrid(),
                                   std::vector<int>{32, 34, 36});
  CHECK(fm->getPmeSplineOrder() == 6);
  CHECK(fm->getVdwType() == VDW_DBEXP);
}

/* *
TEST_CASE("CharmmContextRejectsCoordinatesBeforeBackendInitialization") {
  auto prm = std::make_shared<CharmmParameters>(getDataPath() +
                                                "toppar_water_ions.str");
  auto psf = std::make_shared<CharmmPSF>(getDataPath() + "nacl_pair.psf");
  auto crd = std::make_shared<CharmmCrd>(getDataPath() + "nacl_pair.cor");

  auto ctx = std::make_shared<CharmmContext>(psf, prm);

  CHECK_THROWS_AS(ctx->setCoordinates(crd), std::runtime_error);
}
* */

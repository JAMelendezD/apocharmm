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

  auto ctx = std::make_shared<CharmmContext>(fm);

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

  // Keep this explicit for now. If a later patch auto-initializes in
  // setForceManager(), this remains harmless.
  if (!fm->isInitialized())
    fm->initialize();

  ctx->setCoordinates(crd);

  CHECK(ctx->getNumAtoms() == psf->getNumAtoms());

  auto &coordinatesCharges = ctx->getCoordinatesCharges();
  auto &velocitiesInverseMasses = ctx->getVelocityMass();

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

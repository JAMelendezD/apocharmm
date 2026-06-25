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
#include "apo_test_helpers.h"
#include "catch.hpp"
#include "test_paths.h"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int RANDOM_SEED = 314159;
constexpr double TEMPERATURE = 300.0;

const std::vector<double> BOX_DIMENSIONS = {50.0, 50.0, 50.0};

// std::shared_ptr<CharmmContext>
// CreateContext(const std::shared_ptr<ForceManager> &fm) {
//   auto crd = std::make_shared<CharmmCrd>(getDataPath() + "nacl_pair.cor");
//   auto ctx = std::make_shared<CharmmContext>(fm);
//   ctx->setCoordinates(crd);
//   ctx->useHolonomicConstraints(false);
//   ctx->setRandomSeedForVelocities(RANDOM_SEED);
//   ctx->assignVelocitiesAtTemperature(TEMPERATURE);
//   return ctx;
// }

} // namespace

TEST_CASE("ForceManagerDefaultConstructor") {
  ForceManager fm;

  CHECK(fm.getContext() == nullptr);
  CHECK(fm.hasCharmmContext() == false);
  CHECK(fm.isInitialized() == false);
  CHECK(fm.getBoxDimensions() ==
        std::vector<double>{-9999.9999, -9999.9999, -9999.9999});
  CHECK(fm.getKappa() == Approx(0.34f));
  CHECK(fm.getCutoff() == Approx(14.0f));
  CHECK(fm.getCtonnb() == Approx(12.0f));
  CHECK(fm.getCtofnb() == Approx(10.0f));
  CHECK(fm.getFFTGrid() == std::vector<int>{-1, -1, -1});
  CHECK(fm.getPeriodicBoundaryCondition() == PBC::P1);
  CHECK(fm.isComposite() == false);
  CHECK(fm.getChildren().empty() == true);
}

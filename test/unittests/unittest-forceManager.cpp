// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "CharmmPSF.h"
#include "CharmmParameters.h"
#include "CudaDirectForceTypes.h"
#include "CudaEnergyVirial.h"
#include "Force.h"
#include "ForceManager.h"
#include "PBC.h"
#include "apo_test_helpers.h"
#include "catch.hpp"
#include "test_paths.h"

#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double TOLERANCE = 0.0;

class TestForce {
public:
  static constexpr bool contributesVirial = false;

  TestForce(void)
      : initializeCalls(0), clearCalls(0), calcForceCalls(0),
        setBoxDimensionsCalls(0), initializeNumAtoms(-1), lastXYZQ(nullptr),
        lastCalcEnergy(false), lastCalcVirial(false), lastBoxDimensions(),
        forceValues(std::make_shared<Force<long long int>>()),
        energyVirial(std::make_shared<CudaEnergyVirial>()) {}

  void initialize(const int numAtoms,
                  const std::vector<double> &boxDimensions) {
    initializeCalls++;
    initializeNumAtoms = numAtoms;
    lastBoxDimensions = boxDimensions;
    return;
  }

  void clear(void) {
    clearCalls++;
    return;
  }

  void calcForce(const float4 *xyzq, const bool calcEnergy,
                 const bool calcVirial) {
    calcForceCalls++;
    lastXYZQ = xyzq;
    lastCalcEnergy = calcEnergy;
    lastCalcVirial = calcVirial;
    return;
  }

  void setBoxDimensions(const std::vector<double> &boxDimensions) {
    setBoxDimensionsCalls++;
    lastBoxDimensions = boxDimensions;
    return;
  }

  std::shared_ptr<Force<long long int>> getForce(void) { return forceValues; }

  std::shared_ptr<CudaEnergyVirial> getEnergyVirial(void) {
    return energyVirial;
  }

public:
  int initializeCalls;
  int clearCalls;
  int calcForceCalls;
  int setBoxDimensionsCalls;
  int initializeNumAtoms;
  const float4 *lastXYZQ;
  bool lastCalcEnergy;
  bool lastCalcVirial;
  std::vector<double> lastBoxDimensions;

private:
  std::shared_ptr<Force<long long int>> forceValues;
  std::shared_ptr<CudaEnergyVirial> energyVirial;
};

ForceManager CreateForceManager(void) {
  const std::string dataPath = getDataPath();
  auto prm =
      std::make_shared<CharmmParameters>(dataPath + "toppar_water_ions.str");
  auto psf = std::make_shared<CharmmPSF>(dataPath + "nacl_pair.psf");

  return ForceManager(psf, prm);
}

} // namespace

TEST_CASE("ForceManagerDefaultConstructor") {
  ForceManager fm;

  CHECK(fm.hasCharmmContext() == false);
  CHECK(fm.getContext() == nullptr);
  CHECK(fm.getPsf() == nullptr);
  CHECK(fm.getPrm() == nullptr);
  CHECK(fm.isInitialized() == false);
  CHECK(fm.isComposite() == false);
  CHECK(fm.getChildren().empty() == true);

  apo_test::CheckVectorsClose<double>(
      "default box dimensions", fm.getBoxDimensions(),
      std::vector<double>{-9999.9999, -9999.9999, -9999.9999}, TOLERANCE);

  CHECK(fm.getKappa() == Approx(0.34).margin(TOLERANCE));
  CHECK(fm.getCutoff() == Approx(14.0).margin(TOLERANCE));
  CHECK(fm.getCtonnb() == Approx(12.0).margin(TOLERANCE));
  CHECK(fm.getCtofnb() == Approx(10.0).margin(TOLERANCE));
  apo_test::CheckVectorsEqual<int>("default FFT grid", fm.getFFTGrid(),
                                   std::vector<int>{-1, -1, -1});
  CHECK(fm.getPmeSplineOrder() == 4);
  CHECK(fm.getPeriodicBoundaryCondition() == PBC::P1);
  CHECK(fm.getVdwType() == VDW_VFSW);
}

TEST_CASE("ForceManagerConstructsFromTopologyAndParameters") {
  const std::string dataPath = getDataPath();
  auto prm =
      std::make_shared<CharmmParameters>(dataPath + "toppar_water_ions.str");
  auto psf = std::make_shared<CharmmPSF>(dataPath + "nacl_pair.psf");

  ForceManager fm(psf, prm);

  CHECK(fm.getPsf() == psf);
  CHECK(fm.getPrm() == prm);
  CHECK(fm.getNumAtoms() == 2);
  CHECK(fm.isInitialized() == false);
  CHECK(fm.isComposite() == false);

  CHECK_THROWS_AS(fm.initialize(), std::invalid_argument);
}

TEST_CASE("ForceManagerSettersAndGetters") {
  auto fm = CreateForceManager();

  const std::vector<double> boxDimensions = {20.0, 21.0, 22.0};
  fm.setBoxDimensions(boxDimensions);
  apo_test::CheckVectorsClose<double>("box dimensions", fm.getBoxDimensions(),
                                      boxDimensions, TOLERANCE);

  fm.setKappa(0.45f);
  fm.setCutoff(9.0f);
  fm.setCtonnb(7.5f);
  fm.setCtofnb(8.5f);
  fm.setFFTGrid(32, 34, 36);
  fm.setPmeSplineOrder(6);
  fm.setPeriodicBoundaryCondition(PBC::P21);
  fm.setPrintEnergyDecomposition(true);

  CHECK(fm.getKappa() == Approx(0.45f).margin(TOLERANCE));
  CHECK(fm.getCutoff() == Approx(9.0f).margin(TOLERANCE));
  CHECK(fm.getCtonnb() == Approx(7.5f).margin(TOLERANCE));
  CHECK(fm.getCtofnb() == Approx(8.5f).margin(TOLERANCE));
  apo_test::CheckVectorsEqual<int>("FFT grid", fm.getFFTGrid(),
                                   std::vector<int>{32, 34, 36});
  CHECK(fm.getPmeSplineOrder() == 6);
  CHECK(fm.getPeriodicBoundaryCondition() == PBC::P21);

  const std::vector<int> validVdwTypes = {VDW_VSH,  VDW_VSW, VDW_VFSW,
                                          VDW_VGSH, VDW_CUT, VDW_DBEXP};

  for (const int vdwType : validVdwTypes) {
    fm.setVdwType(vdwType);
    CHECK(fm.getVdwType() == vdwType);
  }
}

TEST_CASE("ForceManagerRejectsInvalidBoxDimensions") {
  auto fm = CreateForceManager();

  CHECK_THROWS_AS(fm.setBoxDimensions({20.0, 20.0}), std::invalid_argument);
  CHECK_THROWS_AS(fm.setBoxDimensions({20.0, 20.0, 20.0, 20.0}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fm.setBoxDimensions({20.0, 0.0, 20.0}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fm.setBoxDimensions({20.0, -1.0, 20.0}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fm.setBoxDimensions(
                      {20.0, std::numeric_limits<double>::quiet_NaN(), 20.0}),
                  std::invalid_argument);
}

TEST_CASE("ForceManagerRejectsInvalidInitializationInputs") {
  auto fm = CreateForceManager();

  SECTION("MissingBoxDimensionsBeforeInitialize") {
    CHECK_THROWS_AS(fm.initialize(), std::invalid_argument);
  }

  SECTION("CutoffLargerThanHalfOfBoxXBeforeInitialize") {
    fm.setBoxDimensions({20.0, 20.0, 20.0});
    fm.setCutoff(11.0f);
    CHECK_THROWS_AS(fm.initialize(), std::invalid_argument);
  }
}

TEST_CASE("ForceManagerCopyConstructorCopiesConfigurationOnly") {
  auto original = CreateForceManager();
  original.setBoxDimensions({30.0, 31.0, 32.0});
  original.setKappa(0.41f);
  original.setCutoff(10.0f);
  original.setCtonnb(8.0f);
  original.setCtofnb(9.0f);
  original.setFFTGrid(30, 32, 34);
  original.setPmeSplineOrder(5);
  original.setPeriodicBoundaryCondition(PBC::P21);
  original.setVdwType(VDW_DBEXP);

  ForceManager copy(original);

  CHECK(copy.getPsf() != original.getPsf());
  CHECK(copy.getPrm() != original.getPrm());
  CHECK(copy.getNumAtoms() == original.getNumAtoms());
  CHECK(copy.hasCharmmContext() == false);
  CHECK(copy.isInitialized() == false);

  apo_test::CheckVectorsClose<double>("copied box dimensions",
                                      copy.getBoxDimensions(),
                                      original.getBoxDimensions(), TOLERANCE);
  CHECK(copy.getKappa() == Approx(original.getKappa()).margin(TOLERANCE));
  CHECK(copy.getCutoff() == Approx(original.getCutoff()).margin(TOLERANCE));
  CHECK(copy.getCtonnb() == Approx(original.getCtonnb()).margin(TOLERANCE));
  CHECK(copy.getCtofnb() == Approx(original.getCtofnb()).margin(TOLERANCE));
  apo_test::CheckVectorsEqual<int>("copied FFT grid", copy.getFFTGrid(),
                                   original.getFFTGrid());
  CHECK(copy.getPmeSplineOrder() == original.getPmeSplineOrder());
  CHECK(copy.getPeriodicBoundaryCondition() ==
        original.getPeriodicBoundaryCondition());
  CHECK(copy.getVdwType() == original.getVdwType());

  original.setBoxDimensions({40.0, 41.0, 42.0});
  original.setKappa(0.11f);
  original.setVdwType(VDW_CUT);
  original.setPmeSplineOrder(7);

  apo_test::CheckVectorsClose<double>(
      "copied box dimensions remain independet", copy.getBoxDimensions(),
      std::vector<double>{30.0, 31.0, 32.0}, TOLERANCE);
  CHECK(copy.getKappa() == Approx(0.41f).margin(TOLERANCE));
  CHECK(copy.getVdwType() == VDW_DBEXP);
  CHECK(copy.getPmeSplineOrder() == 5);
}

TEST_CASE("ForceManagerAddForceManagerRejectsChildren") {
  ForceManager fm;

  CHECK_THROWS_AS(fm.addForceManager(nullptr), std::invalid_argument);

  auto child = std::make_shared<ForceManager>();
  CHECK_THROWS_AS(fm.addForceManager(child), std::invalid_argument);
}

TEST_CASE("ForceManagerSubscribeAndUnsubscribePropagateBoxDimensions") {
  auto fm = CreateForceManager();

  auto force = std::make_shared<TestForce>();
  auto stream = std::make_shared<cudaStream_t>();

  fm.subscribe(force, "test_force", stream, force->getForce(),
               force->getEnergyVirial());

  fm.setBoxDimensions({20.0, 21.0, 22.0});

  CHECK(force->initializeCalls == 0);
  CHECK(force->setBoxDimensionsCalls == 1);
  apo_test::CheckVectorsClose<double>(
      "subscribed force box dimensions", force->lastBoxDimensions,
      std::vector<double>{20.0, 21.0, 22.0}, TOLERANCE);

  fm.unsubscribe(force);
  fm.setBoxDimensions({30.0, 31.0, 32.0});

  CHECK(force->setBoxDimensionsCalls == 1);
  apo_test::CheckVectorsClose<double>(
      "unsubscribed force box dimensions unchanged", force->lastBoxDimensions,
      std::vector<double>{20.0, 21.0, 22.0}, TOLERANCE);
}

TEST_CASE("ForceManagerUnsubscribeByTag") {
  auto fm = CreateForceManager();

  auto force = std::make_shared<TestForce>();
  auto stream = std::make_shared<cudaStream_t>();

  fm.subscribe(force, "test_force", stream, force->getForce(),
               force->getEnergyVirial());

  fm.setBoxDimensions({20.0, 21.0, 22.0});
  CHECK(force->setBoxDimensionsCalls == 1);

  fm.unsubscribe("test_force");
  fm.setBoxDimensions({30.0, 31.0, 32.0});

  CHECK(force->setBoxDimensionsCalls == 1);
}

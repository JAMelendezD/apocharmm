// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author:  James E. Gonzales II
//
// ENDLICENSE

#include "AtomSelection.h"
#include "CudaContainer.h"
#include "CudaEnergyVirial.h"
#include "ForceManager.h"
#include "HarmonicCenterOfMassRestraintForce.h"
#include "catch.hpp"
#include "cuda_utils.h"

#include <array>
#include <cmath>
#include <cuda_runtime.h>
#include <stdexcept>
#include <string>
#include <vector>

using HmcmForce = HarmonicCenterOfMassRestraintForce<long long int, float>;

namespace {

constexpr double forceTolerance = 2.0e-5;
constexpr double energyTolerance = 2.0e-5;
constexpr double virialTolerance = 2.0e-4;

struct HmcmInput {
  int numAtoms = 0;
  std::vector<float4> xyzq;
  std::vector<int> selectedAtomIndices;
  std::vector<double> boxDimensions = {20.0, 20.0, 20.0};
  std::vector<double> referencePosition = {0.0, 0.0, 0.0};
  std::vector<int> referenceMask = {1, 1, 1};
  std::vector<double> masses;
  bool useMassWeighting = true;
  double forceConstant = 0.0;
  double referenceDistance = 0.0;
  bool calcVirial = true;
};

struct HmcmResult {
  double energy = 0.0;
  std::vector<double> fx;
  std::vector<double> fy;
  std::vector<double> fz;
  std::array<double, 9> virial = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
};

AtomSelection make_selection(const int numAtoms,
                             const std::vector<int> &atomIndices) {
  AtomSelection selection(numAtoms, AtomSelection::InitialValue::NONE);

  for (const int atomIndex : atomIndices)
    selection.set(atomIndex);

  return selection;
}

void check_close(const std::string &label, const double observed,
                 const double expected, const double tolerance) {
  INFO(label << ": observed = " << observed << ", expected = " << expected);
  CHECK(std::abs(observed - expected) <= tolerance);
  return;
}

void check_force(const HmcmResult &result, const int atomIndex,
                 const double expectedFx, const double expectedFy,
                 const double expectedFz) {
  check_close("fx[" + std::to_string(atomIndex) + "]", result.fx[atomIndex],
              expectedFx, forceTolerance);
  check_close("fy[" + std::to_string(atomIndex) + "]", result.fy[atomIndex],
              expectedFy, forceTolerance);
  check_close("fz[" + std::to_string(atomIndex) + "]", result.fz[atomIndex],
              expectedFz, forceTolerance);
  return;
}

std::array<double, 9>
expected_virial_from_center_gradient(const std::array<double, 3> &center,
                                     const std::array<double, 3> &gradient) {
  return {-center[0] * gradient[0], -center[0] * gradient[1],
          -center[0] * gradient[2], -center[1] * gradient[0],
          -center[1] * gradient[1], -center[1] * gradient[2],
          -center[2] * gradient[0], -center[2] * gradient[1],
          -center[2] * gradient[2]};
}

void check_virial(const HmcmResult &result,
                  const std::array<double, 9> &expectedVirial) {
  for (int i = 0; i < 9; i++) {
    check_close("virial[" + std::to_string(i) + "]", result.virial[i],
                expectedVirial[i], virialTolerance);
  }
  return;
}

void copy_fixed_point_forces_to_result(HmcmResult &result,
                                       Force<long long int> &force,
                                       const cudaStream_t stream,
                                       const int numAtoms) {
  force.convert<double>(stream);
  const int stride = force.stride();

  std::vector<double> rawForces(static_cast<std::size_t>(3 * stride), 0.0);

  cudaCheck(cudaMemcpy(
      rawForces.data(), reinterpret_cast<const double *>(force.xyz()),
      rawForces.size() * sizeof(double), cudaMemcpyDeviceToHost));

  result.fx.assign(static_cast<std::size_t>(numAtoms), 0.0);
  result.fy.assign(static_cast<std::size_t>(numAtoms), 0.0);
  result.fz.assign(static_cast<std::size_t>(numAtoms), 0.0);

  for (int i = 0; i < numAtoms; i++) {
    result.fx[i] = rawForces[0 * stride + i];
    result.fy[i] = rawForces[1 * stride + i];
    result.fz[i] = rawForces[2 * stride + i];
  }

  return;
}

HmcmResult evaluate_hmcm(const HmcmInput &input) {
  REQUIRE(input.numAtoms > 0);
  REQUIRE(input.xyzq.size() == static_cast<std::size_t>(input.numAtoms));

  HmcmForce hmcm(input.numAtoms);

  hmcm.setBoxDimensions(input.boxDimensions);
  hmcm.setSelection(make_selection(input.numAtoms, input.selectedAtomIndices));
  hmcm.setForceConstant(input.forceConstant);
  hmcm.setReferencePosition(input.referencePosition, input.referenceMask);
  hmcm.setReferenceDistance(input.referenceDistance);

  if (!input.masses.empty()) {
    hmcm.setMasses(input.masses);
    if (!input.useMassWeighting)
      hmcm.setMassWeighting(false);
  }

  CudaContainer<float4> xyzq(input.xyzq);

  hmcm.clear();
  hmcm.calcForce(xyzq.getDeviceArray().data(), true, input.calcVirial);

  std::shared_ptr<cudaStream_t> stream = hmcm.getStream();
  cudaCheck(cudaStreamSynchronize(*stream));

  HmcmResult result;
  copy_fixed_point_forces_to_result(result, *hmcm.getForce(), *stream,
                                    input.numAtoms);

  std::shared_ptr<CudaEnergyVirial> energyVirial = hmcm.getEnergyVirial();
  energyVirial->copyToHost(*stream);
  cudaCheck(cudaStreamSynchronize(*stream));

  result.energy = energyVirial->getEnergy("hmcm");

  if (input.calcVirial) {
    std::shared_ptr<Force<long long int>> force = hmcm.getForce();

    energyVirial->calcVirial(input.numAtoms, xyzq.getDeviceArray().data(),
                             input.boxDimensions[0], input.boxDimensions[1],
                             input.boxDimensions[2], force->stride(),
                             reinterpret_cast<double *>(force->xyz()), *stream);
    cudaCheck(cudaStreamSynchronize(*stream));

    energyVirial->copyToHost(*stream);
    cudaCheck(cudaStreamSynchronize(*stream));

    double virial[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    energyVirial->getVirial(virial);

    for (int i = 0; i < 9; i++)
      result.virial[i] = virial[i];
  }

  return result;
}

} // namespace

TEST_CASE("virial_metadata") {
  static_assert(HarmonicCenterOfMassRestraintForce<long long int,
                                                   float>::contributesVirial);

  HarmonicCenterOfMassRestraintForce<long long int, float> hmcm(2);
  ForceView hmcmView(&hmcm);
  CHECK(hmcmView.contributesVirial());
}

TEST_CASE("validates_input") {
  HarmonicCenterOfMassRestraintForce<long long int, float> hmcm(2);

  CHECK_THROWS_AS(hmcm.setForceConstant(-1.0), std::invalid_argument);
  CHECK_THROWS_AS(hmcm.setReferenceDistance(-1.0), std::invalid_argument);

  CHECK_THROWS_AS(hmcm.setReferencePosition(std::vector<double>{0.0, 0.0},
                                            std::vector<int>{1, 1, 1}),
                  std::invalid_argument);

  CHECK_THROWS_AS(hmcm.setReferencePosition(std::vector<double>{0.0, 0.0, 0.0},
                                            std::vector<int>{1, 1}),
                  std::invalid_argument);

  CHECK_THROWS_AS(hmcm.setReferencePosition(std::vector<double>{0.0, 0.0, 0.0},
                                            std::vector<int>{0, 0, 0}),
                  std::invalid_argument);

  CHECK_THROWS_AS(hmcm.setBoxDimensions(std::vector<double>{10.0, 10.0}),
                  std::invalid_argument);

  CHECK_THROWS_AS(hmcm.setBoxDimensions(std::vector<double>{10.0, -10.0, 10.0}),
                  std::invalid_argument);

  AtomSelection emptySelection(2, AtomSelection::InitialValue::NONE);
  CHECK_THROWS_AS(hmcm.setSelection(emptySelection), std::invalid_argument);

  CHECK_THROWS_AS(hmcm.setMasses(std::vector<double>{0.0, 0.0}),
                  std::invalid_argument);
}

TEST_CASE("unweighted_point_restraint") {
  HmcmInput input;
  input.numAtoms = 3;
  input.xyzq = {make_float4(1.0f, 2.0f, 3.0f, 0.0f),
                make_float4(3.0f, 4.0f, 5.0f, 0.0f),
                make_float4(9.0f, 9.0f, 9.0f, 0.0f)};
  input.selectedAtomIndices = {0, 1};
  input.boxDimensions = {20.0, 20.0, 20.0};
  input.referencePosition = {1.0, 1.0, 1.0};
  input.forceConstant = 2.5;

  const HmcmResult result = evaluate_hmcm(input);

  check_close("energy", result.energy, 35.0, energyTolerance);

  check_force(result, 0, 2.5, 5.0, 7.5);
  check_force(result, 1, 2.5, 5.0, 7.5);
  check_force(result, 2, 0.0, 0.0, 0.0);

  const std::array<double, 9> expectedVirial =
      expected_virial_from_center_gradient({2.0, 3.0, 4.0}, {5.0, 10.0, 15.0});
  check_virial(result, expectedVirial);
}

TEST_CASE("coordinate_mask") {
  HmcmInput input;
  input.numAtoms = 3;
  input.xyzq = {make_float4(1.0f, 2.0f, 3.0f, 0.0f),
                make_float4(3.0f, 4.0f, 5.0f, 0.0f),
                make_float4(9.0f, 9.0f, 9.0f, 0.0f)};
  input.selectedAtomIndices = {0, 1};
  input.boxDimensions = {20.0, 20.0, 20.0};
  input.referencePosition = {1.0, 0.0, 0.0};
  input.referenceMask = {1, 0, 0};
  input.forceConstant = 2.5;

  const HmcmResult result = evaluate_hmcm(input);

  check_close("energy", result.energy, 2.5, energyTolerance);

  check_force(result, 0, 2.5, 0.0, 0.0);
  check_force(result, 1, 2.5, 0.0, 0.0);
  check_force(result, 2, 0.0, 0.0, 0.0);

  const std::array<double, 9> expectedVirial =
      expected_virial_from_center_gradient({2.0, 3.0, 4.0}, {5.0, 0.0, 0.0});
  check_virial(result, expectedVirial);
}

TEST_CASE("mass_weighting") {
  HmcmInput input;
  input.numAtoms = 3;
  input.xyzq = {make_float4(0.0f, 0.0f, 0.0f, 0.0f),
                make_float4(4.0f, 0.0f, 0.0f, 0.0f),
                make_float4(9.0f, 9.0f, 9.0f, 0.0f)};
  input.selectedAtomIndices = {0, 1};
  input.boxDimensions = {20.0, 20.0, 20.0};
  input.referencePosition = {1.0, 0.0, 0.0};
  input.masses = {1.0, 3.0, 100.0};
  input.forceConstant = 2.0;

  const HmcmResult result = evaluate_hmcm(input);

  check_close("energy", result.energy, 8.0, energyTolerance);

  check_force(result, 0, 2.0, 0.0, 0.0);
  check_force(result, 1, 6.0, 0.0, 0.0);
  check_force(result, 2, 0.0, 0.0, 0.0);

  const std::array<double, 9> expectedVirial =
      expected_virial_from_center_gradient({3.0, 0.0, 0.0}, {8.0, 0.0, 0.0});
  check_virial(result, expectedVirial);
}

TEST_CASE("can_disable_mass_weighting") {
  HmcmInput input;
  input.numAtoms = 3;
  input.xyzq = {make_float4(0.0f, 0.0f, 0.0f, 0.0f),
                make_float4(4.0f, 0.0f, 0.0f, 0.0f),
                make_float4(9.0f, 9.0f, 9.0f, 0.0f)};
  input.selectedAtomIndices = {0, 1};
  input.boxDimensions = {20.0, 20.0, 20.0};
  input.referencePosition = {1.0, 0.0, 0.0};
  input.masses = {1.0, 3.0, 100.0};
  input.useMassWeighting = false;
  input.forceConstant = 2.0;

  const HmcmResult result = evaluate_hmcm(input);

  check_close("energy", result.energy, 2.0, energyTolerance);

  check_force(result, 0, 2.0, 0.0, 0.0);
  check_force(result, 1, 2.0, 0.0, 0.0);
  check_force(result, 2, 0.0, 0.0, 0.0);

  const std::array<double, 9> expectedVirial =
      expected_virial_from_center_gradient({2.0, 0.0, 0.0}, {4.0, 0.0, 0.0});
  check_virial(result, expectedVirial);
}

TEST_CASE("reference_distance") {
  HmcmInput input;
  input.numAtoms = 2;
  input.xyzq = {make_float4(3.0f, 0.0f, 0.0f, 0.0f),
                make_float4(5.0f, 0.0f, 0.0f, 0.0f)};
  input.selectedAtomIndices = {0, 1};
  input.boxDimensions = {20.0, 20.0, 20.0};
  input.referencePosition = {0.0, 0.0, 0.0};
  input.referenceDistance = 1.0;
  input.forceConstant = 3.0;

  const HmcmResult result = evaluate_hmcm(input);

  check_close("energy", result.energy, 27.0, energyTolerance);

  check_force(result, 0, 9.0, 0.0, 0.0);
  check_force(result, 1, 9.0, 0.0, 0.0);

  const std::array<double, 9> expectedVirial =
      expected_virial_from_center_gradient({4.0, 0.0, 0.0}, {18.0, 0.0, 0.0});
  check_virial(result, expectedVirial);
}

TEST_CASE("periodic_wrapping_and_virial") {
  HmcmInput input;
  input.numAtoms = 3;
  input.xyzq = {make_float4(9.8f, 0.0f, 0.0f, 0.0f),
                make_float4(0.4f, 0.0f, 0.0f, 0.0f),
                make_float4(5.0f, 5.0f, 5.0f, 0.0f)};
  input.selectedAtomIndices = {0, 1};
  input.boxDimensions = {10.0, 10.0, 10.0};
  input.referencePosition = {1.1, 0.0, 0.0};
  input.forceConstant = 2.0;

  const HmcmResult result = evaluate_hmcm(input);

  // The first selected atom anchors the group near x = 10:
  //   atom 0 effective x before center shift = 9.8
  //   atom 1 effective x before center shift = 10.4
  //   center before center shift = 10.1
  //
  // The center is then minimum-imaged against reference x = 1.1:
  //   center shift = -10.0
  //   displacement = 10.1 - 10.0 - 1.1 = -1.0
  //   center gradient = 2 * k * displacement = -4.0
  //   atom gradients = -2.0 and -2.0
  check_close("energy", result.energy, 2.0, energyTolerance);

  check_force(result, 0, -2.0, 0.0, 0.0);
  check_force(result, 1, -2.0, 0.0, 0.0);
  check_force(result, 2, 0.0, 0.0, 0.0);

  // Effective center after all periodic shifts is x = 0.1.
  // Virial = -center outer centerGradient = -0.1 * -4.0 = 0.4.
  const std::array<double, 9> expectedVirial =
      expected_virial_from_center_gradient({0.1, 0.0, 0.0}, {-4.0, 0.0, 0.0});
  check_virial(result, expectedVirial);
}

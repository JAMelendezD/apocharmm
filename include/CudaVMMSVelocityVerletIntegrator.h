// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: Samarjeet Prasad, James E. Gonzales II
//
// ENDLICENSE

#pragma once

#include "CudaIntegrator.h"

#include "CharmmContext.h"
#include "CudaContainer.h"

#include <memory>
#include <vector>

class CudaVMMSVelocityVerletIntegrator : public CudaIntegrator {
public:
  CudaVMMSVelocityVerletIntegrator(const double timeStep);

public:
  void setCharmmContexts(const std::vector<CharmmContext> &ctxs);
  void setSoluteAtoms(const std::vector<int> &atoms);

protected:
  void initializeImpl(void) override;
  void propagateOneStepImpl(void) override;

private:
  void combineForces(void);

private:
  std::vector<CharmmContext> m_Contexts;
  CudaContainer<int> m_SoluteAtoms;
  std::shared_ptr<Force<float>> m_CombinedForce;

  std::vector<float> m_Weights;
};

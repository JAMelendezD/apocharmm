// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author:  Samarjeet Prasad, James E. Gonzales II
//
// ENDLICENSE

#pragma once

#include "CudaIntegrator.h"

class CudaLeapFrogIntegrator : public CudaIntegrator {
public:
  CudaLeapFrogIntegrator(const double timeStep);

protected:
  void initializeImpl(void) override;
  void propagateOneStepImpl(void) override;
};

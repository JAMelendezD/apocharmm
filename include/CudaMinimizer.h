// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: Samarjeet Prasad
//
// ENDLICENSE

#pragma once

#include "CudaIntegrator.h"

#include <string>

class CudaMinimizer : public CudaIntegrator {
public:
  CudaMinimizer();

  // This should not be a raw pointer
  // void setCharmmContext(std::shared_ptr<CharmmContext> csc);

  void minimize(int numSteps);
  void minimize();

  void setMethod(std::string _method);

  void setVerboseFlag(bool _flag = true);

protected:
  void initializeImpl(void) override;

private:
  int nsteps;
  std::string method;
  bool verboseFlag;
};

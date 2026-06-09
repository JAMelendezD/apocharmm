// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_DETAIL_CUDA_LANGEVIN_THERMOSTAT_INTEGRATOR_HANDLE_H__
#define __APOCHARMM_C_DETAIL_CUDA_LANGEVIN_THERMOSTAT_INTEGRATOR_HANDLE_H__

#include "apocharmm_c/CudaLangevinThermostatIntegrator.h"
#include "apocharmm_c/detail/CudaIntegratorHandle.h"

#include "CudaLangevinThermostatIntegrator.h"

#include <memory>

struct apo_cuda_langevin_thermostat_integrator {
  std::shared_ptr<CudaLangevinThermostatIntegrator> object = nullptr;
  apo_cuda_integrator base;
};

#endif

// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_DETAIL_HARMONIC_CENTER_OF_MASS_RESTRAINT_FORCE_HANDLE_H__
#define __APOCHARMM_C_DETAIL_HARMONIC_CENTER_OF_MASS_RESTRAINT_FORCE_HANDLE_H__

#include "apocharmm_c/HarmonicCenterOfMassRestraintForce.h"

#include "HarmonicCenterOfMassRestraintForce.h"

#include <memory>

struct apo_harmonic_center_of_mass_restraint_force {
  std::shared_ptr<HarmonicCenterOfMassRestraintForce<long long int, float>>
      object = nullptr;
};

#endif

// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_DETAIL_HARMONIC_RESTRAINT_FORCE_HANDLE_H__
#define __APOCHARMM_C_DETAIL_HARMONIC_RESTRAINT_FORCE_HANDLE_H__

#include "apocharmm_c/HarmonicRestraintForce.h"

#include "HarmonicRestraintForce.h"

#include <memory>

/**
 * @brief Owns the native restraint referenced by one opaque C ABI handle.
 *
 * A successfully created handle retains shared ownership of the supported
 * `HarmonicRestraintForce<long long int, float>` specialization. Public C ABI
 * entry points require both the handle pointer and `object` to be non-null.
 *
 * `ForceManager` subscription copies the native shared owner. Deleting this
 * handle therefore does not remove an existing subscription or necessarily
 * destroy the native restraint.
 *
 * This private representation is not part of the stable C ABI.
 */
struct apo_harmonic_restraint_force {
  std::shared_ptr<HarmonicRestraintForce<long long int, float>> object =
      nullptr;
};

#endif

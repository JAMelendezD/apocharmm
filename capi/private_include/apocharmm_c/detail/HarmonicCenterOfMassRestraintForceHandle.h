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

/**
 * @brief Owns one shared native harmonic center-of-mass restraint reference.
 *
 * The public C handle exclusively owns this wrapper object. Its `object` member
 * shares the native `HarmonicCenterOfMassRestraintForce<long long int, float>`
 * instance with any subscribed `ForceManager`. Deleting the C handle releases
 * only this shared reference; manager subscription can keep the native force,
 * stream holder, force storage, and energy-virial storage alive.
 *
 * A non-NULL public handle is valid only while this wrapper remains allocated
 * and `object` is non-NULL. The wrapper provides no host-thread
 * synchronization.
 */
struct apo_harmonic_center_of_mass_restraint_force {
  std::shared_ptr<HarmonicCenterOfMassRestraintForce<long long int, float>>
      object = nullptr;
};

#endif

// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_DETAIL_ATOM_SELECTOR_HANDLE_H__
#define __APOCHARMM_C_DETAIL_ATOM_SELECTOR_HANDLE_H__

#include "apocharmm_c/AtomSelector.h"

#include "AtomSelector.h"

#include <memory>

/**
 * @brief Defines the private owning state behind an atom-selector handle.
 *
 * Successful C construction allocates this wrapper and stores shared ownership
 * of one native @ref AtomSelector. The selector in turn retains shared
 * ownership of its const PSF, so destroying the public PSF handle does not
 * invalidate this wrapper. Public C callers see only the opaque declaration and
 * release the wrapper through @ref apo_atom_selector_destroy.
 */
struct apo_atom_selector {
  /**
   * @brief Retains shared ownership of the native selector and its PSF
   * relationship.
   */
  std::shared_ptr<AtomSelector> object = nullptr;
};

#endif

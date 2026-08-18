// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_DETAIL_ATOM_SELECTION_HANDLE_H__
#define __APOCHARMM_C_DETAIL_ATOM_SELECTION_HANDLE_H__

#include "apocharmm_c/AtomSelection.h"

#include "AtomSelection.h"

#include <memory>

/**
 * @brief Defines the private owning state behind an atom-selection handle.
 *
 * A successful C selection operation allocates this wrapper and stores shared
 * ownership of one immutable native @ref AtomSelection. Public C callers see
 * only the opaque declaration and release the wrapper through
 * @ref apo_atom_selection_destroy.
 */
struct apo_atom_selection {
  /** Retains shared ownership of the native atom-selection result. */
  std::shared_ptr<AtomSelection> object = nullptr;
};

#endif

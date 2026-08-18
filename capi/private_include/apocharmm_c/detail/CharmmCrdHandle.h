// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_DETAIL_CHARMM_CRD_HANDLE_H__
#define __APOCHARMM_C_DETAIL_CHARMM_CRD_HANDLE_H__

#include "apocharmm_c/CharmmCrd.h"

#include "CharmmCrd.h"

#include <memory>

/**
 * @brief Defines the private owning state behind the opaque coordinate handle.
 *
 * Successful C construction allocates this structure and stores one shared
 * owner of the parsed native CharmmCrd. Public getters borrow the structure and
 * never copy or transfer ownership. C callers see only the incomplete public
 * typedef and never access this definition or member directly.
 */
struct apo_charmm_crd {
  /** Retains shared ownership of the native CharmmCrd object. */
  std::shared_ptr<CharmmCrd> object = nullptr;
};

#endif

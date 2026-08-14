// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_DETAIL_FORCE_MANAGER_H__
#define __APOCHARMM_C_DETAIL_FORCE_MANAGER_H__

#include "apocharmm_c/ForceManager.h"

#include "CharmmPSF.h"
#include "CharmmParameters.h"
#include "ForceManager.h"

#include <memory>

/**
 * @brief Defines the private owning state behind the public opaque handle.
 *
 * A successfully created handle retains native shared ownership of the PSF and
 * parameter set independently of the source C handles. It also owns one shared
 * `ForceManager` constructed from those retained collaborators.
 *
 * The public destroy function deletes this structure, releasing all three
 * shared owners. Public C callers never access these members directly.
 */
struct apo_force_manager {
  /** Retained native PSF owner. */
  std::shared_ptr<CharmmPSF> psf = nullptr;

  /** Retained native parameter-set owner. */
  std::shared_ptr<CharmmParameters> parameters = nullptr;

  /** Retained native ForceManager owner exposed through the C ABI. */
  std::shared_ptr<ForceManager> object = nullptr;
};

#endif

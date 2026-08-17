// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_DETAIL_CHARMM_PSF_HANDLE_H__
#define __APOCHARMM_C_DETAIL_CHARMM_PSF_HANDLE_H__

#include "apocharmm_c/CharmmPsf.h"

#include "CharmmPSF.h"

#include <memory>

/**
 * @brief Defines the private owning state behind the opaque PSF handle.
 *
 * Successful C construction allocates this structure and stores one shared
 * owner of the parsed native PSF. ForceManager and CharmmContext handle
 * implementations may copy `object` so the native topology remains alive
 * independently of this public handle. C callers never access this definition
 * or member directly.
 */
struct apo_charmm_psf {
  /** Retains shared ownership of the native CharmmPSF object. */
  std::shared_ptr<CharmmPSF> object = nullptr;
};

#endif

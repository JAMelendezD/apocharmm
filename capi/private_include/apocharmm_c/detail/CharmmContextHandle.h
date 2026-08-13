// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_DETAIL_CHARMM_CONTEXT_HANDLE_H__
#define __APOCHARMM_C_DETAIL_CHARMM_CONTEXT_HANDLE_H__

#include "apocharmm_c/CharmmContext.h"

#include "CharmmContext.h"
#include "CharmmPSF.h"
#include "CharmmParameters.h"
#include "ForceManager.h"

#include <memory>

/**
 * @brief Defines the private ownership payload of an apo_charmm_context handle.
 *
 * `object` is the native context represented by the public opaque handle.
 * The collaborator members retain the same native PSF, parameter set, and
 * force manager independently of the source C handles used to construct or
 * update the context.
 *
 * Members are declared so normal reverse destruction releases `object` before
 * the separately retained force manager, parameters, and PSF. The public
 * destroy function deletes this struct and accepts a null pointer.
 *
 * @note This definition is private to the C ABI implementation and does not
 * expose C++ representation through the stable public ABI.
 */
struct apo_charmm_context {
  /** Shared native PSF retained for C-handle lifetime. */
  std::shared_ptr<CharmmPSF> psf = nullptr;
  /** Shared native parameter set retained for C-handle lifetime. */
  std::shared_ptr<CharmmParameters> parameters = nullptr;
  /** Shared native force manager retained for C-handle lifetime. */
  std::shared_ptr<ForceManager> force_manager = nullptr;
  /** Shared native context represented by this handle. */
  std::shared_ptr<CharmmContext> object = nullptr;
};

#endif

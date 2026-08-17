// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_DETAIL_CHARMM_PARAMETERS_H__
#define __APOCHARMM_C_DETAIL_CHARMM_PARAMETERS_H__

#include "apocharmm_c/CharmmParameters.h"

#include "CharmmParameters.h"

#include <memory>

/**
 * @brief Defines the private owning state behind the public opaque handle.
 *
 * A successful C construction allocates this structure and stores one shared
 * owner of the parsed native parameter set. Other native C handle structures
 * may copy `object` to preserve the parameter lifetime independently of this
 * public handle. C callers never access this member directly.
 */
struct apo_charmm_parameters {
  std::shared_ptr<CharmmParameters> object = nullptr;
};

#endif

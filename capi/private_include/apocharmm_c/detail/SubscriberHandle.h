// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_DETAIL_SUBSCRIBER_HANDLE_H__
#define __APOCHARMM_C_DETAIL_SUBSCRIBER_HANDLE_H__

#include "apocharmm_c/Subscriber.h"

#include "Subscriber.h"

#include <memory>

/**
 * @brief Implements a C base-subscriber view with shared native ownership.
 *
 * Concrete C handles embed this structure and assign @ref object to the same
 * `std::shared_ptr` stored by their concrete member. The public pointer to the
 * embedded structure is borrowed and becomes invalid when the concrete C handle
 * is deleted, even when another native owner keeps the C++ object alive.
 */
struct apo_subscriber {
  /** @brief Shares ownership of the native polymorphic subscriber object. */
  std::shared_ptr<Subscriber> object = nullptr;
};

#endif

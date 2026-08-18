// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_DETAIL_RESTART_SUBSCRIBER_HANDLE_H__
#define __APOCHARMM_C_DETAIL_RESTART_SUBSCRIBER_HANDLE_H__

#include "apocharmm_c/RestartSubscriber.h"
#include "apocharmm_c/detail/SubscriberHandle.h"

#include "RestartSubscriber.h"

#include <memory>

/**
 * @brief Implements the owned C handle for a native restart subscriber.
 *
 * @ref object and `base.object` share ownership of the same native instance.
 * The address of @ref base is returned as a borrowed polymorphic view and is
 * valid only while this enclosing handle remains allocated.
 */
struct apo_restart_subscriber {
  /**
   * @brief Owns the concrete native restart subscriber through shared
   * ownership.
   */
  std::shared_ptr<RestartSubscriber> object = nullptr;

  /** @brief Stores the embedded borrowed base view of @ref object. */
  apo_subscriber base;
};

#endif

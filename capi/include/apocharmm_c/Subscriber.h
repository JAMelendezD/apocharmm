// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_SUBSCRIBER_H__
#define __APOCHARMM_C_SUBSCRIBER_H__

#include "apocharmm_c/Export.h"
#include "apocharmm_c/Status.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Represents a borrowed base-subscriber view in the C ABI.
 *
 * Concrete subscriber conversion functions, such as
 * @ref apo_dcd_subscriber_as_subscriber, return this view. The view does not
 * own the concrete handle and has no destroy function. It remains valid only
 * while the concrete C handle that contains it remains alive.
 *
 * All functions in this header borrow the view for one call and do not retain
 * the C pointer. Native integrator subscription can retain the underlying C++
 * subscriber object, but it does not extend the lifetime of this embedded C
 * view. Preserve the concrete C handle whenever the view will be used again.
 *
 * Passing `NULL`, a dangling view, a view whose concrete owner was destroyed,
 * or a foreign pointer is invalid. Do not overlap mutation, destruction of the
 * concrete owner, or integrator attachment operations on the same object from
 * multiple threads.
 *
 * Status-returning calls clear the calling thread's previous diagnostic at
 * entry. On failure, call @ref apo_last_error immediately on the same thread.
 * The returned diagnostic pointer is borrowed and is invalidated by the next
 * diagnostic-changing C ABI call on that thread.
 *
 * @see subscriber
 */
typedef struct apo_subscriber apo_subscriber;

/**
 * @brief Sets a subscriber's reporting interval.
 *
 * @param[in] subscriber Borrowed live base view. The pointer may not be `NULL`
 * and is not retained by this call.
 * @param[in] report_frequency Positive, dimensionless number of propagated
 * steps between scheduled updates.
 *
 * @retval APO_STATUS_OK The native reporting interval was updated.
 * @retval APO_STATUS_INVALID_ARGUMENT `subscriber` is `NULL`, its native object
 * is absent, or `report_frequency` is not positive.
 * @retval APO_STATUS_RUNTIME_ERROR An unexpected standard or nonstandard C++
 * exception crossed the C ABI boundary.
 *
 * @post On a rejected frequency, the previous native value is unchanged.
 * @note The function clears the previous thread-local diagnostic at entry.
 * Success leaves it empty; failure leaves text available through
 * @ref apo_last_error.
 * @warning An integrator caches the frequency during subscription. Changing the
 * native subscriber afterward does not update the scheduling cache.
 */
APOCHARMM_C_API apo_status apo_subscriber_set_report_frequency(
    apo_subscriber *subscriber, const int report_frequency);

/**
 * @brief Returns a subscriber's reporting interval.
 *
 * @param[out] report_frequency Non-`NULL` output pointer receiving a positive,
 * dimensionless step count on success. The function does not initialize this
 * location before validating the handle and pointer.
 * @param[in] subscriber Borrowed live base view. The pointer may not be `NULL`
 * and is not retained by this call.
 *
 * @retval APO_STATUS_OK The current interval was written to
 * `*report_frequency`.
 * @retval APO_STATUS_INVALID_ARGUMENT `subscriber` is `NULL`, its native object
 * is absent, or `report_frequency` is `NULL`.
 * @retval APO_STATUS_RUNTIME_ERROR An unexpected standard or nonstandard C++
 * exception crossed the C ABI boundary.
 *
 * @post On failure, a valid `report_frequency` location retains its previous
 * contents.
 * @note The function clears the previous thread-local diagnostic at entry.
 * Success leaves it empty; failure leaves text available through
 * @ref apo_last_error.
 */
APOCHARMM_C_API apo_status apo_subscriber_get_report_frequency(
    int *report_frequency, const apo_subscriber *subscriber);

#ifdef __cplusplus
}
#endif

#endif

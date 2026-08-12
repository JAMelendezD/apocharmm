// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_STATUS_H__
#define __APOCHARMM_C_STATUS_H__

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Identifies success or a categorized failure returned by the C ABI.
 *
 * The explicit numeric values are part of the public C ABI. After a nonzero
 * status, call @ref apo_last_error on the same thread before another
 * status-returning apoCHARMM C API call clears or replaces the diagnostic.
 *
 * @note Every declared nonzero value represents failure.
 * @see ApoCharmmErrorCode
 */
typedef enum apo_status {
  /** Indicates that the operation completed successfully. */
  APO_STATUS_OK = 0,
  /** Indicates that an argument, pointer, size, or value was rejected. */
  APO_STATUS_INVALID_ARGUMENT = 1,
  /** Indicates an uncategorized runtime or translated exception failure. */
  APO_STATUS_RUNTIME_ERROR = 2,
  /** Indicates that a CUDA operation failed. */
  APO_STATUS_CUDA_ERROR = 3,
  /** Indicates that required native state has not been initialized. */
  APO_STATUS_NOT_INITIALIZED = 4,
  /** Indicates that a recognized operation is not implemented. */
  APO_STATUS_NOT_IMPLEMENTED = 5,
} apo_status;

#ifdef __cplusplus
}
#endif

#endif

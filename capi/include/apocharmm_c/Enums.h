// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#ifndef __APOCHARMM_C_ENUMS_H__
#define __APOCHARMM_C_ENUMS_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef enum apo_pbc {
  APO_PBC_NONE = 0,
  APO_PBC_P1 = 1,
  APO_PBC_P21 = 2
} apo_pbc;

typedef enum apo_crystal_type {
  APO_CRYSTAL_TYPE_NONE = 0,
  APO_CRYSTAL_TYPE_CUBIC = 1,
  APO_CRYSTAL_TYPE_TETRAGONAL = 2,
  APO_CRYSTAL_TYPE_ORTHORHOMBIC = 3,
} apo_crystal_type;

#ifdef __cplusplus
}
#endif

#endif

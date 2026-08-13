// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: Antti-Pekka Hynninen, Samarjeet Prasad
//
// ENDLICENSE

#pragma once

/**
 * @brief Identifies the periodic boundary condition used by apoCHARMM.
 *
 * The value is stored by @ref CharmmContext and forwarded to its associated
 * @ref ForceManager. `P1` is the default imported from a newly constructed
 * `ForceManager`.
 */
enum class PBC {
  /** Selects no periodic boundary condition. */
  NONE,
  /** Selects conventional three-dimensional translational periodicity. */
  P1,
  /** Selects the P2_1 screw-symmetry periodic boundary condition. */
  P21
};

enum class CRYSTAL { NONE, CUBIC, TETRAGONAL, ORTHORHOMBIC };

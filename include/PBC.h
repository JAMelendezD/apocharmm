// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: Antti-Pekka Hynninen, Samarjeet Prasad
//
// ENDLICENSE

/**
 * @file
 * @brief Declares apoCHARMM periodic-boundary and crystal identifiers.
 */

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

/**
 * @brief Selects the crystal symmetry used by the Langevin-piston integrator.
 *
 * The selected symmetry determines the number and ordering of active piston
 * degrees of freedom. @ref CudaLangevinPistonIntegrator accepts only `CUBIC`,
 * `TETRAGONAL`, and `ORTHORHOMBIC`; `NONE` represents an unconfigured state.
 *
 * @see cuda_integrators
 */
enum class CRYSTAL {
  /** Represents an unconfigured crystal with no piston degree of freedom. */
  NONE,
  /** Couples X, Y, and Z to one isotropic piston degree of freedom. */
  CUBIC,
  /** Couples X and Y and assigns Z a second piston degree of freedom. */
  TETRAGONAL,
  /** Assigns independent piston degrees of freedom to X, Y, and Z. */
  ORTHORHOMBIC
};

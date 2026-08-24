// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author:  James E. Gonzales II
//
// ENDLICENSE

#pragma once

#include "AtomSelection.h"
#include "CudaContainer.h"
#include "CudaEnergyVirial.h"
#include "Force.h"

#include <cuda_runtime.h>
#include <memory>
#include <vector>

/**
 * @brief Computes a periodic harmonic restraint on a selected group center.
 *
 * The selected atoms define either a uniformly weighted geometric center or a
 * mass-weighted center. The first selected atom is used as an anchor while
 * selected coordinates are unwrapped through an orthorhombic periodic box.
 * Let `w_i` be one when mass weighting is disabled, or the configured mass
 * when it is enabled. The resulting center is:
 *
 * @code
 * W = sum over selected atoms i of w_i
 * C = (1 / W) * sum over selected atoms i of w_i * r_i_unwrapped
 * @endcode
 *
 * A three-element mask selects the Cartesian components used in the
 * center-to-reference displacement. Active components are minimum-imaged
 * relative to the reference position; inactive components are zero. For
 * displacement magnitude `r`, force constant `k`, and reference distance
 * `r0`, the implemented energy is:
 *
 * @code
 * U_hmcm = k * r^2             when r0 == 0
 * U_hmcm = k * (r - r0)^2      when r0 > 0
 * @endcode
 *
 * The per-atom force array receives each atom's weight fraction times the
 * potential gradient with respect to the center. apoCHARMM's force aggregation
 * subtracts this stored gradient to obtain the physical restoring force.
 * Energy is accumulated in the named term `hmcm`.
 *
 * Configuration is split between host scalar state and synchronized
 * host/device arrays. Selection or weight changes rebuild device-resident
 * selected-index and selected-weight arrays and can reallocate their storage.
 * `calcForce()` enqueues a three-kernel pipeline on a private CUDA stream and
 * does not synchronize before returning.
 *
 * @tparam AT Device force-accumulator representation. The library explicitly
 * instantiates `long long int`, using apoCHARMM's fixed-point force scale.
 * @tparam CT Arithmetic representation used while distributing the center
 * gradient. The library explicitly instantiates `float` and `double`.
 *
 * @warning The class provides no internal host-thread synchronization. Callers
 * must serialize configuration, clearing, force evaluation, and destruction.
 * @warning Do not copy or copy-assign this class. Its compiler-generated copy
 * operations share stream and output owners without defining an independent
 * CUDA-resource lifecycle.
 * @warning Call `clear()` before the first calculation and before every
 * calculation that requires fresh rather than accumulated output.
 *
 * @see harmonic_center_of_mass_restraint_force
 * @see AtomSelection
 * @see ForceManager
 */
template <typename AT, typename CT> class HarmonicCenterOfMassRestraintForce {
public:
  /**
   * @brief Reports that this restraint contributes periodic virial state.
   *
   * `ForceManager` passes `calcVirial == true` to this force only when global
   * virial calculation is requested. The restraint then accumulates periodic
   * shift-force records into its `CudaEnergyVirial` object.
   */
  static constexpr bool contributesVirial = true;

public:
  /**
   * @brief Prevents construction without an explicit atom count.
   */
  HarmonicCenterOfMassRestraintForce(void) = delete;

  /**
   * @brief Constructs a restraint for a fixed positive atom count.
   *
   * Every atom is initially selected. The force constant and reference
   * distance are zero, mass weighting is disabled, all stored masses are one,
   * the reference position is `(0, 0, 0)`, every Cartesian mask component is
   * active, and box dimensions are unset at zero.
   *
   * Construction creates coherent host/device mass, selected-index,
   * selected-weight, partial-sum, and restraint-state storage. It also creates
   * the `hmcm` energy term, a `Force<AT>` sized for `numAtoms`, and a private
   * CUDA stream. Output force and energy storage must be cleared before use.
   *
   * @param[in] numAtoms Dimensionless atom count. The value must be greater
   * than zero and remains fixed for the lifetime of the restraint.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `numAtoms` is not positive.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if CUDA
   * allocation, transfer, synchronization, force allocation, or stream
   * creation fails.
   * @throws std::bad_alloc If host storage, shared ownership state, or an error
   * diagnostic cannot be allocated.
   * @throws std::length_error If requested storage or an error diagnostic
   * exceeds an implementation-defined limit.
   *
   * @post `getForce()->size() == numAtoms`.
   * @post Every atom is selected with uniform weight one.
   */
  HarmonicCenterOfMassRestraintForce(const int numAtoms);

  /**
   * @brief Destroys the private CUDA stream and releases owned resources.
   *
   * CUDA stream destruction is attempted without propagating cleanup failure.
   * A shared stream holder previously returned by `getStream()` can outlive
   * this object, but its stored stream value is reset during destruction.
   * Shared force and energy objects can also outlive this restraint while
   * another owner exists, but no restraint remains to update them.
   *
   * @post The private CUDA stream is no longer valid.
   */
  ~HarmonicCenterOfMassRestraintForce(void) noexcept;

public:
  /**
   * @brief Sets the nonempty atom selection used to define the group center.
   *
   * The selection is copied. Its selected indices are materialized in
   * ascending order, and selected weights are rebuilt from either one or the
   * configured masses. Both arrays are copied to device memory. Partial-sum
   * storage is resized to one element per 256-thread calculation block and
   * internal state is cleared.
   *
   * @param[in] selection Selection with exactly the restraint's fixed atom
   * count and at least one selected atom. The object is borrowed for the call,
   * copied, and not retained by reference.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the atom count differs, no atom is
   * selected, or mass weighting is active and the selected masses sum to zero.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if selected
   * arrays cannot be allocated, transferred, or synchronized.
   * @throws std::bad_alloc If selection, index, weight, or diagnostic storage
   * cannot be allocated.
   * @throws std::length_error If copied storage or a diagnostic exceeds an
   * implementation-defined limit.
   *
   * @post On success, the restraint owns an independent selection copy and
   * selected device arrays describe exactly that selection.
   * @warning Failure after selection assignment begins can leave the copied
   * selection or selected count newer than the device index and weight arrays.
   * @warning Successful rebuilding can invalidate device pointers previously
   * borrowed from the selected arrays or partial-sum storage.
   */
  void setSelection(const AtomSelection &selection);

  /**
   * @brief Sets the scalar harmonic force constant.
   *
   * The value is retained in host scalar state and passed by value to later
   * kernel launches. This method performs no CUDA transfer or synchronization.
   *
   * @param[in] forceConstant Finite, non-negative force constant in
   * kilocalories per mole per square angstrom.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `forceConstant` is negative, NaN,
   * or infinite.
   * @throws std::bad_alloc If an error diagnostic cannot be allocated.
   * @throws std::length_error If an error diagnostic exceeds an
   * implementation-defined limit.
   *
   * @post On success, subsequent evaluations use `forceConstant`.
   * @note Zero disables force evaluation and permits `calcForce()` to return
   * without configured box dimensions.
   */
  void setForceConstant(const double forceConstant);

  /**
   * @brief Sets a Cartesian reference position with all components active.
   *
   * This overload delegates to the masked overload with `{1, 1, 1}`. The
   * values are copied into host scalar state; no caller storage is retained and
   * no CUDA transfer occurs.
   *
   * @param[in] referencePosition Exactly three finite Cartesian coordinates in
   * `[x, y, z]` order and angstrom units.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the vector length is not three or
   * any component is NaN or infinite.
   * @throws std::bad_alloc If an error diagnostic cannot be allocated.
   * @throws std::length_error If an error diagnostic exceeds an
   * implementation-defined limit.
   *
   * @post On success, every Cartesian component participates in the restraint.
   */
  void setReferencePosition(const std::vector<double> &referencePosition);

  /**
   * @brief Sets the Cartesian reference position and active-component mask.
   *
   * Active components use minimum-image center-to-reference displacement.
   * Inactive components contribute neither energy nor gradient. The position
   * and mask are copied into host scalar state; no CUDA transfer occurs.
   *
   * @param[in] referencePosition Exactly three finite coordinates in
   * `[x, y, z]` order and angstrom units.
   * @param[in] referenceMask Exactly three dimensionless values in
   * `[x, y, z]` order. Each value must be zero or one, and at least one value
   * must be one.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if either vector length is not three,
   * a position is NaN or infinite, a mask value is not zero or one, or every
   * mask value is zero.
   * @throws std::bad_alloc If an error diagnostic cannot be allocated.
   * @throws std::length_error If an error diagnostic exceeds an
   * implementation-defined limit.
   *
   * @post On success, subsequent evaluations use independent scalar copies of
   * the supplied position and mask.
   */
  void setReferencePosition(const std::vector<double> &referencePosition,
                            const std::vector<int> &referenceMask);

  /**
   * @brief Sets the target radial distance from the reference position.
   *
   * A zero target selects `k * r^2`. A positive target selects
   * `k * (r - referenceDistance)^2`. The scalar is stored on the host and no
   * CUDA transfer occurs.
   *
   * @param[in] referenceDistance Finite, non-negative target distance in
   * angstroms.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `referenceDistance` is negative,
   * NaN, or infinite.
   * @throws std::bad_alloc If an error diagnostic cannot be allocated.
   * @throws std::length_error If an error diagnostic exceeds an
   * implementation-defined limit.
   *
   * @post On success, subsequent evaluations use `referenceDistance`.
   * @note For a positive target and an exactly zero current displacement, the
   * implementation assigns a zero gradient at the nondifferentiable point.
   */
  void setReferenceDistance(const double referenceDistance);

  /**
   * @brief Sets per-atom masses and enables mass weighting.
   *
   * The full atom-order array is copied into owned host/device storage. The
   * selected weights are then rebuilt from the copied masses and transferred to
   * device memory. Only relative selected masses affect the normalized center;
   * multiplying every selected mass by the same positive factor has no effect.
   *
   * @param[in] masses Exactly one finite, non-negative value per atom, in atom
   * order and normally atomic mass units. The vector is borrowed for the call,
   * copied, and not retained by reference. Selected masses must have positive
   * total weight.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the vector length differs from the
   * atom count, a mass is negative, NaN, or infinite, or selected masses sum to
   * zero.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if mass or
   * selected-weight storage cannot be allocated, transferred, or synchronized.
   * @throws std::bad_alloc If host, device-metadata, or diagnostic storage
   * cannot be allocated.
   * @throws std::length_error If copied storage or a diagnostic exceeds an
   * implementation-defined limit.
   *
   * @post On success, mass weighting is enabled and host/device selected
   * weights correspond to `masses`.
   * @warning Failure after mass assignment begins can leave new masses or a new
   * weighting mode paired with older selected device arrays.
   * @warning Successful rebuilding can invalidate previously borrowed selected
   * array or partial-sum device pointers.
   */
  void setMasses(const std::vector<double> &masses);

  /**
   * @brief Enables mass weighting or restores uniform weighting.
   *
   * The selected-weight array is rebuilt and copied to device memory. `true`
   * uses the current per-atom masses; `false` assigns every selected atom
   * weight one. The stored mass array is retained when weighting is disabled.
   *
   * @param[in] useMassWeighting `true` to use configured masses or `false` to
   * use uniform dimensionless weights.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if weighting is enabled while the
   * selected masses have zero total weight.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if selected
   * arrays cannot be allocated, transferred, or synchronized.
   * @throws std::bad_alloc If selected-index, selected-weight, or diagnostic
   * storage cannot be allocated.
   * @throws std::length_error If rebuilt storage or a diagnostic exceeds an
   * implementation-defined limit.
   *
   * @post On success, selected device weights match the requested mode.
   * @warning Failure after the mode flag changes can leave the mode and
   * selected device arrays inconsistent.
   * @warning Successful rebuilding can invalidate previously borrowed selected
   * array or partial-sum device pointers.
   */
  void setMassWeighting(const bool useMassWeighting);

  /**
   * @brief Sets the orthorhombic periodic box lengths.
   *
   * The three values are copied into host scalar state. They are used to unwrap
   * selected atoms relative to the first selected atom and to minimum-image the
   * active center-to-reference displacement. No CUDA transfer or
   * synchronization occurs, and changing the box does not rescale the reference
   * position.
   *
   * @param[in] boxDimensions Exactly three finite positive lengths in
   * `[x, y, z]` order and angstrom units.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the vector length is not three or
   * any length is non-positive, NaN, or infinite.
   * @throws std::bad_alloc If an error diagnostic cannot be allocated.
   * @throws std::length_error If an error diagnostic exceeds an
   * implementation-defined limit.
   *
   * @post On success, subsequent nonzero-force evaluations use the supplied
   * orthorhombic box.
   * @warning Triclinic box vectors and angles are not represented by this API.
   */
  void setBoxDimensions(const std::vector<double> &boxDimensions);

public:
  /**
   * @brief Returns shared ownership of energy and virial storage.
   *
   * No transfer or synchronization occurs. The object contains the named
   * energy term `hmcm` and periodic virial work storage.
   *
   * @return A copied shared owner of the restraint's `CudaEnergyVirial` object.
   * The object can outlive this restraint, but no restraint remains to update
   * it after destruction.
   */
  std::shared_ptr<CudaEnergyVirial> getEnergyVirial(void);

  /**
   * @brief Returns shared ownership of device force storage.
   *
   * The force array contains `m_NumAtoms` vectors in structure-of-arrays
   * component layout. For the public fixed-point specialization, multiply raw
   * `long long int` values by `INV_FORCE_SCALE` to recover potential-gradient
   * values in kilocalories per mole per angstrom.
   *
   * @return A copied shared owner of the restraint's `Force<AT>` object. No
   * transfer or synchronization occurs.
   */
  std::shared_ptr<Force<AT>> getForce(void);

  /**
   * @brief Returns shared access to the private CUDA stream holder.
   *
   * The returned shared pointer owns only the host scalar containing the stream
   * handle. This restraint remains responsible for destroying the underlying
   * CUDA stream and resets the scalar during destruction.
   *
   * @return A copied shared owner of the stream holder. The CUDA stream remains
   * valid only while the restraint remains alive.
   * @warning Do not destroy the CUDA stream through the returned holder.
   */
  std::shared_ptr<cudaStream_t> getStream(void);

public:
  /**
   * @brief Initializes the restraint for an atom count and periodic box.
   *
   * This is the `ForceManager` initialization hook. It verifies that the
   * manager atom count matches construction and then delegates to
   * `setBoxDimensions()`. It does not allocate a new force array, rebuild the
   * selection, clear outputs, or set a separate initialized flag.
   *
   * @param[in] numAtoms Dimensionless atom count supplied by the manager. It
   * must equal the fixed construction count.
   * @param[in] boxDimensions Exactly three finite positive orthorhombic box
   * lengths in `[x, y, z]` order and angstrom units.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the atom count differs or the box
   * contract is violated.
   * @throws std::bad_alloc If an error diagnostic cannot be allocated.
   * @throws std::length_error If an error diagnostic exceeds an
   * implementation-defined limit.
   *
   * @post On success, subsequent nonzero-force calculations have valid box
   * dimensions.
   */
  void initialize(const int numAtoms, const std::vector<double> &boxDimensions);

  /**
   * @brief Enqueues clearing of force, energy, and virial outputs.
   *
   * Both clears are issued on the private CUDA stream. The method checks
   * immediate CUDA failures but does not synchronize before returning.
   *
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if output
   * clearing or an immediate CUDA operation fails.
   *
   * @post After the private stream reaches the queued clears, force, energy,
   * and virial storage is zero.
   */
  void clear(void);

  /**
   * @brief Enqueues center calculation, state calculation, and force updates.
   *
   * The method reads `m_NumAtoms` device `float4` coordinate-charge records,
   * computes the selected center under orthorhombic periodic imaging, computes
   * the energy and center gradient, and distributes that gradient to selected
   * atoms in proportion to their weights. All three kernels use the private
   * stream and return after immediate launch checks.
   *
   * Force output is always accumulated. Energy is accumulated only when
   * `calcEnergy` is true. Periodic virial shift-force records are accumulated
   * only when `calcVirial` is true. Disabled outputs retain their previous
   * values and can therefore be stale.
   *
   * @param[in] xyzq Borrowed non-NULL CUDA-device pointer to at least
   * `m_NumAtoms` contiguous `float4` records in `[x, y, z, charge]` order.
   * Coordinates use angstroms and the charge component is ignored. The pointer
   * must remain valid until the private stream completes all queued work.
   * @param[in] calcEnergy Whether to add the current restraint energy to the
   * device `hmcm` term.
   * @param[in] calcVirial Whether to add periodic shift-force records used by
   * later virial aggregation.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::NotInitialized` if the force constant is nonzero and
   * positive box dimensions have not been set.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if any kernel
   * launch reports an immediate CUDA failure.
   *
   * @pre `xyzq` is device-accessible, correctly sized, and does not overlap an
   * incompatible concurrent write.
   * @pre `clear()` has been ordered before this call when fresh output is
   * required.
   * @post Successful return means all requested kernels were enqueued; it does
   * not mean device execution has completed.
   * @note If the force constant is zero, the method returns without reading
   * `xyzq`, validating the box, or modifying outputs.
   * @warning A later asynchronous CUDA failure is reported only by a subsequent
   * CUDA synchronization or CUDA API operation.
   * @warning A launch failure after earlier kernels were enqueued can leave a
   * partially updated calculation once the stream completes.
   * @warning Virial center-shift bookkeeping represents only image counts from
   * `-1` through `1` on each axis. A larger center-to-reference image count is
   * omitted from that bookkeeping even though force and energy are evaluated.
   */
  void calcForce(const float4 *xyzq, const bool calcEnergy,
                 const bool calcVirial);

private:
  /**
   * @brief Rebuilds selected indices, selected weights, and reduction storage.
   *
   * The method materializes ascending selected atom indices, chooses mass or
   * uniform weights, verifies positive total selected weight, copies indices
   * and weights to device memory, resizes block partial sums, and clears
   * partial and restraint state. The `CudaContainer` operations perform
   * device-wide synchronization.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if selected total weight is not
   * positive.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if allocation,
   * transfer, clearing, or synchronization fails.
   * @throws std::bad_alloc If temporary or container host storage cannot be
   * allocated.
   * @throws std::length_error If temporary or container storage exceeds an
   * implementation-defined limit.
   *
   * @post On success, selected arrays and reduction storage match
   * `m_Selection` and `m_UseMassWeighting`.
   * @warning The method is sequential rather than transactional and can
   * invalidate previously borrowed device pointers through reallocation.
   */
  void updateSelectedAtoms(void);

private:
  /** Fixed dimensionless atom count established during construction. */
  int m_NumAtoms;
  /** Current selected count; updated before selected-array replacement. */
  int m_NumSelected;

  /** Harmonic coefficient in kilocalories per mole per square angstrom. */
  double m_ForceConstant;
  /** Non-negative target center distance in angstroms. */
  double m_ReferenceDistance;
  /** Selects configured masses or uniform weights. */
  bool m_UseMassWeighting;

  /** Owned host/device selected indices in ascending atom order. */
  CudaContainer<int> m_AtomIndices;
  /** Owned host/device selected weights corresponding to `m_AtomIndices`. */
  CudaContainer<double> m_AtomWeights;
  /** Owned host/device per-atom masses in atom order. */
  CudaContainer<double> m_Masses;

  /** One device reduction record per 256-thread selected-atom block. */
  CudaContainer<double4> m_PartialSums;

  /**
   * Two device state records. Element zero stores
   * `{dU/dCx, dU/dCy, dU/dCz, 1 / totalWeight}`. Element one stores the
   * center-to-reference periodic image counts as `{x, y, z, 0}`.
   */
  CudaContainer<double4> m_RestraintState;

  /** Host-only orthorhombic X box length in angstroms. */
  double m_BoxDimX;
  /** Host-only orthorhombic Y box length in angstroms. */
  double m_BoxDimY;
  /** Host-only orthorhombic Z box length in angstroms. */
  double m_BoxDimZ;

  /** Host-only reference position in angstroms; component `w` is unused. */
  double4 m_ReferencePosition;
  /** Dimensionless X-component activation flag, either zero or one. */
  int m_ReferenceMaskX;
  /** Dimensionless Y-component activation flag, either zero or one. */
  int m_ReferenceMaskY;
  /** Dimensionless Z-component activation flag, either zero or one. */
  int m_ReferenceMaskZ;

  /** Owned host selection used when rebuilding selected arrays. */
  AtomSelection m_Selection;

  /** Shared owner of the `hmcm` energy and periodic virial state. */
  std::shared_ptr<CudaEnergyVirial> m_EnergyVirial;
  /** Shared owner of device structure-of-arrays force storage. */
  std::shared_ptr<Force<AT>> m_Forces;
  /** Shared host holder for the privately owned CUDA stream handle. */
  std::shared_ptr<cudaStream_t> m_Stream;
};

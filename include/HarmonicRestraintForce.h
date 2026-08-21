// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author:  James E. Gonzales II, Samarjeet Prasad
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
 * @brief Computes absolute mass-weighted harmonic restraints for atoms.
 *
 * For atom `i`, the restraint stores an absolute reference position
 * `(x0_i, y0_i, z0_i)`, a non-negative force constant `k_i`, and a
 * non-negative mass or weight `m_i`. Define the Cartesian displacement as:
 *
 * @code
 * dx_i = x_i - x0_i
 * dy_i = y_i - y0_i
 * dz_i = z_i - z0_i
 * @endcode
 *
 * The selected atoms contribute the following energy:
 *
 * @code
 * U_harm = sum over i of
 *          k_i * m_i * (dx_i^2 + dy_i^2 + dz_i^2)
 * @endcode
 *
 * The device force array receives the following potential-gradient
 * components:
 *
 * @code
 * gradient_x_i = 2 * k_i * m_i * dx_i
 * gradient_y_i = 2 * k_i * m_i * dy_i
 * gradient_z_i = 2 * k_i * m_i * dz_i
 * @endcode
 *
 * apoCHARMM integration code subtracts this stored potential gradient, so the
 * corresponding physical restoring force has the opposite sign.
 * Contributions are accumulated atomically into the force and optional
 * `harm` energy buffers.
 *
 * Selection is implemented by setting the force constant of every unselected
 * atom to zero. Current coordinates are read from a borrowed device `float4`
 * array; its `w` charge component is ignored. Reference coordinates and
 * current coordinates are compared directly without minimum-image wrapping.
 *
 * Construction allocates per-atom host/device parameter mirrors, one CUDA
 * stream, one fixed-point force array, and one energy-virial object. Parameter
 * setters copy caller data and perform synchronous host-to-device transfers.
 * Force calculation and output clearing enqueue work on the private stream
 * without synchronizing it before return.
 *
 * @tparam AT Device force-accumulator representation. The library explicitly
 * instantiates only `long long int`, using apoCHARMM's fixed-point force scale.
 * @tparam CT Arithmetic representation used inside the restraint kernel. The
 * library explicitly instantiates `float` and `double`.
 *
 * @warning The class provides no internal host-thread synchronization. Do not
 * overlap configuration, force calculation, clearing, or destruction.
 * @warning Do not copy or copy-assign this class. The current
 * compiler-generated operations combine copied parameter containers with
 * shared stream, force, and energy owners and do not define a safe independent
 * object lifecycle.
 * @warning Call `clear()` before the first calculation and before any
 * calculation that requires fresh rather than accumulated output.
 *
 * @see harmonic_restraint_force
 * @see AtomSelection
 * @see ForceManager
 */
template <typename AT, typename CT> class HarmonicRestraintForce {
public:
  /**
   * @brief Reports that this restraint does not contribute to virial
   * aggregation.
   */
  static constexpr bool contributesVirial = false;

public:
  /**
   * @brief Prevents construction without an explicit positive atom count.
   */
  HarmonicRestraintForce(void) = delete;

  /**
   * @brief Constructs a harmonic restraint for a fixed number of atoms.
   *
   * All atoms are initially selected. Force constants and reference
   * coordinates are initialized to zero, per-atom weights are initialized to
   * one, and the stored box dimensions are initialized to zero. Host and
   * device parameter mirrors are coherent when construction succeeds.
   *
   * The constructor also creates a private CUDA stream, allocates a
   * `Force<AT>` with `numAtoms` active vectors, and creates a
   * `CudaEnergyVirial` containing the named energy term `harm`. Output force
   * and energy storage must be cleared before its first use.
   *
   * @param[in] numAtoms Dimensionless atom count. The value must be greater
   * than zero and remains fixed for the lifetime of the restraint.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `numAtoms` is not positive.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if CUDA
   * allocation, parameter initialization, force allocation, transfer,
   * synchronization, or stream creation fails.
   * @throws std::bad_alloc If host storage, shared ownership state, or an
   * error diagnostic cannot be allocated.
   * @throws std::length_error If requested host storage or an error diagnostic
   * exceeds an implementation-defined limit.
   *
   * @post `getForce()->size() == numAtoms`.
   * @post Every atom is selected, every force constant is zero, every
   * reference position is `(0, 0, 0)`, and every weight is one.
   */
  HarmonicRestraintForce(const int numAtoms);

  /**
   * @brief Destroys the private CUDA stream and releases owned resources.
   *
   * The underlying CUDA stream is destroyed without propagating cleanup
   * failures. A `std::shared_ptr<cudaStream_t>` previously returned by
   * `getStream()` may keep the host-side stream holder alive, but the stored
   * CUDA stream value is set to `nullptr` during destruction.
   *
   * Shared force and energy objects returned by `getForce()` and
   * `getEnergyVirial()` may remain allocated while another shared owner
   * exists. They no longer have an owning restraint capable of calculating
   * into them.
   *
   * @post The private CUDA stream is no longer valid.
   */
  ~HarmonicRestraintForce(void) noexcept;

public:
  /**
   * @brief Sets the atoms that may receive harmonic restraint constants.
   *
   * The selection is copied and is not retained. Every atom that is
   * unselected by the new value has its stored force constant set to zero.
   * Selecting an atom does not restore a force constant that was zeroed by an
   * earlier selection change; call `setForceConstant()` or
   * `setForceConstants()` after expanding a selection.
   *
   * Updated force constants are transferred from host to device, followed by
   * a device-wide CUDA synchronization.
   *
   * @param[in] selection Selection with exactly the restraint's fixed atom
   * count. The object is borrowed for the call and remains unchanged.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the represented atom count
   * differs from the restraint atom count.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if the
   * force-constant transfer or synchronization fails.
   * @throws std::bad_alloc If selection copying or diagnostic construction
   * cannot allocate storage.
   * @throws std::length_error If copied selection storage or a diagnostic
   * exceeds an implementation-defined limit.
   *
   * @post On success, the restraint owns an independent copy of `selection`
   * and every unselected atom has a zero force constant.
   * @warning The operation is not transactional after selection copying
   * begins. An allocation or CUDA failure can leave the copied selection,
   * host constants, and device constants representing different stages.
   */
  void setSelection(const AtomSelection &selection);

  /**
   * @brief Sets one force constant for every currently selected atom.
   *
   * The method first sets every per-atom force constant to zero and then
   * assigns `forceConstant` to the selected atoms. Unselected atoms therefore
   * remain disabled.
   *
   * With weights supplied in atomic mass units, `forceConstant` has units of
   * kilocalories per mole per square angstrom per atomic mass unit. When every
   * weight is one, its numerical effect is the ordinary per-atom harmonic
   * constant in kilocalories per mole per square angstrom.
   *
   * @param[in] forceConstant Finite, non-negative force constant applied to
   * every selected atom.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `forceConstant` is negative, NaN,
   * or infinite.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if clearing,
   * transferring, or synchronizing force constants fails.
   * @throws std::bad_alloc If selected-index storage or an error diagnostic
   * cannot be allocated.
   * @throws std::length_error If selected-index storage or a diagnostic
   * exceeds an implementation-defined limit.
   *
   * @post On success, selected atoms have `forceConstant` and all other atoms
   * have zero force constants.
   * @note The current implementation performs a synchronized zero fill and a
   * second synchronized host-to-device transfer after assigning selected
   * values.
   * @warning A failure after the zero fill can leave all constants zero or can
   * leave the host and device mirrors divergent.
   */
  void setForceConstant(const double forceConstant);

  /**
   * @brief Sets one force constant per atom.
   *
   * The vector is indexed in atom order and must contain exactly one value for
   * every atom. Only values belonging to currently selected atoms are
   * validated and used. Values at unselected indices are ignored, and the
   * corresponding stored constants are set to zero.
   *
   * With weights supplied in atomic mass units, each selected value has units
   * of kilocalories per mole per square angstrom per atomic mass unit.
   *
   * @param[in] forceConstants Contiguous host vector in atom order. Its length
   * must equal the fixed atom count. Every selected value must be finite and
   * non-negative. The vector is borrowed for the call and is not retained.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the vector length differs from
   * the atom count or a selected value is negative, NaN, or infinite.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if clearing,
   * transferring, or synchronizing force constants fails.
   * @throws std::bad_alloc If selected-index storage or diagnostic
   * construction cannot allocate storage.
   * @throws std::length_error If selected-index storage or a diagnostic
   * exceeds an implementation-defined limit.
   *
   * @post On success, every selected atom has the corresponding input value
   * and every unselected atom has a zero force constant.
   * @post A validation failure before mutation leaves the stored constants
   * unchanged.
   * @note The current implementation performs a synchronized zero fill and a
   * second synchronized host-to-device transfer after assigning selected
   * values.
   * @warning A CUDA or allocation failure after mutation starts can leave the
   * host and device mirrors at different stages.
   */
  void setForceConstants(const std::vector<double> &forceConstants);

  /**
   * @brief Sets one Cartesian reference position per atom.
   *
   * Components `x`, `y`, and `z` are copied into the host reference mirror in
   * atom order. The existing per-atom mass or weight stored in component `w`
   * is preserved. The complete `double4` reference array is then transferred
   * to device memory, followed by a device-wide CUDA synchronization.
   *
   * Reference positions are absolute. They are not wrapped, imaged, or
   * rescaled when box dimensions change.
   *
   * @param[in] referenceCoordinates Contiguous host vector containing exactly
   * one finite `double3` position per atom, in angstroms. The vector is
   * borrowed for the call and is not retained.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the vector length differs from
   * the atom count or any Cartesian component is NaN or infinite.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if the
   * host-to-device transfer or synchronization fails.
   * @throws std::bad_alloc If an error diagnostic cannot be allocated.
   * @throws std::length_error If an error diagnostic exceeds an
   * implementation-defined limit.
   *
   * @post On success, host and device reference positions equal the supplied
   * values and all existing weights are preserved.
   * @warning Validation and host mutation are interleaved by atom. If a later
   * coordinate is invalid, an earlier host prefix may already contain new
   * values while the device mirror still contains the previous array.
   */
  void
  setReferenceCoordinates(const std::vector<double3> &referenceCoordinates);

  /**
   * @brief Sets reference positions from nested Cartesian rows.
   *
   * The outer vector must contain one row per atom. Every row must contain
   * exactly three values in `[x, y, z]` order and angstrom units. The method
   * copies the rows into a temporary `double3` vector and delegates to the
   * primary reference-coordinate setter.
   *
   * @param[in] referenceCoordinates Host rows in atom order. The outer vector
   * must have the fixed atom count, and every row must have exactly three
   * finite values. No caller storage is retained.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the row count differs from the
   * atom count, a row does not have exactly three values, or a Cartesian value
   * is NaN or infinite.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if the
   * delegated host-to-device transfer or synchronization fails.
   * @throws std::bad_alloc If the temporary coordinate vector or an error
   * diagnostic cannot be allocated.
   * @throws std::length_error If temporary storage or an error diagnostic
   * exceeds an implementation-defined limit.
   *
   * @post On success, the result is identical to calling the `double3`
   * overload with the converted rows.
   * @post An outer-size or row-size validation failure leaves restraint state
   * unchanged.
   * @warning A non-finite value detected by the delegated primary setter can
   * leave an earlier host prefix updated while the device mirror remains old.
   */
  void setReferenceCoordinates(
      const std::vector<std::vector<double>> &referenceCoordinates);

  /**
   * @brief Sets the per-atom masses or restraint weights.
   *
   * Each value is stored in component `w` of the corresponding reference
   * `double4`. Existing reference-position components are preserved. The
   * complete reference array is transferred to device memory, followed by a
   * device-wide CUDA synchronization.
   *
   * A zero value is accepted and suppresses the harmonic energy and
   * force-array contribution for that atom, even when the atom is selected and
   * has a nonzero force constant.
   *
   * @param[in] masses Contiguous host vector containing exactly one finite,
   * non-negative value per atom, normally in atomic mass units. The vector is
   * borrowed for the call and is not retained.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the vector length differs from
   * the atom count or any value is negative, NaN, or infinite.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if the
   * host-to-device transfer or synchronization fails.
   * @throws std::bad_alloc If an error diagnostic cannot be allocated.
   * @throws std::length_error If an error diagnostic exceeds an
   * implementation-defined limit.
   *
   * @post On success, host and device weight values equal `masses` and all
   * reference-position components are preserved.
   * @warning Validation and host mutation are interleaved by atom. If a later
   * mass is invalid, an earlier host prefix may already contain new values
   * while the device mirror still contains the previous array.
   */
  void setMasses(const std::vector<double> &masses);

public:
  /**
   * @brief Returns shared access to the private CUDA stream holder.
   *
   * The returned `std::shared_ptr` is a copied shared owner of the host scalar
   * containing the CUDA stream handle. No stream synchronization is performed.
   *
   * The restraint remains responsible for destroying the underlying CUDA
   * stream. Its destructor writes `nullptr` into the shared scalar before
   * releasing its own shared owner.
   *
   * @return A copied shared owner of the stream holder. The pointer is
   * non-null for a successfully constructed live restraint.
   * @warning Retaining the returned `std::shared_ptr` does not extend the
   * lifetime of the underlying CUDA stream beyond restraint destruction.
   */
  std::shared_ptr<cudaStream_t> getStream(void);

public:
  /**
   * @brief Validates manager state and stores current box dimensions.
   *
   * This method is the initialization hook used by `ForceManager`. It verifies
   * that the manager atom count equals the immutable construction count and
   * then delegates to `setBoxDimensions()`.
   *
   * The call does not allocate new per-atom storage, clear outputs, mark a
   * separate initialized state, or alter reference coordinates. Direct native
   * force calculation is not gated by a successful prior call, although
   * manager-driven use calls this hook as part of subscription initialization.
   *
   * @param[in] numAtoms Dimensionless atom count supplied by the manager. It
   * must equal the construction count.
   * @param[in] boxDimensions Three finite, positive box lengths in
   * `[x, y, z]` order and angstrom units.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `numAtoms` differs from the fixed
   * atom count or `boxDimensions` violates the box contract.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if a changed
   * box cannot be transferred or synchronized.
   * @throws std::bad_alloc If box copying or diagnostic construction cannot
   * allocate storage.
   * @throws std::length_error If box storage or an error diagnostic exceeds an
   * implementation-defined limit.
   *
   * @post On success, the stored box dimensions equal `boxDimensions`.
   */
  void initialize(const int numAtoms, const std::vector<double> &boxDimensions);

  /**
   * @brief Enqueues clearing of the force, energy, and virial output storage.
   *
   * Device clearing is enqueued on the restraint's private CUDA stream. The
   * method does not synchronize that stream before return. The host buffer
   * owned by `CudaEnergyVirial` is zeroed during the call.
   *
   * Call this method before the first calculation and before every calculation
   * that requires fresh output, because `calcForce()` atomically adds rather
   * than replaces force and energy contributions.
   *
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if output
   * allocation, clearing, or immediate CUDA launch checking fails.
   * @throws std::bad_alloc If output or diagnostic storage cannot be
   * allocated.
   * @throws std::length_error If output or diagnostic storage exceeds an
   * implementation-defined limit.
   *
   * @post After all previously ordered work on `getStream()` completes, the
   * restraint force array and energy-virial device buffer are zero.
   * @warning Clearing the energy buffer and force buffer is sequential rather
   * than transactional. A failure can leave only one output object cleared.
   */
  void clear(void);

  /**
   * @brief Enqueues harmonic force-array and optional energy accumulation.
   *
   * One logical input record is read for each atom. Components `x`, `y`, and
   * `z` are Cartesian coordinates in angstroms; `w` is ignored. Force-array
   * contributions are always accumulated. The `harm` energy is accumulated
   * only when `calcEnergy` is `true`.
   *
   * The method launches one of four compile-time kernel variants on the
   * private stream and performs only the immediate kernel-launch check. It
   * does not synchronize the stream. Asynchronous execution errors can surface
   * from a later CUDA synchronization or dependent operation.
   *
   * @param[in] xyzq Borrowed CUDA-device pointer to at least the fixed atom
   * count of contiguous `float4` records. The pointer may not be null and must
   * remain valid until the enqueued kernel completes.
   * @param[in] calcEnergy `true` to atomically add the harmonic energy to the
   * named `harm` device term; `false` to leave that energy term unchanged.
   * @param[in] calcVirial Accepted for the force-manager interface but ignored
   * because this restraint does not currently calculate a virial contribution.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if CUDA reports
   * an immediate kernel-launch failure.
   * @throws std::bad_alloc If construction of a CUDA failure diagnostic cannot
   * allocate storage.
   * @throws std::length_error If a CUDA failure diagnostic exceeds an
   * implementation-defined limit.
   *
   * @pre The restraint parameters required by the calculation have completed
   * their host-to-device transfers.
   * @pre `xyzq` addresses at least the fixed atom count of live device records.
   * @pre Call `clear()` first when fresh rather than accumulated output is
   * required.
   * @post On successful immediate launch checking, calculation is enqueued on
   * `getStream()` but need not be complete.
   * @warning This performance-critical method does not validate `xyzq`, output
   * storage, parameter sizes, or stream state.
   */
  void calcForce(const float4 *xyzq, const bool calcEnergy,
                 const bool calcVirial);

  /**
   * @brief Stores the current orthogonal box lengths.
   *
   * The vector must contain `[x, y, z]` box lengths in angstroms. When all
   * three values exactly equal the stored values, the method returns without a
   * transfer. Otherwise, the vector is copied into the host mirror and then
   * transferred to the device with a device-wide CUDA synchronization.
   *
   * The current force kernel does not read these dimensions, apply
   * minimum-image displacement, or rescale stored reference positions.
   *
   * @param[in] boxDimensions Host vector containing exactly three finite,
   * strictly positive lengths in `[x, y, z]` order and angstrom units. The
   * vector is borrowed for the call and is not retained.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the vector does not contain
   * exactly three values or any value is non-positive, NaN, or infinite.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if copying or
   * synchronizing a changed box fails.
   * @throws std::bad_alloc If host copying or diagnostic construction cannot
   * allocate storage.
   * @throws std::length_error If host storage or an error diagnostic exceeds
   * an implementation-defined limit.
   *
   * @post On success, host and device box mirrors contain `boxDimensions`.
   * @post Input validation failure leaves the stored box unchanged.
   * @warning A CUDA failure after host assignment can leave the host and
   * device box mirrors divergent.
   */
  void setBoxDimensions(const std::vector<double> &boxDimensions);

  /**
   * @brief Returns shared ownership of the device force-array storage.
   *
   * The force object contains the fixed atom count of vectors in
   * structure-of-arrays form: X components begin at offset zero, Y components
   * begin at `stride()`, and Z components begin at `2 * stride()`.
   *
   * For the supported `AT = long long int` instantiations, values use
   * apoCHARMM's fixed-point scale of `2^40`. Divide each stored integer by
   * `2^40` to recover the force-array quantity in kilocalories per mole per
   * angstrom.
   *
   * No transfer, clearing, or stream synchronization is performed.
   *
   * @return A copied shared owner of the force object. The pointer is non-null
   * for a successfully constructed restraint and can outlive the restraint.
   * @warning A retained force object does not retain a valid calculation
   * stream or a restraint capable of updating it.
   * @warning Contents are additive and can be stale or uninitialized unless
   * the caller orders `clear()`, `calcForce()`, and stream synchronization
   * correctly.
   */
  std::shared_ptr<Force<AT>> getForce(void);

  /**
   * @brief Returns shared ownership of the energy-virial output object.
   *
   * The object contains the named scalar energy term `harm`, in kilocalories
   * per mole. The restraint does not populate a virial contribution. Device
   * energy is accumulated by `calcForce()` only when energy calculation is
   * requested.
   *
   * No device-to-host transfer or stream synchronization is performed.
   * `CudaEnergyVirial::getEnergy("harm")` reads the host mirror and therefore
   * remains stale until the caller orders an appropriate `copyToHost()` and
   * synchronization, or until a manager operation performs them.
   *
   * @return A copied shared owner of the energy-virial object. The pointer is
   * non-null for a successfully constructed restraint and can outlive the
   * restraint.
   * @warning A retained object does not retain a valid calculation stream or a
   * restraint capable of updating it.
   */
  std::shared_ptr<CudaEnergyVirial> getEnergyVirial(void);

private:
  int m_NumAtoms;

  CudaContainer<double> m_ForceConstants;

  CudaContainer<double4> m_ReferenceCoordinates;

  CudaContainer<double> m_BoxDimensions;

  AtomSelection m_Selection;

  std::shared_ptr<CudaEnergyVirial> m_EnergyVirial;

  std::shared_ptr<Force<AT>> m_Forces;

  std::shared_ptr<cudaStream_t> m_Stream;
};

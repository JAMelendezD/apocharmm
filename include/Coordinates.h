// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: Samarjeet Prasad, James E. Gonzales II
//
// ENDLICENSE

#pragma once

#include <vector>
#include <vector_types.h>

/**
 * @brief Owns host-resident Cartesian coordinates in two precisions.
 *
 * `Coordinates` stores one `double3` record and one `float3` record per atom in
 * atom order. Each record contains `x`, `y`, and `z` Cartesian coordinates in
 * angstroms. Both `std::vector` instances are owned by this object. The class
 * does not allocate CUDA device memory, transfer data, or synchronize a CUDA
 * stream.
 *
 * Constructors initialize both precision representations from one input.
 * Mutable access through `getCoordinatesDP()` or `getCoordinatesSP()` does not
 * update the other representation or the stored atom count. Callers using
 * those escape hatches must preserve the normal invariant that
 * `getNumAtoms()`, `getCoordinatesDP().size()`, and
 * `getCoordinatesSP().size()` describe the same number of atoms.
 *
 * Default construction represents unset state. Construction from an empty
 * input or a successful `setNumAtoms(0)` represents an explicitly initialized
 * zero-atom object. Compiler-generated copy operations deep-copy both vectors;
 * compiler-generated move operations transfer their storage. Assignment or
 * resizing can change the vector contents and invalidate element references,
 * pointers, and iterators.
 *
 * @warning The destructor is not virtual. Do not delete a derived
 * @ref CharmmCrd or @ref PDB object through a `Coordinates` pointer.
 * @warning The class provides no internal synchronization. Externally
 * serialize overlapping access whenever an operation can mutate the object or
 * either returned vector.
 * @see coordinates
 * @see CharmmContext
 */
class Coordinates {
public:
  /**
   * @brief Constructs an unset coordinate object.
   *
   * @post `getNumAtoms() == -1`.
   * @post Both coordinate vectors are empty.
   * @note The `-1` count is a dimensionless sentinel that distinguishes unset
   * state from an explicitly initialized zero-atom object.
   */
  Coordinates(void);

  /**
   * @brief Constructs both representations from double-precision coordinates.
   *
   * The double-precision records are copied unchanged. Each component is also
   * converted with `static_cast<float>` for the single-precision
   * representation. No finiteness or physical-range validation is performed.
   *
   * @param[in] coords Borrowed host vector containing one `double3` record per
   * atom in `[x, y, z]` order and angstrom units. The values are copied and the
   * vector is not retained.
   * @throws std::bad_alloc If either owned vector cannot allocate storage.
   * @throws std::length_error If an owned vector rejects the requested length.
   *
   * @pre `coords.size()` is at most `INT_MAX` because the atom count is stored
   * as `int`.
   * @pre Every component is representable as `float`.
   * @post Both owned vectors contain `coords.size()` records and
   * `getNumAtoms() == static_cast<int>(coords.size())`.
   */
  Coordinates(const std::vector<double3> &coords);

  /**
   * @brief Constructs both representations from single-precision coordinates.
   *
   * The single-precision records are copied unchanged. Each component is also
   * converted to `double` for the double-precision representation. No
   * finiteness or physical-range validation is performed.
   *
   * @param[in] coords Borrowed host vector containing one `float3` record per
   * atom in `[x, y, z]` order and angstrom units. The values are copied and the
   * vector is not retained.
   * @throws std::bad_alloc If either owned vector cannot allocate storage.
   * @throws std::length_error If an owned vector rejects the requested length.
   *
   * @pre `coords.size()` is at most `INT_MAX` because the atom count is stored
   * as `int`.
   * @post Both owned vectors contain `coords.size()` records and
   * `getNumAtoms() == static_cast<int>(coords.size())`.
   */
  Coordinates(const std::vector<float3> &coords);

  /**
   * @brief Constructs both representations from nested double rows.
   *
   * Each outer element represents one atom and must contain exactly three
   * values in `[x, y, z]` order. Values use angstroms. Double values are copied
   * to `double3` records and converted with `static_cast<float>` for the
   * `float3` records.
   *
   * @param[in] coords Borrowed host rows in atom order. All values are copied
   * and no row or outer-vector storage is retained.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if any row does not contain exactly
   * three values.
   * @throws std::bad_alloc If vector or error-diagnostic storage cannot be
   * allocated.
   * @throws std::length_error If a vector or diagnostic exceeds an
   * implementation limit.
   *
   * @pre `coords.size()` is at most `INT_MAX` because the atom count is stored
   * as `int`.
   * @pre Every component is representable as `float`.
   * @post Both owned vectors contain `coords.size()` records and
   * `getNumAtoms() == static_cast<int>(coords.size())`.
   */
  Coordinates(const std::vector<std::vector<double>> &coords);

  /**
   * @brief Constructs both representations from nested float rows.
   *
   * Each outer element represents one atom and must contain exactly three
   * values in `[x, y, z]` order. Values use angstroms. Float values are copied
   * to `float3` records and converted to `double` for the `double3` records.
   *
   * @param[in] coords Borrowed host rows in atom order. All values are copied
   * and no row or outer-vector storage is retained.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if any row does not contain exactly
   * three values.
   * @throws std::bad_alloc If vector or error-diagnostic storage cannot be
   * allocated.
   * @throws std::length_error If a vector or diagnostic exceeds an
   * implementation limit.
   *
   * @pre `coords.size()` is at most `INT_MAX` because the atom count is stored
   * as `int`.
   * @post Both owned vectors contain `coords.size()` records and
   * `getNumAtoms() == static_cast<int>(coords.size())`.
   */
  Coordinates(const std::vector<std::vector<float>> &coords);

public:
  /**
   * @brief Sets the atom count and resizes both coordinate vectors.
   *
   * Existing elements in the retained prefix are preserved. New elements are
   * value-initialized by `std::vector::resize`. This host-only operation does
   * not transfer coordinates to a CUDA device or update another object.
   *
   * @param[in] numAtoms Dimensionless atom count. The value must be
   * non-negative.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `numAtoms` is negative.
   * @throws std::bad_alloc If either vector cannot allocate storage.
   * @throws std::length_error If either vector rejects the requested length.
   *
   * @post On success, `getNumAtoms() == numAtoms` and both vector sizes equal
   * `static_cast<std::size_t>(numAtoms)`.
   * @post If `numAtoms` is negative, the object is unchanged.
   * @warning Successful reallocation invalidates element references, pointers,
   * and iterators into the affected vector. References to the vector objects
   * themselves remain valid.
   * @warning The two resizes are sequential rather than transactional. An
   * allocation failure can leave the stored count and vector sizes
   * inconsistent.
   */
  void setNumAtoms(const int numAtoms);

public:
  /**
   * @brief Returns the stored atom count.
   *
   * @return `-1` for the default-constructed unset state, or a non-negative
   * dimensionless count after successful initialization.
   * @note The value is stored independently and is not recomputed from either
   * vector after mutable access.
   */
  int getNumAtoms(void) const;

  /**
   * @brief Returns the double-precision coordinate vector.
   *
   * @return A borrowed read-only reference to the owned host vector. Each
   * `double3` is one atom in `[x, y, z]` order and angstrom units. The vector
   * reference remains valid for the lifetime of this object. Element
   * references, pointers, and iterators follow `std::vector` invalidation
   * rules.
   * @note This call performs no copy, precision conversion, CUDA transfer, or
   * synchronization.
   */
  const std::vector<double3> &getCoordinatesDP(void) const;

  /**
   * @brief Returns mutable double-precision coordinate storage.
   *
   * @return A borrowed mutable reference to the owned host vector. Each
   * `double3` is one atom in `[x, y, z]` order and angstrom units. The vector
   * reference remains valid for the lifetime of this object. Element
   * references, pointers, and iterators follow `std::vector` invalidation
   * rules.
   * @warning Mutation is unchecked. Changing values does not update the
   * single-precision vector, and changing the vector length does not update
   * `m_NumAtoms` or the single-precision vector.
   */
  std::vector<double3> &getCoordinatesDP(void);

  /**
   * @brief Returns the single-precision coordinate vector.
   *
   * @return A borrowed read-only reference to the owned host vector. Each
   * `float3` is one atom in `[x, y, z]` order and angstrom units. The vector
   * reference remains valid for the lifetime of this object. Element
   * references, pointers, and iterators follow `std::vector` invalidation
   * rules.
   * @note This call performs no copy, precision conversion, CUDA transfer, or
   * synchronization.
   */
  const std::vector<float3> &getCoordinatesSP(void) const;

  /**
   * @brief Returns mutable single-precision coordinate storage.
   *
   * @return A borrowed mutable reference to the owned host vector. Each
   * `float3` is one atom in `[x, y, z]` order and angstrom units. The vector
   * reference remains valid for the lifetime of this object. Element
   * references, pointers, and iterators follow `std::vector` invalidation
   * rules.
   * @warning Mutation is unchecked. Changing values does not update the
   * double-precision vector, and changing the vector length does not update
   * `m_NumAtoms` or the double-precision vector.
   */
  std::vector<float3> &getCoordinatesSP(void);

protected:
  /** Stores the dimensionless atom count, or `-1` for unset state. */
  int m_NumAtoms;
  /** Owns double-precision host coordinates in atom-major `[x, y, z]` order. */
  std::vector<double3> m_CoordinatesDP;
  /** Owns single-precision host coordinates in atom-major `[x, y, z]` order. */
  std::vector<float3> m_CoordinatesSP;
};

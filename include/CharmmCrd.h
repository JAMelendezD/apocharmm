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

#include "Coordinates.h"

#include <string>

/**
 * @brief Owns Cartesian coordinates parsed from a CHARMM CRD or COR file.
 *
 * Construction reads and parses the file synchronously into the inherited
 * host-resident double- and single-precision coordinate vectors. Only Cartesian
 * coordinates are retained. Atom identifiers, residue metadata, segment data,
 * weights, and any other trailing fields are ignored. No file path is retained
 * after construction, and no CUDA allocation, transfer, or synchronization is
 * performed.
 *
 * Compiler-generated copy and move operations have the same owned-vector
 * semantics as @ref Coordinates. The inherited mutable accessors remain
 * unchecked escape hatches.
 *
 * @warning The inherited @ref Coordinates destructor is not virtual. Do not
 * delete this object through a `Coordinates` pointer.
 * @see coordinates
 * @see apo_charmm_crd
 * @see python_charmm_crd
 */
class CharmmCrd : public Coordinates {
public:
  /**
   * @brief Constructs coordinates by parsing one CHARMM CRD or COR file.
   *
   * Consecutive nonempty lines beginning with `*` in column one are treated as
   * title records. The first following line supplies the atom count in its
   * first space-delimited token. An exact uppercase `EXT` second token
   * selects extended records; otherwise standard records are assumed.
   *
   * Standard records read X, Y, and Z from one-based columns 21-30, 31-40, and
   * 41-50. Extended records read them from columns 41-60, 61-80, and 81-100.
   * The parsed values must be finite doubles and use angstroms. The file-name
   * extension is not validated.
   *
   * @param[in] fileName Borrowed path string. It must be nonempty. The path and
   * file contents are copied as needed during construction and are not retained
   * after return.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `fileName` is empty.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Runtime` if the file
   * cannot be opened, sized, sought, or read; the atom count is missing,
   * malformed, negative, or greater than `INT_MAX`; a required coordinate
   * record or field is missing; or a coordinate is malformed, out of range, or
   * non-finite.
   * @throws std::bad_alloc If file, vector, or diagnostic storage cannot be
   * allocated.
   * @throws std::length_error If file, vector, or diagnostic storage exceeds an
   * implementation limit.
   *
   * @post On success, `getNumAtoms()` is the parsed non-negative count and both
   * inherited coordinate vectors contain that many atom records.
   * @warning The parser verifies that each coordinate is a finite `double`, but
   * does not verify representability as `float` before constructing the
   * single-precision representation.
   */
  CharmmCrd(const std::string &fileName);

private:
  /**
   * @brief Parses the file and initializes both inherited coordinate vectors.
   *
   * This construction-only helper loads the complete file into host memory,
   * selects the standard or extended fixed-width layout, resizes both vectors,
   * and writes coordinates in file-record order.
   *
   * @param[in] fileName Borrowed nonempty path string used only for this call.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Runtime` for verified
   * file I/O, count, record-shape, and coordinate-value failures.
   * @throws std::bad_alloc If host or diagnostic storage cannot be allocated.
   * @throws std::length_error If host or diagnostic storage exceeds an
   * implementation limit.
   *
   * @pre The constructor has already rejected an empty path.
   * @post On success, the inherited count and both coordinate-vector lengths
   * equal the parsed atom count.
   * @note The helper performs no CUDA work and retains no file metadata.
   */
  void readCharmmCrdFile(const std::string &fileName);
};

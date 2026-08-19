// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

/**
 * @file
 * @brief Declares lexical token types and token records for atom selections.
 */

#pragma once

#include <string>

/**
 * @brief Identifies one lexical element in an atom-selection expression.
 *
 * The values are internal parser categories rather than a serialized ABI. Text
 * matching for recognized keywords and dotted operators is case-insensitive.
 * See @ref atom_selection for the user-facing language.
 */
enum class SelectionTokenType {
  /** Marks the byte position immediately after the complete expression. */
  End, // 0

  /** Represents an unclassified bare value. */
  Identifier, // 1
  /** Represents a bare signed decimal integer. */
  Integer, // 2
  /** Represents a bare decimal value containing one decimal point. */
  Real, // 3

  /** Represents `(`. */
  LeftParenthesis, // 4
  /** Represents `)`. */
  RightParenthesis, // 5
  /** Represents the inclusive range separator `:`. */
  Colon, // 6

  /** Represents the primary selection keyword `ALL`. */
  All, // 7
  /** Represents the primary selection keyword `NONE`. */
  None, // 8

  /** Represents the binary intersection operator `.AND.`. */
  And, // 9
  /** Represents the binary union operator `.OR.`. */
  Or, // 10
  /** Represents the prefix complement operator `.NOT.`. */
  Not, // 11

  /** Represents the prefix residue-expansion operator `.BYRES.`. */
  ByResidue, // 12
  /** Represents the prefix group-expansion operator `.BYGROUP.`. */
  ByGroup, // 13
  /** Represents the prefix direct-neighbor operator `.BONDED.`. */
  Bonded, // 14

  /** Represents selection by the PSF atom-name field. */
  Type, // 15
  /** Represents selection by the PSF atom-type field. */
  Chemical, // 16
  /** Represents selection by the PSF segment identifier. */
  SegmentIdentifier, // 17
  /** Represents selection by the PSF residue identifier. */
  ResidueIdentifier, // 18
  /** Represents selection by the PSF residue name. */
  ResidueName, // 19
  /** Represents a three-field segment, residue, and atom-name selection. */
  Atom, // 20
  /** Represents selection by one-based atom number. */
  ByNumber // 21
};

/**
 * @brief Owns one token and its source position.
 *
 * Token text is copied from the source expression. `pos` is a zero-based byte
 * offset, not a Unicode code-point index. The terminal token owns an empty
 * string and records the source byte length.
 */
struct SelectionToken {
  /** Stores the parser category assigned by @ref SelectionTokenizer. */
  SelectionTokenType type;
  /** Owns the exact source bytes for this token. */
  std::string text;
  /** Stores the zero-based source byte offset as an `int`. */
  int pos;
};

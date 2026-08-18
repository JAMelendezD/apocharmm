// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#pragma once

#include "AtomSelection.h"
#include "CharmmPSF.h"
#include "SelectionToken.h"

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

/**
 * @brief Evaluates tokenized CHARMM-style atom-selection expressions.
 *
 * The parser is instantiated for one call to `parse()`. It retains the supplied
 * PSF and token vector while building a result, constructs atom-to-residue and
 * atom-to-group lookup arrays, and evaluates the expression with separate
 * operator and selection stacks. `.AND.` binds more tightly than `.OR.`;
 * operators of equal precedence are reduced left to right. Prefix operators are
 * applied after their primary selection or parenthesized expression is read.
 *
 * Field values come from host-resident @ref CharmmPSF arrays. `TYPE` maps to
 * PSF atom names, `CHEM` to atom types, `SEGI` to segment identifiers, `RESI`
 * to decimal residue identifiers, `RESN` to residue names, and `BYNU` to
 * one-based atom numbers. The parser performs no CUDA allocation, transfer,
 * stream work, or synchronization.
 *
 * The parser provides no internal locking. The PSF must remain immutable and
 * internally consistent throughout parsing. The public entry point is intended
 * to consume the complete vector returned by @ref SelectionTokenizer, including
 * its terminal @ref SelectionTokenType::End token.
 *
 * @see atom_selection
 */
class SelectionParser {
public:
  /**
   * @brief Parses a complete token vector against a PSF.
   *
   * The shared pointer and token vector are passed by value, then moved into a
   * temporary parser. The returned @ref AtomSelection owns independent host
   * storage and retains neither the parser, the token vector, nor the PSF.
   *
   * @param[in] psf Shared pointer to the topology whose host metadata and
   * derived tables are read. The pointer is retained only until parsing
   * finishes and may not be null.
   * @param[in] tokens Complete token sequence, normally produced by
   * @ref SelectionTokenizer::tokenize. The vector must contain one reachable
   * terminal token after the expression.
   * @return A newly owned selection whose atom count equals
   * `psf->getNumAtoms()`.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::InvalidArgument` if `psf`
   * is null, the tokens encode an invalid expression or range, or a PSF
   * bonded-neighbor index lies outside the atom range.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::NotInitialized` if the PSF
   * atom count is negative.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::Runtime` if token-stack,
   * operator-stack, residue, group, bonded-connectivity, or field-token
   * invariants fail.
   * @throws std::bad_alloc If retained ownership, lookup arrays, stacks,
   * selections, field text, or diagnostics cannot allocate storage.
   * @throws std::length_error If a parser container, string, selection, or
   * diagnostic exceeds an implementation-defined limit.
   *
   * @pre Every per-atom PSF metadata array used by the requested fields has at
   * least `psf->getNumAtoms()` elements.
   * @pre Residue and group records use inclusive zero-based atom ranges, and
   * direct bonded-connectivity entries contain valid zero-based indices.
   * @post On success, the supplied PSF object is unchanged.
   * @warning A token vector with no reachable terminal token is an internal-use
   * contract violation and can report `ApoCharmmErrorCode::Runtime` rather than
   * a user syntax error.
   */
  static AtomSelection parse(std::shared_ptr<const CharmmPSF> psf,
                             std::vector<SelectionToken> tokens);

private:
  /**
   * @brief Constructs one parser and builds topology lookup arrays.
   *
   * @param[in] psf Shared pointer copied into the parser. The pointer may not
   * be null and its atom count must be initialized.
   * @param[in] tokens Token vector moved into parser-owned storage.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::InvalidArgument` if `psf`
   * is null.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::NotInitialized` if the PSF
   * atom count is negative.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::Runtime` if a residue or
   * group range is outside the atom range.
   */
  SelectionParser(std::shared_ptr<const CharmmPSF> psf,
                  std::vector<SelectionToken> tokens);

  /**
   * @brief Evaluates the parser-owned token sequence.
   *
   * @return A newly owned selection moved from the final selection-stack entry.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::InvalidArgument` for an
   * invalid expression.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::Runtime` for an internal
   * parser invariant failure.
   *
   * @post On success, exactly one selection result has been produced.
   */
  AtomSelection parse_impl(void);

private:
  /**
   * @brief Reads one primary selection beginning at the current token.
   *
   * @return An owned selection for `ALL`, `NONE`, `ATOM`, or one field token.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::InvalidArgument` when the
   * current token cannot begin a primary selection.
   */
  AtomSelection readPrimarySelection(void);

  /**
   * @brief Reads `ATOM` segment, residue, and atom-name values.
   *
   * @return The intersection of the corresponding `SEGI`, `RESI`, and `TYPE`
   * field selections.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::InvalidArgument` if any of
   * the three required values is missing.
   */
  AtomSelection readAtomSelection(void);

  /**
   * @brief Reads one field value or inclusive range.
   *
   * @param[in] fieldTokenType Field whose PSF values are tested.
   * @param[in] fieldName Borrowed diagnostic name for the field. The view is
   * not retained.
   * @return An owned field or range selection.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::InvalidArgument` if a
   * value or range endpoint is missing or a `BYNU` range is not integral.
   */
  AtomSelection readFieldSelection(const SelectionTokenType fieldTokenType,
                                   const std::string_view fieldName);

  /**
   * @brief Selects field values matching one wildcard pattern.
   *
   * @param[in] fieldTokenType Field whose PSF values are tested.
   * @param[in] fieldName Borrowed case-insensitive pattern. The view is not
   * retained.
   * @return An owned selection containing every matching zero-based atom index.
   */
  AtomSelection makeFieldSelection(const SelectionTokenType fieldTokenType,
                                   const std::string_view fieldName) const;

  /**
   * @brief Selects field values inside an inclusive range.
   *
   * Reversed endpoints are normalized. `BYNU` endpoints must be integers, are
   * clamped to the one-based atom-number range, and are converted to zero-based
   * indices. Other integer endpoints use numeric comparison; all other ranges
   * use case-insensitive lexicographic comparison.
   *
   * @param[in] fieldTokenType Field whose PSF values are tested.
   * @param[in] firstText Borrowed first endpoint. The view is not retained.
   * @param[in] lastText Borrowed second endpoint. The view is not retained.
   * @return An owned selection containing every atom inside the normalized
   * range.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::InvalidArgument` if a
   * `BYNU` endpoint is not an integer.
   */
  AtomSelection makeRangeSelection(const SelectionTokenType fieldTokenType,
                                   const std::string_view firstText,
                                   const std::string_view lastText) const;

private:
  /**
   * @brief Removes and returns the top selection-stack entry.
   *
   * @return The top selection moved into an owned result.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::Runtime` if the stack is
   * empty.
   */
  AtomSelection popSelection(void);

  /**
   * @brief Complements a selection across the complete PSF atom range.
   *
   * @param[in] selection Selection borrowed for the operation.
   * @return An owned selection containing exactly the previously unselected
   * atoms.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::Runtime` if the selection
   * atom count differs from the PSF atom count.
   */
  AtomSelection invertSelection(const AtomSelection &selection) const;

  /**
   * @brief Expands selected atoms to complete PSF residue intervals.
   *
   * @param[in] selection Selection borrowed for the operation.
   * @return An owned selection containing every atom in each represented
   * residue.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::Runtime` if a selected
   * atom has no valid residue mapping.
   */
  AtomSelection expandByResidue(const AtomSelection &selection) const;

  /**
   * @brief Expands selected atoms to complete PSF group intervals.
   *
   * @param[in] selection Selection borrowed for the operation.
   * @return An owned selection containing every atom in each represented group.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::Runtime` if a selected
   * atom has no valid group mapping.
   */
  AtomSelection expandByGroup(const AtomSelection &selection) const;

  /**
   * @brief Replaces selected atoms with their direct bonded neighbors.
   *
   * @param[in] selection Selection borrowed for the operation.
   * @return An owned selection containing the union of direct 1-2 neighbors.
   * The original atoms are included only when present in their own neighbor
   * sets.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::Runtime` if the PSF
   * connectivity-vector length differs from the atom count.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::InvalidArgument` if a
   * stored neighbor index is outside `[0, atom_count)`.
   */
  AtomSelection expandBonded(const AtomSelection &selection) const;

  /**
   * @brief Applies and removes the top operator-stack entry.
   *
   * Binary operators pop two selections. Prefix operators pop one selection.
   * The computed selection is pushed back onto the selection stack.
   *
   * @throws ApoCharmmError With `ApoCharmmErrorCode::Runtime` if either stack
   * lacks a required entry or the token is not applicable.
   */
  void applyTopOperator(void);

  /**
   * @brief Applies consecutive pending prefix operators.
   *
   * Operators are popped from the stack, so nested prefixes are applied from
   * the primary selection outward.
   */
  void applyPendingPrefixOperators(void);

  /**
   * @brief Reduces binary operators that bind at least as tightly as one token.
   *
   * @param[in] currOp Incoming binary operator used for precedence comparison.
   */
  void reduceBinaryOperatorsFor(const SelectionTokenType currOp);

  /**
   * @brief Reduces operators through the matching left parenthesis.
   *
   * @throws ApoCharmmError With `ApoCharmmErrorCode::InvalidArgument` if no
   * matching left parenthesis exists.
   */
  void reduceUntilLeftParenthesis(void);

  /**
   * @brief Reduces every operator remaining at end of input.
   *
   * @throws ApoCharmmError With `ApoCharmmErrorCode::InvalidArgument` if an
   * unmatched left parenthesis remains.
   */
  void reduceRemainingOperators(void);

private:
  /**
   * @brief Returns the current token without advancing.
   *
   * @return A borrowed const reference into `m_Tokens`. The reference remains
   * valid until the parser or token vector is destroyed or modified.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::Runtime` if the cursor is
   * outside the token vector.
   */
  const SelectionToken &peek(void) const;

  /**
   * @brief Returns the current token and advances the cursor.
   *
   * @return A borrowed const reference into `m_Tokens`. The reference remains
   * valid until the parser or token vector is destroyed or modified.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::Runtime` if the cursor is
   * outside the token vector.
   */
  const SelectionToken &consume(void);

  /**
   * @brief Consumes one token accepted as a selection value.
   *
   * @param[in] context Borrowed field or construct name used in diagnostics.
   * The view is not retained.
   * @return A borrowed const reference into `m_Tokens`.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::InvalidArgument` if the
   * current token is not an accepted value token.
   */
  const SelectionToken &consumeSelectionValue(const std::string_view context);

  /**
   * @brief Consumes the current token when it has a requested type.
   *
   * @param[in] type Token category to match.
   * @return `true` after consuming a matching token; otherwise `false` without
   * advancing.
   */
  bool match(const SelectionTokenType type);

  /**
   * @brief Returns one atom's text value for a selection field.
   *
   * @param[in] fieldTokenType Field to read.
   * @param[in] atomIndex Zero-based atom index into PSF per-atom arrays.
   * @return An owned string containing the mapped field value. `RESI` and
   * `BYNU` values are formatted as decimal integers.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::Runtime` if the token is
   * not a field category.
   *
   * @pre Every mapped PSF array contains `atomIndex`.
   */
  std::string getFieldValue(const SelectionTokenType fieldTokenType,
                            const int atomIndex) const;

private:
  /**
   * @brief Builds the zero-based atom-to-residue lookup array.
   *
   * Every PSF residue record is interpreted as an inclusive `int2{x, y}` host
   * range. Atoms absent from all records retain the sentinel `-1`.
   *
   * @throws ApoCharmmError With `ApoCharmmErrorCode::Runtime` if a range is
   * negative, reversed, or outside the atom count.
   */
  void buildResidueIndex(void);

  /**
   * @brief Builds the zero-based atom-to-group lookup array.
   *
   * Every PSF group record is interpreted as an inclusive `int2{x, y}` host
   * range. Atoms absent from all records retain the sentinel `-1`.
   *
   * @throws ApoCharmmError With `ApoCharmmErrorCode::Runtime` if a range is
   * negative, reversed, or outside the atom count.
   */
  void buildGroupIndex(void);

private:
  /**
   * @brief Parses a complete signed decimal integer.
   *
   * @param[out] value Parsed value when the function returns `true`. Its value
   * is unspecified when the function returns `false`.
   * @param[in] str Borrowed text to parse. The view is not retained.
   * @return `true` only when `std::stoll` consumes every source byte without
   * throwing; otherwise `false`.
   */
  static bool parseInteger(long long int &value, const std::string_view str);

  /**
   * @brief Tests one PSF field value against a wildcard pattern.
   *
   * Matching is case-insensitive. `*` matches zero or more arbitrary bytes, `#`
   * matches zero or more decimal digits, `%` matches exactly one arbitrary
   * byte, and `+` matches exactly one decimal digit.
   *
   * @param[in] value Borrowed field value. The view is not retained.
   * @param[in] pattern Borrowed pattern. The view is not retained.
   * @return `true` when the complete value matches the complete pattern.
   */
  static bool doSelectionPatternsMatch(const std::string_view value,
                                       const std::string_view pattern);

  /**
   * @brief Tests whether a token category may be used as a field value.
   *
   * @param[in] type Token category to classify.
   * @return `true` for identifiers, numeric tokens, and keyword tokens accepted
   * literally after a field token.
   */
  static bool isSelectionValueToken(const SelectionTokenType type);

  /**
   * @brief Tests whether a token category begins a field selection.
   *
   * @param[in] type Token category to classify.
   * @return `true` for `TYPE`, `CHEM`, `SEGI`, `RESI`, `RESN`, or `BYNU`.
   */
  static bool isFieldToken(const SelectionTokenType type);

  /**
   * @brief Tests whether a token category begins a primary selection.
   *
   * @param[in] type Token category to classify.
   * @return `true` for `ALL`, `NONE`, `ATOM`, or a field token.
   */
  static bool isPrimarySelectionStart(const SelectionTokenType type);

  /**
   * @brief Tests whether a token category is a prefix operator.
   *
   * @param[in] type Token category to classify.
   * @return `true` for `.NOT.`, `.BYRES.`, `.BYGROUP.`, or `.BONDED.`.
   */
  static bool isPrefixOperator(const SelectionTokenType type);

  /**
   * @brief Tests whether a token category is a binary operator.
   *
   * @param[in] type Token category to classify.
   * @return `true` for `.AND.` or `.OR.`.
   */
  static bool isBinaryOperator(const SelectionTokenType type);

  /**
   * @brief Returns the binary precedence of a token category.
   *
   * @param[in] type Token category to inspect.
   * @return `2` for `.AND.`, `1` for `.OR.`, and `0` otherwise.
   */
  static int getOperatorPrecedence(const SelectionTokenType type);

  /**
   * @brief Returns the canonical diagnostic name of a field token.
   *
   * @param[in] type Field token category to name.
   * @return A borrowed view of program-lifetime static text.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::Runtime` if `type` is not
   * a field token.
   */
  static const std::string_view getFieldName(const SelectionTokenType type);

private:
  /**
   * @brief Throws an invalid-argument error at the current token position.
   *
   * @param[in] message Borrowed diagnostic prefix. The text is copied into the
   * exception diagnostic and is not retained.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::InvalidArgument` always,
   * after appending the current token's zero-based byte position.
   */
  [[noreturn]] void throwErrorAtCurrent(const std::string_view message) const;

private:
  /** Retains the topology throughout one parse. */
  std::shared_ptr<const CharmmPSF> m_Psf;
  /** Owns the complete token sequence, including the terminal token. */
  std::vector<SelectionToken> m_Tokens;
  /** Stores the zero-based index of the current token. */
  std::size_t m_Position;

  /** Owns pending prefix operators, binary operators, and left parentheses. */
  std::vector<SelectionTokenType> m_OperatorStack;
  /** Owns partial selection values during expression reduction. */
  std::vector<AtomSelection> m_SelectionStack;

  /** Maps each atom to a residue index, or `-1` when unmapped. */
  std::vector<int> m_ResidueIndex;
  /** Maps each atom to a group index, or `-1` when unmapped. */
  std::vector<int> m_GroupIndex;
};

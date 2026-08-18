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

#include "SelectionToken.h"

#include <string_view>
#include <vector>

/**
 * @brief Converts atom-selection expression bytes into owned parser tokens.
 *
 * The tokenizer is stateless and operates synchronously on the calling host
 * thread. It skips whitespace, recognizes parentheses and range separators,
 * classifies supported dotted operators and bare keywords without regard to
 * ASCII letter case, and appends exactly one @ref SelectionTokenType::End
 * token. Every nonterminal token owns a copy of its source spelling.
 *
 * The tokenizer performs no PSF access, CUDA operation, transfer, stream work,
 * or synchronization. It provides no internal locking, but independent calls
 * share no mutable state.
 *
 * @see atom_selection
 */
class SelectionTokenizer {
public:
  /**
   * @brief Tokenizes one atom-selection expression.
   *
   * Source positions are zero-based byte offsets. The current implementation
   * stores them in `int`, so expressions longer than `INT_MAX` bytes cannot
   * preserve positions reliably.
   *
   * @param[in] selectionString Borrowed expression view. The referenced bytes
   * need not be null-terminated and are not retained. Embedded control bytes,
   * including `\0`, are rejected unless classified as whitespace by the C
   * locale.
   * @return A newly owned token vector in source order followed by one terminal
   * token whose position is `selectionString.size()` after narrowing to `int`.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::InvalidArgument` for an
   * unterminated or unknown dotted operator, or for an unexpected source byte.
   * @throws std::bad_alloc If token, text, or diagnostic allocation fails.
   * @throws std::length_error If a token vector, token string, or diagnostic
   * exceeds an implementation-defined limit.
   *
   * @post The input bytes are unchanged and no reference to them is retained.
   * @warning The implementation does not reject source lengths greater than
   * `INT_MAX` before narrowing token positions.
   */
  static std::vector<SelectionToken>
  tokenize(const std::string_view selectionString);

private:
  /**
   * @brief Tests whether one byte can continue a bare token.
   *
   * @param[in] c Source byte to classify.
   * @return `false` for control characters, whitespace, parentheses, or `:`;
   * otherwise `true`.
   */
  static bool isBareTokenCharacter(const char c);

  /**
   * @brief Tests whether a token is a signed decimal integer.
   *
   * @param[in] str Borrowed token view. The text is not retained.
   * @return `true` for one or more decimal digits with an optional leading
   * minus sign; otherwise `false`.
   */
  static bool isInteger(const std::string_view str);

  /**
   * @brief Tests whether a token is a simple decimal value.
   *
   * @param[in] str Borrowed token view. The text is not retained.
   * @return `true` when the text contains at least one decimal digit, exactly
   * one decimal point, an optional leading minus sign, and no other bytes.
   *
   * @note Exponent notation and a leading plus sign are not recognized.
   */
  static bool isReal(const std::string_view str);

  /**
   * @brief Classifies a bare token as a keyword, number, or identifier.
   *
   * @param[in] str Borrowed nonempty token view. The text is not retained.
   * @return The recognized keyword or numeric category, or
   * @ref SelectionTokenType::Identifier.
   * @throws std::bad_alloc If uppercasing or trimming cannot allocate storage.
   * @throws std::length_error If temporary text exceeds a string limit.
   *
   * @note Keywords of four or more characters are currently classified from
   * their first four uppercase bytes. Prefer the canonical spellings documented
   * on @ref atom_selection.
   */
  static SelectionTokenType getBareTokenType(const std::string_view str);

  /**
   * @brief Classifies one complete dotted operator.
   *
   * @param[in] str Borrowed token view including both dots. The text is not
   * retained.
   * @return The matching case-insensitive dotted-operator category.
   * @throws ApoCharmmError With `ApoCharmmErrorCode::InvalidArgument` if the
   * spelling is not a supported dotted operator.
   * @throws std::bad_alloc If uppercase or diagnostic allocation fails.
   * @throws std::length_error If temporary or diagnostic text exceeds a string
   * limit.
   */
  static SelectionTokenType getDottedTokenType(const std::string_view str);
};

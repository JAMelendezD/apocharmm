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

#include <stdexcept>
#include <string>
#include <string_view>

/**
 * @brief Identifies the category of an apoCHARMM native failure.
 *
 * The category is stored independently from the human-readable diagnostic.
 * C ABI boundaries map each declared value to the corresponding @ref apo_status
 * value without relying on this enum's underlying integer representation.
 */
enum class ApoCharmmErrorCode {
  /** Indicates that the caller supplied an invalid argument. */
  InvalidArgument,
  /** Indicates a runtime failure without a more specific category. */
  Runtime,
  /** Indicates that a CUDA operation failed. */
  Cuda,
  /** Indicates that an operation requires state that is not initialized. */
  NotInitialized,
  /** Indicates that a recognized operation is not implemented. */
  NotImplemented
};

/**
 * @brief Returns the diagnostic name of an apoCHARMM error category.
 *
 * @param[in] code Category to name.
 * @return A borrowed view of program-lifetime static storage. Declared
 * categories return their exact enumerator names; any other value returns
 * `"Unknown"`.
 */
std::string_view
GetApoCharmmErrorCodeName(const ApoCharmmErrorCode code) noexcept;

/**
 * @brief Represents a categorized native error with source metadata.
 *
 * The exception owns copies of its human message, normalized source-file
 * spelling, source function, and source line. The inherited `what()` diagnostic
 * begins with the apoCHARMM ASCII banner and then reports the category,
 * message, source location, and function. The exception remains catchable as
 * `std::runtime_error` and `std::exception`.
 *
 * Views returned by the string accessors borrow storage from this object.
 * Copying the exception creates an independent logical copy of that metadata.
 * Moving transfers the metadata, and may invalidate views into the moved-from
 * exception.
 *
 * @warning Do not allow this C++ exception to cross a C ABI boundary. Use the
 * apoCHARMM C ABI guards, which translate it to an @ref apo_status and a
 * thread-local diagnostic.
 * @see apocharmm_error
 */
class ApoCharmmError : public std::runtime_error {
public:
  /**
   * @brief Constructs a categorized exception and renders its diagnostic.
   *
   * All character data is copied before the constructor returns. A source path
   * beneath the configured apoCHARMM source root is stored relative to that
   * root. Any other path containing `/` or `\` is reduced to its final path
   * component.
   *
   * @param[in] code Failure category to retain.
   * @param[in] message Human-readable diagnostic. The text may be empty or
   * multiline and is copied without normalization.
   * @param[in] sourceFile Source-file spelling to normalize and copy.
   * @param[in] sourceFunction Source-function spelling to copy unchanged.
   * @param[in] sourceLine Source-line value to store unchanged.
   * @throws std::bad_alloc If owned diagnostic storage cannot be allocated.
   * @throws std::length_error If an input or rendered diagnostic exceeds an
   * implementation-defined string limit.
   *
   * @note Prefer @ref APOCHARMM_THROW or @ref APOCHARMM_REQUIRE so the source
   * metadata is captured at the failure site.
   */
  ApoCharmmError(const ApoCharmmErrorCode code, const std::string_view message,
                 const std::string_view sourceFile,
                 const std::string_view sourceFunction, const int sourceLine);

public:
  /**
   * @brief Returns the retained failure category.
   *
   * @return The category supplied at construction.
   */
  ApoCharmmErrorCode getCode(void) const noexcept;

  /**
   * @brief Returns the retained human-readable message.
   *
   * @return A borrowed read-only view into storage owned by this exception. The
   * view remains valid until the exception is destroyed or assigned a new
   * value. Moving from the exception may also invalidate the view.
   */
  std::string_view getMessage(void) const noexcept;

  /**
   * @brief Returns the normalized source-file spelling.
   *
   * @return A borrowed read-only view into storage owned by this exception. The
   * view remains valid until the exception is destroyed or assigned a new
   * value. Moving from the exception may also invalidate the view.
   */
  std::string_view getSourceFile(void) const noexcept;

  /**
   * @brief Returns the retained source-function spelling.
   *
   * @return A borrowed read-only view into storage owned by this exception. The
   * view remains valid until the exception is destroyed or assigned a new
   * value. Moving from the exception may also invalidate the view.
   */
  std::string_view getSourceFunction(void) const noexcept;

  /**
   * @brief Returns the retained source-line value.
   *
   * @return The source-line value supplied at construction. The
   * @ref APOCHARMM_THROW macro supplies the compiler-provided `__LINE__` value.
   */
  int getSourceLine(void) const noexcept;

private:
  ApoCharmmErrorCode m_Code;
  std::string m_Message;
  std::string m_SourceFile;
  std::string m_SourceFunction;
  int m_SourceLine;
};

/**
 * @brief Throws a categorized apoCHARMM error at the current source location.
 *
 * The macro evaluates each argument exactly once and captures `__FILE__`,
 * `__func__`, and `__LINE__` from the expansion site.
 *
 * @param[in] errorCode Expression yielding the retained
 * @ref ApoCharmmErrorCode.
 * @param[in] message Expression convertible to `std::string_view`; its text is
 * copied by the exception constructor.
 * @throws ApoCharmmError Always when exception construction succeeds.
 * @throws std::bad_alloc If exception storage cannot be allocated.
 * @throws std::length_error If the diagnostic exceeds a string limit.
 *
 * @pre The macro is expanded in function scope so `__func__` is available.
 * @warning Do not use this macro from a `noexcept` cleanup path or destructor.
 */
#define APOCHARMM_THROW(errorCode, message)                                    \
  throw ::ApoCharmmError((errorCode), (message), __FILE__, __func__, __LINE__)

/**
 * @brief Throws a categorized apoCHARMM error unless a condition is true.
 *
 * The macro evaluates `condition` exactly once. It evaluates `errorCode` and
 * `message` only when the condition is false, then delegates source capture and
 * construction to @ref APOCHARMM_THROW.
 *
 * @param[in] condition Expression whose contextual Boolean value must be true.
 * @param[in] errorCode Expression yielding the retained
 * @ref ApoCharmmErrorCode when `condition` is false.
 * @param[in] message Expression convertible to `std::string_view`; evaluated
 * and copied only when `condition` is false.
 * @throws ApoCharmmError If `condition` is false and construction succeeds.
 * @throws std::bad_alloc If failure-path storage cannot be allocated.
 * @throws std::length_error If the failure diagnostic exceeds a string limit.
 *
 * @pre The macro is expanded in function scope so `__func__` is available.
 * @post If `condition` is true, execution continues and neither `errorCode` nor
 * `message` is evaluated.
 * @warning Do not use this macro from a `noexcept` cleanup path or destructor.
 */
#define APOCHARMM_REQUIRE(condition, errorCode, message)                       \
  do {                                                                         \
    if (!(condition)) {                                                        \
      APOCHARMM_THROW((errorCode), (message));                                 \
    }                                                                          \
  } while (false)

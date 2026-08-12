# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

"""
@brief Provides Python status constants and native-error translation.

The public constants and `ApoCharmmError` are re-exported by the top-level
`apocharmm` package. Native wrapper functions use `configure_status_function()`
to convert nonzero C ABI statuses into one Python exception type while
retaining the numeric status, operation context, and copied native diagnostic.

@anchor python_apocharmm_error_module
"""

from collections.abc import Callable, Sequence
import ctypes
from typing import Protocol

## @brief Indicates that a C ABI operation completed successfully.
APO_STATUS_OK: int = 0

## @brief Indicates that a C ABI argument, pointer, size, or value was rejected.
APO_STATUS_INVALID_ARGUMENT: int = 1

## @brief Indicates an uncategorized native runtime or exception failure.
APO_STATUS_RUNTIME_ERROR: int = 2

## @brief Indicates that a CUDA operation failed.
APO_STATUS_CUDA_ERROR: int = 3

## @brief Indicates that required native state has not been initialized.
APO_STATUS_NOT_INITIALIZED: int = 4

## @brief Indicates that a recognized native operation is not implemented.
APO_STATUS_NOT_IMPLEMENTED: int = 5


_STATUS_NAMES: dict[int, str] = {
    APO_STATUS_OK: "APO_STATUS_OK",
    APO_STATUS_INVALID_ARGUMENT: "APO_STATUS_INVALID_ARGUMENT",
    APO_STATUS_RUNTIME_ERROR: "APO_STATUS_RUNTIME_ERROR",
    APO_STATUS_CUDA_ERROR: "APO_STATUS_CUDA_ERROR",
    APO_STATUS_NOT_INITIALIZED: "APO_STATUS_NOT_INITIALIZED",
    APO_STATUS_NOT_IMPLEMENTED: "APO_STATUS_NOT_IMPLEMENTED",
}
_UNKNOWN_STATUS_NAME: str = "APO_STATUS_UNKNOWN"
_NATIVE_DIAGNOSTIC_FALLBACK: str = "Unknown apoCHARMM C API error"


class _StatusFunction(Protocol):
    """
    @brief Describes a mutable `ctypes` status-returning function object.

    `configure_status_function()` uses this structural protocol to assign the
    function's argument types, integer return type, and error callback without
    depending on a concrete private `ctypes` function class.
    """

    ## @brief Holds the declared `ctypes` argument conversions.
    argtypes: Sequence[object] | None

    ## @brief Holds the declared `ctypes` return conversion.
    restype: object

    ## @brief Holds the callback that validates each converted return value.
    errcheck: Callable[[int, "_StatusFunction", tuple[object, ...]], int] | None


class ApoCharmmError(RuntimeError):
    """
    @brief Represents one failed apoCHARMM native operation in Python.

    The exception retains the exact numeric C ABI status, a known exported
    status name or `APO_STATUS_UNKNOWN`, the Python operation context, the
    decoded native diagnostic, and the final rendered message. It is also a
    `RuntimeError`, so callers may catch either class.

    `check_status()` copies and decodes the native thread-local C string before
    constructing this object. The exception therefore owns Python references
    that remain valid after later native calls. The rendered message normalizes
    CRLF and lone carriage returns to line feeds, while `native_diagnostic`
    retains the original decoded text.

    @anchor python_apocharmm_error_class
    """

    def __init__(self, status: int, context: str, native_diagnostic: str) -> None:
        """
        @brief Constructs a Python exception from status and diagnostic fields.

        `status` is converted with `int()`. Unknown integer values are retained
        and use `APO_STATUS_UNKNOWN`. A single-line diagnostic is separated from
        the context header by `": "`. A multiline diagnostic begins on the next
        line unless it already starts with a newline, which preserves the native
        apoCHARMM banner layout.

        @param[in] status C ABI status or another object accepted by `int()`.
        @param[in] context Python-facing operation description retained in the
        rendered header.
        @param[in] native_diagnostic Decoded native diagnostic retained without
        newline normalization.
        @throws TypeError If `status` cannot be converted with `int()`.
        @throws ValueError If `int()` rejects the supplied status value.
        @throws OverflowError If `int()` reports an overflow while converting
        `status`.
        """
        self._status: int = int(status)
        self._status_name: str = _STATUS_NAMES.get(self._status, _UNKNOWN_STATUS_NAME)
        self._context: str = context
        self._native_diagnostic: str = native_diagnostic

        rendered_native_diagnostic: str = native_diagnostic.replace(
            "\r\n", "\n"
        ).replace("\r", "\n")
        header: str = "{} [{} ({})]".format(
            self._context, self._status_name, self._status
        )

        separator: str = ": "
        if "\n" in rendered_native_diagnostic:
            separator = "" if rendered_native_diagnostic.startswith("\n") else "\n"

        self._message: str = "{}{}{}".format(
            header, separator, rendered_native_diagnostic
        )

        super().__init__(self._message)
        return

    @property
    def status(self) -> int:
        """
        @brief Returns the exact numeric C ABI status.

        @return The integer supplied at construction, including an unknown
        nonzero value.
        """
        return self._status

    @property
    def status_name(self) -> str:
        """
        @brief Returns the symbolic name associated with the numeric status.

        @return The matching exported `APO_STATUS_*` name, or
        `"APO_STATUS_UNKNOWN"` when the numeric value is not recognized.
        """
        return self._status_name

    @property
    def context(self) -> str:
        """
        @brief Returns the Python-facing operation context.

        @return The context supplied at construction.
        """
        return self._context

    @property
    def native_diagnostic(self) -> str:
        """
        @brief Returns the copied and decoded native diagnostic.

        @return The diagnostic supplied at construction, before
        rendered-message newline normalization.
        """
        return self._native_diagnostic

    @property
    def message(self) -> str:
        """
        @brief Returns the complete rendered Python exception message.

        @return The context, status name, numeric status, and
        newline-normalized native diagnostic in the same text returned by
        `str(error)`.
        """
        return self._message


def check_status(status: int, context: str) -> None:
    """
    @brief Raises `ApoCharmmError` for a failing C ABI status.

    `APO_STATUS_OK` returns immediately without reading the native diagnostic.
    Every other integer causes the current thread's `apo_last_error()` pointer
    to be copied immediately. A NULL or empty C string uses the defined
    fallback; other bytes are decoded as UTF-8 with replacement for invalid
    sequences.

    @param[in] status Numeric C ABI status to inspect. Unknown nonzero values
    are preserved in the raised exception.
    @param[in] context Python-facing operation description retained by the
    exception.
    @return `None` when `status` is `APO_STATUS_OK`.
    @throws ApoCharmmError If `status` is any nonzero value.

    @post A raised exception owns Python text independent of the native
    diagnostic pointer's lifetime.
    @anchor python_apocharmm_check_status
    """
    if status == APO_STATUS_OK:
        return

    from ._lib import lib

    message_address: int | None = lib().apo_last_error()
    if message_address is None:
        native_diagnostic: str = _NATIVE_DIAGNOSTIC_FALLBACK
    else:
        message_bytes: bytes = ctypes.string_at(message_address)
        native_diagnostic = message_bytes.decode("utf-8", errors="replace")
        if native_diagnostic == "":
            native_diagnostic = _NATIVE_DIAGNOSTIC_FALLBACK

    raise ApoCharmmError(status, context, native_diagnostic)


def configure_status_function(
    function: _StatusFunction, argtypes: Sequence[object], context: str
) -> None:
    """
    @brief Configures a `ctypes` function to translate nonzero C ABI statuses.

    The function's `argtypes` are replaced with a new list copied from
    `argtypes`, its `restype` is set to `ctypes.c_int`, and its `errcheck`
    callback is replaced. The callback returns a successful integer unchanged
    and delegates every nonzero result to `check_status()` using the retained
    operation context.

    @param[in,out] function Borrowed mutable `ctypes` function object. The
    object is configured in place and is not otherwise retained.
    @param[in] argtypes Sequence of `ctypes` argument declarations. The
    sequence is copied to a list; the declaration objects themselves are
    retained by that list.
    @param[in] context Python-facing operation description captured by the
    installed callback for as long as the function retains that callback.
    @return `None`.

    @post `function` has integer status conversion and automatic
    `ApoCharmmError` translation installed.
    @anchor python_apocharmm_configure_status_function
    """

    def raise_on_error(
        result: int, _function: _StatusFunction, _arguments: tuple[object, ...]
    ) -> int:
        check_status(result, context)
        return result

    function.argtypes = list(argtypes)
    function.restype = ctypes.c_int
    function.errcheck = raise_on_error

    return

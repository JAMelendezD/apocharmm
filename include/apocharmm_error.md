# ApoCharmmError {#apocharmm_error}

## Purpose

The ApoCharmmError subsystem carries categorized native failures from C++
through the C ABI and into Python without allowing a C++ exception to cross the
C boundary. Use `ApoCharmmError` in native code when a failure needs an
apoCHARMM category, a human-readable message, and source metadata. C and Python
callers consume the corresponding `apo_status` and copied diagnostic instead of
a native exception object.

## Quick Start

Python callers normally catch the top-level `apocharmm.ApoCharmmError`. This
complete example deliberately supplies an empty parameter-file sequence so the
current API produces a verified invalid-argument error:

```python
import apocharmm as apo

try:
    apo.CharmmParameters([])
except apo.ApoCharmmError as error:
    print(error.status_name)
    print(error.context)
    print(error.native_diagnostic)
```

The three lines are:

```text
APO_STATUS_INVALID_ARGUMENT
CharmmParameters construction
apo_charmm_parameters_create_from_files: paths must contain at least one parameter file
```

Native C++ code should use the source-capturing macros rather than spelling
`__FILE__`, `__func__`, and `__LINE__` at each throw site:

```cpp
#include "ApoCharmmError.h"

#include <iostream>

int main() {
  try {
    APOCHARMM_REQUIRE(false, ApoCharmmErrorCode::InvalidArgument,
                      "demonstration failure");
  } catch (const ApoCharmmError &error) {
    std::cerr << GetApoCharmmErrorCodeName(error.getCode()) << ": "
              << error.getMessage() << '\n';
  }
}
```

This example writes `InvalidArgument: demonstration failure`.

## Construction and Required State

Constructing `ApoCharmmError` requires no CUDA context, simulation object, or
other initialized apoCHARMM state. The constructor copies the category, human
message, normalized source-file spelling, source function, and source line. The
`APOCHARMM_THROW` and `APOCHARMM_REQUIRE` macros must be expanded in function
scope because they capture `__func__`.

A status-returning C ABI entry point runs its implementation inside
`apocharmm_c::guard`. The guard clears the current thread's stale diagnostic
before invoking the operation. A returned nonzero status is completed with a
diagnostic when necessary, and every C++ exception is translated before control
returns to C. Call `apo_last_error()` immediately after observing a nonzero
status and before making another status-returning C ABI call on the same thread.

Python wrapper modules configure status-returning `ctypes` functions with
`configure_status_function()`. The installed `errcheck` callback calls
`check_status()` on the same thread immediately after the C function returns.
No separate Python initialization is required beyond loading the shared library
through the normal wrapper API.

## Ownership and Lifetime

A native `ApoCharmmError` owns its message and source strings. Views returned by
`getMessage()`, `getSourceFile()`, and `getSourceFunction()` borrow storage from
the exception and become invalid when that exception is destroyed or assigned
a new value. `GetApoCharmmErrorCodeName()` instead returns a view of static
storage that remains valid for the process lifetime.

The C ABI owns one diagnostic slot per thread. `apo_last_error()` returns a
borrowed, immutable pointer into that slot. Copy the bytes before another
status-returning call on the same thread clears or replaces the slot. Repeated
calls to `apo_last_error()` do not alter it. Successful destruction preserves
an existing diagnostic so a language-runtime finalizer cannot erase the error
that triggered cleanup; a destruction failure may replace the slot.

`check_status()` copies the C diagnostic bytes and decodes them into a Python
`str` before constructing the Python exception. Consequently,
`ApoCharmmError.status`, `status_name`, `context`, `native_diagnostic`, and
`message` remain valid independently of later C ABI calls and native object
lifetimes.

## Data, Shapes, and AKMA Units

This subsystem does not represent coordinates, forces, energies, time, or any
other physical quantity. It performs no host/device transfer and has no array
shape, stride, or CUDA-residence contract. Error categories, C status values,
and source-line numbers are dimensionless. Native and Python diagnostics are
text; the C ABI exposes them as a null-terminated byte string.

## Errors

The native categories and their C/Python mappings are:

| Native category | C status | Python status name |
| --- | --- | --- |
| `ApoCharmmErrorCode::InvalidArgument` | `APO_STATUS_INVALID_ARGUMENT` | `APO_STATUS_INVALID_ARGUMENT` |
| `ApoCharmmErrorCode::Runtime` | `APO_STATUS_RUNTIME_ERROR` | `APO_STATUS_RUNTIME_ERROR` |
| `ApoCharmmErrorCode::Cuda` | `APO_STATUS_CUDA_ERROR` | `APO_STATUS_CUDA_ERROR` |
| `ApoCharmmErrorCode::NotInitialized` | `APO_STATUS_NOT_INITIALIZED` | `APO_STATUS_NOT_INITIALIZED` |
| `ApoCharmmErrorCode::NotImplemented` | `APO_STATUS_NOT_IMPLEMENTED` | `APO_STATUS_NOT_IMPLEMENTED` |

At the C boundary, a legacy `std::invalid_argument` maps to
`APO_STATUS_INVALID_ARGUMENT`. Any other `std::exception`, an unknown C++
exception, or an unrecognized `ApoCharmmErrorCode` maps to
`APO_STATUS_RUNTIME_ERROR`. A callable that directly returns a nonzero
`apo_status` retains that exact status; the guard only ensures that a nonempty,
function-prefixed diagnostic is available.

A successful guarded C ABI call leaves the thread-local diagnostic empty. A
failing call returns one exact nonzero status and records a diagnostic
retrievable with `apo_last_error()`. The diagnostic does not determine whether
another function's output parameters were modified; consult that function's own
contract for output state after failure.

Python-side type or value validation performed before the C call continues to
raise the documented Python exception, such as `TypeError` or `ValueError`. A
nonzero native status raises the single Python `ApoCharmmError` class. Unknown
numeric statuses are preserved in `status` and use `APO_STATUS_UNKNOWN` as
`status_name`. A null or empty native diagnostic becomes
`Unknown apoCHARMM C API error`, and invalid UTF-8 bytes are decoded with the
Unicode replacement character.

## Important Behavior and Limitations

The native `what()` text begins with two newline characters and the apoCHARMM
ASCII banner, followed by the category, human message, normalized source path,
source line, and source function. The C ABI preserves that multiline text when
it catches an `ApoCharmmError`. Python normalizes CRLF and lone carriage returns
to line feeds only in the rendered `message`; `native_diagnostic` retains the
original decoded text.

The C diagnostic store has one slot per thread rather than a queue or stack.
Later failures overwrite earlier diagnostics, and successful status-returning
calls clear them. Diagnostics from different threads do not overwrite one
another, but a borrowed C pointer must not outlive its originating thread or the
next invalidating call on that thread.

Error formatting is an error-path operation and may allocate. The C ABI catches
formatting/storage failures and substitutes a fallback diagnostic. Its
function-prefixed emergency buffer is fixed at 1024 bytes, including the null
terminator, so an emergency diagnostic may be truncated. No C++ exception is
allowed to escape either the status-returning or destruction boundary.

Choosing `ApoCharmmErrorCode::Cuda` records a category; this subsystem does not
query CUDA state, synchronize a stream or device, or copy device data. The code
that detects the CUDA failure remains responsible for supplying the precise
human diagnostic.

## Related Subsystems

- @ref CharmmParameters supplies the minimal Python failure example above.
- @ref CharmmContext and @ref ForceManager report state and runtime failures
  through this subsystem.
- @ref CudaIntegrator and its derived integrators use the same native-to-Python
  error boundary for integration failures.

## Developer Architecture

The native layer consists of `ApoCharmmErrorCode`, `ApoCharmmError`,
`GetApoCharmmErrorCodeName()`, `APOCHARMM_THROW`, and `APOCHARMM_REQUIRE`.
`src/ApoCharmmError.cpp` owns the diagnostic banner and source-path
normalization. The native library receives its configured source root through
the `APOCHARMM_SOURCE_ROOT` compile definition, allowing macro-captured paths to
avoid exposing an absolute developer-machine prefix.

The C ABI layer has three roles. `apo_status` is the stable C-visible category.
The private `status_from_error_code()` switch performs the explicit native-to-C
mapping. `guard()` and `guard_destroy()` form the nonthrowing ABI boundary.
`LastErrorState` provides normal string storage, a no-throw emergency buffer,
and the borrowed pointer returned by `apo_last_error()`. Validation helpers set
a specific status and diagnostic, while `APOCHARMM_C_RETURN_IF_ERROR` propagates
that status from an entry-point body without throwing.

The Python layer configures `apo_last_error()` to return a raw address, not an
automatically converted `ctypes.c_char_p`. `check_status()` can therefore
distinguish a null pointer from an empty C string, copy the bytes while the C
lifetime is valid, decode with replacement, and construct one Python exception.
`configure_status_function()` installs this conversion as the `errcheck`
callback and retains the operation context in its closure.

The ownership graph is linear on failure: the native exception owns its C++
strings; the C boundary copies `what()` into thread-local C ABI storage; the
Python boundary copies and decodes that storage into Python strings. No layer
retains a pointer or view owned by the preceding layer after its documented
lifetime.

To add a new stable error category, update the native enum and name switch, the
public C `apo_status` enum, `status_from_error_code()`, the exported Python
constant and `_STATUS_NAMES`, and the top-level Python re-export. Extend the
native mapping tests, direct C ABI mapping tests, and Python field/rendering
tests in the same change. A new public C status value is an ABI change and its
numeric value must not reuse an existing value.

The focused regression coverage is in
`test/unittests/unittest-apoCharmmError.cpp`,
`test/unittests/unittest-capiError.cpp`, and
`test/pytest/python_api_error.py`. These tests fix the category names and
mappings, source capture, multiline rendering, C pointer lifetime, thread-local
isolation, exception translation, destroy behavior, unknown Python statuses,
and decoding fallbacks.

The principal architectural constraints are the single diagnostic slot per
thread and the void C destruction API. A destroy failure cannot return a
status; it can only replace the thread-local diagnostic. Error-path storage is
designed to remain nonthrowing at the C boundary rather than to preserve an
arbitrarily large diagnostic under allocation failure.

## API Reference

- Native C++: @ref ApoCharmmError, @ref ApoCharmmErrorCode,
  @ref GetApoCharmmErrorCodeName(), @ref APOCHARMM_THROW, and
  @ref APOCHARMM_REQUIRE.
- C ABI: @ref apo_status and @ref apo_last_error().
- Python: @ref python_apocharmm_error_module "error module",
  @ref python_apocharmm_error_class "ApoCharmmError",
  @ref python_apocharmm_check_status "check_status()", and
  @ref python_apocharmm_configure_status_function
  "configure_status_function()".

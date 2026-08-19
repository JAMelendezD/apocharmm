# Coordinates {#coordinates}

## Purpose

The Coordinates subsystem represents Cartesian atom positions before and at the
boundary of an apoCHARMM simulation context. The native
[Coordinates](@ref Coordinates) class owns parallel host vectors of `double3`
and `float3` records. [CharmmCrd](@ref CharmmCrd) constructs that representation
from a CHARMM coordinate file.

Use [CharmmCrd](@ref CharmmCrd) when coordinates are supplied in standard or
extended CHARMM CRD/COR format. Use [Coordinates](@ref Coordinates) directly
when native C++ code already has one `[x, y, z]` record per atom. The public C
ABI and Python API expose CharmmCrd, but do not expose direct construction of
the Coordinates base class.

These objects are input-side host representations. Runtime coordinates and
charges are owned separately by [CharmmContext](@ref charmm_context) in
[CudaContainer](@ref cuda_container) instances after `setCoordinates()` copies
those values.

## Quick Start

The highest-level supported workflow is the Python wrapper:

```python
from pathlib import Path

import apocharmm as apo

with apo.CharmmCrd(Path("test/data/nacl_pair.cor")) as coordinates:
    assert coordinates.getNumAtoms() == 2
    assert coordinates.getCoordinates() == [
        [0.0, 0.0, 0.0],
        [2.82, 2.82, 2.82],
    ]
```

Run the example from the repository root. Construction finishes parsing before
return, the nested result is an independent Python copy, and context-manager
exit releases the native handle.

Direct native use follows the same parse-before-return model:

```cpp
#include "CharmmCrd.h"

#include <vector>

int main() {
  CharmmCrd coordinates("test/data/nacl_pair.cor");
  const std::vector<double3> &xyz = coordinates.getCoordinatesDP();

  return coordinates.getNumAtoms() == 2 && xyz.size() == 2 ? 0 : 1;
}
```

Link the example against the native apoCHARMM library. Reading coordinates does
not itself allocate or synchronize CUDA storage.

## Construction and Required State

Default construction of [Coordinates](@ref Coordinates) creates unset state:
`getNumAtoms() == -1` and both vectors are empty. The four value constructors
accept either `std::vector<double3>`, `std::vector<float3>`, nested double rows,
or nested float rows. An empty input constructs an explicit zero-atom object.
Nested rows must contain exactly three values.

`setNumAtoms()` accepts a non-negative dimensionless count and resizes both
vectors. Existing values in the retained prefix survive a successful resize;
new records are value-initialized. The normal native invariant is:

```text
getNumAtoms() == getCoordinatesDP().size()
              == getCoordinatesSP().size()
```

The public mutable vector accessors can violate that invariant and do not
perform precision conversion. Code that mutates one representation must update
the other representation and atom count explicitly when later consumers require
coherence.

[CharmmCrd](@ref CharmmCrd) requires a nonempty path. It reads the complete file
into host memory, skips consecutive nonempty title records beginning with `*` in
column one, and parses the atom count from the next line's first space-delimited
token. An exact uppercase `EXT` second token selects extended fixed-width
records; all other count lines select standard records. Counts from zero through
`INT_MAX` are accepted. The parser retains only X, Y, and Z in file-record order
and does not validate the file-name extension.

A [CharmmContext](@ref CharmmContext) must already have a positive atom count
before coordinates are assigned. The coordinate count must exactly match the
context atom count. For a context created from a PSF and parameters, the PSF
establishes that count. Coordinates may be assigned before or after valid box
dimensions; final context initialization occurs when all required context and
force-manager state is available.

## Ownership and Lifetime

A native Coordinates object exclusively owns its two host vectors and stored
count. Constructor inputs are copied and never retained. Compiler-generated
copy operations create independent vectors. Compiler-generated move operations
transfer vector storage and leave the source in a valid but otherwise
unspecified moved-from state.

Const and mutable coordinate accessors return borrowed references to the vector
members. The vector objects remain alive with their owning Coordinates object.
Element references, pointers, and iterators can be invalidated by resizing,
vector assignment, or other structural mutation. Object destruction ends every
borrowed reference. Mutating returned storage is never synchronized with the
other precision representation.

Coordinates is used as a public base class but does not have a virtual
destructor. Derived objects must be destroyed through their actual type rather
than through a Coordinates pointer.

A public [apo_charmm_crd](@ref apo_charmm_crd) handle owns a private
`std::shared_ptr<CharmmCrd>`. Release it with
[apo_charmm_crd_destroy](@ref apo_charmm_crd_destroy). C ABI getter inputs are
borrowed. Passing a handle to
[apo_charmm_context_set_coordinates](@ref apo_charmm_context_set_coordinates)
copies its coordinates into the context and does not retain the coordinate
handle.

The Python [python_charmm_crd](@ref python_charmm_crd) wrapper owns one C
handle. `close()`, `destroy()`, and context-manager exit release it
idempotently. The wrapper does not retain the path object. `getCoordinates()`
returns new Python lists and no returned value aliases native storage.
`CharmmContext.setCoordinates()` borrows a CharmmCrd wrapper for the native call
and does not retain it.

No layer provides internal synchronization. Externally serialize concurrent
mutation, destruction, and any read that overlaps mutation.

## Data, Shapes, and AKMA Units

Cartesian coordinates use angstroms throughout this subsystem. Atom counts,
array lengths, and indices are dimensionless.

| Interface or storage | Shape | Element order | Residence and ownership |
| --- | --- | --- | --- |
| `Coordinates::m_CoordinatesDP` | `N` `double3` records | atom-major `x, y, z` | owned host vector |
| `Coordinates::m_CoordinatesSP` | `N` `float3` records | atom-major `x, y, z` | owned host vector |
| `apo_charmm_crd_get_coordinates()` | flat `3 * N` doubles | `x_i, y_i, z_i` in increasing atom order | caller-owned host buffer |
| `CharmmCrd.getCoordinates()` | nested `(N, 3)` Python list | one `[x, y, z]` row per atom | newly owned Python copy |

There is no fourth coordinate component and this subsystem does not store
charge. Context coordinate-and-charge arrays use separate `float4` and
`double4` representations documented by [CharmmContext](@ref charmm_context).

Double-input native constructors and CharmmCrd preserve parsed values in the
double-precision vector and convert each component for the single-precision
vector. Float-input constructors preserve the single-precision values and widen
them for the double-precision vector.

Coordinates and its file readers perform no device allocation or transfer.
`CharmmContext::setCoordinates()` copies the double-precision records into the
context's coordinate/charge mirrors, preserves existing charges, transfers both
context mirrors to device, and can image coordinates and rebuild the neighbor
list when the force manager is already initialized.

## Errors

### Native C++

Nested Coordinates constructors throw [ApoCharmmError](@ref ApoCharmmError) with
`ApoCharmmErrorCode::InvalidArgument` when a row does not contain exactly three
values. `setNumAtoms()` uses the same code for a negative count and rejects that
value before mutating the object.

CharmmCrd construction uses `ApoCharmmErrorCode::InvalidArgument` for an empty
path. File-open, sizing, seeking, reading, end-of-file, count, fixed-width
record, and coordinate-value failures use `ApoCharmmErrorCode::Runtime`.
Malformed, out-of-range, and non-finite coordinate fields are rejected.

Native vector, string, file-buffer, and diagnostic allocation can propagate
`std::bad_alloc`. Requests exceeding implementation container limits can
propagate `std::length_error`. `setNumAtoms()` is not transactional: an
allocation failure after an earlier mutation can leave the stored count and two
vector lengths inconsistent.

PDB does not use the unified error utility for its parser. A fixed-width
field whose starting column is beyond the line or an out-of-range conversion
can throw `std::out_of_range`; unconvertible coordinate text can throw
`std::invalid_argument`. Failure to open the file prints a message and calls
`exit(0)` instead of throwing.

### C ABI

[apo_charmm_crd_create](@ref apo_charmm_crd_create) returns
`APO_STATUS_INVALID_ARGUMENT` for a NULL output slot or a NULL or empty path.
Verified file, parse, count-range, host allocation, and translated
standard-exception failures return `APO_STATUS_RUNTIME_ERROR`. Coordinate
construction performs no CUDA work and has no verified `APO_STATUS_CUDA_ERROR`
path.

[apo_charmm_crd_get_num_atoms](@ref apo_charmm_crd_get_num_atoms) sets a
validated output slot to zero before handle validation. A later invalid-handle
failure leaves zero visible. Its explicit negative-native-count guard returns
`APO_STATUS_RUNTIME_ERROR`.

[apo_charmm_crd_get_coordinates](@ref apo_charmm_crd_get_coordinates) validates
the handle and complete required capacity before writing. A rejected handle,
NULL nonempty buffer, or short buffer returns `APO_STATUS_INVALID_ARGUMENT` and
leaves all caller bytes unchanged. For a zero-length coordinate vector, a NULL
output buffer is accepted. On success, trailing buffer elements beyond `3 * N`
remain unchanged.

Each status-returning call clears the calling thread's previous diagnostic.
Failure leaves borrowed thread-local text available through
[apo_last_error](@ref apo_last_error) until the next guarded call on that thread
or thread exit. Normal destruction accepts NULL, does not throw across the C
boundary, and preserves the previous diagnostic.

### Python

The Python constructor accepts `str`, `bytes`, `os.PathLike[str]`, and
`os.PathLike[bytes]`. `os.fsencode` conversion can raise `TypeError`. Missing
library configuration raises `RuntimeError`, loading failure raises `OSError`,
and native nonzero statuses raise
[python_apocharmm_error_class](@ref python_apocharmm_error_class).

Using a closed wrapper raises `RuntimeError`. `getCoordinates()` can raise
`OverflowError` if the native atom count cannot be represented as a ctypes
array length, or `MemoryError` while allocating its temporary C buffer or
returned lists. Native file and parse failures are raised as `ApoCharmmError`
with the numeric status, operation context, and a copied diagnostic.

## Important Behavior and Limitations

The two precision vectors are mirrors only by convention. Mutable access does
not set a dirty flag, resize the other vector, update the atom count, or perform
any conversion. Consumers that assume coherent lengths or values can observe
inconsistent state after unchecked mutation.

`setNumAtoms()` assigns the new count before resizing the double vector and then
the float vector. Its successful path preserves the normal invariant, but its
allocation-failure path provides only partial exception safety. Element
addresses can change after either resize.

The native value constructors convert `std::vector::size()` to `int` without
runtime validation. Inputs larger than `INT_MAX` are outside the supported
contract. Double-to-float construction also lacks a runtime representability
check. Keep every double component within the representable float range when a
single-precision mirror must be produced.

CharmmCrd reads the entire file before parsing. It recognizes only an exact
uppercase `EXT` marker, uses fixed field positions, ignores all metadata and
trailing fields, and does not retain the source path. Finite double fields that
are outside the representable float range are not rejected before the
single-precision conversion.

The subsystem performs no CUDA operation despite the `.cu` implementation-file
suffix. CUDA failures become relevant only after coordinates are copied into a
CharmmContext or another GPU-owning subsystem.

The class hierarchy is not safe for polymorphic deletion because Coordinates
has a non-virtual destructor. The subsystem also has no virtual parsing or
storage extension interface; derived readers write the protected representation
directly.

## Related Subsystems

- [CharmmContext](@ref charmm_context) copies coordinate input into simulation
  state and controls host/device transfers, imaging, and neighbor-list
  rebuilding.
- [CharmmPSF](@ref charmm_psf) supplies the topology whose atom order and count
  must match the coordinate records.
- [CudaContainer](@ref cuda_container) owns explicit host/device mirrors used
  after coordinates enter a context.
- [ApoCharmmError](@ref apocharmm_error) documents native categories, C ABI
  translation, and Python exceptions.

## Developer Architecture

The native base declaration is `include/Coordinates.h`; host-only definitions
are in `src/Coordinates.cu`. The `.cu` suffix does not imply device storage or
kernel execution. `include/CharmmCrd.h` and `src/CharmmCrd.cu` implement the
supported CHARMM file reader. Both documented sources are registered in
`src/CMakeLists.txt`.

The public C ABI is limited to CharmmCrd. Its opaque declaration is in
`capi/include/apocharmm_c/CharmmCrd.h`; the private shared-ownership payload is
in `capi/private_include/apocharmm_c/detail/CharmmCrdHandle.h`; and guarded
entry points are in `capi/src/apocharmm_c/CharmmCrd.cpp`. The source is
registered in `capi/CMakeLists.txt`.

The Python wrapper is `python/apocharmm/charmm_crd.py` and is re-exported as
`apocharmm.CharmmCrd`. It configures ctypes prototypes lazily, converts paths
with `os.fsencode`, owns one opaque handle, and copies native doubles through a
flat temporary C buffer into a nested Python list.

The normal native-to-runtime data flow is:

```text
host vector or coordinate file
        |
        v
Coordinates
  |-- std::vector<double3>  host DP representation
  `-- std::vector<float3>   host SP representation
        |
        | CharmmContext::setCoordinates() reads DP and copies values
        v
CharmmContext host float4/double4 coordinate-charge mirrors
        |
        | synchronous CudaContainer host-to-device transfers
        v
CharmmContext device coordinate-charge mirrors
```

The C ABI boundary validates pointers and lengths, catches native exceptions,
maps [ApoCharmmError](@ref ApoCharmmError) categories to exact
[apo_status](@ref apo_status) values, and stores thread-local diagnostics. The
Python ctypes `errcheck` boundary copies that diagnostic and raises one owned
Python exception.

Normal representation invariants are a non-negative initialized count, equal
DP and SP vector lengths, atom-major record order, and corresponding components
across precisions subject to float rounding. The protected members and mutable
accessors permit derived readers and low-level callers to bypass every
invariant; no repair or validation pass is available.

Extension work should preserve these boundaries:

- public native contracts belong in `include/Coordinates.h` and reader headers;
- file-format algorithms and non-obvious parsing constraints belong in the
  corresponding source;
- C ABI functions must validate complete output capacity before writing;
- Python outputs should remain owned copies unless a new lifetime-safe view type
  is designed;
- context assignment should continue to copy input rather than retain mutable
  coordinate objects.

Performance-sensitive work is linear in the atom count. Native construction
allocates two vectors and converts every record. CharmmCrd additionally reads
the complete file into one string. C ABI extraction copies `3 * N` doubles, and
Python extraction performs a second conversion into `N` inner lists. Passing
coordinates to a context performs further host copies and two device transfers.

Focused native coverage is in `test/unittests/unittest-coordinates.cpp` and
`test/unittests/unittest-charmmCrd.cpp`. Direct C ABI coverage is in
`test/unittests/unittest-capiCharmmCrd.cpp`. Python coverage is in
`test/pytest/python_api_charmm_crd.py`. These files are registered by the
corresponding test CMake lists. Stable minimal data is
`test/data/nacl_pair.cor`.

Known architectural constraints visible in the documented implementation include
the non-virtual base destructor, unchecked mutable representation access,
`int`-sized atom counts, non-transactional resizing, unchecked double-to-float
range conversion, and whole-file CharmmCrd loading.

## API Reference

- [Coordinates](@ref Coordinates) is the native host representation and mutable
  input base.
- [CharmmCrd](@ref CharmmCrd) is the supported native CHARMM coordinate-file
  reader.
- [apo_charmm_crd](@ref apo_charmm_crd) is the public opaque C ABI handle.
- [apo_charmm_crd_create](@ref apo_charmm_crd_create) and
  [apo_charmm_crd_destroy](@ref apo_charmm_crd_destroy) manage C handle
  lifetime.
- [apo_charmm_crd_get_num_atoms](@ref apo_charmm_crd_get_num_atoms) and
  [apo_charmm_crd_get_coordinates](@ref apo_charmm_crd_get_coordinates) copy C
  ABI outputs.
- [python_charmm_crd](@ref python_charmm_crd) is the owning Python wrapper
  reference.
- [CharmmContext](@ref charmm_context) documents the native consumer of
  coordinate input.

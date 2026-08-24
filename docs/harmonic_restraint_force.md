# HarmonicRestraintForce {#harmonic_restraint_force}

## Purpose

`HarmonicRestraintForce` applies absolute Cartesian harmonic restraints to
individual atoms. It is appropriate when selected atoms should remain near
fixed reference positions rather than near another atom, a group center, or a
periodically imaged location.

For atom `i`, let `(x_i, y_i, z_i)` be the current Cartesian position and
`(x0_i, y0_i, z0_i)` be the absolute reference position. Define the
displacement components as:

```text
dx_i = x_i - x0_i
dy_i = y_i - y0_i
dz_i = z_i - z0_i
```

The implemented energy is:

```text
U_harm = sum over all atoms i of
         k_i * m_i * (dx_i^2 + dy_i^2 + dz_i^2)
```

Here, `k_i` is the per-atom force constant and `m_i` is the configured mass or
weight.

The values accumulated into apoCHARMM's force array are the
potential-gradient components:

```text
gradient_x_i = 2 * k_i * m_i * dx_i
gradient_y_i = 2 * k_i * m_i * dy_i
gradient_z_i = 2 * k_i * m_i * dz_i
```

apoCHARMM integration code subtracts this stored potential gradient. The
physical restoring-force components therefore have the opposite signs.

Selection is implemented through the per-atom force constants. Unselected
atoms have `k_i=0`.

## Quick Start

The Python API is the highest-level supported interface. This example uses
repository test data and assumes it is run from the repository root.

```python
import apocharmm as apo

parameters = apo.CharmmParameters("test/data/toppar_water_ions.str")
psf = apo.CharmmPsf("test/data/nacl_pair.psf")
coordinates = apo.CharmmCrd("test/data/nacl_pair.cor")

context = apo.CharmmContext(psf, parameters)
context.setBoxDimensions((50.0, 50.0, 50.0))
context.setCoordinates(coordinates)

selector = apo.AtomSelector(psf)
selection = selector.select("all")

restraint = apo.HarmonicRestraintForce(psf.getNumAtoms())
restraint.setSelection(selection)
restraint.setReferenceCoordinates(coordinates.getCoordinates())
restraint.setMasses(psf.getMasses())
restraint.setForceConstant(10.0)

force_manager = context.getForceManager()
force_manager.subscribe(restraint)

context.calculatePotentialEnergy()

force_manager.unsubscribe(restraint)
restraint.close()
```

Because the current and reference coordinates are identical in this example,
the harmonic contribution is zero. The total potential energy can still
contain the ordinary force-field terms.

Direct C++ use supplies the same resources explicitly during subscription:

```cpp
#include "AtomSelection.h"
#include "CharmmContext.h"
#include "CharmmCrd.h"
#include "CharmmPSF.h"
#include "CharmmParameters.h"
#include "HarmonicRestraintForce.h"

#include <memory>
#include <string>
#include <vector>

const std::string dataPath = "test/data/";

auto parameters =
    std::make_shared<CharmmParameters>(dataPath + "toppar_water_ions.str");
auto psf = std::make_shared<CharmmPSF>(dataPath + "nacl_pair.psf");
auto coordinates =
    std::make_shared<CharmmCrd>(dataPath + "nacl_pair.cor");

auto context = std::make_shared<CharmmContext>(psf, parameters);
context->setBoxDimensions({50.0, 50.0, 50.0});
context->setCoordinates(coordinates);

AtomSelection selection(psf->getNumAtoms(),
                        AtomSelection::InitialValue::ALL);

auto restraint =
    std::make_shared<HarmonicRestraintForce<long long int, float>>(
        psf->getNumAtoms());

restraint->setSelection(selection);
restraint->setReferenceCoordinates(coordinates->getCoordinatesDP());
restraint->setMasses(psf->getMasses());
restraint->setForceConstant(10.0);

auto forceManager = context->getForceManager();
forceManager->subscribe(restraint, "harm", restraint->getStream(),
                        restraint->getForce(),
                        restraint->getEnergyVirial());

context->calculatePotentialEnergy();

forceManager->unsubscribe(restraint);
```

## Construction and Required State

Construction requires a strictly positive atom count. The count is fixed for the
object's lifetime. The constructor establishes these defaults:

- every atom is selected;
- every force constant is zero;
- every reference position is `(0, 0, 0)`;
- every per-atom weight is one;
- stored box dimensions are `(0, 0, 0)`;
- the output energy term is named `harm`;
- a private CUDA stream and fixed-point force array are allocated.

Output accumulators are not a construction-time result. Direct native callers
must call `clear()` before the first `calcForce()` and before every calculation
that requires fresh output.

`initialize(numAtoms, boxDimensions)` is the force-manager initialization hook.
It checks that `numAtoms` matches the construction count and stores the box.
There is no separate initialized flag, no per-atom reallocation, and no
reference-coordinate rescaling. A manager that is already initialized invokes
this hook immediately when the restraint is subscribed. Otherwise, manager
initialization invokes it later.

Normal manager-driven order is:

1. construct the restraint with the manager atom count;
2. set selection, references, masses, and force constants;
3. subscribe the restraint;
4. let `CharmmContext` or an integrator drive force calculation;
5. unsubscribe before closing or destroying the public handle.

Native direct calculation additionally requires a valid device coordinate
pointer, `clear()`, `calcForce()`, and synchronization of `getStream()` before
host code consumes the result.

## Ownership and Lifetime

The native object owns or retains:

- host/device mirrors of force constants, references, weights, and box lengths;
- an `AtomSelection` copy;
- shared ownership of one `Force<AT>`;
- shared ownership of one `CudaEnergyVirial`;
- shared ownership of a host scalar containing the CUDA stream handle.

`getForce()`, `getEnergyVirial()`, and `getStream()` return copied shared owners.
A returned force or energy object can outlive the restraint. A returned stream
holder can also outlive it, but the restraint destructor destroys the
underlying CUDA stream and writes `nullptr` into the shared holder.

The public C handle owns a shared pointer to the native
`HarmonicRestraintForce<long long int, float>`. Successful manager subscription
copies that native owner. Destroying the C handle therefore does not remove an
existing subscription.

The Python wrapper owns its C handle. `ForceManager.subscribe()` retains the
Python wrapper after native subscription succeeds, and the native manager
retains the underlying restraint independently. Unsubscribe before calling
`close()`. Closing a subscribed wrapper destroys the identifying C handle but
does not remove the native manager subscription.

Selections and all scalar arrays passed to setters are copied. No C++, C, or
Python caller buffer is retained.

Do not copy or copy-assign the native C++ class. Its current implicit copy
operations do not establish safe independent ownership of the stream and output
resources.

## Data, Shapes, and AKMA Units

For a restraint containing `N` atoms:

| Quantity | Native C++ | C ABI | Python | Meaning and units |
|---|---|---|---|---|
| Selection | `AtomSelection` with `N` atoms | `apo_atom_selection` with `N` atoms | `AtomSelection` | Zero-based, dimensionless atom membership |
| Force constants | `std::vector<double>` of length `N` | `double[N]` | Iterable of `N` float-convertible values | Per-atom `k_i` |
| References | `std::vector<double3>` or `N` rows of length 3 | Flat `double[3N]` | `N` iterables of length 3 | `[x, y, z]` in angstroms |
| Masses or weights | `std::vector<double>` of length `N` | `double[N]` | Iterable of `N` float-convertible values | Normally atomic mass units |
| Box | `std::vector<double>` of length 3 | `double[3]` | Iterable of 3 float-convertible values | `[x, y, z]` lengths in angstroms |
| Current coordinates | Device `float4[N]` | Not exposed directly | Not exposed directly | `x`, `y`, `z` in angstroms; `w` ignored |
| Force output | `Force<long long int>` | Managed through `ForceManager` | Managed through `ForceManager` | Structure-of-arrays potential gradient |
| Energy output | `CudaEnergyVirial` term `harm` | Managed through `ForceManager` | Managed through `ForceManager` | Kilocalories per mole |

The implemented energy contains the product `k_i * m_i`. When `m_i` is in
atomic mass units, use `k_i` in kilocalories per mole per square angstrom
per atomic mass unit. With every weight set to one, the numerical behavior is
the ordinary per-atom harmonic expression with `k_i` in kilocalories per
mole per square angstrom.

Masses and force constants must be finite and non-negative. Zero is accepted
for either quantity and disables that atom's contribution. Reference
coordinates must be finite. Box lengths must be finite and strictly positive.

For the supported fixed-point native instantiations, force-array integers use a
scale of `2^40`. Multiplying by `1 / 2^40` recovers the stored
potential-gradient value in kilocalories per mole per angstrom. The active
force vectors use structure-of-arrays layout with implementation-selected
stride:

```text
x: force[0 * stride + atom]
y: force[1 * stride + atom]
z: force[2 * stride + atom]
```

The restraint compares current and reference coordinates directly. It does not
apply minimum-image wrapping.

## Errors

### Native C++

Verified argument failures raise
[ApoCharmmError](@ref ApoCharmmError) with
`ApoCharmmErrorCode::InvalidArgument`:

- construction with a non-positive atom count;
- a selection or initialization atom-count mismatch;
- a force-constant vector with the wrong length;
- a negative or non-finite selected force constant;
- a reference-coordinate collection with the wrong shape;
- a non-finite reference-coordinate component;
- a mass vector with the wrong length;
- a negative or non-finite mass;
- a box vector whose length is not three;
- a non-positive or non-finite box length.

CUDA allocation, transfer, synchronization, stream creation, output clearing,
and immediate kernel-launch failures raise `ApoCharmmError` with
`ApoCharmmErrorCode::Cuda`.

Host allocation and implementation-limit failures can raise standard
`std::bad_alloc` or `std::length_error`. Some setters are not transactional;
consult their symbol contracts before assuming failure preserves all host and
device state.

`calcForce()` performs no pointer, size, stream, or output-state validation.
Invalid direct-native inputs are outside its checked contract.

### C ABI

Every status-returning C function clears the calling thread's previous
diagnostic at entry. Success returns `APO_STATUS_OK` and leaves the diagnostic
empty. Failure leaves nonempty text available through
[apo_last_error](@ref apo_last_error) on the same thread.

`APO_STATUS_INVALID_ARGUMENT` covers invalid handles, null required pointers,
shape or length mismatches, invalid numeric values, duplicate subscription, and
unsubscription of a restraint that is not subscribed.

`APO_STATUS_CUDA_ERROR` covers categorized native CUDA failures.

`APO_STATUS_RUNTIME_ERROR` covers host allocation, diagnostic construction,
and uncategorized standard or nonstandard C++ exceptions.

Creation sets a valid output slot to `NULL` before native construction and
leaves it `NULL` on failure. Setter functions have no output parameters.
Destruction accepts `NULL`, does not throw across the ABI boundary, and normally
preserves a stale thread-local diagnostic.

### Python

Python conversion and shape checks occur before the C call:

- incompatible objects or non-iterable inputs can raise `TypeError`;
- failed scalar conversion or a reference row whose length is not three can
  raise `ValueError`;
- conversion of an excessively large numeric value can raise `OverflowError`;
- temporary list or `ctypes` allocation can raise `MemoryError`;
- using a closed wrapper raises `RuntimeError`;
- UTF-8 tag encoding can raise `UnicodeEncodeError`.

A nonzero native status raises `ApoCharmmError` with the numeric C status,
symbolic status name, Python operation context, copied native diagnostic, and
rendered message.

## Important Behavior and Limitations

`calcForce()` and `clear()` enqueue work on the restraint's private CUDA stream.
Neither method synchronizes the stream before returning. The manager orders and
synchronizes subscribed streams during aggregate force and energy processing.
Direct native callers must synchronize explicitly.

Parameter setters use `CudaContainer` transfers that call
`cudaDeviceSynchronize()`. They therefore wait for all previously requested
work on the current CUDA device, not only restraint work. The scalar and
per-atom force-constant setters currently perform a synchronized zero fill
followed by a second synchronized transfer.

Force and energy updates are additive. When `calcEnergy` is false, the existing
`harm` energy remains unchanged and can be stale. `calcVirial` is accepted for
interface compatibility but has no effect. The restraint declares
`contributesVirial == false`.

Changing the selection zeros constants for newly unselected atoms. Expanding
the selection does not restore old constants; set the scalar or per-atom
constants again.

Reference-coordinate and mass validation is interleaved with host mutation. A
later invalid value can leave an earlier host prefix updated while the device
mirror remains old.

Box dimensions are stored and transferred but are not read by the current
kernel. Box changes do not rescale references, and displacements are not
minimum-image wrapped.

The linked library provides explicit native instantiations only for:

```cpp
HarmonicRestraintForce<long long int, float>
HarmonicRestraintForce<long long int, double>
```

The public C ABI and Python API use the first specialization.

The native class is not safely copyable with its current compiler-generated
copy operations. Treat each instance as a unique configuration object managed
through `std::shared_ptr`.

No layer provides internal host-thread synchronization.

## Related Subsystems

- [AtomSelection](@ref atom_selection) defines the selected atom set.
- [ForceManager](@ref force_manager) owns subscribed calculation resources and
  aggregates their outputs.
- [CharmmContext](@ref charmm_context) owns molecular state and drives manager
  initialization and force calculation.
- [CudaContainer](@ref cuda_container) defines the parameter mirrors and their
  synchronous transfer behavior.
- `Force` defines fixed-point structure-of-arrays force storage.
- `CudaEnergyVirial` stores the named `harm` energy.
- [HarmonicCenterOfMassRestraintForce](@ref harmonic_center_of_mass_restraint_force)
  restrains a selected group center rather than each atom independently.
- [ApoCharmmError](@ref apocharmm_error) describes native, C ABI, and Python
  error propagation.

## Developer Architecture

### Public and private layers

The public native template is declared in
`include/HarmonicRestraintForce.h` and implemented in
`src/HarmonicRestraintForce.cu`.

The stable C ABI is declared in
`capi/include/apocharmm_c/HarmonicRestraintForce.h`. Its private shared-owner
handle is defined in
`capi/private_include/apocharmm_c/detail/HarmonicRestraintForceHandle.h`, and
guarded entry points are implemented in
`capi/src/apocharmm_c/HarmonicRestraintForce.cpp`.

The Python wrapper is implemented in
`python/apocharmm/harmonic_restraint_force.py`. It configures `ctypes`
prototypes, owns one C handle, converts Python inputs, and delegates nonzero
status handling to the shared Python error layer.

### Ownership graph

After successful Python subscription, the relevant ownership graph is:

```text
Python HarmonicRestraintForce --owns--> C restraint handle
C restraint handle --shared--> native restraint
Python ForceManager --retains--> Python HarmonicRestraintForce
native ForceManager --shared--> native restraint
native ForceManager --shared--> restraint stream holder
native ForceManager --shared--> restraint Force
native ForceManager --shared--> restraint CudaEnergyVirial
native restraint --shared--> stream holder, Force, CudaEnergyVirial
```

Unsubscription removes the manager-held native owners and then the Python
manager-held wrapper reference. Destruction of the restraint C handle removes
only the handle's native owner.

### State and invariants

The immutable atom count determines every per-atom array length. Normal
representation requires:

```text
forceConstants.size() == atomCount
referenceCoordinates.size() == atomCount
selection.getNumAtoms() == atomCount
boxDimensions.size() == 3
force.size() == atomCount
referenceCoordinates[i].w == configured mass or weight
unselected forceConstants[i] == 0
```

`m_ReferenceCoordinates` stores `[x0, y0, z0, mass]` in one `double4` per atom.
Selection is not read by the kernel; it affects calculation only through
zeroed force constants.

### Host and device data flow

Configuration begins in host memory. `CudaContainer` copies or updates the host
mirror and then transfers the complete active range to device memory.
Reference-position and mass setters share one `double4` array, so each transfer
moves both kinds of data.

The kernel reads:

```text
float4 current coordinate records
double per-atom force constants
double4 reference-position and mass records
```

It writes fixed-point X, Y, and Z force arrays atomically. When energy is
requested, each block reduces a double-precision partial sum and one thread
atomically adds it to the `harm` device scalar.

The current implementation launches 256 threads per block and selects one of
four template variants for the two Boolean calculation flags. The virial
template flag currently produces no virial write.

### Error boundaries

Native setters perform categorized validation before or during mutation.
`calcForce()` intentionally omits validation in the hot path and checks only
the immediate CUDA launch result.

C ABI entry points validate opaque handles and pointers before dereferencing
them. The common guard translates `ApoCharmmErrorCode` to `apo_status`, maps
legacy `std::invalid_argument` to `APO_STATUS_INVALID_ARGUMENT`, maps other
exceptions to `APO_STATUS_RUNTIME_ERROR`, and prevents exceptions from crossing
the ABI boundary.

Python wrappers perform Python type, conversion, and row-shape work before the
C call. Configured status callbacks convert every nonzero C result into
`ApoCharmmError`.

### Extension points and constraints

`ForceManager::subscribe()` accepts the restraint through its generic
force-view machinery together with the stream, force storage, energy storage,
and static `contributesVirial` value. A new compatible force type must provide
the same resource and lifecycle interface.

Adding periodic restraint behavior would require an explicit displacement
policy rather than merely reading the currently stored box vector. Adding a
virial requires defining and implementing the external-virial convention that
the current kernel comments defer.

Adding another `AT` or `CT` combination requires a valid device write
implementation and a corresponding explicit template instantiation in
`src/HarmonicRestraintForce.cu`.

### Performance-sensitive behavior

Force and energy writes are atomic. The force representation avoids floating
point accumulation in the shared manager output by using fixed-point integers.
Energy accumulation uses a block reduction followed by one atomic addition per
block.

Configuration setters are not hot-path operations and currently use
device-wide synchronization. Do not move those transfers into the per-step
force path.

### Tests and example

Relevant validation and behavior coverage is located in:

- `test/unittests/unittest-harmonicRestraintForce.cpp`;
- `test/unittests/unittest-capiHarmonicRestraintForce.cpp`;
- `test/pytest/python_api_harmonic_restraint_force.py`;
- `example/cons_harm.py`.

The native tests verify constructor defaults, validation, selection behavior,
force and energy values, clearing, and force-manager subscription. The direct
C ABI tests verify handles, arrays, statuses, diagnostics, subscription, and
diagnostic clearing. The Python tests verify conversion, validation,
subscription, error translation, and short propagation.

## API Reference

### Native C++

- [HarmonicRestraintForce](@ref HarmonicRestraintForce)
- [AtomSelection](@ref AtomSelection)
- [ForceManager](@ref ForceManager)

### C ABI

- [apo_harmonic_restraint_force](@ref apo_harmonic_restraint_force)
- [apo_harmonic_restraint_force_create](@ref apo_harmonic_restraint_force_create)
- [apo_harmonic_restraint_force_destroy](@ref apo_harmonic_restraint_force_destroy)
- [apo_harmonic_restraint_force_set_selection](@ref apo_harmonic_restraint_force_set_selection)
- [apo_harmonic_restraint_force_set_force_constant](@ref apo_harmonic_restraint_force_set_force_constant)
- [apo_harmonic_restraint_force_set_force_constants](@ref apo_harmonic_restraint_force_set_force_constants)
- [apo_harmonic_restraint_force_set_reference_coordinates](@ref apo_harmonic_restraint_force_set_reference_coordinates)
- [apo_harmonic_restraint_force_set_masses](@ref apo_harmonic_restraint_force_set_masses)
- [apo_harmonic_restraint_force_set_box_dimensions](@ref apo_harmonic_restraint_force_set_box_dimensions)
- [apo_force_manager_subscribe_harmonic_restraint_force](@ref apo_force_manager_subscribe_harmonic_restraint_force)
- [apo_force_manager_unsubscribe_harmonic_restraint_force](@ref apo_force_manager_unsubscribe_harmonic_restraint_force)

### Python

- [apocharmm.HarmonicRestraintForce](@ref python_harmonic_restraint_force)
- [apocharmm.ForceManager](@ref python_force_manager)
- [apocharmm.AtomSelection](@ref python_atom_selection)

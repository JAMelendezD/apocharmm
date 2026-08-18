# Subscriber {#subscriber}

## Purpose

The Subscriber subsystem schedules side effects during dynamics propagation.
A concrete @ref Subscriber writes selected simulation state when an attached
@ref CudaIntegrator reaches the subscriber's reporting interval. The interval
is a positive, dimensionless number of propagated steps.

Use @ref DcdSubscriber for coordinate trajectories and @ref RestartSubscriber
for restart state. Those two reporters are available through C++, the C ABI,
and Python. The repository also contains native-only text reporters for state,
dynamics, composite-force, alchemical, MBAR, and indexed-coordinate output.
Several native-only reporters expose current implementation limitations and
should be selected only after reading their class contracts.

Subscribers are synchronous observers, not background tasks. An update runs on
the propagation caller's thread after the corresponding dynamics step and may
perform force evaluation, CUDA transfers, device-wide synchronization, and file
I/O before propagation continues.

## Quick Start

The Python example below writes two DCD frames while propagating four steps. It
uses stable repository test data and unsubscribes before the context manager
closes the DCD handle.

```python
import apocharmm as apo

prm = apo.CharmmParameters("test/data/toppar_water_ions.str")
psf = apo.CharmmPsf("test/data/nacl_pair.psf")
crd = apo.CharmmCrd("test/data/nacl_pair.cor")

ctx = apo.CharmmContext(psf, prm)
ctx.setBoxDimensions([50.0, 50.0, 50.0])
ctx.setCoordinates(crd)
ctx.useHolonomicConstraints(False)
ctx.assignVelocitiesAtTemperature(300.0)

integrator = apo.CudaNoseHooverIntegrator(0.001)
integrator.setReferenceTemperature(300.0)
integrator.setNoseHooverPistonMass(300.0)
integrator.useOldTemperature(False)
integrator.setCharmmContext(ctx)

with apo.DcdSubscriber("trajectory.dcd", 2) as dcd:
    integrator.subscribe(dcd)
    integrator.propagate(4)
    integrator.unsubscribe(dcd)
```

The equivalent direct C++ workflow uses shared ownership because native
subscription obtains a shared self-reference from the integrator.

```cpp
#include "CharmmContext.h"
#include "CharmmCrd.h"
#include "CharmmPSF.h"
#include "CharmmParameters.h"
#include "CudaNoseHooverIntegrator.h"
#include "DcdSubscriber.h"

#include <memory>
#include <vector>

int main() {
  const auto prm = std::make_shared<CharmmParameters>(
      "test/data/toppar_water_ions.str");
  const auto psf =
      std::make_shared<CharmmPSF>("test/data/nacl_pair.psf");
  const auto crd =
      std::make_shared<CharmmCrd>("test/data/nacl_pair.cor");

  const auto ctx = std::make_shared<CharmmContext>(psf, prm);
  ctx->setBoxDimensions(std::vector<double>{50.0, 50.0, 50.0});
  ctx->setCoordinates(crd);
  ctx->useHolonomicConstraints(false);
  ctx->assignVelocitiesAtTemperature(300.0);

  const auto integrator =
      std::make_shared<CudaNoseHooverIntegrator>(0.001);
  integrator->setReferenceTemperature(300.0);
  integrator->setNoseHooverPistonMass(300.0);
  integrator->useOldTemperature(false);
  integrator->setCharmmContext(ctx);

  const auto dcd = std::make_shared<DcdSubscriber>("trajectory.dcd", 2);
  integrator->subscribe(dcd);
  integrator->propagate(4);
  integrator->unsubscribe(dcd);

  return 0;
}
```

## Construction and Required State

A concrete file-writing subscriber normally validates its positive report
frequency, copies its path, checks the nonempty parent prefix with `stat()`, and
creates or truncates the output during construction. @ref DcdSubscriber opens
binary output. The base @ref Subscriber and the remaining current reporters use
text output unless they override @ref Subscriber::openFile.

Establish state in this order:

1. Construct and initialize the @ref CharmmContext.
2. Construct the concrete integrator and attach the context.
3. Construct the subscriber with its final path and report frequency.
4. Subscribe it to the integrator.
5. Propagate dynamics.
6. Unsubscribe before closing or destroying either wrapper.

Native `CudaIntegrator::subscribe()` requires the integrator itself to be owned
by `std::shared_ptr`. It calls `Subscriber::setCharmmContext()` with the
integrator's current context, then `Subscriber::setIntegrator()`, then appends
the subscriber and a copy of its current frequency to parallel arrays. Set the
integrator context first. Subscribing while that context is null can succeed,
and a later context setter does not update the already attached subscriber.

The scheduling interval is a snapshot. Changing a subscriber's native
frequency after subscription changes the property read by some file metadata,
but it does not change the integrator's cached callback interval. Each separate
`propagate(num_steps)` call numbers steps locally from one. A subscriber with
frequency 100 is therefore called at local steps 100, 200, and so on in every
call, rather than from an absolute lifetime-step modulus.

@ref DcdSubscriber requires a live binary stream, a context with a positive
atom count, exactly three positive box lengths, and a valid single-precision
coordinate/charge device container. @ref RestartSubscriber additionally
requires the subscriber and integrator to retain the same context and supports
only @ref CudaNoseHooverIntegrator, @ref CudaLangevinPistonIntegrator, and
@ref CudaLangevinThermostatIntegrator.

Native-only reporters have additional preconditions documented on their class
and update methods. Several do not validate null attachments, force-manager
type, container shape, or stream state before dereferencing.

## Ownership and Lifetime

At the native layer, an integrator retains every subscribed object through
`std::shared_ptr<Subscriber>`. The subscriber retains the context and integrator
through shared pointers. This produces the following strong ownership graph:

```text
CudaIntegrator -> Subscriber -> CudaIntegrator
                       |
                       +------> CharmmContext
```

`CudaIntegrator::unsubscribe()` erases the integrator-to-subscriber edge, but it
does not clear the subscriber's context or integrator pointers. Unsubscribe
before releasing final owners. The retained backlinks also mean that the same
native subscriber cannot currently be resubscribed through the normal
attachment path.

The public C ABI owns concrete handles such as @ref apo_dcd_subscriber and
@ref apo_restart_subscriber. Each concrete handle embeds an @ref apo_subscriber
base view. Conversion returns a borrowed pointer into that concrete allocation;
it is not a separately owned handle and has no destroy function. Destroying the
concrete C handle invalidates the base pointer even when a native integrator
continues to retain the underlying C++ object. Keep the concrete handle alive
until after C ABI unsubscription.

Python `DcdSubscriber` and `RestartSubscriber` wrappers own their concrete C
handles and store the borrowed base pointer. A Python `CudaIntegrator` appends a
successfully subscribed wrapper to its `_subscribers` list, preserving the C
handle and base-view lifetime. `unsubscribe()` removes that retained Python
reference after native removal. Closing a subscriber clears its borrowed base
view, so close only after unsubscription.

The subsystem performs no internal host-thread locking. Serialize propagation,
updates, frequency changes, attachment, unsubscription, stream operations, and
destruction involving shared objects.

## Data, Shapes, and AKMA Units

The reporting interval and frame counters are dimensionless step counts.
Integrator constructor time steps and public `getTimeStep()` values use
picoseconds.

@ref DcdSubscriber writes:

- one six-`double` unit-cell record per frame, with X, Y, and Z box lengths in
  elements 0, 2, and 5;
- three Fortran-style coordinate records containing all X values, all Y values,
  and all Z values;
- one IEEE-like native `float` per coordinate component in atom order;
- coordinates and box lengths in angstroms;
- no charge, velocity, mass, force, or energy values.

The current DCD implementation writes native integer and floating-point
representations and does not perform byte-order conversion.

@ref RestartSubscriber writes a version-50 CHARMM-style text restart. Verified
quantities include:

| Section or value | Shape and order | Unit |
| --- | --- | --- |
| Box lengths | Three diagonal lengths in X, Y, Z order | Angstrom |
| `XOLD, YOLD, ZOLD` | One double-precision XYZ row per atom | Angstrom |
| `VX, VY, VZ` | One double-precision velocity row per atom | Angstrom per AKMA time unit |
| `X, Y, Z` | One previous coordinate-delta XYZ row per atom | Stored integrator representation |
| Average temperature | One scalar | Kelvin |
| Atom, step, degree-of-freedom, and seed fields | Scalar integer fields | Dimensionless |

Restart Nose-Hoover and Langevin-piston state is serialized directly from the
integrator. The current repository does not establish every one of those
field-specific units clearly enough to make them a stable subscriber contract.

Native text energy reporters use the force-manager convention of kilocalories
per mole. @ref StateSubscriber also reports time in picoseconds, temperature in
kelvin, box lengths in angstroms, volume in cubic angstroms, and density in grams
per cubic centimeter. Its pressure columns are labeled in atmospheres but are
currently populated with zeros. @ref XYZSubscriber writes zero-based atom
indices and coordinates in angstroms; despite its name, its output is not the
standard XYZ file format.

## Errors

C++ constructors and base mutators use @ref ApoCharmmError. Empty paths,
nonexistent checked parent paths, and nonpositive frequencies use
`ApoCharmmErrorCode::InvalidArgument`. Missing required update state uses
`ApoCharmmErrorCode::NotInitialized`. Unsupported restart-integrator types use
`ApoCharmmErrorCode::NotImplemented`. File open/write and state-consistency
failures use `ApoCharmmErrorCode::Runtime`. CUDA transfers, kernels, and
synchronization failures use `ApoCharmmErrorCode::Cuda`.

Legacy native-only reporters are not uniform. @ref DynaSubscriber and
@ref XYZSubscriber still expose verified standard-exception paths, and several
specialized reporters rely on unchecked preconditions that can fail before an
apoCHARMM exception is produced. Their individual API contracts identify those
boundaries.

C ABI status functions clear the calling thread's previous diagnostic at entry.
Native error categories map to the matching `APO_STATUS_*` value. Unexpected
standard or nonstandard exceptions map to `APO_STATUS_RUNTIME_ERROR`. On
failure, copy @ref apo_last_error immediately on the same thread; its pointer is
borrowed and changes after the next diagnostic-changing call on that thread.
Concrete destroy functions are void, accept `NULL`, and do not normally clear a
previous diagnostic.

Python configures each status-returning function with a shared error callback.
Every nonzero status raises `ApoCharmmError` containing the numeric status,
symbolic status name, operation context, and copied native diagnostic. Python
also raises `TypeError` for wrapper/type or `ctypes` conversion failures,
`ValueError` for values outside a signed 32-bit C `int`, and `RuntimeError` for
closed wrappers or impossible successful-null-handle results.

A propagation failure is not transactional at any interface. The integrator
increments its total requested-step count before the loop. Completed dynamics
steps, changed host/device state, and file bytes written before a later callback
failure remain observable.

## Important Behavior and Limitations

- Subscriber callbacks run synchronously and can substantially increase the
  wall time of the propagation step that triggers them.
- DCD coordinate output and restart output perform device-to-host transfers that
  call `cudaDeviceSynchronize()` on the current device. They do not use a
  subscriber-specific CUDA stream.
- `Subscriber::setFileName()` changes only the stored logical name. It does not
  reopen the stream. Mutable `getFileName()` access bypasses all validation.
- DCD files are native-endian and frame writes can be partial. DCD `NSTEP` is
  derived from frame count times the subscriber's current frequency, not the
  integrator's absolute step.
- Restart updates truncate before writing and currently fill many CHARMM energy
  and statistics fields with zero. The recorded `NSAVC` and `NSAVV` values do
  not represent independent DCD coordinate and velocity saving frequencies.
- @ref StateSubscriber does not currently emit box-size or volume values through
  its public flags, does not compute pressure values, and uses a kinetic-energy
  conversion inconsistent with the current scalar context API.
- @ref DynaSubscriber is implemented only for the Langevin-piston integrator and
  emits several placeholder values.
- @ref CompositeSubscriber and @ref EDSSubscriber report stored energies rather
  than requesting a fresh calculation.
- @ref BEDSSubscriber, @ref FEPSubscriber, and @ref MBARSubscriber contain
  unchecked force-manager casts. Their documented force-manager-type
  preconditions are mandatory.
- @ref MBARSubscriber's auxiliary force-manager list is not used by its current
  update path, and its derived context setter is not invoked polymorphically by
  integrator subscription.
- @ref DualTopologySubscriber and @ref XYZSubscriber do not write explicit frame
  delimiters suitable for their conventional format names.

## Related Subsystems

- @subpage charmm_context supplies coordinates, velocities, box state, energies,
  and the force manager queried by subscribers.
- @subpage force_manager supplies potential energies and specialized composite
  state.
- @subpage cuda_container defines the explicit device-to-host transfers and
  synchronization used by output writers.
- @subpage apocharmm_error defines native errors, C ABI diagnostics, and Python
  error translation.
- @subpage coordinates describes the coordinate inputs used to initialize the
  context before propagation.

## Developer Architecture

The native public layer begins at @ref Subscriber. Its extension point is the
pure virtual `update()` method. A concrete reporter normally delegates path,
stream, frequency, context, and integrator storage to the base and adds only the
format-specific state required to write one update. Override `openFile()` only
when the format requires non-default stream flags or format state reset, as
@ref DcdSubscriber does for binary output.

The scheduling collaborator is @ref CudaIntegrator. It stores parallel
`m_Subscribers` and `m_ReportFreqList` arrays. Entry `i` in one must always
correspond to entry `i` in the other. Subscription appends both; unsubscription
erases both; `reportIfNeeded()` indexes both and invokes `update()`. Mutable
accessors expose these arrays, so callers can violate the invariant. A future
reorganization should replace the parallel representation with one attachment
record and define explicit backlink clearing and resubscription semantics.

The native implementation layer under `src/` performs format generation and
CUDA transfers. DCD writes header, unit-cell, and coordinate records directly to
`std::fstream`. Restart assembles one monolithic CHARMM-style file and reads
integrator-specific containers through dynamic casts. Native-only reporters are
compiled into `apoCHARMMlib` by `src/CMakeLists.txt`; disabled checkpoint and
NetCDF subscriber sources are not part of the current target.

The C ABI exposes only the base operations, DCD writer, restart writer, and the
integrator attachment/propagation operations. Private handle structs store a
concrete `std::shared_ptr` and an embedded base struct that shares the same
native owner. C functions validate handles, translate exceptions with
`apocharmm_c::guard`, and never expose the private C++ layout as a stable public
contract.

The Python layer mirrors that mapping. Concrete wrappers own their C handles,
base methods operate on borrowed embedded views, and the integrator retains
Python wrappers to preserve those views. Status translation is centralized in
`configure_status_function()` rather than duplicated by subscriber wrappers.

Host/device data flow is explicit. DCD copies `CudaContainer<float4>` device
coordinates to host before splitting components. Restart copies each required
context or integrator container to host before text formatting. These transfers
are blocking at the device level. Other native reporters vary: some transfer
stored energy scalars, some invoke force calculations, and some assume helper
getters already expose host values.

The error boundary is layered. Public C++ contracts name native categories.
`apocharmm_c::guard` maps those categories and captures diagnostics. Python
`errcheck` immediately copies the thread-local C diagnostic into an owned
exception string. Destroy functions use a separate non-throwing boundary because
they cannot return status.

Relevant tests are registered in `test/unittests/CMakeLists.txt`. Native coverage
is concentrated in `unittest-subscriber.cpp`, `unittest-dcdSubscriber.cpp`,
`unittest-restartSubscriber.cpp`, and `unittest-cudaIntegrator.cpp`. Direct C
ABI coverage uses the corresponding `unittest-capi*.cpp` files. Python coverage
uses the subscriber, DCD-subscriber, restart-subscriber, and CUDA-integrator
pytest modules. The native-only specialized reporters currently lack dedicated
registered tests.

## API Reference

C++:

- @ref Subscriber
- @ref DcdSubscriber
- @ref RestartSubscriber
- @ref CudaIntegrator::subscribe
- @ref CudaIntegrator::unsubscribe
- @ref CudaIntegrator::propagate

C ABI:

- @ref apo_subscriber
- @ref apo_subscriber_set_report_frequency
- @ref apo_subscriber_get_report_frequency
- @ref apo_dcd_subscriber
- @ref apo_dcd_subscriber_create
- @ref apo_dcd_subscriber_create_with_report_frequency
- @ref apo_dcd_subscriber_destroy
- @ref apo_dcd_subscriber_as_subscriber
- @ref apo_restart_subscriber
- @ref apo_restart_subscriber_create
- @ref apo_restart_subscriber_create_with_report_frequency
- @ref apo_restart_subscriber_destroy
- @ref apo_restart_subscriber_as_subscriber
- @ref apo_cuda_integrator_subscribe
- @ref apo_cuda_integrator_unsubscribe
- @ref apo_cuda_integrator_propagate

Python:

- @ref python_subscriber
- @ref python_dcd_subscriber
- @ref python_restart_subscriber

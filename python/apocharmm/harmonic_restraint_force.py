# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

"""
@brief Exposes the owning Python wrapper for absolute harmonic restraints.

`HarmonicRestraintForce` copies Python configuration through the C ABI and can
be subscribed to `ForceManager`. Force calculation and native output access are
not exposed directly by this module.

@anchor python_harmonic_restraint_force_module
@see harmonic_restraint_force
"""

from collections.abc import Sequence
import ctypes

from ._base import _ApoObject
from ._lib import lib
from .atom_selection import AtomSelection
from .error import configure_status_function
from .force_manager import ForceManager

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    configure_status_function(
        lib().apo_harmonic_restraint_force_create,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_int],
        "HarmonicRestraintForce construction",
    )

    lib().apo_harmonic_restraint_force_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_harmonic_restraint_force_destroy.restype = None

    configure_status_function(
        lib().apo_harmonic_restraint_force_set_selection,
        [ctypes.c_void_p, ctypes.c_void_p],
        "HarmonicRestraintForce.setSelection(selection)",
    )

    configure_status_function(
        lib().apo_harmonic_restraint_force_set_force_constant,
        [ctypes.c_void_p, ctypes.c_double],
        "HarmonicRestraintForce.setForceConstant(force_constant)",
    )

    configure_status_function(
        lib().apo_harmonic_restraint_force_set_force_constants,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "HarmonicRestraintForce.setForceConstants(force_constants)",
    )

    configure_status_function(
        lib().apo_harmonic_restraint_force_set_reference_coordinates,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "HarmonicRestraintForce.setReferenceCoordinates(reference_coordinates)",
    )

    configure_status_function(
        lib().apo_harmonic_restraint_force_set_masses,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "HarmonicRestraintForce.setMasses(masses)",
    )

    configure_status_function(
        lib().apo_harmonic_restraint_force_set_box_dimensions,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "HarmonicRestraintForce.setBoxDimensions(box_dimensions)",
    )

    configure_status_function(
        lib().apo_force_manager_subscribe_harmonic_restraint_force,
        [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_char_p],
        "ForceManager.subscribe(HarmonicRestraintForce)",
    )

    configure_status_function(
        lib().apo_force_manager_unsubscribe_harmonic_restraint_force,
        [ctypes.c_void_p, ctypes.c_void_p],
        "ForceManager.unsubscribe(HarmonicRestraintForce)",
    )

    _prototypes_initialized = True

    return


class HarmonicRestraintForce(_ApoObject):
    """
    @brief Configures an owned native absolute harmonic restraint.

    The restraint applies one absolute Cartesian reference position, force
    constant, and mass or weight per atom. See the subsystem page for the
    energy equation, force-array sign convention, data layout, and units.

    This wrapper owns one C handle for the native
    `HarmonicRestraintForce<long long int, float>` specialization. Python
    sequences are converted to `float`, copied into temporary contiguous
    `ctypes` buffers, and copied again into native storage. No caller sequence
    or buffer is retained.

    `ForceManager.subscribe()` retains this wrapper after native subscription
    succeeds. The native manager independently retains the restraint and its
    calculation resources. Unsubscribe before calling `close()`: closing a
    subscribed wrapper destroys the C handle but does not remove the native
    subscription, and the closed wrapper can no longer identify that
    subscription for normal unsubscription.

    `close()`, `destroy()`, context-manager exit, and finalization release the
    owned C handle. Explicit closure is idempotent. Methods that access
    `self.handle` after closure raise `RuntimeError`.

    The wrapper and native object provide no internal synchronization. Do not
    overlap configuration, manager subscription changes, calculation, or
    closure across threads.

    @anchor python_harmonic_restraint_force
    @see harmonic_restraint_force
    """

    _destroy_function_name = "apo_harmonic_restraint_force_destroy"

    def __init__(self, num_atoms: int) -> None:
        """
        @brief Constructs a restraint for a fixed atom count.

        `num_atoms` is checked as a Python `int`, converted to signed C `int`,
        and passed to the native constructor. The wrapper owns the newly
        returned handle and sets `default_force_tag` to `"harm"`.

        @param[in] num_atoms Dimensionless positive atom count representable by
        signed 32-bit C `int`. The current `isinstance(value, int)` check also
        accepts Python `bool`; `True` reaches native construction as one and
        `False` is rejected natively as zero.
        @throws TypeError If `num_atoms` is not a Python `int`.
        @throws ValueError If `num_atoms` is outside the signed 32-bit C `int`
        range.
        @throws ApoCharmmError If native validation rejects a non-positive
        value, native CUDA setup fails, or native allocation fails.
        @throws RuntimeError If native construction reports success but
        produces a NULL handle.
        @throws MemoryError If Python or `ctypes` bookkeeping cannot allocate
        storage.

        @post On success, this wrapper owns a live native handle and
        `default_force_tag == "harm"`.
        """
        super().__init__()

        if not isinstance(num_atoms, int):
            raise TypeError("num_atoms must be an int")

        if num_atoms < -(2**31) or num_atoms > 2**31 - 1:
            raise ValueError("num_atoms must fit in int")

        _initialize_prototypes()

        handle: ctypes.c_void_p = ctypes.c_void_p()
        c_num_atoms: ctypes.c_int = ctypes.c_int(num_atoms)

        lib().apo_harmonic_restraint_force_create(ctypes.byref(handle), c_num_atoms)

        if handle.value is None:
            raise RuntimeError(
                "apo_harmonic_restraint_force_create returned success but produced a NULL handle"
            )

        self._handle = handle
        self._default_force_tag: str = "harm"

        return

    @property
    def default_force_tag(self) -> str:
        """
        @brief Returns the default manager energy tag.

        This property reads Python-owned text and does not access the native
        handle, so it remains available after `close()`.

        @return The Python string `"harm"`.
        """
        return self._default_force_tag

    def setSelection(self, selection: AtomSelection) -> None:
        """
        @brief Copies the atom selection used by the restraint.

        The native selection is copied. This wrapper does not retain
        `selection`. Unselected atoms have their force constants set to zero;
        expanding the selection later does not restore constants that were
        previously zeroed.

        @param[in] selection Live `AtomSelection` with the same atom count as
        this restraint.
        @throws TypeError If `selection` is not an `AtomSelection`.
        @throws RuntimeError If this wrapper or `selection` has been closed.
        @throws ApoCharmmError If native handle validation or atom-count
        validation fails, a native allocation fails, or the native
        host-to-device update fails.
        """
        _initialize_prototypes()

        if not isinstance(selection, AtomSelection):
            raise TypeError(
                "HarmonicRestraintForce.setSelection expects an AtomSelection"
            )

        lib().apo_harmonic_restraint_force_set_selection(self.handle, selection.handle)

        return

    def setForceConstant(self, force_constant: float) -> None:
        """
        @brief Sets one force constant for all selected atoms.

        `force_constant` is converted with `float()` and then to C `double`.
        Unselected atoms receive zero. With masses in atomic mass units, the
        value has units of kilocalories per mole per square angstrom per atomic
        mass unit.

        @param[in] force_constant Any Python object accepted by `float()`. The
        resulting value must be finite and non-negative.
        @throws TypeError If `float(force_constant)` rejects the object.
        @throws ValueError If Python conversion rejects the value.
        @throws OverflowError If Python conversion overflows.
        @throws RuntimeError If this wrapper has been closed.
        @throws ApoCharmmError If native validation rejects the converted
        value or the native CUDA update fails.
        """
        _initialize_prototypes()

        c_force_constant: ctypes.c_double = ctypes.c_double(float(force_constant))

        lib().apo_harmonic_restraint_force_set_force_constant(
            self.handle, c_force_constant
        )

        return

    def setForceConstants(self, force_constants: Sequence[float]) -> None:
        """
        @brief Copies one force constant per atom.

        The method consumes `force_constants` as an iterable, converts every
        element with `float()`, materializes a new Python list, and copies that
        list into a temporary contiguous C `double` array. Native code then
        copies the array again.

        The iterable must produce exactly one value per atom. Only currently
        selected values are validated and used; unselected values are ignored
        and their native constants are zero.

        @param[in] force_constants Iterable of values accepted by `float()`, in
        atom order. Selected values must be finite and non-negative.
        @throws TypeError If the object is not iterable or an element cannot be
        converted with `float()`.
        @throws ValueError If element conversion fails.
        @throws OverflowError If element conversion overflows.
        @throws MemoryError If the temporary list or C buffer cannot be
        allocated.
        @throws RuntimeError If this wrapper has been closed.
        @throws ApoCharmmError If the native length or selected-value contract
        is violated or the native CUDA update fails.
        """
        _initialize_prototypes()

        force_constant_values: list[float] = [float(value) for value in force_constants]

        c_buffer_type = ctypes.c_double * len(force_constant_values)
        c_buffer = c_buffer_type(*force_constant_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(force_constant_values))

        lib().apo_harmonic_restraint_force_set_force_constants(
            self.handle, c_buffer, c_buffer_len
        )

        return

    def setReferenceCoordinates(
        self, reference_coordinates: Sequence[Sequence[float]]
    ) -> None:
        """
        @brief Copies absolute Cartesian reference coordinates.

        The method consumes an outer iterable of coordinate rows. Every row is
        consumed as an iterable, every element is converted with `float()`, and
        each converted row must contain exactly three values in `[x, y, z]`
        order. Rows are flattened into a newly allocated row-major C `double`
        buffer.

        The outer iterable must produce exactly one row per atom. Coordinates
        use angstroms and are not minimum-image wrapped or rescaled when the box
        changes.

        @param[in] reference_coordinates Iterable of three-element coordinate
        iterables in atom order. Every converted value must be finite.
        @throws TypeError If an outer or inner value is not iterable or a
        scalar cannot be converted with `float()`.
        @throws ValueError If scalar conversion fails or a converted row does
        not contain exactly three elements.
        @throws OverflowError If scalar conversion overflows.
        @throws MemoryError If the flattened list or C buffer cannot be
        allocated.
        @throws RuntimeError If this wrapper has been closed.
        @throws ApoCharmmError If the native atom count is wrong, a converted
        coordinate is non-finite, or the native CUDA transfer fails.

        @warning A native non-finite-value failure can leave an earlier host
        prefix updated while the native device mirror remains old.
        """
        _initialize_prototypes()

        flattened_coordinates: list[float] = []

        for i, coordinate in enumerate(reference_coordinates):
            coordinate_values: list[float] = [float(value) for value in coordinate]

            if len(coordinate_values) != 3:
                raise ValueError(
                    f"reference_coordinates[{i}] must contain exactly 3 elements"
                )

            flattened_coordinates.extend(coordinate_values)

        c_buffer_type = ctypes.c_double * len(flattened_coordinates)
        c_buffer = c_buffer_type(*flattened_coordinates)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(flattened_coordinates))

        lib().apo_harmonic_restraint_force_set_reference_coordinates(
            self.handle, c_buffer, c_buffer_len
        )

        return

    def setMasses(self, masses: Sequence[float]) -> None:
        """
        @brief Copies one mass or restraint weight per atom.

        The method consumes `masses` as an iterable, converts every element with
        `float()`, materializes a new Python list, and copies it into a
        temporary contiguous C `double` buffer. Native code copies the values
        again.

        Values normally use atomic mass units. Every value must be finite and
        non-negative. Zero is accepted and suppresses that atom's harmonic
        contribution.

        @param[in] masses Iterable producing exactly one mass or weight per
        atom, in atom order.
        @throws TypeError If `masses` is not iterable or an element cannot be
        converted with `float()`.
        @throws ValueError If element conversion fails.
        @throws OverflowError If element conversion overflows.
        @throws MemoryError If the temporary list or C buffer cannot be
        allocated.
        @throws RuntimeError If this wrapper has been closed.
        @throws ApoCharmmError If the native length or value contract is
        violated or the native CUDA transfer fails.

        @warning A native invalid-value failure can leave an earlier host
        prefix updated while the native device mirror remains old.
        """
        _initialize_prototypes()

        mass_values: list[float] = [float(value) for value in masses]

        c_buffer_type = ctypes.c_double * len(mass_values)
        c_buffer = c_buffer_type(*mass_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(mass_values))

        lib().apo_harmonic_restraint_force_set_masses(
            self.handle, c_buffer, c_buffer_len
        )

        return

    def setBoxDimensions(self, box_dimensions: Sequence[float]) -> None:
        """
        @brief Copies three orthogonal box lengths.

        The method consumes `box_dimensions` as an iterable, converts every
        value with `float()`, and copies the resulting list through a temporary
        C `double` buffer. Exactly three finite, positive values are required in
        `[x, y, z]` order and angstrom units.

        The current restraint kernel stores but does not use these dimensions.
        Changing them does not wrap displacements or rescale references.

        @param[in] box_dimensions Iterable of exactly three box lengths.
        @throws TypeError If the object is not iterable or a value cannot be
        converted with `float()`.
        @throws ValueError If value conversion fails.
        @throws OverflowError If value conversion overflows.
        @throws MemoryError If the temporary list or C buffer cannot be
        allocated.
        @throws RuntimeError If this wrapper has been closed.
        @throws ApoCharmmError If native length or value validation fails or a
        changed box cannot be transferred to the device.
        """
        _initialize_prototypes()

        box_values: list[float] = [float(value) for value in box_dimensions]

        c_buffer_type = ctypes.c_double * len(box_values)
        c_buffer = c_buffer_type(*box_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(box_values))

        lib().apo_harmonic_restraint_force_set_box_dimensions(
            self.handle, c_buffer, c_buffer_len
        )

        return

    def _subscribe_to_force_manager(
        self, force_manager: ForceManager, force_tag: str | None = None
    ) -> None:
        """
        @brief Implements the `ForceManager.subscribe()` callback.

        `None` selects `default_force_tag`; otherwise `force_tag` must be a
        Python `str`. The selected text is encoded as UTF-8 and passed as a
        borrowed null-terminated C string. Native code copies the tag and
        retains the restraint and calculation resources on success.

        This callback performs the native subscription. The calling
        `ForceManager.subscribe()` method retains this Python wrapper only after
        the callback succeeds.

        @param[in] force_manager Live `ForceManager` receiving the restraint.
        @param[in] force_tag Python `str` tag, or `None` to use `"harm"`.
        @throws TypeError If `force_manager` is not a `ForceManager` or a
        non-`None` tag is not a `str`.
        @throws UnicodeEncodeError If `force_tag` cannot be encoded as UTF-8.
        @throws RuntimeError If this wrapper or `force_manager` has been closed.
        @throws ApoCharmmError If the tag is empty, either native handle is
        invalid, the restraint is already subscribed, initialization detects a
        mismatch, or native allocation or CUDA setup fails.

        @warning An embedded NUL character terminates the C string and causes
        native code to observe only the preceding tag prefix.
        """
        _initialize_prototypes()

        if not isinstance(force_manager, ForceManager):
            raise TypeError("force_manager must be a ForceManager")

        if force_tag is None:
            force_tag_value: str = self._default_force_tag
        else:
            if not isinstance(force_tag, str):
                raise TypeError("force_tag must be a str")

            force_tag_value = force_tag

        encoded_force_tag: bytes = force_tag_value.encode("utf-8")
        c_force_tag: ctypes.c_char_p = ctypes.c_char_p(encoded_force_tag)

        lib().apo_force_manager_subscribe_harmonic_restraint_force(
            force_manager.handle, self.handle, c_force_tag
        )

        return

    def _unsubscribe_from_force_manager(self, force_manager: ForceManager) -> None:
        """
        @brief Implements the `ForceManager.unsubscribe()` callback.

        Native unsubscription releases the manager's retained restraint,
        stream, force-array, and energy owners. After this callback succeeds,
        the calling `ForceManager.unsubscribe()` method releases its retained
        Python reference.

        @param[in] force_manager Live `ForceManager` currently containing this
        restraint.
        @throws RuntimeError If this wrapper or `force_manager` has been closed.
        @throws ApoCharmmError If either native handle is invalid, this
        restraint is not subscribed to the manager, or an unexpected native
        runtime failure occurs.
        """
        _initialize_prototypes()

        lib().apo_force_manager_unsubscribe_harmonic_restraint_force(
            force_manager.handle, self.handle
        )

        return

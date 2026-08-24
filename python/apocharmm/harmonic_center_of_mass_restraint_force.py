# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

"""
@brief Exposes the owning Python wrapper for harmonic group-center restraints.

`HarmonicCenterOfMassRestraintForce` copies Python configuration through the C
ABI and can be subscribed to `ForceManager`. Force calculation and native
output access are not exposed directly by this module.

@anchor python_harmonic_center_of_mass_restraint_force_module
@see harmonic_center_of_mass_restraint_force
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
        lib().apo_harmonic_center_of_mass_restraint_force_create,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_int],
        "HarmonicCenterOfMassRestraintForce construction",
    )

    lib().apo_harmonic_center_of_mass_restraint_force_destroy.argtypes = [
        ctypes.c_void_p
    ]
    lib().apo_harmonic_center_of_mass_restraint_force_destroy.restype = None

    configure_status_function(
        lib().apo_harmonic_center_of_mass_restraint_force_set_selection,
        [ctypes.c_void_p, ctypes.c_void_p],
        "HarmonicCenterOfMassRestraintForce.setSelection(selection)",
    )

    configure_status_function(
        lib().apo_harmonic_center_of_mass_restraint_force_set_force_constant,
        [ctypes.c_void_p, ctypes.c_double],
        "HarmonicCenterOfMassRestraintForce.setForceConstant(force_constant)",
    )

    configure_status_function(
        lib().apo_harmonic_center_of_mass_restraint_force_set_reference_position,
        [
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_int),
            ctypes.c_size_t,
        ],
        "HarmonicCenterOfMassRestraintForce.setReferencePosition(reference_position)",
    )

    configure_status_function(
        lib().apo_harmonic_center_of_mass_restraint_force_set_reference_distance,
        [ctypes.c_void_p, ctypes.c_double],
        "HarmonicCenterOfMassRestraintForce.setReferenceDistance(reference_distance)",
    )

    configure_status_function(
        lib().apo_harmonic_center_of_mass_restraint_force_set_masses,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "HarmonicCenterOfMassRestraintForce.setMasses(masses)",
    )

    configure_status_function(
        lib().apo_harmonic_center_of_mass_restraint_force_set_mass_weighting,
        [ctypes.c_void_p, ctypes.c_bool],
        "HarmonicCenterOfMassRestraintForce.setMassWeighting(flag)",
    )

    configure_status_function(
        lib().apo_force_manager_subscribe_harmonic_center_of_mass_restraint_force,
        [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_char_p],
        "ForceManager.subscribe(HarmonicCenterOfMassRestraintForce)",
    )

    configure_status_function(
        lib().apo_force_manager_unsubscribe_harmonic_center_of_mass_restraint_force,
        [ctypes.c_void_p, ctypes.c_void_p],
        "ForceManager.unsubscribe(HarmonicCenterOfMassRestraintForce)",
    )

    _prototypes_initialized = True

    return


class HarmonicCenterOfMassRestraintForce(_ApoObject):
    """
    @brief Configures an owned native harmonic group-center restraint.

    The selected atoms define a uniformly weighted geometric center by default
    or a mass-weighted center after masses are supplied or mass weighting is
    enabled. The native force applies periodic CONS HMCM-style energy and
    gradient behavior described on the subsystem page.

    This wrapper owns one C handle for the native
    `HarmonicCenterOfMassRestraintForce<long long int, float>` specialization.
    Python iterables are materialized, converted to native scalar types, copied
    into temporary contiguous `ctypes` arrays, and copied again into native
    storage. No caller iterable or buffer is retained.

    `ForceManager.subscribe()` retains this wrapper only after native
    subscription succeeds. The native manager independently retains the
    restraint and its CUDA resources. Unsubscribe before calling `close()`:
    closing a subscribed wrapper destroys the identifying C handle but does not
    remove the native subscription, and the closed wrapper can no longer be
    used for normal unsubscription.

    `close()`, `destroy()`, context-manager exit, and finalization release the
    owned C handle. Explicit closure is idempotent. Methods that access
    `self.handle` after closure raise `RuntimeError`; `default_force_tag` remains
    available because it is Python-owned text.

    The wrapper and native object provide no internal synchronization. Do not
    overlap configuration, manager operations, calculation, or closure across
    threads.

    @anchor python_harmonic_center_of_mass_restraint_force
    @see harmonic_center_of_mass_restraint_force
    """

    _destroy_function_name = "apo_harmonic_center_of_mass_restraint_force_destroy"

    def __init__(self, num_atoms: int) -> None:
        """
        @brief Constructs a restraint for a fixed atom count.

        `num_atoms` is checked with `isinstance(value, int)`, range-checked for
        signed 32-bit C `int`, and passed to native construction. The wrapper
        owns the returned handle and sets `default_force_tag` to `"hmcm"`.

        @param[in] num_atoms Positive dimensionless atom count representable by
        signed 32-bit C `int`. Because Python `bool` is an `int` subclass,
        `True` reaches native construction as one and `False` is rejected
        natively as zero.
        @throws TypeError If `num_atoms` is not a Python `int`.
        @throws ValueError If `num_atoms` is outside the signed 32-bit C `int`
        range.
        @throws ApoCharmmError If native validation rejects a non-positive
        value, CUDA setup fails, or native allocation fails.
        @throws RuntimeError If native construction reports success but returns
        a NULL handle.
        @throws MemoryError If Python or `ctypes` bookkeeping cannot allocate
        storage.

        @post On success, this wrapper owns a live native handle and
        `default_force_tag == "hmcm"`.
        """
        super().__init__()

        if not isinstance(num_atoms, int):
            raise TypeError("num_atoms must be an int")

        if num_atoms < -(2**31) or num_atoms > 2**31 - 1:
            raise ValueError("num_atoms must fit in int")

        _initialize_prototypes()

        handle: ctypes.c_void_p = ctypes.c_void_p()
        c_num_atoms: ctypes.c_int = ctypes.c_int(num_atoms)

        lib().apo_harmonic_center_of_mass_restraint_force_create(
            ctypes.byref(handle), c_num_atoms
        )

        if handle.value is None:
            raise RuntimeError(
                "apo_harmonic_center_of_mass_restraint_force_create returned success but produced a NULL handle"
            )

        self._handle = handle
        self._default_force_tag: str = "hmcm"

        return

    @property
    def default_force_tag(self) -> str:
        """
        @brief Returns the default force-manager energy tag.

        This property reads Python-owned text and does not access the native
        handle, so it remains available after `close()`.

        @return The Python string `"hmcm"`.
        """
        return self._default_force_tag

    def setSelection(self, selection: AtomSelection) -> None:
        """
        @brief Copies the nonempty selection used to define the group center.

        The native selection is copied; this wrapper retains no Python
        reference to `selection`. Native selected-index and selected-weight
        arrays are rebuilt and copied to device memory.

        @param[in] selection Live `AtomSelection` with the same atom count as
        this restraint and at least one selected atom.
        @throws TypeError If `selection` is not an `AtomSelection`.
        @throws RuntimeError If this wrapper or `selection` has been closed.
        @throws ApoCharmmError If native handle, atom-count, nonempty-selection,
        or positive-total-weight validation fails, native allocation fails, or
        the native CUDA update fails.

        @warning A native failure after mutation begins can leave newer
        selection state paired with older selected device arrays.
        """
        _initialize_prototypes()

        if not isinstance(selection, AtomSelection):
            raise TypeError(
                "HarmonicCenterOfMassRestraintForce.setSelection expects an AtomSelection"
            )

        lib().apo_harmonic_center_of_mass_restraint_force_set_selection(
            self.handle, selection.handle
        )

        return

    def setForceConstant(self, force_constant: float) -> None:
        """
        @brief Sets the scalar harmonic force constant.

        `force_constant` is converted with `float()` and then to C `double`.
        The value uses kilocalories per mole per square angstrom. Zero disables
        evaluation.

        @param[in] force_constant Any Python object accepted by `float()`. The
        converted value must be finite and non-negative.
        @throws TypeError If `float(force_constant)` rejects the object.
        @throws ValueError If Python floating-point conversion rejects the
        value.
        @throws OverflowError If conversion to Python `float` overflows.
        @throws RuntimeError If this wrapper has been closed.
        @throws ApoCharmmError If native validation rejects the converted value.
        """
        _initialize_prototypes()

        c_force_constant: ctypes.c_double = ctypes.c_double(float(force_constant))

        lib().apo_harmonic_center_of_mass_restraint_force_set_force_constant(
            self.handle, c_force_constant
        )

        return

    def setReferencePosition(
        self,
        reference_position: Sequence[float],
        reference_mask: Sequence[int] | None = None,
    ) -> None:
        """
        @brief Copies the Cartesian reference position and component mask.

        `reference_position` is consumed as an iterable and every value is
        converted with `float()`. It must produce exactly three values in
        `[x, y, z]` order and angstrom units.

        When `reference_mask` is `None`, all three components are active. A
        supplied mask is consumed as an iterable and every value is converted
        with `int()`, then narrowed through `ctypes.c_int`, before native
        validation. Integer conversion can truncate non-integral numeric
        objects, and narrowing can change an out-of-range integer; the resulting
        C values must each be zero or one and at least one must be one.

        Both converted lists are copied into temporary contiguous `ctypes`
        arrays and copied again into native host scalar state. No caller object
        or buffer is retained. For example,
        `setReferencePosition((0.0, 0.0, 0.0), (1, 1, 0))` restrains only the X
        and Y components.

        @param[in] reference_position Iterable of exactly three values accepted
        by `float()`, in `[x, y, z]` order. Converted values must be finite.
        @param[in] reference_mask `None` or an iterable of exactly three values
        accepted by `int()`. Converted values must form an active zero-or-one
        mask in `[x, y, z]` order.
        @throws TypeError If either supplied object is not iterable or an
        element cannot be converted by `float()` or `int()`.
        @throws ValueError If scalar conversion rejects a value.
        @throws OverflowError If scalar conversion overflows.
        @throws MemoryError If a temporary list or C array cannot be allocated.
        @throws RuntimeError If this wrapper has been closed.
        @throws ApoCharmmError If native length, finite-value, mask-value, or
        active-component validation fails.
        """
        _initialize_prototypes()

        reference_position_values: list[float] = [
            float(value) for value in reference_position
        ]

        if reference_mask is None:
            reference_mask_values: list[int] = [1, 1, 1]
        else:
            reference_mask_values: list[int] = [int(value) for value in reference_mask]

        c_buffer_type1 = ctypes.c_double * len(reference_position_values)
        c_buffer1 = c_buffer_type1(*reference_position_values)
        c_buffer_len1: ctypes.c_size_t = ctypes.c_size_t(len(reference_position_values))

        c_buffer_type2 = ctypes.c_int * len(reference_mask_values)
        c_buffer2 = c_buffer_type2(*reference_mask_values)
        c_buffer_len2: ctypes.c_size_t = ctypes.c_size_t(len(reference_mask_values))

        lib().apo_harmonic_center_of_mass_restraint_force_set_reference_position(
            self.handle, c_buffer1, c_buffer_len1, c_buffer2, c_buffer_len2
        )

        return

    def setReferenceDistance(self, reference_distance: float) -> None:
        """
        @brief Sets the target radial distance from the reference position.

        `reference_distance` is converted with `float()` and then to C `double`.
        Zero selects `k * r^2`; a positive value selects
        `k * (r - reference_distance)^2`.

        @param[in] reference_distance Any Python object accepted by `float()`.
        The converted distance must be finite, non-negative, and is measured in
        angstroms.
        @throws TypeError If `float(reference_distance)` rejects the object.
        @throws ValueError If Python floating-point conversion rejects the
        value.
        @throws OverflowError If conversion to Python `float` overflows.
        @throws RuntimeError If this wrapper has been closed.
        @throws ApoCharmmError If native validation rejects the converted value.

        @note For a positive target and exactly zero current displacement, the
        native implementation assigns a zero gradient at the
        nondifferentiable point.
        """
        _initialize_prototypes()

        c_reference_distance: ctypes.c_double = ctypes.c_double(
            float(reference_distance)
        )

        lib().apo_harmonic_center_of_mass_restraint_force_set_reference_distance(
            self.handle, c_reference_distance
        )

        return

    def setMasses(self, masses: Sequence[float]) -> None:
        """
        @brief Copies per-atom masses and enables mass weighting.

        `masses` is consumed as an iterable, every element is converted with
        `float()`, and the result is copied through a temporary contiguous C
        `double` array. Native code copies the complete atom-order array and
        rebuilds selected device weights. No caller iterable or buffer is
        retained.

        Values normally use atomic mass units. Only relative selected masses
        affect the normalized center. Every mass must be finite and
        non-negative, and selected masses must have positive total weight.

        @param[in] masses Iterable producing exactly one value per atom, in atom
        order. Every element must be accepted by `float()`.
        @throws TypeError If `masses` is not iterable or an element cannot be
        converted with `float()`.
        @throws ValueError If element conversion rejects a value.
        @throws OverflowError If conversion to Python `float` overflows.
        @throws MemoryError If the temporary list or C array cannot be
        allocated.
        @throws RuntimeError If this wrapper has been closed.
        @throws ApoCharmmError If native length, finite-value, non-negativity,
        or positive-total-weight validation fails, or a native CUDA update
        fails.

        @post On success, mass weighting is enabled.
        @warning A native failure after mutation begins can leave new masses or
        a new weighting mode paired with older selected device arrays.
        """
        _initialize_prototypes()

        mass_values: list[float] = [float(value) for value in masses]

        c_buffer_type = ctypes.c_double * len(mass_values)
        c_buffer = c_buffer_type(*mass_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(mass_values))

        lib().apo_harmonic_center_of_mass_restraint_force_set_masses(
            self.handle, c_buffer, c_buffer_len
        )

        return

    def setMassWeighting(self, flag: bool) -> None:
        """
        @brief Enables configured mass weights or restores uniform weights.

        The implementation passes `flag` to `ctypes.c_bool` without an explicit
        Python type check. False-valued objects disable mass weighting and
        true-valued objects enable it. Native selected weights are rebuilt and
        copied to device storage; disabling does not discard stored masses.

        @param[in] flag Object whose truth value can initialize
        `ctypes.c_bool`. Use an actual Python `bool` for an unambiguous API call.
        @throws Exception If custom truth-value conversion raises an exception.
        @throws RuntimeError If this wrapper has been closed.
        @throws ApoCharmmError If enabling produces zero selected total mass, a
        native allocation fails, or the native CUDA update fails.

        @warning A native failure after the mode flag changes can leave the
        weighting mode and selected device arrays inconsistent.
        """
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool(flag)

        lib().apo_harmonic_center_of_mass_restraint_force_set_mass_weighting(
            self.handle, c_flag
        )

        return

    def _subscribe_to_force_manager(
        self, force_manager: ForceManager, force_tag: str | None = None
    ) -> None:
        """
        @brief Implements the `ForceManager.subscribe()` callback.

        `None` selects `default_force_tag`; otherwise `force_tag` must be a
        Python `str`. The selected text is encoded as UTF-8 and passed as a
        temporary null-terminated C string. Native subscription retains the
        restraint and its CUDA resources, not this Python wrapper or its C
        handle. The public manager method retains this wrapper only after this
        callback succeeds.

        @param[in] force_manager Live `ForceManager` receiving the restraint.
        @param[in] force_tag `None` or a Python string. The encoded string must
        be nonempty after C-string interpretation.
        @throws TypeError If `force_manager` is not a `ForceManager` or
        `force_tag` is neither `None` nor `str`.
        @throws UnicodeEncodeError If UTF-8 encoding rejects the string.
        @throws RuntimeError If either native wrapper is closed.
        @throws ApoCharmmError If native handle validation, tag validation,
        duplicate detection, immediate initialization, or manager allocation
        fails.

        @post On success, the native manager owns independent shared references
        to the restraint and its resources.
        @warning An embedded NUL character truncates the tag at the C ABI
        boundary.
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

        lib().apo_force_manager_subscribe_harmonic_center_of_mass_restraint_force(
            force_manager.handle, self.handle, c_force_tag
        )

        return

    def _unsubscribe_from_force_manager(self, force_manager: ForceManager) -> None:
        """
        @brief Implements the `ForceManager.unsubscribe()` callback.

        Native unsubscription removes the matching restraint by object identity
        and releases manager-owned shared references. The public manager method
        removes its retained Python reference only after this callback succeeds.

        @param[in] force_manager Live `ForceManager` that currently owns the
        native subscription.
        @throws RuntimeError If this wrapper or `force_manager` is closed.
        @throws ApoCharmmError If native handle validation fails or the
        restraint is not subscribed to that manager.

        @post On success, the restraint no longer participates in manager
        calculation or box propagation.
        """
        _initialize_prototypes()

        lib().apo_force_manager_unsubscribe_harmonic_center_of_mass_restraint_force(
            force_manager.handle, self.handle
        )

        return

# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

from collections.abc import Sequence
import ctypes

from ._lib import lib
from .enums import CrystalType
from .error import check_status
from .cuda_integrator import CudaIntegrator
import math

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    lib().apo_cuda_langevin_piston_integrator_create.argtypes = [
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_double,
    ]
    lib().apo_cuda_langevin_piston_integrator_create.restype = ctypes.c_int

    lib().apo_cuda_langevin_piston_integrator_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_cuda_langevin_piston_integrator_destroy.restype = None

    lib().apo_cuda_langevin_piston_integrator_use_nose_hoover_thermostat.argtypes = [
        ctypes.c_void_p,
        ctypes.c_bool,
    ]
    lib().apo_cuda_langevin_piston_integrator_use_nose_hoover_thermostat.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_set_reference_temperature.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
    ]
    lib().apo_cuda_langevin_piston_integrator_set_reference_temperature.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_set_nose_hoover_piston_mass.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
    ]
    lib().apo_cuda_langevin_piston_integrator_set_nose_hoover_piston_mass.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_use_old_temperature.argtypes = [
        ctypes.c_void_p,
        ctypes.c_bool,
    ]
    lib().apo_cuda_langevin_piston_integrator_use_old_temperature.restype = ctypes.c_int

    lib().apo_cuda_langevin_piston_integrator_set_reference_pressure.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
    ]
    lib().apo_cuda_langevin_piston_integrator_set_reference_pressure.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_set_constant_surface_tension.argtypes = [
        ctypes.c_void_p,
        ctypes.c_bool,
    ]
    lib().apo_cuda_langevin_piston_integrator_set_constant_surface_tension.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_set_crystal_type.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
    ]
    lib().apo_cuda_langevin_piston_integrator_set_crystal_type.restype = ctypes.c_int

    lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_mass.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
    ]
    lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_mass.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_friction_seed.argtypes = [
        ctypes.c_void_p,
        ctypes.c_uint64,
    ]
    lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_friction_seed.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_friction.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
    ]
    lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_friction.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_reset_averages.argtypes = [
        ctypes.c_void_p
    ]
    lib().apo_cuda_langevin_piston_integrator_reset_averages.restype = ctypes.c_int

    lib().apo_cuda_langevin_piston_integrator_get_reference_temperature.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_piston_integrator_get_reference_temperature.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_get_nose_hoover_piston_mass.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_piston_integrator_get_nose_hoover_piston_mass.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_get_average_temperature.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_piston_integrator_get_average_temperature.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_get_reference_pressure_tensor.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_piston_integrator_get_reference_pressure_tensor.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_get_crystal_type.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_piston_integrator_get_crystal_type.restype = ctypes.c_int

    lib().apo_cuda_langevin_piston_integrator_get_langevin_piston_mass.argtypes = [
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_piston_integrator_get_langevin_piston_mass.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_tensor.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_tensor.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_scalar.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_scalar.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_get_average_pressure_tensor.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_piston_integrator_get_average_pressure_tensor.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_get_average_pressure_scalar.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_piston_integrator_get_average_pressure_scalar.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_get_langevin_piston_friction_seed.argtypes = [
        ctypes.POINTER(ctypes.c_uint64),
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_piston_integrator_get_langevin_piston_friction_seed.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_get_instantaneous_temperature.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_piston_integrator_get_instantaneous_temperature.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_piston_integrator_as_cuda_integrator.argtypes = [
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_piston_integrator_as_cuda_integrator.restype = ctypes.c_int

    _prototypes_initialized = True

    return


class CudaLangevinPistonIntegrator(CudaIntegrator):
    _destroy_function_name = "apo_cuda_langevin_piston_integrator_destroy"

    def __init__(self, time_step: float) -> None:
        _initialize_prototypes()
        super().__init__()

        handle: ctypes.c_void_p = ctypes.c_void_p()
        c_time_step: ctypes.c_double = ctypes.c_double(time_step)

        status = lib().apo_cuda_langevin_piston_integrator_create(
            ctypes.byref(handle), c_time_step
        )

        check_status(status, "CudaLangevinPistonIntegrator construction failed")

        if handle.value is None:
            raise RuntimeError(
                "apo_cuda_langevin_piston_integrator_create returned success but produced a NULL handle"
            )

        self._handle = handle

        integrator_handle: ctypes.c_void_p = ctypes.c_void_p()

        status = lib().apo_cuda_langevin_piston_integrator_as_cuda_integrator(
            ctypes.byref(integrator_handle), self.handle
        )

        check_status(
            status, "CudaLangevinPistonIntegrator base-integrator conversion failed"
        )

        if integrator_handle.value is None:
            raise RuntimeError(
                "apo_cuda_langevin_piston_integrator_as_cuda_integrator returned success but produced a NULL handle"
            )

        self._integrator_handle = integrator_handle

        return

    def useNoseHooverThermostat(self, flag: bool = True) -> None:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool(flag)

        status = lib().apo_cuda_langevin_piston_integrator_use_nose_hoover_thermostat(
            self.handle, c_flag
        )

        check_status(
            status, "CudaLangevinPistonIntegrator.useNoseHooverThermostat(flag) failed"
        )

        return

    def setReferenceTemperature(self, temperature: float) -> None:
        _initialize_prototypes()

        c_temperature: ctypes.c_double = ctypes.c_double(temperature)

        status = lib().apo_cuda_langevin_piston_integrator_set_reference_temperature(
            self.handle, c_temperature
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.setReferenceTemperature(temperature) failed",
        )

        return

    def setNoseHooverPistonMass(self, mass: float) -> None:
        _initialize_prototypes()

        c_mass: ctypes.c_double = ctypes.c_double(mass)

        status = lib().apo_cuda_langevin_piston_integrator_set_nose_hoover_piston_mass(
            self.handle, c_mass
        )

        check_status(
            status, "CudaLangevinPistonIntegrator.setNoseHooverPistonMass(mass) failed"
        )

        return

    def useOldTemperature(self, flag: bool = True) -> None:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool(flag)

        status = lib().apo_cuda_langevin_piston_integrator_use_old_temperature(
            self.handle, c_flag
        )

        check_status(
            status, "CudaLangevinPistonIntegrator.useOldTemperature(flag) failed"
        )

        return

    def setReferencePressure(self, pressure_tensor: Sequence[float]) -> None:
        _initialize_prototypes()

        if isinstance(pressure_tensor, (str, bytes, bytearray)):
            raise TypeError("pressure_tensor must be a 9-element sequence of floats")

        pressure_values: list[float] = []

        try:
            iterator = iter(pressure_tensor)
        except TypeError as exc:
            raise TypeError("pressure_tensor must be a 9-element sequence of floats")

        for index, value in enumerate(iterator):
            if isinstance(value, bool):
                raise TypeError(f"pressure_tensor[{index}] must be a float, not bool")

            pressure_value: float = float(value)

            if not math.isfinite(pressure_value):
                raise ValueError(f"pressure_tensor[{index}] must be finite")

            pressure_values.append(pressure_value)

        if len(pressure_values) != 9:
            raise ValueError(
                "pressure_tensor must contain exactly 9 elements (XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ)"
            )

        pressure_tensor_array = ctypes.c_double * 9
        c_pressure_tensor = pressure_tensor_array(*pressure_values)
        c_pressure_tensor_len: ctypes.c_size_t = ctypes.c_size_t(len(pressure_values))

        status = lib().apo_cuda_langevin_piston_integrator_set_reference_pressure(
            self.handle, c_pressure_tensor, c_pressure_tensor_len
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.setReferencePressure(pressure_tensor) failed",
        )

        return

    def setConstantSurfaceTension(self, flag: bool = True) -> None:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool(flag)

        status = lib().apo_cuda_langevin_piston_integrator_set_constant_surface_tension(
            self.handle, c_flag
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.setConstantSurfaceTension(flag) failed",
        )

        return

    def setCrystalType(self, crystal_type: CrystalType | int) -> None:
        _initialize_prototypes()

        if isinstance(crystal_type, bool):
            raise TypeError("crystal_type must be a CrystalType")

        try:
            crystal_type_value: CrystalType = CrystalType(crystal_type)
        except ValueError as exc:
            raise ValueError(f"invalid crystal_type: {crystal_type!r}") from exc

        if crystal_type_value == CrystalType.NONE:
            raise ValueError("crystal_type must be CUBIC, TETRAGONAL, or ORTHORHOMBIC")

        c_crystal_type: ctypes.c_int = ctypes.c_int(int(crystal_type_value))

        status = lib().apo_cuda_langevin_piston_integrator_set_crystal_type(
            self.handle, c_crystal_type
        )

        check_status(
            status, "CudaLangevinPistonIntegrator.setCrystalType(crystal_type) failed"
        )

        return

    def setLangevinPistonMass(self, mass: list[float]) -> None:
        _initialize_prototypes()

        if isinstance(mass, (str, bytes, bytearray)):
            raise TypeError("mass must be a sequence of 1, 2, or 3 non-negative floats")

        mass_values: list[float] = []

        try:
            iterator = iter(mass)
        except TypeError as exc:
            raise TypeError("mass must be a sequence of 1, 2, or 3 non-negative floats")

        for index, value in enumerate(iterator):
            if isinstance(value, bool):
                raise TypeError(f"mass[{index}] must be a float, not bool")

            mass_value: float = float(value)

            if not math.isfinite(mass_value):
                raise ValueError(f"mass[{index}] must be finite")

            if mass_value < 0.0:
                raise ValueError(f"mass[{index}] must be non-negative")

            mass_values.append(mass_value)

        if len(mass_values) == 0 or len(mass_values) > 3:
            raise ValueError(
                "mass must contain 1, 2, or 3 elements depending on crystal type"
            )

        mass_array = ctypes.c_double * len(mass_values)
        c_mass = mass_array(*mass_values)
        c_mass_len: ctypes.c_size_t = ctypes.c_size_t(len(mass_values))

        status = lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_mass(
            self.handle, c_mass, c_mass_len
        )

        check_status(
            status, "CudaLangevinPistonIntegrator.setLangevinPistonMass(mass) failed"
        )

        return

    def setLangevinPistonFrictionSeed(self, seed: int) -> None:
        _initialize_prototypes()

        if isinstance(seed, bool) or seed < 0 or seed > 2**64 - 1:
            raise ValueError("seed must fit in uint64_t")

        c_seed: ctypes.c_uint64 = ctypes.c_uint64(seed)

        status = (
            lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_friction_seed(
                self.handle, c_seed
            )
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.setLangevinPistonFrictionSeed(seed) failed",
        )

        return

    def setLangevinPistonFriction(self, friction: float) -> None:
        _initialize_prototypes()

        c_friction: ctypes.c_double = ctypes.c_double(friction)

        status = lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_friction(
            self.handle, c_friction
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.setLangevinPistonFriction(friction) failed",
        )

        return

    def resetAverages(self) -> None:
        _initialize_prototypes()

        status = lib().apo_cuda_langevin_piston_integrator_reset_averages(self.handle)

        check_status(status, "CudaLangevinPistonIntegrator.resetAverages() failed")

        return

    def getReferenceTemperature(self) -> float:
        _initialize_prototypes()

        reference_temperature = ctypes.c_double()

        status = lib().apo_cuda_langevin_piston_integrator_get_reference_temperature(
            ctypes.byref(reference_temperature), self.handle
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.getReferenceTemperature() failed",
        )

        return float(reference_temperature.value)

    def getNoseHooverPistonMass(self) -> float:
        _initialize_prototypes()

        nose_hoover_piston_mass = ctypes.c_double()

        status = lib().apo_cuda_langevin_piston_integrator_get_nose_hoover_piston_mass(
            ctypes.byref(nose_hoover_piston_mass), self.handle
        )

        check_status(
            status, "CudaLangevinPistonIntegrator.getNoseHooverPistonMass() failed"
        )

        return float(nose_hoover_piston_mass.value)

    def getAverageTemperature(self) -> float:
        _initialize_prototypes()

        average_temperature = ctypes.c_double()

        status = lib().apo_cuda_langevin_piston_integrator_get_average_temperature(
            ctypes.byref(average_temperature), self.handle
        )

        check_status(
            status, "CudaLangevinPistonIntegrator.getAverageTemperature() failed"
        )

        return float(average_temperature.value)

    def getReferencePressureTensor(self) -> tuple[float, ...]:
        _initialize_prototypes()

        pressure_tensor_array = ctypes.c_double * 9
        c_pressure_tensor = pressure_tensor_array()
        c_pressure_tensor_len: ctypes.c_size_t = ctypes.c_size_t(9)

        status = (
            lib().apo_cuda_langevin_piston_integrator_get_reference_pressure_tensor(
                c_pressure_tensor, c_pressure_tensor_len, self.handle
            )
        )

        check_status(
            status, "CudaLangevinPistonIntegrator.getReferencePressureTensor() failed"
        )

        return tuple(float(c_pressure_tensor[index]) for index in range(9))

    def getCrystalType(self) -> CrystalType:
        _initialize_prototypes()

        crystal_type = ctypes.c_int()

        status = lib().apo_cuda_langevin_piston_integrator_get_crystal_type(
            ctypes.byref(crystal_type), self.handle
        )

        check_status(status, "CudaLangevinPistonIntegrator.getCrystalType() failed")

        return CrystalType(crystal_type.value)

    def getLangevinPistonMass(self) -> tuple[float, ...]:
        _initialize_prototypes()

        max_num_mass: int = 3
        mass_array = ctypes.c_double * max_num_mass
        c_mass = mass_array()
        c_mass_len: ctypes.c_size_t = ctypes.c_size_t(max_num_mass)
        c_num_mass: ctypes.c_size_t = ctypes.c_size_t()

        status = lib().apo_cuda_langevin_piston_integrator_get_langevin_piston_mass(
            ctypes.byref(c_num_mass), c_mass, c_mass_len, self.handle
        )

        check_status(
            status, "CudaLangevinPistonIntegrator.getLangevinPistonMass() failed"
        )

        num_mass: int = int(c_num_mass.value)

        if num_mass > max_num_mass:
            raise RuntimeError(
                "apo_cuda_langevin_piston_integrator_get_langevin_piston_mass returned more than 3 masses"
            )

        return tuple(float(c_mass[index]) for index in range(num_mass))

    def getInstantaneousPressureTensor(self) -> tuple[float, ...]:
        _initialize_prototypes()

        pressure_tensor_array = ctypes.c_double * 9
        c_pressure_tensor = pressure_tensor_array()
        c_pressure_tensor_len: ctypes.c_size_t = ctypes.c_size_t(9)

        status = (
            lib().apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_tensor(
                c_pressure_tensor, c_pressure_tensor_len, self.handle
            )
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.getInstantaneousPressureTensor() failed",
        )

        return tuple(float(c_pressure_tensor[index]) for index in range(9))

    def getInstantaneousPressureScalar(self) -> float:
        _initialize_prototypes()

        pressure_scalar = ctypes.c_double()

        status = (
            lib().apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_scalar(
                ctypes.byref(pressure_scalar), self.handle
            )
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.getInstantaneousPressureScalar() failed",
        )

        return float(pressure_scalar.value)

    def getAveragePressureTensor(self) -> tuple[float, ...]:
        _initialize_prototypes()

        pressure_tensor_array = ctypes.c_double * 9
        c_pressure_tensor = pressure_tensor_array()
        c_pressure_tensor_len: ctypes.c_size_t = ctypes.c_size_t(9)

        status = lib().apo_cuda_langevin_piston_integrator_get_average_pressure_tensor(
            c_pressure_tensor, c_pressure_tensor_len, self.handle
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.getAveragePressureTensor() failed",
        )

        return tuple(float(c_pressure_tensor[index]) for index in range(9))

    def getAveragePressureScalar(self) -> float:
        _initialize_prototypes()

        pressure_scalar = ctypes.c_double()

        status = lib().apo_cuda_langevin_piston_integrator_get_average_pressure_scalar(
            ctypes.byref(pressure_scalar), self.handle
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.getAveragePressureScalar() failed",
        )

        return float(pressure_scalar.value)

    def getLangevinPistonFrictionSeed(self) -> int:
        _initialize_prototypes()

        seed = ctypes.c_uint64()

        status = (
            lib().apo_cuda_langevin_piston_integrator_get_langevin_piston_friction_seed(
                ctypes.byref(seed), self.handle
            )
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.getLangevinPistonFrictionSeed() failed",
        )

        return int(seed.value)

    def getInstantaneousTemperature(self) -> float:
        _initialize_prototypes()

        instantaneous_temperature = ctypes.c_double()

        status = (
            lib().apo_cuda_langevin_piston_integrator_get_instantaneous_temperature(
                ctypes.byref(instantaneous_temperature), self.handle
            )
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.getInstantaneousTemperature() failed",
        )

        return float(instantaneous_temperature.value)

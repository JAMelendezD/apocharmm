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

    def setReferencePressure(self, pressure_tensor: Sequence[Sequence[float]]) -> None:
        _initialize_prototypes()

        flattened_pressure_tensor: list[float] = []

        for i, pressure in enumerate(pressure_tensor):
            pressure_values: list[float] = [float(value) for value in pressure]

            if len(pressure_values) != 3:
                raise ValueError(
                    f"pressure_tensor[{i}] must contain exactly 3 elements"
                )

            flattened_pressure_tensor.extend(pressure_values)

        c_buffer_type = ctypes.c_double * len(flattened_pressure_tensor)
        c_buffer = c_buffer_type(*flattened_pressure_tensor)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(flattened_pressure_tensor))

        status = lib().apo_cuda_langevin_piston_integrator_set_reference_pressure(
            self.handle, c_buffer, c_buffer_len
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

        try:
            crystal_type_value: CrystalType = CrystalType(crystal_type)
        except ValueError as exc:
            raise ValueError(f"invalid crystal_type: {crystal_type!r}") from exc

        c_crystal_type: ctypes.c_int = ctypes.c_int(int(crystal_type_value))

        status = lib().apo_cuda_langevin_piston_integrator_set_crystal_type(
            self.handle, c_crystal_type
        )

        check_status(
            status, "CudaLangevinPistonIntegrator.setCrystalType(crystal_type) failed"
        )

        return

    def setLangevinPistonMass(self, masses: Sequence[float]) -> None:
        _initialize_prototypes()

        mass_values: list[float] = [float(value) for value in masses]

        c_buffer_type = ctypes.c_double * len(mass_values)
        c_buffer = c_buffer_type(*mass_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(mass_values))

        status = lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_mass(
            self.handle, c_buffer, c_buffer_len
        )

        check_status(
            status, "CudaLangevinPistonIntegrator.setLangevinPistonMass(mass) failed"
        )

        return

    def setLangevinPistonFrictionSeed(self, seed: int) -> None:
        _initialize_prototypes()

        if seed < 0 or seed > 2**61 - 1:
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

        c_temperature = ctypes.c_double()

        status = lib().apo_cuda_langevin_piston_integrator_get_reference_temperature(
            ctypes.byref(c_temperature), self.handle
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.getReferenceTemperature() failed",
        )

        return float(c_temperature.value)

    def getNoseHooverPistonMass(self) -> float:
        _initialize_prototypes()

        c_mass = ctypes.c_double()

        status = lib().apo_cuda_langevin_piston_integrator_get_nose_hoover_piston_mass(
            ctypes.byref(c_mass), self.handle
        )

        check_status(
            status, "CudaLangevinPistonIntegrator.getNoseHooverPistonMass() failed"
        )

        return float(c_mass.value)

    def getAverageTemperature(self) -> float:
        _initialize_prototypes()

        c_temperature = ctypes.c_double()

        status = lib().apo_cuda_langevin_piston_integrator_get_average_temperature(
            ctypes.byref(c_temperature), self.handle
        )

        check_status(
            status, "CudaLangevinPistonIntegrator.getAverageTemperature() failed"
        )

        return float(c_temperature.value)

    def getReferencePressureTensor(self) -> list[list[float]]:
        _initialize_prototypes()

        c_buffer_type = ctypes.c_double * 9
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(9)

        status = (
            lib().apo_cuda_langevin_piston_integrator_get_reference_pressure_tensor(
                c_buffer, c_buffer_len, self.handle
            )
        )

        check_status(
            status, "CudaLangevinPistonIntegrator.getReferencePressureTensor() failed"
        )

        pressure_tensor: list[list[float]] = []
        for i in range(3):
            pressure_tensor.append(
                [
                    float(c_buffer[i * 3 + 0]),
                    float(c_buffer[i * 3 + 1]),
                    float(c_buffer[i * 3 + 2]),
                ]
            )

        return pressure_tensor

    def getCrystalType(self) -> CrystalType:
        _initialize_prototypes()

        c_crystal_type = ctypes.c_int()

        status = lib().apo_cuda_langevin_piston_integrator_get_crystal_type(
            ctypes.byref(c_crystal_type), self.handle
        )

        check_status(status, "CudaLangevinPistonIntegrator.getCrystalType() failed")

        return CrystalType(c_crystal_type.value)

    def getLangevinPistonMass(self) -> list[float]:
        _initialize_prototypes()

        c_num_mass: ctypes.c_size_t = ctypes.c_size_t()

        c_buffer_type = ctypes.c_double * 3
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(3)

        status = lib().apo_cuda_langevin_piston_integrator_get_langevin_piston_mass(
            ctypes.byref(c_num_mass), c_buffer, c_buffer_len, self.handle
        )

        check_status(
            status, "CudaLangevinPistonIntegrator.getLangevinPistonMass() failed"
        )

        num_mass: int = int(c_num_mass.value)

        if num_mass > 3:
            raise RuntimeError(
                "apo_cuda_langevin_piston_integrator_get_langevin_piston_mass returned more than 3 masses"
            )

        masses: list[float] = []
        for i in range(num_mass):
            masses.append(float(c_buffer[i]))

        return masses

    def getInstantaneousPressureTensor(self) -> list[list[float]]:
        _initialize_prototypes()

        c_buffer_type = ctypes.c_double * 9
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(9)

        status = (
            lib().apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_tensor(
                c_buffer, c_buffer_len, self.handle
            )
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.getInstantaneousPressureTensor() failed",
        )

        pressure_tensor: list[list[float]] = []
        for i in range(3):
            pressure_tensor.append(
                [
                    float(c_buffer[i * 3 + 0]),
                    float(c_buffer[i * 3 + 1]),
                    float(c_buffer[i * 3 + 2]),
                ]
            )

        return pressure_tensor

    def getInstantaneousPressureScalar(self) -> float:
        _initialize_prototypes()

        c_pressure_scalar = ctypes.c_double()

        status = (
            lib().apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_scalar(
                ctypes.byref(c_pressure_scalar), self.handle
            )
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.getInstantaneousPressureScalar() failed",
        )

        return float(c_pressure_scalar.value)

    def getAveragePressureTensor(self) -> list[list[float]]:
        _initialize_prototypes()

        c_buffer_type = ctypes.c_double * 9
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(9)

        status = lib().apo_cuda_langevin_piston_integrator_get_average_pressure_tensor(
            c_buffer, c_buffer_len, self.handle
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.getAveragePressureTensor() failed",
        )

        pressure_tensor: list[list[float]] = []
        for i in range(3):
            pressure_tensor.append(
                [
                    float(c_buffer[i * 3 + 0]),
                    float(c_buffer[i * 3 + 1]),
                    float(c_buffer[i * 3 + 2]),
                ]
            )

        return pressure_tensor

    def getAveragePressureScalar(self) -> float:
        _initialize_prototypes()

        c_pressure_scalar = ctypes.c_double()

        status = lib().apo_cuda_langevin_piston_integrator_get_average_pressure_scalar(
            ctypes.byref(c_pressure_scalar), self.handle
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.getAveragePressureScalar() failed",
        )

        return float(c_pressure_scalar.value)

    def getLangevinPistonFrictionSeed(self) -> int:
        _initialize_prototypes()

        c_seed = ctypes.c_uint64()

        status = (
            lib().apo_cuda_langevin_piston_integrator_get_langevin_piston_friction_seed(
                ctypes.byref(c_seed), self.handle
            )
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.getLangevinPistonFrictionSeed() failed",
        )

        return int(c_seed.value)

    def getInstantaneousTemperature(self) -> float:
        _initialize_prototypes()

        c_temperature = ctypes.c_double()

        status = (
            lib().apo_cuda_langevin_piston_integrator_get_instantaneous_temperature(
                ctypes.byref(c_temperature), self.handle
            )
        )

        check_status(
            status,
            "CudaLangevinPistonIntegrator.getInstantaneousTemperature() failed",
        )

        return float(c_temperature.value)

# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

import ctypes

from ._lib import lib
from .error import check_status
from .cuda_integrator import CudaIntegrator

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    lib().apo_cuda_langevin_thermostat_integrator_create.argtypes = [
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_double,
    ]
    lib().apo_cuda_langevin_thermostat_integrator_create.restype = ctypes.c_int

    lib().apo_cuda_langevin_thermostat_integrator_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_cuda_langevin_thermostat_integrator_destroy.restype = None

    lib().apo_cuda_langevin_thermostat_integrator_set_reference_temperature.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
    ]
    lib().apo_cuda_langevin_thermostat_integrator_set_reference_temperature.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_thermostat_integrator_set_thermostat_friction.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
    ]
    lib().apo_cuda_langevin_thermostat_integrator_set_thermostat_friction.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_thermostat_integrator_set_thermostat_rng_seed.argtypes = [
        ctypes.c_void_p,
        ctypes.c_uint64,
    ]
    lib().apo_cuda_langevin_thermostat_integrator_set_thermostat_rng_seed.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_thermostat_integrator_reset_average_temperature.argtypes = [
        ctypes.c_void_p
    ]
    lib().apo_cuda_langevin_thermostat_integrator_reset_average_temperature.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_thermostat_integrator_get_reference_temperature.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_thermostat_integrator_get_reference_temperature.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_thermostat_integrator_get_thermostat_friction.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_thermostat_integrator_get_thermostat_friction.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_thermostat_integrator_get_thermostat_rng_seed.argtypes = [
        ctypes.POINTER(ctypes.c_uint64),
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_thermostat_integrator_get_thermostat_rng_seed.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_thermostat_integrator_get_average_temperature.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_thermostat_integrator_get_average_temperature.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_thermostat_integrator_get_instantaneous_temperature.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_thermostat_integrator_get_instantaneous_temperature.restype = (
        ctypes.c_int
    )

    lib().apo_cuda_langevin_thermostat_integrator_as_cuda_integrator.argtypes = [
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_void_p,
    ]
    lib().apo_cuda_langevin_thermostat_integrator_as_cuda_integrator.restype = (
        ctypes.c_int
    )

    _prototypes_initialized = True

    return


class CudaLangevinThermostatIntegrator(CudaIntegrator):
    _destroy_function_name = "apo_cuda_langevin_thermostat_integrator_destroy"

    def __init__(self, time_step: float) -> None:
        _initialize_prototypes()
        super().__init__()

        handle: ctypes.c_void_p = ctypes.c_void_p()
        c_time_step: ctypes.c_double = ctypes.c_double(time_step)

        status = lib().apo_cuda_langevin_thermostat_integrator_create(
            ctypes.byref(handle), c_time_step
        )

        check_status(status, "CudaLangevinThermostatIntegrator construction failed")

        if handle.value is None:
            raise RuntimeError(
                "apo_cuda_langevin_thermostat_integrator_create returned success but produced a NULL handle"
            )

        self._handle = handle

        integrator_handle: ctypes.c_void_p = ctypes.c_void_p()

        status = lib().apo_cuda_langevin_thermostat_integrator_as_cuda_integrator(
            ctypes.byref(integrator_handle), self.handle
        )

        check_status(
            status, "CudaLangevinThermostatIntegrator base-integrator conversion failed"
        )

        if integrator_handle.value is None:
            raise RuntimeError(
                "apo_cuda_langevin_thermostat_integrator_as_cuda_integrator returned success but produced a NULL handle"
            )

        self._integrator_handle = integrator_handle

        return

    def setReferenceTemperature(self, temperature: float) -> None:
        _initialize_prototypes()

        c_temperature: ctypes.c_double = ctypes.c_double(temperature)

        status = (
            lib().apo_cuda_langevin_thermostat_integrator_set_reference_temperature(
                self.handle, c_temperature
            )
        )

        check_status(
            status,
            "CudaLangevinThermostatIntegrator.setReferenceTemperature(temperature) failed",
        )

        return

    def setThermostatFriction(self, friction: float) -> None:
        _initialize_prototypes()

        c_friction: ctypes.c_double = ctypes.c_double(friction)

        status = lib().apo_cuda_langevin_thermostat_integrator_set_thermostat_friction(
            self.handle, c_friction
        )

        check_status(
            status,
            "CudaLangevinThermostatIntegrator.setThermostatFriction(friction) failed",
        )

        return

    def setThermostatRngSeed(self, seed: int) -> None:
        _initialize_prototypes()

        if isinstance(seed, bool) or seed < 0 or seed > 2**64 - 1:
            raise ValueError("seed must fit in uint64_t")

        c_seed: ctyps.c_uint64 = ctypes.c_uint64(seed)

        status = lib().apo_cuda_langevin_thermostat_integrator_set_thermostat_rng_seed(
            self.handle, c_seed
        )

        check_status(
            status, "CudaLangevinThermostatIntegrator.setThermostatRngSeed(seed) faield"
        )

        return

    def resetAverageTemperature(self) -> None:
        _initialize_prototypes()

        status = (
            lib().apo_cuda_langevin_thermostat_integrator_reset_average_temperature(
                self.handle
            )
        )

        check_status(
            status, "CudaLangevinThermostatIntegrator.resetAverageTemperature() failed"
        )

        return

    def getReferenceTemperature(self) -> float:
        _initialize_prototypes()

        reference_temperature = ctypes.c_double()

        status = (
            lib().apo_cuda_langevin_thermostat_integrator_get_reference_temperature(
                ctypes.byref(reference_temperature), self.handle
            )
        )

        check_status(
            status,
            "CudaLangevinThermostatIntegrator.getReferenceTemperature() failed",
        )

        return float(reference_temperature.value)

    def getThermostatFriction(self) -> float:
        _initialize_prototypes()

        thermostat_friction = ctypes.c_double()

        status = lib().apo_cuda_langevin_thermostat_integrator_get_thermostat_friction(
            ctypes.byref(thermostat_friction), self.handle
        )

        check_status(
            status,
            "CudaLangevinThermostatIntegrator.getThermostatFriction() failed",
        )

        return float(thermostat_friction.value)

    def getThermostatRngSeed(self) -> int:
        _initialize_prototypes()

        seed = ctypes.c_uint64()

        status = lib().apo_cuda_langevin_thermostat_integrator_get_thermostat_rng_seed(
            ctypes.byref(seed), self.handle
        )

        check_status(
            status,
            "CudaLangevinThermostatIntegrator.getThermostatRngSeed() failed",
        )

        return float(seed.value)

    def getAverageTemperature(self) -> float:
        _initialize_prototypes()

        average_temperature = ctypes.c_double()

        status = lib().apo_cuda_langevin_thermostat_integrator_get_average_temperature(
            ctypes.byref(average_temperature), self.handle
        )

        check_status(
            status, "CudaLangevinThermostatIntegrator.getAverageTemperature() failed"
        )

        return float(average_temperature.value)

    def getInstantaneousTemperature(self) -> float:
        _initialize_prototypes()

        instantaneous_temperature = ctypes.c_double()

        status = (
            lib().apo_cuda_langevin_thermostat_integrator_get_instantaneous_temperature(
                ctypes.byref(instantaneous_temperature), self.handle
            )
        )

        check_status(
            status,
            "CudaLangevinThermostatIntegrator.getInstantaneousTemperature() failed",
        )

        return float(instantaneous_temperature.value)

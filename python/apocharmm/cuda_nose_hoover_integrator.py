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
from .error import configure_status_function
from .cuda_integrator import CudaIntegrator

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    configure_status_function(
        lib().apo_cuda_nose_hoover_integrator_create,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_double],
        "CudaNoseHooverIntegrator construction",
    )

    lib().apo_cuda_nose_hoover_integrator_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_cuda_nose_hoover_integrator_destroy.restype = None

    configure_status_function(
        lib().apo_cuda_nose_hoover_integrator_set_reference_temperature,
        [ctypes.c_void_p, ctypes.c_double],
        "CudaNoseHooverIntegrator.setReferenceTemperature(temperature)",
    )

    configure_status_function(
        lib().apo_cuda_nose_hoover_integrator_set_nose_hoover_piston_mass,
        [ctypes.c_void_p, ctypes.c_double],
        "CudaNoseHooverIntegrator.setNoseHooverPistonMass(mass)",
    )

    configure_status_function(
        lib().apo_cuda_nose_hoover_integrator_use_old_temperature,
        [ctypes.c_void_p, ctypes.c_bool],
        "CudaNoseHooverIntegrator.useOldTemperature(flag)",
    )

    configure_status_function(
        lib().apo_cuda_nose_hoover_integrator_reset_average_temperature,
        [ctypes.c_void_p],
        "CudaNoseHooverIntegrator.resetAverageTemperature()",
    )

    configure_status_function(
        lib().apo_cuda_nose_hoover_integrator_get_reference_temperature,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CudaNoseHooverIntegrator.getReferenceTemperature()",
    )

    configure_status_function(
        lib().apo_cuda_nose_hoover_integrator_get_nose_hoover_piston_mass,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CudaNoseHooverIntegrator.getNoseHooverPistonMass()",
    )

    configure_status_function(
        lib().apo_cuda_nose_hoover_integrator_get_average_temperature,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CudaNoseHooverIntegrator.getAverageTemperature()",
    )

    configure_status_function(
        lib().apo_cuda_nose_hoover_integrator_get_instantaneous_temperature,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CudaNoseHooverIntegrator.getInstantaneousTemperature()",
    )

    configure_status_function(
        lib().apo_cuda_nose_hoover_integrator_as_cuda_integrator,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p],
        "CudaNoseHooverIntegrator base-integrator coversion",
    )

    _prototypes_initialized = True

    return


class CudaNoseHooverIntegrator(CudaIntegrator):
    _destroy_function_name = "apo_cuda_nose_hoover_integrator_destroy"

    def __init__(self, time_step: float) -> None:
        _initialize_prototypes()
        super().__init__()

        handle: ctypes.c_void_p = ctypes.c_void_p()
        c_time_step: ctypes.c_double = ctypes.c_double(time_step)

        lib().apo_cuda_nose_hoover_integrator_create(ctypes.byref(handle), c_time_step)

        if handle.value is None:
            raise RuntimeError(
                "apo_cuda_nose_hoover_integrator_create returned success but produced a NULL handle"
            )

        self._handle = handle

        integrator_handle: ctypes.c_void_p = ctypes.c_void_p()

        lib().apo_cuda_nose_hoover_integrator_as_cuda_integrator(
            ctypes.byref(integrator_handle), self.handle
        )

        if integrator_handle.value is None:
            raise RuntimeError(
                "apo_cuda_nose_hoover_integrator_as_cuda_integrator returned success but produced a NULL handle"
            )

        self._integrator_handle = integrator_handle

        return

    def setReferenceTemperature(self, temperature: float) -> None:
        _initialize_prototypes()

        c_temperature: ctypes.c_double = ctypes.c_double(temperature)

        lib().apo_cuda_nose_hoover_integrator_set_reference_temperature(
            self.handle, c_temperature
        )

        return

    def setNoseHooverPistonMass(self, mass: float) -> None:
        _initialize_prototypes()

        c_mass: ctypes.c_double = ctypes.c_double(mass)

        lib().apo_cuda_nose_hoover_integrator_set_nose_hoover_piston_mass(
            self.handle, c_mass
        )

        return

    def useOldTemperature(self, flag: bool = True) -> None:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool(flag)

        lib().apo_cuda_nose_hoover_integrator_use_old_temperature(self.handle, c_flag)

        return

    def resetAverageTemperature(self) -> None:
        _initialize_prototypes()

        lib().apo_cuda_nose_hoover_integrator_reset_average_temperature(self.handle)

        return

    def getReferenceTemperature(self) -> float:
        _initialize_prototypes()

        c_temperature = ctypes.c_double()

        lib().apo_cuda_nose_hoover_integrator_get_reference_temperature(
            ctypes.byref(c_temperature), self.handle
        )

        return float(c_temperature.value)

    def getNoseHooverPistonMass(self) -> float:
        _initialize_prototypes()

        c_mass = ctypes.c_double()

        lib().apo_cuda_nose_hoover_integrator_get_nose_hoover_piston_mass(
            ctypes.byref(c_mass), self.handle
        )

        return float(c_mass.value)

    def getAverageTemperature(self) -> float:
        _initialize_prototypes()

        c_temperature = ctypes.c_double()

        lib().apo_cuda_nose_hoover_integrator_get_average_temperature(
            ctypes.byref(c_temperature), self.handle
        )

        return float(c_temperature.value)

    def getInstantaneousTemperature(self) -> float:
        _initialize_prototypes()

        c_temperature = ctypes.c_double()

        lib().apo_cuda_nose_hoover_integrator_get_instantaneous_temperature(
            ctypes.byref(c_temperature), self.handle
        )

        return float(c_temperature.value)

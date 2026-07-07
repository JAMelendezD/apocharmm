# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

from __future__ import annotations

from pathlib import Path
import struct
import sys

import apocharmm as apo

from python_api_test_helpers import (
    get_repo_root,
    get_data_path,
    require_file,
    assert_file_created,
    remove_if_exists,
    assert_equal,
    assert_close,
    assert_finite_temperature,
    assert_nested_sequence_close,
    expect_exception,
)

BOX_DIMENSIONS: list[float] = [50.0, 50.0, 50.0]
RANDOM_SEED: int = 314159
TEMPERATURE: float = 300.0
TIME_STEP: float = 0.001
THERMOSTAT_FRICTION: float = 0.0
REPORT_FREQUENCY: int = 1
NUM_STEPS: int = 2
TOLERANCE: float = 1.0e-6


class DcdFrame:
    def __init__(
        self, xtal: list[float], x: list[float], y: list[float], z: list[float]
    ) -> None:
        self.xtal: list[float] = xtal
        self.x: list[float] = x
        self.y: list[float] = y
        self.z: list[float] = z
        return


class DcdFile:
    def __init__(self, data: bytes) -> None:
        self.cord: int = 0
        self.num_frames: int = 0
        self.start: int = 0
        self.report_frequency: int = 0
        self.num_steps: int = 0
        self.num_degrees_of_freedom: int = 0
        self.time_step: float = 0.0
        self.has_xtal_data: int = 0
        self.version: int = 0
        self.num_atoms: int = 0
        self.frames: list[DcdFrame] = []

        self._parse(data)

        return

    def _read_int(self, data: bytes, offset: int) -> tuple[int, int]:
        return struct.unpack_from("<i", data, offset)[0], offset + 4

    def _read_float_block(
        self, data: bytes, offset: int, num_atoms: int
    ) -> tuple[list[float], int]:
        block_size, offset = self._read_int(data, offset)
        expected_block_size: int = num_atoms * 4
        assert_equal("DCD coordinate block size", block_size, expected_block_size)

        values: list[float] = list(struct.unpack_from(f"<{num_atoms}f", data, offset))
        offset += block_size

        end_block_size, offset = self._read_int(data, offset)
        assert_equal("DCD coordinate block size", end_block_size, block_size)

        return values, offset

    def _parse(self, data: bytes) -> None:
        offset: int = 0

        control_size, offset = self._read_int(data, offset)
        assert_equal("DCD control block size", control_size, 84)

        control_offset: int = offset
        self.cord = struct.unpack_from("<i", data, control_offset + 0)[0]
        self.num_frames = struct.unpack_from("<i", data, control_offset + 4)[0]
        self.start = struct.unpack_from("<i", data, control_offset + 8)[0]
        self.report_frequency = struct.unpack_from("<i", data, control_offset + 12)[0]
        self.num_steps = struct.unpack_from("<i", data, control_offset + 16)[0]
        self.num_degrees_of_freedom = struct.unpack_from(
            "<i", data, control_offset + 32
        )[0]
        self.time_step = struct.unpack_from("<f", data, control_offset + 40)[0]
        self.has_xtal_data = struct.unpack_from("<i", data, control_offset + 44)[0]
        self.version = struct.unpack_from("<i", data, control_offset + 80)[0]
        offset += control_size

        control_end_size, offset = self._read_int(data, offset)
        assert_equal("DCD control end block size", control_end_size, control_size)

        title_size, offset = self._read_int(data, offset)
        assert_equal("DCD title block size", title_size, 164)
        offset += title_size
        title_end_size, offset = self._read_int(data, offset)
        assert_equal("DCD title end block size", title_end_size, title_size)

        atom_size, offset = self._read_int(data, offset)
        assert_equal("DCD atom block size", atom_size, 4)
        self.num_atoms, offset = self._read_int(data, offset)
        atom_end_size, offset = self._read_int(data, offset)
        assert_equal("DCD atom end block size", atom_end_size, atom_size)

        while offset < len(data):
            xtal_size, offset = self._read_int(data, offset)
            assert_equal("DCD xtal block size", xtal_size, 48)
            xtal: list[float] = list(struct.unpack_from("<6d", data, offset))
            offset += xtal_size
            xtal_end_size, offset = self._read_int(data, offset)
            assert_equal("DCD xtal end block size", xtal_end_size, xtal_size)

            x, offset = self._read_float_block(data, offset, self.num_atoms)
            y, offset = self._read_float_block(data, offset, self.num_atoms)
            z, offset = self._read_float_block(data, offset, self.num_atoms)

            self.frames.append(DcdFrame(xtal, x, y, z))

        return


def read_dcd_file(path: Path) -> DcdFile:
    return DcdFile(path.read_bytes())


def create_context() -> apo.CharmmContext:
    prm_path: str = require_file(get_data_path() / "toppar_water_ions.str")
    psf_path: str = require_file(get_data_path() / "nacl_pair.psf")
    crd_path: str = require_file(get_data_path() / "nacl_pair.cor")

    prm = apo.CharmmParameters(prm_path)
    psf = apo.CharmmPsf(psf_path)
    crd = apo.CharmmCrd(crd_path)

    ctx = apo.CharmmContext(psf, prm)
    ctx.setBoxDimensions(BOX_DIMENSIONS)
    ctx.setCoordinates(crd)
    ctx.useHolonomicConstraints(False)
    ctx.setRandomSeed(RANDOM_SEED)
    ctx.assignVelocitiesAtTemperature(TEMPERATURE)

    return ctx


def create_integrator(ctx: apo.CharmmContext) -> apo.CudaLangevinThermostatIntegrator:
    integrator = apo.CudaLangevinThermostatIntegrator(TIME_STEP)
    integrator.setReferenceTemperature(TEMPERATURE)
    integrator.setThermostatFriction(THERMOSTAT_FRICTION)
    integrator.setThermostatRngSeed(RANDOM_SEED)
    integrator.setCharmmContext(ctx)
    return integrator


def check_constructor_validation(output_dir: Path) -> None:
    print("Checking DcdSubscriber construction and validation...")

    default_path: Path = output_dir / "tmp_python_api_dcd_default.dcd"
    frequency_path: Path = output_dir / "tmp_python_api_dcd_frequency.dcd"

    remove_if_exists(default_path)
    remove_if_exists(frequency_path)

    default_dcd = apo.DcdSubscriber(default_path)
    assert_equal(
        "default DcdSubscriber report frequency", default_dcd.getReportFrequency(), 1000
    )
    default_dcd.close()

    dcd = apo.DcdSubscriber(frequency_path, REPORT_FREQUENCY)
    assert_equal(
        "DcdSubscriber report frequency", dcd.getReportFrequency(), REPORT_FREQUENCY
    )
    dcd.setReportFrequency(2)
    assert_equal("DcdSubscriber updated report frequency", dcd.getReportFrequency(), 2)
    dcd.close()

    expect_exception(
        "DcdSubscriber rejects zero report frequency",
        ValueError,
        lambda: apo.DcdSubscriber(output_dir / "tmp_zero.dcd", 0),
    )
    expect_exception(
        "DcdSubscriber rejects negative report frequency",
        ValueError,
        lambda: apo.DcdSubscriber(output_dir / "tmp_negative.dcd", -1),
    )
    expect_exception(
        "DcdSubscriber rejects too-large report frequency",
        ValueError,
        lambda: apo.DcdSubscriber(output_dir / "tmp_large.dcd", 2**31),
    )
    expect_exception(
        "DcdSubscriber rejects missing output directory",
        apo.ApoCharmmError,
        lambda: apo.DcdSubscriber(output_dir / "missing_dir" / "tmp.dcd"),
    )

    remove_if_exists(default_path)
    remove_if_exists(frequency_path)

    return


def check_dcd_output(dcd_path: Path) -> None:
    print("Checking DcdSubscriber DCD output...")

    remove_if_exists(dcd_path)

    ctx = create_context()
    integrator = create_integrator(ctx)

    dcd = apo.DcdSubscriber(dcd_path, REPORT_FREQUENCY)
    integrator.subscribe(dcd)

    integrator.propagate(NUM_STEPS)
    assert_finite_temperature("post-DCD propagation", ctx.computeTemperature())

    final_coordinates_charges: list[list[float]] = ctx.getCoordinatesCharges()

    integrator.close()
    dcd.close()

    assert_file_created(dcd_path)

    dcd_file = read_dcd_file(dcd_path)

    assert_equal("DCD CORD marker", dcd_file.cord, 1146244931)
    assert_equal("DCD num frames", dcd_file.num_frames, NUM_STEPS)
    assert_equal("DCD start", dcd_file.start, 0)
    assert_equal("DCD report frequency", dcd_file.report_frequency, REPORT_FREQUENCY)
    assert_equal("DCD num steps", dcd_file.num_steps, NUM_STEPS * REPORT_FREQUENCY)
    assert_close("DCD time step", dcd_file.time_step, TIME_STEP, TOLERANCE)
    assert_equal("DCD has xtal data", dcd_file.has_xtal_data, 1)
    assert_equal("DCD version", dcd_file.version, 49)
    assert_equal("DCD num atoms", dcd_file.num_atoms, 2)
    assert_equal("DCD parsed frame count", len(dcd_file.frames), NUM_STEPS)

    last_frame: DcdFrame = dcd_file.frames[-1]

    assert_close("DCD xtal a", last_frame.xtal[0], BOX_DIMENSIONS[0], TOLERANCE)
    assert_close("DCD xtal b", last_frame.xtal[2], BOX_DIMENSIONS[1], TOLERANCE)
    assert_close("DCD xtal c", last_frame.xtal[5], BOX_DIMENSIONS[2], TOLERANCE)

    final_xyz: list[list[float]] = [
        [coord[0], coord[1], coord[2]] for coord in final_coordinates_charges
    ]
    dcd_xyz: list[list[float]] = [
        [last_frame.x[i], last_frame.y[i], last_frame.z[i]]
        for i in range(dcd_file.num_atoms)
    ]

    assert_nested_sequence_close(
        f"DCD final frame coordinate", dcd_xyz, final_xyz, TOLERANCE
    )

    return


def check_close_invalidates_handle(output_dir: Path) -> None:
    print("Checking DcdSubscriber close invalidates handle...")

    dcd_path: Path = output_dir / "tmp_python_api_dcd_closed.dcd"
    remove_if_exists(dcd_path)

    dcd = apo.DcdSubscriber(dcd_path, REPORT_FREQUENCY)
    dcd.close()

    expect_exception(
        "DcdSubscriber.getReportFrequency rejects closed handle",
        RuntimeError,
        lambda: dcd.getReportFrequency(),
    )

    remove_if_exists(dcd_path)

    return


def main(argc: int, argv: list[str]) -> int:
    output_dir: Path = get_repo_root() / "test/pytest"
    output_dir.mkdir(parents=True, exist_ok=True)

    dcd_path: Path = output_dir / "tmp_python_api_dcd_subscriber.dcd"

    try:
        check_constructor_validation(output_dir)
        check_dcd_output(dcd_path)
        check_close_invalidates_handle(output_dir)
    finally:
        print("Cleaning up DcdSubscriber Python API test files...")
        remove_if_exists(dcd_path)

    print("\033[32m" + "PASS: DcdSubscriber Python API tests completed." + "\033[0m")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))

# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

from __future__ import annotations

import math
from pathlib import Path
import sys

import apocharmm as apo

from python_api_test_helpers import (
    require_file,
    assert_equal,
    assert_close,
    assert_sequence_close,
    expect_exception,
)

TOLERANCE: float = 1.0e-7
BOX_DIMENSIONS: tuple[float, float, float] = (20.0, 21.0, 22.0)
FFT_GRID: tuple[int, int, int] = (32, 34, 36)


class TestForce:
    def __init__(self) -> None:
        self.subscribe_calls: list[tuple[apo.ForceManager, str | None]] = []
        self.unsubscribe_calls: list[apo.ForceManager] = []
        return

    def _subscribe_to_force_manager(
        self, force_manager: apo.ForceManager, force_tag: str | None = None
    ) -> None:
        self.subscribe_calls.append((force_manager, force_tag))
        return

    def _unsubscribe_from_force_manager(self, force_manager: apo.ForceManager) -> None:
        self.unsubscribe_calls.append(force_manager)
        return


def create_force_manager(repo_root: Path) -> apo.ForceManager:
    prm_path: str = require_file(repo_root / "test/data/toppar_water_ions.str")
    psf_path: str = require_file(repo_root / "test/data/nacl_pair.psf")

    prm = apo.CharmmParameters(prm_path)
    psf = apo.CharmmPsf(psf_path)

    return apo.ForceManager(psf, prm)


def check_default_state(repo_root: Path) -> None:
    print("Checking ForceManager default state...")

    fm = create_force_manager(repo_root)

    assert_equal("ForceManager.getNumAtoms", fm.getNumAtoms(), 2)
    assert_equal("ForceManager.isInitialized", fm.isInitialized(), False)
    assert_equal("ForceManager.isComposite", fm.isComposite(), False)

    assert_sequence_close(
        "ForceManager.getBoxDimensions default",
        fm.getBoxDimensions(),
        (-9999.9999, -9999.9999, -9999.9999),
        TOLERANCE,
    )

    assert_close("ForceManager.getKappa default", fm.getKappa(), 0.34, TOLERANCE)
    assert_close("ForceManager.getCutoff default", fm.getCutoff(), 14.0, TOLERANCE)
    assert_close("ForceManager.getCtonnb default", fm.getCtonnb(), 12.0, TOLERANCE)
    assert_close("ForceManager.getCtofnb default", fm.getCtofnb(), 10.0, TOLERANCE)
    assert_equal("ForceManager.getFFTGrid default", fm.getFFTGrid(), (-1, -1, -1))
    assert_equal("ForceManager.getPmeSplineOrder default", fm.getPmeSplineOrder(), 4)
    assert_equal(
        "ForceManager.getPeriodicBoundaryCondition default",
        fm.getPeriodicBoundaryCondition(),
        apo.PeriodicBoundaryCondition.P1,
    )
    assert_equal("ForceManager.getVdwType default", fm.getVdwType(), apo.VdwType.VFSW)

    return


def check_setters_and_getters(repo_root: Path) -> None:
    print("Checking ForceManager setter/getter round trips...")

    fm = create_force_manager(repo_root)

    fm.setBoxDimensions(BOX_DIMENSIONS)
    fm.setKappa(0.45)
    fm.setCutoff(9.0)
    fm.setCtonnb(8.5)
    fm.setCtofnb(7.5)
    fm.setFFTGrid(FFT_GRID)
    fm.setPmeSplineOrder(6)
    fm.setPeriodicBoundaryCondition(apo.PeriodicBoundaryCondition.P21)
    fm.setVdwType(apo.VdwType.DBEXP)
    fm.setPrintEnergyDecomposition(True)
    fm.setPrintEnergyDecomposition(False)

    assert_sequence_close(
        "ForceManager.getBoxDimensions",
        fm.getBoxDimensions(),
        BOX_DIMENSIONS,
        TOLERANCE,
    )
    assert_close("ForceManager.getKappa", fm.getKappa(), 0.45, TOLERANCE)
    assert_close("ForceManager.getCutoff", fm.getCutoff(), 9.0, TOLERANCE)
    assert_close("ForceManager.getCtonnb", fm.getCtonnb(), 8.5, TOLERANCE)
    assert_close("ForceManager.getCtofnb", fm.getCtofnb(), 7.5, TOLERANCE)
    assert_equal("ForceManager.getFFTGRid", fm.getFFTGrid(), FFT_GRID)
    assert_equal("ForceManager.getPmeSplineOrder", fm.getPmeSplineOrder(), 6)
    assert_equal(
        "ForceManager.getPeriodicBoundaryCondition",
        fm.getPeriodicBoundaryCondition(),
        apo.PeriodicBoundaryCondition.P21,
    )
    assert_equal("ForceManager.getVdwType", fm.getVdwType(), apo.VdwType.DBEXP)

    for vdw_type in (
        apo.VdwType.VSH,
        apo.VdwType.VSW,
        apo.VdwType.VFSW,
        apo.VdwType.VGSH,
        apo.VdwType.CUT,
        apo.VdwType.DBEXP,
    ):
        fm.setVdwType(vdw_type)
        assert_equal(f"ForceManager.getVdwType {vdw_type}", fm.getVdwType(), vdw_type)

    return


def check_validation(repo_root: Path) -> None:
    print("Checking ForceManager validation...")

    prm_path: str = require_file(repo_root / "test/data/toppar_water_ions.str")
    psf_path: str = require_file(repo_root / "test/data/nacl_pair.psf")

    prm = apo.CharmmParameters(prm_path)
    psf = apo.CharmmPsf(psf_path)

    fm = apo.ForceManager(psf, prm)

    expect_exception(
        "ForceManager rejects non-CharmmPsf psf",
        TypeError,
        lambda: apo.ForceManager(object(), prm),
    )
    expect_exception(
        "ForceManager rejects non-CharmmParameters parameters",
        TypeError,
        lambda: apo.ForceManager(psf, object()),
    )

    expect_exception(
        "ForceManager.setBoxDimensions rejects short input",
        apo.ApoCharmmError,
        lambda: fm.setBoxDimensions([1.0, 2.0]),
    )
    expect_exception(
        "ForceManager.setBoxDimensions rejects long input",
        apo.ApoCharmmError,
        lambda: fm.setBoxDimensions([1.0, 2.0, 3.0, 4.0]),
    )
    expect_exception(
        "ForceManager.setBoxDimensions rejects negative dimensions",
        apo.ApoCharmmError,
        lambda: fm.setBoxDimensions([20.0, -1.0, 20.0]),
    )
    expect_exception(
        "ForceManager.setBoxDimensions rejects zero dimensions",
        apo.ApoCharmmError,
        lambda: fm.setBoxDimensions([20.0, 0.0, 20.0]),
    )
    expect_exception(
        "ForceManager.setBoxDimensions rejects NaN",
        apo.ApoCharmmError,
        lambda: fm.setBoxDimensions([20.0, math.nan, 20.0]),
    )
    expect_exception(
        "ForceManager.setKappa rejects NaN",
        apo.ApoCharmmError,
        lambda: fm.setKappa(math.nan),
    )
    expect_exception(
        "ForceManager.setKappa rejects negative values",
        apo.ApoCharmmError,
        lambda: fm.setKappa(-1.0),
    )
    expect_exception(
        "ForceManager.setCutoff rejects zero",
        apo.ApoCharmmError,
        lambda: fm.setCutoff(0.0),
    )
    expect_exception(
        "ForceManager.setCtonnb rejects zero",
        apo.ApoCharmmError,
        lambda: fm.setCtonnb(0.0),
    )
    expect_exception(
        "ForceManager.setCtofnb rejects zero",
        apo.ApoCharmmError,
        lambda: fm.setCtofnb(0.0),
    )
    expect_exception(
        "ForceManager.setFFTGrid rejects short input",
        apo.ApoCharmmError,
        lambda: fm.setFFTGrid([32, 34]),
    )
    expect_exception(
        "ForceManager.setPmeSplineOrder rejects zero",
        apo.ApoCharmmError,
        lambda: fm.setPmeSplineOrder(0),
    )
    expect_exception(
        "ForceManager.setPeriodicBoundaryCondition rejects invalid enum value",
        ValueError,
        lambda: fm.setPeriodicBoundaryCondition(99),
    )
    expect_exception(
        "ForceManager.setVdwType rejects invalid enum value",
        ValueError,
        lambda: fm.setVdwType(99),
    )
    expect_exception(
        "ForceManager.setVdwType rejects VdwType.NONE through C ABI",
        apo.ApoCharmmError,
        lambda: fm.setVdwType(apo.VdwType.NONE),
    )

    return


def check_subscription_wrapper(repo_root: Path) -> None:
    print("Checking ForceManager Python subscription bookkeeping...")

    fm = create_force_manager(repo_root)

    force = TestForce()
    fm.subscribe(force)

    assert_equal("ForceManager.subscribe call count", len(force.subscribe_calls), 1)
    assert_equal(
        "ForceManager.subscribe manager identity",
        force.subscribe_calls[0][0] is fm,
        True,
    )
    assert_equal(
        "ForceManager.subscribe default tag", force.subscribe_calls[0][1], None
    )

    expect_exception(
        "ForceManager.subscribe rejects duplicate force",
        ValueError,
        lambda: fm.subscribe(force),
    )

    fm.unsubscribe(force)

    assert_equal("ForceManager.unsubscribe call count", len(force.unsubscribe_calls), 1)
    assert_equal(
        "ForceManager.unsubscribe manager identity",
        force.unsubscribe_calls[0] is fm,
        True,
    )

    tagged_force = TestForce()
    fm.subscribe(tagged_force, "custom_tag")

    assert_equal(
        "ForceManager.subscribe custom tag",
        tagged_force.subscribe_calls[0][1],
        "custom_tag",
    )
    fm.unsubscribe(tagged_force)

    expect_exception(
        "ForceManager.subscribe rejects non-subscribable object",
        TypeError,
        lambda: fm.subscribe(object()),
    )
    expect_exception(
        "ForceManager.subscribe rejects non-string force_tag",
        TypeError,
        lambda: fm.subscribe(TestForce(), 1),
    )
    expect_exception(
        "ForceManager.subscribe rejects empty force_tag",
        ValueError,
        lambda: fm.subscribe(TestForce(), ""),
    )
    expect_exception(
        "ForceManager.unsubscribe rejects non-subscribable object",
        TypeError,
        lambda: fm.unsubscribe(object()),
    )
    expect_exception(
        "ForceManager.unsubscribe rejects unregistered force",
        ValueError,
        lambda: fm.unsubscribe(TestForce()),
    )

    return


def main(argc: int, argv: list[str]) -> int:
    repo_root: Path = Path(argv[1]) if argc > 1 else Path(".")

    check_default_state(repo_root)
    check_setters_and_getters(repo_root)
    check_validation(repo_root)
    check_subscription_wrapper(repo_root)

    print("\033[32m" + "PASS: ForceManager Python API tests completed." + "\033[0m")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))

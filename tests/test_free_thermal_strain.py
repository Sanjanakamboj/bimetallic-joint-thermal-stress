"""Tests H-L: free (unrestrained) thermal strain."""

from __future__ import annotations

import math

import pytest

from thermal_joint import AxialMember, ThermoelasticMaterial

ALPHA = 23.0e-6
MATERIAL = ThermoelasticMaterial("M", 70.0e9, ALPHA, 270.0e6)
MEMBER = AxialMember(MATERIAL, 100.0e-6)


def test_h_heating_hand_calculation() -> None:
    """H. alpha = 23e-6 /K over dT = +100 K gives eps = 2.3e-3, positive."""
    strain = MATERIAL.free_thermal_strain(100.0)
    assert strain == pytest.approx(2.3e-3, rel=1e-15)
    assert strain > 0.0


def test_i_cooling_hand_calculation() -> None:
    """I. The same material over dT = -140 K gives eps = -3.22e-3, negative."""
    strain = MATERIAL.free_thermal_strain(-140.0)
    assert strain == pytest.approx(-3.22e-3, rel=1e-15)
    assert strain < 0.0


def test_j_zero_delta_temperature_gives_zero_strain() -> None:
    """J. dT = 0 gives exactly zero free thermal strain."""
    assert MATERIAL.free_thermal_strain(0.0) == 0.0


@pytest.mark.parametrize("factor", [0.5, 2.0, 3.0, -1.0])
def test_k_linear_in_alpha(factor: float) -> None:
    """K. Free thermal strain scales linearly with alpha."""
    scaled = ThermoelasticMaterial("M", 70.0e9, ALPHA * factor, 270.0e6)
    assert scaled.free_thermal_strain(80.0) == pytest.approx(
        factor * MATERIAL.free_thermal_strain(80.0), rel=1e-15
    )


@pytest.mark.parametrize("factor", [0.5, 2.0, 3.0, -1.0])
def test_l_linear_in_delta_temperature(factor: float) -> None:
    """L. Free thermal strain scales linearly with dT."""
    base = MATERIAL.free_thermal_strain(80.0)
    assert MATERIAL.free_thermal_strain(80.0 * factor) == pytest.approx(
        factor * base, rel=1e-15
    )


def test_member_delegates_free_thermal_strain_to_its_material() -> None:
    """A member reports the same free thermal strain as its material."""
    assert MEMBER.free_thermal_strain(-60.0) == MATERIAL.free_thermal_strain(-60.0)


@pytest.mark.parametrize("bad_dt", [math.nan, math.inf, -math.inf])
def test_non_finite_delta_temperature_rejected(bad_dt: float) -> None:
    """Free thermal strain rejects non-finite temperature changes."""
    with pytest.raises(ValueError):
        MATERIAL.free_thermal_strain(bad_dt)

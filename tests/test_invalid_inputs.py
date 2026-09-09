"""Section 37: a systematic sweep of rejected inputs across the public API."""

from __future__ import annotations

import math

import pytest

from thermal_joint import (
    AxialMember,
    ThermalEnvironment,
    ThermoelasticMaterial,
    YieldBasis,
    common_joint_strain,
    solve_bimetallic_joint,
)

NON_FINITE = [math.nan, math.inf, -math.inf]
NON_POSITIVE = [0.0, -1.0, -1.0e9]

GOOD_MATERIAL = ThermoelasticMaterial("good", 70.0e9, 23.0e-6, 270.0e6)
GOOD_MEMBER_1 = AxialMember(GOOD_MATERIAL, 100.0e-6)
GOOD_MEMBER_2 = AxialMember(
    ThermoelasticMaterial("good-2", 110.0e9, 8.5e-6, 830.0e6), 100.0e-6
)


@pytest.mark.parametrize("bad", NON_FINITE + NON_POSITIVE)
def test_modulus_must_be_finite_and_positive(bad: float) -> None:
    with pytest.raises(ValueError):
        ThermoelasticMaterial("M", bad, 23.0e-6, 270.0e6)


@pytest.mark.parametrize("bad", NON_FINITE)
def test_cte_must_be_finite(bad: float) -> None:
    with pytest.raises(ValueError):
        ThermoelasticMaterial("M", 70.0e9, bad, 270.0e6)


@pytest.mark.parametrize("bad", NON_FINITE + NON_POSITIVE)
def test_yield_strength_must_be_finite_and_positive(bad: float) -> None:
    with pytest.raises(ValueError):
        ThermoelasticMaterial("M", 70.0e9, 23.0e-6, bad)


@pytest.mark.parametrize("bad", NON_FINITE + NON_POSITIVE)
def test_density_must_be_finite_and_positive(bad: float) -> None:
    with pytest.raises(ValueError):
        ThermoelasticMaterial("M", 70.0e9, 23.0e-6, 270.0e6, density=bad)


@pytest.mark.parametrize("bad", NON_FINITE + NON_POSITIVE)
def test_area_must_be_finite_and_positive(bad: float) -> None:
    with pytest.raises(ValueError):
        AxialMember(GOOD_MATERIAL, bad)


@pytest.mark.parametrize("bad", NON_FINITE)
def test_delta_temperature_must_be_finite(bad: float) -> None:
    with pytest.raises(ValueError):
        GOOD_MATERIAL.free_thermal_strain(bad)
    with pytest.raises(ValueError):
        common_joint_strain(GOOD_MEMBER_1, GOOD_MEMBER_2, bad)
    with pytest.raises(ValueError):
        solve_bimetallic_joint(GOOD_MEMBER_1, GOOD_MEMBER_2, bad)


@pytest.mark.parametrize("bad", NON_FINITE + [0.0, 0.5, 0.999, -2.0])
def test_design_factor_must_be_finite_and_at_least_one(bad: float) -> None:
    with pytest.raises(ValueError):
        YieldBasis(bad)


@pytest.mark.parametrize("bad", NON_FINITE)
def test_margin_rejects_non_finite_stress(bad: float) -> None:
    with pytest.raises(ValueError):
        YieldBasis().margin_of_safety(GOOD_MATERIAL, bad)


@pytest.mark.parametrize("bad", NON_FINITE)
def test_environment_rejects_non_finite_temperatures(bad: float) -> None:
    with pytest.raises(ValueError):
        ThermalEnvironment(bad, -120.0, 120.0)
    with pytest.raises(ValueError):
        ThermalEnvironment(20.0, bad, 120.0)
    with pytest.raises(ValueError):
        ThermalEnvironment(20.0, -120.0, bad)


@pytest.mark.parametrize(
    "reference, cold, hot",
    [(20.0, 30.0, 120.0), (20.0, -120.0, 5.0), (20.0, 120.0, -120.0)],
)
def test_environment_rejects_malformed_ordering(reference: float, cold: float, hot: float) -> None:
    with pytest.raises(ValueError):
        ThermalEnvironment(reference, cold, hot)


@pytest.mark.parametrize("bad_name", ["", "   ", "\t\n"])
def test_names_must_be_non_empty(bad_name: str) -> None:
    with pytest.raises(ValueError):
        ThermoelasticMaterial(bad_name, 70.0e9, 23.0e-6, 270.0e6)
    with pytest.raises(ValueError):
        AxialMember(GOOD_MATERIAL, 100.0e-6, label=bad_name)


@pytest.mark.parametrize("bad_name", [None, 3, 4.5, []])
def test_names_must_be_strings(bad_name: object) -> None:
    with pytest.raises(TypeError):
        ThermoelasticMaterial(bad_name, 70.0e9, 23.0e-6, 270.0e6)  # type: ignore[arg-type]


def test_wrong_types_are_rejected_across_the_api() -> None:
    with pytest.raises(TypeError):
        AxialMember(object(), 100.0e-6)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        YieldBasis().allowable_stress(GOOD_MEMBER_1)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        solve_bimetallic_joint(GOOD_MATERIAL, GOOD_MEMBER_2, 100.0)  # type: ignore[arg-type]

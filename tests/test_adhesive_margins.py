"""Tests AI-AO: adhesive shear allowable, margins and governing extreme."""

from __future__ import annotations

import math

import pytest
from shear_lag_cases import HAND_ADHESIVE, HAND_GEOMETRY, HAND_MEMBER_1, HAND_MEMBER_2

from thermal_joint import (
    AdhesiveMaterial,
    AdhesiveShearBasis,
    AxialMember,
    ThermalEnvironment,
    ThermoelasticMaterial,
    assess_shear_lag_extremes,
    solve_shear_lag,
)
from thermal_joint.illustrative import (
    ADHESIVE_LIKE,
    aluminium_like_member,
    radiator_joint_environment,
    radiator_joint_overlap,
    titanium_like_member,
)

MEMBER_AL = aluminium_like_member(100.0)
MEMBER_TI = titanium_like_member(100.0)
ENVIRONMENT = radiator_joint_environment()
GEOMETRY = radiator_joint_overlap()
BASIS = AdhesiveShearBasis(1.25)


def test_ai_allowable_hand_calculation() -> None:
    """AI. tau_allow = 25 MPa / 1.25 = 20 MPa."""
    assert BASIS.allowable_shear_stress(ADHESIVE_LIKE) == pytest.approx(20.0e6, rel=1e-15)
    assert AdhesiveShearBasis().allowable_shear_stress(ADHESIVE_LIKE) == (
        ADHESIVE_LIKE.shear_strength
    )


@pytest.mark.parametrize("factor", [1.0, 1.25, 1.5, 2.0])
def test_aj_valid_design_factors_accepted(factor: float) -> None:
    """AJ. Any finite factor >= 1 is a valid basis."""
    assert AdhesiveShearBasis(factor).design_factor == factor


@pytest.mark.parametrize("bad", [0.0, 0.5, 0.999, -1.0, math.nan, math.inf])
def test_aj_invalid_design_factors_rejected(bad: float) -> None:
    """AJ. A factor below 1 would raise the allowable above the strength."""
    with pytest.raises(ValueError):
        AdhesiveShearBasis(bad)


@pytest.mark.parametrize("shear", [-40.0e6, 40.0e6])
def test_ak_margin_hand_calculation(shear: float) -> None:
    """AK. MS = (25/1.25)/40 - 1 = -0.5, for either sign of shear."""
    assert BASIS.margin_of_safety(ADHESIVE_LIKE, shear) == pytest.approx(-0.5, rel=1e-14)


def test_ak_margin_uses_absolute_shear() -> None:
    """AK. Sign of the shear does not change the margin."""
    assert BASIS.margin_of_safety(ADHESIVE_LIKE, 10.0e6) == BASIS.margin_of_safety(
        ADHESIVE_LIKE, -10.0e6
    )


def test_al_exact_boundary_passes() -> None:
    """AL. A demand exactly at the allowable gives MS = 0 exactly, and passes."""
    basis = AdhesiveShearBasis()
    adhesive = AdhesiveMaterial("boundary", 1.0e9, 33.0e6, "illustrative")
    margin = basis.margin_of_safety(adhesive, -33.0e6)
    assert margin == 0.0
    assert margin >= 0.0


def test_al_near_boundary_passes_through_a_solved_demand() -> None:
    """AL. Tuning the strength to the computed peak gives a zero margin."""
    demand = solve_shear_lag(HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, HAND_ADHESIVE, 50.0)
    tuned = AdhesiveMaterial(
        "tuned", HAND_ADHESIVE.shear_modulus, demand.peak_shear_stress, "illustrative"
    )
    assert AdhesiveShearBasis().margin_of_safety(tuned, demand.peak_shear_stress) == 0.0


def test_am_overstressed_bondline_fails() -> None:
    """AM. The canonical illustrative joint fails the adhesive screen."""
    assessment = assess_shear_lag_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE, BASIS
    )
    assert assessment.minimum_margin < 0.0
    assert not assessment.passes
    assert assessment.cold_margin == pytest.approx(-0.6987, rel=1e-3)
    assert assessment.hot_margin == pytest.approx(-0.5782, rel=1e-3)


def test_an_zero_demand_gives_an_infinite_non_governing_margin() -> None:
    """AN. Zero shear never divides by zero; it yields math.inf."""
    assert math.isinf(BASIS.margin_of_safety(ADHESIVE_LIKE, 0.0))

    alpha = 15.0e-6
    member_1 = AxialMember(ThermoelasticMaterial("A", 70.0e9, alpha, 270.0e6), 100.0e-6)
    member_2 = AxialMember(ThermoelasticMaterial("B", 110.0e9, alpha, 830.0e6), 100.0e-6)
    assessment = assess_shear_lag_extremes(
        member_1, member_2, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE, BASIS
    )
    assert math.isinf(assessment.minimum_margin)
    assert assessment.passes


def test_an_zero_excursion_environment_is_non_governing() -> None:
    """AN. A zero-dT environment leaves the bondline unloaded."""
    assessment = assess_shear_lag_extremes(
        MEMBER_AL,
        MEMBER_TI,
        ThermalEnvironment(20.0, 20.0, 20.0),
        GEOMETRY,
        ADHESIVE_LIKE,
        BASIS,
    )
    assert math.isinf(assessment.minimum_margin)
    assert assessment.passes


def test_ao_governing_extreme_is_computed_from_the_margins() -> None:
    """AO. Cold governs canonically because |dT_cold| is larger - computed, not fixed."""
    assessment = assess_shear_lag_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE, BASIS
    )
    assert assessment.governing_extreme == "cold"
    assert assessment.minimum_margin == min(assessment.cold_margin, assessment.hot_margin)
    assert assessment.governing_demand is assessment.cold_demand
    assert assessment.governing_peak_shear_stress == assessment.cold_demand.peak_shear_stress


def test_ao_hot_governs_when_the_hot_excursion_is_larger() -> None:
    """AO. Reversing the asymmetry makes hot govern."""
    assessment = assess_shear_lag_extremes(
        MEMBER_AL,
        MEMBER_TI,
        ThermalEnvironment(20.0, -30.0, 320.0),
        GEOMETRY,
        ADHESIVE_LIKE,
        BASIS,
    )
    assert assessment.governing_extreme == "hot"
    assert assessment.minimum_margin == assessment.hot_margin


def test_ao_symmetric_excursions_tie_and_resolve_to_cold() -> None:
    """AO. Equal |dT| gives equal peaks and margins; the tie resolves to cold."""
    assessment = assess_shear_lag_extremes(
        MEMBER_AL,
        MEMBER_TI,
        ThermalEnvironment(20.0, -80.0, 120.0),
        GEOMETRY,
        ADHESIVE_LIKE,
        BASIS,
    )
    assert assessment.cold_demand.peak_shear_stress == pytest.approx(
        assessment.hot_demand.peak_shear_stress, rel=1e-13
    )
    assert assessment.cold_margin == pytest.approx(assessment.hot_margin, rel=1e-13)
    assert assessment.governing_extreme == "cold"


def test_assessment_is_deterministic_and_type_checked() -> None:
    """Repeat runs are identical; bad arguments are rejected."""
    first = assess_shear_lag_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE, BASIS
    )
    second = assess_shear_lag_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE, BASIS
    )
    assert first == second
    assert first.allowable_shear_stress == pytest.approx(20.0e6, rel=1e-15)
    with pytest.raises(TypeError):
        assess_shear_lag_extremes(
            MEMBER_AL, MEMBER_TI, "nope", GEOMETRY, ADHESIVE_LIKE, BASIS
        )  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        assess_shear_lag_extremes(
            MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE, "nope"
        )  # type: ignore[arg-type]


def test_assessment_defaults_to_a_unity_basis() -> None:
    """Omitting the basis is the same as AdhesiveShearBasis()."""
    assert assess_shear_lag_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE
    ) == assess_shear_lag_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE, AdhesiveShearBasis(1.0)
    )

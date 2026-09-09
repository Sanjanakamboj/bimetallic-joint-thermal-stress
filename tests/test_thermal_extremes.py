"""Tests AL-AR: thermal environment validation and hot/cold extreme assessment."""

from __future__ import annotations

import math

import pytest

from thermal_joint import (
    AxialMember,
    ThermalEnvironment,
    ThermoelasticMaterial,
    YieldBasis,
    assess_temperature_extremes,
)
from thermal_joint.illustrative import aluminium_like_member, titanium_like_member

MEMBER_1 = aluminium_like_member(100.0)
MEMBER_2 = titanium_like_member(100.0)


def test_al_valid_environment_is_accepted() -> None:
    """AL. cold <= reference <= hot is the admissible ordering."""
    environment = ThermalEnvironment(20.0, -120.0, 120.0)
    assert environment.reference_temperature == 20.0
    assert environment.cold_temperature == -120.0
    assert environment.hot_temperature == 120.0


def test_al_equal_temperatures_are_allowed() -> None:
    """AL. Degenerate environments with zero excursions are permitted."""
    environment = ThermalEnvironment(20.0, 20.0, 20.0)
    assert environment.delta_temperature_cold == 0.0
    assert environment.delta_temperature_hot == 0.0


@pytest.mark.parametrize(
    "reference, cold, hot",
    [
        (20.0, 50.0, 120.0),   # cold above reference
        (20.0, -120.0, 10.0),  # hot below reference
        (20.0, 120.0, -120.0), # cold and hot swapped
        (200.0, -120.0, 120.0),# reference above hot
    ],
)
def test_al_malformed_ordering_rejected(reference: float, cold: float, hot: float) -> None:
    """AL. Any violation of cold <= reference <= hot is rejected."""
    with pytest.raises(ValueError):
        ThermalEnvironment(reference, cold, hot)


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_al_non_finite_temperatures_rejected(bad: float) -> None:
    """AL. Every temperature must be finite."""
    with pytest.raises(ValueError):
        ThermalEnvironment(bad, -120.0, 120.0)
    with pytest.raises(ValueError):
        ThermalEnvironment(20.0, bad, 120.0)
    with pytest.raises(ValueError):
        ThermalEnvironment(20.0, -120.0, bad)


def test_am_delta_temperature_calculation() -> None:
    """AM. dT_cold = -140 K and dT_hot = +100 K for the illustrative excursion."""
    environment = ThermalEnvironment(20.0, -120.0, 120.0)
    assert environment.delta_temperature_cold == pytest.approx(-140.0, rel=1e-15)
    assert environment.delta_temperature_hot == pytest.approx(100.0, rel=1e-15)


def test_am_celsius_and_kelvin_give_the_same_deltas() -> None:
    """AM. Only differences matter, so a 273.15 K offset changes nothing."""
    celsius = ThermalEnvironment(20.0, -120.0, 120.0)
    kelvin = ThermalEnvironment(20.0 + 273.15, -120.0 + 273.15, 120.0 + 273.15)
    assert kelvin.delta_temperature_cold == pytest.approx(
        celsius.delta_temperature_cold, rel=1e-12
    )
    assert kelvin.delta_temperature_hot == pytest.approx(
        celsius.delta_temperature_hot, rel=1e-12
    )


def test_an_symmetric_excursions_give_equal_margins() -> None:
    """AN. |dT_hot| = |dT_cold| gives identical stress magnitudes and margins."""
    environment = ThermalEnvironment(20.0, -80.0, 120.0)
    assessment = assess_temperature_extremes(MEMBER_1, MEMBER_2, environment, YieldBasis(1.25))

    assert abs(environment.delta_temperature_cold) == abs(environment.delta_temperature_hot)
    assert abs(assessment.cold_result.member_1_stress) == pytest.approx(
        abs(assessment.hot_result.member_1_stress), rel=1e-12
    )
    assert assessment.cold_assessment.minimum_margin == pytest.approx(
        assessment.hot_assessment.minimum_margin, rel=1e-12
    )
    assert assessment.cold_result.member_1_stress * assessment.hot_result.member_1_stress < 0.0


def test_ao_larger_cold_excursion_governs() -> None:
    """AO. With |dT_cold| > |dT_hot| the cold case governs."""
    environment = ThermalEnvironment(20.0, -120.0, 120.0)
    assessment = assess_temperature_extremes(MEMBER_1, MEMBER_2, environment, YieldBasis(1.25))
    assert abs(environment.delta_temperature_cold) > abs(environment.delta_temperature_hot)
    assert assessment.governing_extreme == "cold"
    assert assessment.minimum_yield_margin == assessment.cold_assessment.minimum_margin


def test_ao_larger_hot_excursion_governs() -> None:
    """AO. Reversing the asymmetry makes the hot case govern instead."""
    environment = ThermalEnvironment(20.0, -30.0, 220.0)
    assessment = assess_temperature_extremes(MEMBER_1, MEMBER_2, environment, YieldBasis(1.25))
    assert abs(environment.delta_temperature_hot) > abs(environment.delta_temperature_cold)
    assert assessment.governing_extreme == "hot"
    assert assessment.minimum_yield_margin == assessment.hot_assessment.minimum_margin


def test_ap_hot_and_cold_stress_signs_reverse() -> None:
    """AP. The high-CTE member is compressive when hot and tensile when cold."""
    environment = ThermalEnvironment(20.0, -120.0, 120.0)
    assessment = assess_temperature_extremes(MEMBER_1, MEMBER_2, environment)

    assert assessment.hot_result.member_1_stress < 0.0
    assert assessment.hot_result.member_2_stress > 0.0
    assert assessment.cold_result.member_1_stress > 0.0
    assert assessment.cold_result.member_2_stress < 0.0


def test_aq_governing_extreme_is_computed_not_hard_coded() -> None:
    """AQ. The same members give either extreme as governing, depending on dT."""
    cold_biased = assess_temperature_extremes(
        MEMBER_1, MEMBER_2, ThermalEnvironment(20.0, -200.0, 40.0)
    )
    hot_biased = assess_temperature_extremes(
        MEMBER_1, MEMBER_2, ThermalEnvironment(20.0, 0.0, 300.0)
    )
    assert cold_biased.governing_extreme == "cold"
    assert hot_biased.governing_extreme == "hot"
    assert cold_biased.minimum_yield_margin != hot_biased.minimum_yield_margin


def test_aq_governing_extreme_matches_the_smaller_minimum_margin() -> None:
    """AQ. The reported extreme is always the one with the smaller margin."""
    for environment in (
        ThermalEnvironment(20.0, -120.0, 120.0),
        ThermalEnvironment(20.0, -30.0, 220.0),
        ThermalEnvironment(20.0, -80.0, 120.0),
    ):
        assessment = assess_temperature_extremes(MEMBER_1, MEMBER_2, environment)
        expected = min(
            assessment.cold_assessment.minimum_margin,
            assessment.hot_assessment.minimum_margin,
        )
        assert assessment.minimum_yield_margin == expected
        assert assessment.governing_assessment.minimum_margin == expected
        assert (
            assessment.governing_member
            == assessment.governing_assessment.governing_member
        )


def test_aq_governing_extreme_ties_resolve_to_cold() -> None:
    """AQ. Symmetric excursions tie; the tie resolves deterministically to cold."""
    assessment = assess_temperature_extremes(
        MEMBER_1, MEMBER_2, ThermalEnvironment(20.0, -80.0, 120.0)
    )
    assert assessment.governing_extreme == "cold"


def test_ar_repeated_assessment_is_deterministic() -> None:
    """AR. Repeating the assessment reproduces every field exactly."""
    environment = ThermalEnvironment(20.0, -120.0, 120.0)
    basis = YieldBasis(1.25)
    first = assess_temperature_extremes(MEMBER_1, MEMBER_2, environment, basis)
    second = assess_temperature_extremes(MEMBER_1, MEMBER_2, environment, basis)
    assert first == second
    assert first.cold_result == second.cold_result
    assert first.minimum_yield_margin == second.minimum_yield_margin


def test_extreme_assessment_defaults_and_type_checking() -> None:
    """The basis defaults to unity and a non-environment argument is rejected."""
    environment = ThermalEnvironment(20.0, -120.0, 120.0)
    assert assess_temperature_extremes(MEMBER_1, MEMBER_2, environment) == (
        assess_temperature_extremes(MEMBER_1, MEMBER_2, environment, YieldBasis(1.0))
    )
    with pytest.raises(TypeError):
        assess_temperature_extremes(MEMBER_1, MEMBER_2, "not-an-environment")  # type: ignore[arg-type]


def test_zero_excursion_environment_is_non_governing() -> None:
    """A zero-excursion environment produces zero stress and infinite margins."""
    assessment = assess_temperature_extremes(
        MEMBER_1, MEMBER_2, ThermalEnvironment(20.0, 20.0, 20.0)
    )
    assert assessment.cold_result.member_1_stress == 0.0
    assert assessment.hot_result.member_2_stress == 0.0
    assert assessment.minimum_yield_margin == math.inf
    assert assessment.passes


def test_governing_member_can_differ_from_the_lower_yield_material() -> None:
    """The governing member at an extreme is computed, not taken from yield alone."""
    low_yield_large_area = AxialMember(
        ThermoelasticMaterial("low yield", 70.0e9, 23.0e-6, 100.0e6), 1000.0e-6
    )
    high_yield_small_area = AxialMember(
        ThermoelasticMaterial("high yield", 110.0e9, 8.5e-6, 500.0e6), 10.0e-6
    )
    assessment = assess_temperature_extremes(
        low_yield_large_area,
        high_yield_small_area,
        ThermalEnvironment(20.0, -120.0, 120.0),
    )
    assert assessment.governing_member == 2

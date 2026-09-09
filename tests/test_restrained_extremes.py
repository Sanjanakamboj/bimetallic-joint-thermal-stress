"""Sections 18-21: yield margins and governing selection under external restraint."""

from __future__ import annotations

import math

import pytest

from thermal_joint import (
    AxialMember,
    AxialRestraint,
    ThermalEnvironment,
    ThermoelasticMaterial,
    YieldBasis,
    assess_restrained_temperature_extremes,
    assess_temperature_extremes,
    assess_yield,
    solve_restrained_joint,
)
from thermal_joint.illustrative import (
    aluminium_like_member,
    radiator_joint_environment,
    titanium_like_member,
)

MEMBER_AL = aluminium_like_member(100.0)
MEMBER_TI = titanium_like_member(100.0)
ENVIRONMENT = radiator_joint_environment()
BASIS = YieldBasis(1.25)
RATIOS = (0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0)


def assess_at(ratio: float, environment: ThermalEnvironment = ENVIRONMENT, basis: YieldBasis = BASIS):
    restraint = AxialRestraint.from_stiffness_ratio(ratio, MEMBER_AL, MEMBER_TI)
    return assess_restrained_temperature_extremes(
        MEMBER_AL, MEMBER_TI, environment, restraint, basis
    )


def test_yield_semantics_are_reused_not_reimplemented() -> None:
    """Section 18. The restrained assessment is the Milestone 1 assess_yield."""
    restraint = AxialRestraint.from_stiffness_ratio(1.0, MEMBER_AL, MEMBER_TI)
    result = solve_restrained_joint(MEMBER_AL, MEMBER_TI, -140.0, restraint)
    expected = assess_yield(result, MEMBER_AL, MEMBER_TI, BASIS)
    assert assess_at(1.0).cold_assessment == expected


def test_margins_use_absolute_stress_exactly_as_milestone_1() -> None:
    """Section 18. Tension and compression of equal magnitude give equal margins."""
    assessment = assess_at(1.0)
    margin = assessment.cold_assessment.member_1_margin
    assert margin.margin_of_safety == pytest.approx(
        BASIS.allowable_stress(MEMBER_AL.material) / abs(margin.stress) - 1.0, rel=1e-14
    )


def test_zero_restraint_reproduces_the_milestone_1_extreme_assessment() -> None:
    """Section 19. eta_r = 0 gives the Milestone 1 governing result exactly."""
    free = assess_temperature_extremes(MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS)
    restrained = assess_at(0.0)
    assert restrained.governing_extreme == free.governing_extreme
    assert restrained.governing_member == free.governing_member
    assert restrained.minimum_yield_margin == pytest.approx(
        free.minimum_yield_margin, rel=1e-12
    )
    assert restrained.minimum_yield_margin == pytest.approx(1.4874, rel=1e-3)


def test_restraint_reduces_the_minimum_margin_for_this_joint() -> None:
    """Section 19. Adding restraint eats margin, monotonically here."""
    margins = [assess_at(r).minimum_yield_margin for r in RATIOS]
    assert margins == sorted(margins, reverse=True)
    assert margins[0] > 0.0


def test_canonical_restrained_case_values() -> None:
    """Section 14. The canonical eta_r = 1.0 case is useful but not pathological."""
    assessment = assess_at(1.0)
    assert assessment.hot_result.member_1_stress == pytest.approx(-111.5139e6, rel=1e-5)
    assert assessment.hot_result.member_2_stress == pytest.approx(-15.7361e6, rel=1e-5)
    assert assessment.cold_result.member_1_stress == pytest.approx(+156.1194e6, rel=1e-5)
    assert assessment.cold_result.member_2_stress == pytest.approx(+22.0306e6, rel=1e-5)
    assert assessment.governing_extreme == "cold"
    assert assessment.governing_member == 1
    assert assessment.minimum_yield_margin == pytest.approx(0.3836, rel=1e-3)
    assert assessment.passes


def test_symmetric_excursions_give_equal_margins_under_restraint() -> None:
    """Section 20. |dT_hot| = |dT_cold| with eps_ref = 0 keeps hot/cold symmetric."""
    symmetric = ThermalEnvironment(20.0, -80.0, 120.0)
    for ratio in RATIOS:
        assessment = assess_at(ratio, environment=symmetric)
        assert abs(assessment.cold_result.member_1_stress) == pytest.approx(
            abs(assessment.hot_result.member_1_stress), rel=1e-12
        )
        assert abs(assessment.cold_result.member_2_stress) == pytest.approx(
            abs(assessment.hot_result.member_2_stress), rel=1e-12
        )
        assert assessment.cold_assessment.minimum_margin == pytest.approx(
            assessment.hot_assessment.minimum_margin, rel=1e-12
        )


def test_symmetric_excursions_reverse_every_sign_under_restraint() -> None:
    """Section 20. Only the signs differ between the two symmetric extremes."""
    symmetric = ThermalEnvironment(20.0, -80.0, 120.0)
    assessment = assess_at(1.0, environment=symmetric)
    assert assessment.cold_result.member_1_stress == pytest.approx(
        -assessment.hot_result.member_1_stress, rel=1e-12
    )
    assert assessment.cold_result.restraint_force == pytest.approx(
        -assessment.hot_result.restraint_force, rel=1e-12
    )


def test_cold_governs_at_every_restraint_level_for_the_canonical_case() -> None:
    """Section 21. Computed, not assumed: cold governs across the whole sweep."""
    for ratio in RATIOS:
        assessment = assess_at(ratio)
        expected = min(
            assessment.cold_assessment.minimum_margin,
            assessment.hot_assessment.minimum_margin,
        )
        assert assessment.minimum_yield_margin == expected
        assert assessment.governing_extreme == "cold"


def test_hot_can_govern_under_restraint_when_the_excursion_is_reversed() -> None:
    """Section 21. The governing extreme is not hard-coded to cold."""
    hot_biased = ThermalEnvironment(20.0, -30.0, 320.0)
    for ratio in (0.0, 1.0, 10.0):
        restraint = AxialRestraint.from_stiffness_ratio(ratio, MEMBER_AL, MEMBER_TI)
        assessment = assess_restrained_temperature_extremes(
            MEMBER_AL, MEMBER_TI, hot_biased, restraint, BASIS
        )
        assert assessment.governing_extreme == "hot"


def test_governing_member_is_computed_under_restraint() -> None:
    """Section 19. The governing member follows the smaller margin, not the yield."""
    low_yield_large_area = AxialMember(
        ThermoelasticMaterial("low yield", 70.0e9, 23.0e-6, 100.0e6), 1000.0e-6
    )
    high_yield_small_area = AxialMember(
        ThermoelasticMaterial("high yield", 110.0e9, 8.5e-6, 500.0e6), 10.0e-6
    )
    restraint = AxialRestraint.from_stiffness_ratio(
        0.5, low_yield_large_area, high_yield_small_area
    )
    assessment = assess_restrained_temperature_extremes(
        low_yield_large_area, high_yield_small_area, ENVIRONMENT, restraint
    )
    assert assessment.governing_member == assessment.governing_assessment.governing_member


def test_governing_extreme_tie_resolves_to_cold() -> None:
    """Section 19. Symmetric excursions tie; the tie resolves deterministically."""
    assert assess_at(1.0, environment=ThermalEnvironment(20.0, -80.0, 120.0)).governing_extreme == "cold"


def test_repeated_restrained_assessment_is_deterministic() -> None:
    """Repeating the assessment reproduces every field exactly."""
    assert assess_at(1.0) == assess_at(1.0)


def test_defaults_and_type_checking() -> None:
    """Omitting the restraint or basis falls back to free joint / unity factor."""
    assert assess_restrained_temperature_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT
    ) == assess_restrained_temperature_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, AxialRestraint(0.0), YieldBasis(1.0)
    )
    with pytest.raises(TypeError):
        assess_restrained_temperature_extremes(MEMBER_AL, MEMBER_TI, "nope")  # type: ignore[arg-type]


def test_near_zero_stress_member_is_effectively_non_governing() -> None:
    """At the crossing the titanium-like member has a vast margin, so member 1 governs.

    The crossing stress is a floating-point residue (order 1e-8 Pa) rather than
    a bit-exact zero, so the margin is enormous but finite. The exactly-zero
    case returning ``math.inf`` is covered by the Milestone 1 margin tests.
    """
    from thermal_joint import zero_stress_restraint_stiffness

    crossing = zero_stress_restraint_stiffness(MEMBER_AL, MEMBER_TI, -140.0, 2)
    restraint = AxialRestraint(crossing)
    result = solve_restrained_joint(MEMBER_AL, MEMBER_TI, -140.0, restraint)
    assessment = assess_yield(result, MEMBER_AL, MEMBER_TI, BASIS)

    assert abs(result.member_2_stress) < 1.0e-3
    assert assessment.member_2_margin.margin_of_safety > 1.0e12
    assert assessment.governing_member == 1
    assert math.isfinite(assessment.minimum_margin)


def test_exactly_zero_stress_under_restraint_still_gives_an_infinite_margin() -> None:
    """A bit-exact zero stress under restraint returns math.inf, as in Milestone 1."""
    athermal = ThermoelasticMaterial("athermal", 100.0e9, 0.0, 300.0e6)
    member_1 = AxialMember(athermal, 100.0e-6)
    member_2 = AxialMember(athermal, 200.0e-6)
    result = solve_restrained_joint(
        member_1, member_2, 100.0, AxialRestraint.from_stiffness_ratio(2.0, member_1, member_2)
    )
    assessment = assess_yield(result, member_1, member_2, BASIS)
    assert result.member_1_stress == 0.0
    assert math.isinf(assessment.minimum_margin)

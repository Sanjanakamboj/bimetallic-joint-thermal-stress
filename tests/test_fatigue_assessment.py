"""The integrated static + fatigue assessment and its governing-component logic."""

from __future__ import annotations

import math

import pytest
from fatigue_cases import (
    ADHESIVE,
    ADHESIVE_BASIS,
    BASELINE_GEOMETRY,
    BASELINE_TAU_A,
    CURVES,
    ENVIRONMENT,
    FREE_ALTERNATING,
    LIFE_BASELINE_ADHESIVE,
    LIFE_MEMBER_1,
    LIFE_MEMBER_2,
    LIFE_SELECTED_ADHESIVE,
    MEMBER_AL,
    MEMBER_TI,
    REQUIREMENT,
    SELECTED_GEOMETRY,
    SELECTED_TAU_A,
    YIELD_BASIS,
)

from thermal_joint import (
    AdhesiveShearBasis,
    BasquinFatigueCurve,
    FatigueCurveSet,
    ThermalCycleRequirement,
    assess_shear_lag_extremes,
    assess_temperature_extremes,
    assess_thermal_cycle_fatigue,
)
from thermal_joint.illustrative import illustrative_cycle_requirement as _req


def assess(geometry, requirement=REQUIREMENT):
    return assess_thermal_cycle_fatigue(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, geometry, ADHESIVE, CURVES, requirement,
        YIELD_BASIS, ADHESIVE_BASIS,
    )


def test_baseline_assessment_canonical_values() -> None:
    """M3 baseline: metals hugely safe, adhesive fails statically and in fatigue."""
    result = assess(BASELINE_GEOMETRY)
    assert result.member_1_fatigue.alternating_stress == pytest.approx(FREE_ALTERNATING, rel=1e-5)
    assert result.member_1_fatigue.predicted_cycles_to_failure == pytest.approx(
        LIFE_MEMBER_1, rel=1e-3
    )
    assert result.member_2_fatigue.predicted_cycles_to_failure == pytest.approx(
        LIFE_MEMBER_2, rel=1e-3
    )
    assert result.adhesive_fatigue.alternating_stress == pytest.approx(BASELINE_TAU_A, rel=1e-5)
    assert result.adhesive_fatigue.predicted_cycles_to_failure == pytest.approx(
        LIFE_BASELINE_ADHESIVE, rel=1e-3
    )
    assert not result.static_adhesive_feasible
    assert result.static_yield_feasible
    assert not result.fatigue_feasible
    assert not result.overall_feasible


def test_selected_point_passes_statically_but_fails_fatigue() -> None:
    """The headline Milestone 5 result for the M4 selected bounded-grid point."""
    result = assess(SELECTED_GEOMETRY)
    assert result.static_yield_feasible
    assert result.static_adhesive_feasible
    assert result.minimum_adhesive_margin == pytest.approx(0.0648, abs=1e-3)

    assert result.adhesive_fatigue.alternating_stress == pytest.approx(SELECTED_TAU_A, rel=1e-5)
    assert result.adhesive_fatigue.predicted_cycles_to_failure == pytest.approx(
        LIFE_SELECTED_ADHESIVE, rel=1e-3
    )
    assert result.adhesive_fatigue.life_ratio == pytest.approx(0.6439, rel=1e-3)
    assert not result.adhesive_fatigue.passes
    assert not result.fatigue_feasible
    assert not result.overall_feasible
    assert result.governing_fatigue_component == "adhesive shear"


def test_metal_fatigue_is_not_governing_for_the_canonical_joint() -> None:
    """Both metals clear 1e4 cycles by many orders of magnitude - stated explicitly."""
    result = assess(SELECTED_GEOMETRY)
    assert result.member_1_fatigue.passes
    assert result.member_2_fatigue.passes
    assert result.member_1_fatigue.life_ratio > 1.0e4
    assert result.member_2_fatigue.life_ratio > 1.0e9


def test_governing_component_is_computed_not_assumed() -> None:
    """A curve set that makes member 1 weak moves the governing component."""
    weak_metal = FatigueCurveSet(
        member_1=BasquinFatigueCurve("weak Al", 90.0e6, -0.12, "synthetic test fixture"),
        member_2=CURVES.member_2,
        adhesive_shear=CURVES.adhesive_shear,
    )
    result = assess_thermal_cycle_fatigue(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, SELECTED_GEOMETRY, ADHESIVE, weak_metal,
        REQUIREMENT, YIELD_BASIS, ADHESIVE_BASIS,
    )
    assert result.governing_fatigue_component == "member 1"
    assert result.minimum_life_ratio == result.member_1_fatigue.life_ratio
    assert result.governing_fatigue_result is result.member_1_fatigue


def test_governing_component_matches_the_smallest_life_ratio() -> None:
    for geometry in (BASELINE_GEOMETRY, SELECTED_GEOMETRY):
        result = assess(geometry)
        ratios = [r.life_ratio for r in result.fatigue_results]
        assert result.minimum_life_ratio == min(ratios)
        assert result.governing_fatigue_result.life_ratio == min(ratios)


def test_overall_feasibility_is_a_three_way_boolean_and() -> None:
    """Static yield AND static adhesive AND every fatigue component."""
    for geometry in (BASELINE_GEOMETRY, SELECTED_GEOMETRY):
        result = assess(geometry)
        assert result.overall_feasible == (
            result.static_yield_feasible
            and result.static_adhesive_feasible
            and result.fatigue_feasible
        )
        assert result.fatigue_feasible == all(r.passes for r in result.fatigue_results)


def test_static_and_fatigue_margins_stay_separate() -> None:
    """No blended margin is exposed; the static values are the M1/M3 ones."""
    result = assess(SELECTED_GEOMETRY)
    yield_assessment = assess_temperature_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, YIELD_BASIS
    )
    shear_assessment = assess_shear_lag_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, SELECTED_GEOMETRY, ADHESIVE, ADHESIVE_BASIS
    )
    assert result.minimum_yield_margin == yield_assessment.minimum_yield_margin
    assert result.minimum_adhesive_margin == shear_assessment.minimum_margin
    assert result.minimum_yield_margin != result.minimum_adhesive_margin
    assert result.minimum_life_ratio != result.minimum_adhesive_margin

    attributes = set(dir(result))
    assert not {"combined_margin", "overall_margin", "blended_margin"} & attributes


def test_requirement_sensitivity_for_the_selected_point() -> None:
    """The 1e4 choice is transparent: 1e3 passes, 1e4 and above fail."""
    ratios = {}
    for cycles in (1.0e3, 1.0e4, 1.0e5, 1.0e6):
        result = assess(SELECTED_GEOMETRY, ThermalCycleRequirement(cycles, "sweep"))
        ratios[cycles] = result.adhesive_fatigue.life_ratio
        assert result.governing_fatigue_component == "adhesive shear"

    assert ratios[1.0e3] == pytest.approx(6.4393, rel=1e-3)
    assert ratios[1.0e4] == pytest.approx(0.6439, rel=1e-3)
    assert ratios[1.0e3] > 1.0
    assert ratios[1.0e4] < 1.0
    # The ratio is exactly inversely proportional to the requirement.
    assert ratios[1.0e3] / ratios[1.0e4] == pytest.approx(10.0, rel=1e-12)


def test_assessment_is_deterministic() -> None:
    assert assess(SELECTED_GEOMETRY) == assess(SELECTED_GEOMETRY)


def test_assessment_defaults_and_type_checking() -> None:
    from thermal_joint import YieldBasis

    assert assess_thermal_cycle_fatigue(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, SELECTED_GEOMETRY, ADHESIVE, CURVES, REQUIREMENT
    ) == assess_thermal_cycle_fatigue(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, SELECTED_GEOMETRY, ADHESIVE, CURVES, REQUIREMENT,
        YieldBasis(1.0), AdhesiveShearBasis(1.0),
    )
    with pytest.raises(TypeError):
        assess_thermal_cycle_fatigue(
            MEMBER_AL, MEMBER_TI, ENVIRONMENT, SELECTED_GEOMETRY, ADHESIVE,
            "not-a-curve-set", REQUIREMENT,
        )  # type: ignore[arg-type]


def test_curve_set_type_checking() -> None:
    with pytest.raises(TypeError):
        FatigueCurveSet(CURVES.member_1, CURVES.member_2, "nope")  # type: ignore[arg-type]


def test_illustrative_requirement_helper_accepts_an_override() -> None:
    assert _req().required_cycles == 1.0e4
    assert _req(1.0e5).required_cycles == 1.0e5
    assert "not a qualification requirement" in _req(1.0e5).label


def test_zero_excursion_environment_has_no_fatigue_demand() -> None:
    """A joint that never cycles has infinite life on every component."""
    from thermal_joint import ThermalEnvironment

    result = assess_thermal_cycle_fatigue(
        MEMBER_AL, MEMBER_TI, ThermalEnvironment(20.0, 20.0, 20.0), SELECTED_GEOMETRY,
        ADHESIVE, CURVES, REQUIREMENT, YIELD_BASIS, ADHESIVE_BASIS,
    )
    for component in result.fatigue_results:
        assert math.isinf(component.predicted_cycles_to_failure)
        assert component.passes
    assert result.fatigue_feasible
    assert result.overall_feasible

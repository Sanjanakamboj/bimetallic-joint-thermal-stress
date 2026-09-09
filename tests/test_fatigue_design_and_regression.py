"""Allowable cycle scale, the bounded fatigue design map, and M1-M4 regression."""

from __future__ import annotations

import pytest
from fatigue_cases import (
    ADHESIVE,
    ADHESIVE_BASIS,
    BASELINE_GEOMETRY,
    CURVES,
    ENVIRONMENT,
    MEMBER_AL,
    MEMBER_TI,
    MM,
    REQUIREMENT,
    SELECTED_GEOMETRY,
    YIELD_BASIS,
)

from thermal_joint import (
    AxialRestraint,
    ThermalCycleRequirement,
    allowable_cycle_scale_for_fatigue,
    assess_restrained_temperature_extremes,
    assess_shear_lag_extremes,
    assess_temperature_extremes,
    assess_thermal_cycle_fatigue,
    fatigue_design_map,
    scale_environment,
    select_preliminary_fatigue_design,
    select_preliminary_joint_design,
    solve_bimetallic_joint,
    solve_shear_lag,
    width_thickness_design_map,
)

MAP_WIDTHS_MM = (20.0, 30.0, 40.0, 50.0, 75.0, 100.0)
MAP_THICKNESSES_MM = (0.5, 0.75, 1.0, 1.5, 2.0)


def canonical_map():
    return fatigue_design_map(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASELINE_GEOMETRY, ADHESIVE, CURVES,
        REQUIREMENT, [w * MM for w in MAP_WIDTHS_MM], [t * MM for t in MAP_THICKNESSES_MM],
        YIELD_BASIS, ADHESIVE_BASIS,
    )


# ------------------------------------------------- allowable cycle scale

def test_allowable_cycle_scale_canonical_values() -> None:
    """The adhesive governs; the canonical excursion is ~6% too large."""
    result = allowable_cycle_scale_for_fatigue(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, SELECTED_GEOMETRY, ADHESIVE, CURVES, REQUIREMENT
    )
    assert result.governing_component == "adhesive shear"
    assert result.allowable_scale == pytest.approx(0.936107, rel=1e-5)
    assert result.adhesive_scale == pytest.approx(0.936107, rel=1e-5)
    assert result.member_1_scale == pytest.approx(4.0038, rel=1e-4)
    assert result.member_2_scale == pytest.approx(10.697, rel=1e-4)
    assert not result.canonical_cycle_passes


def test_allowable_cycle_scale_round_trip_is_exact() -> None:
    """At the returned scale the minimum life ratio is 1 - exact, not iterative."""
    result = allowable_cycle_scale_for_fatigue(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, SELECTED_GEOMETRY, ADHESIVE, CURVES, REQUIREMENT
    )
    scaled = assess_thermal_cycle_fatigue(
        MEMBER_AL, MEMBER_TI, scale_environment(ENVIRONMENT, result.allowable_scale),
        SELECTED_GEOMETRY, ADHESIVE, CURVES, REQUIREMENT, YIELD_BASIS, ADHESIVE_BASIS,
    )
    assert scaled.minimum_life_ratio == pytest.approx(1.0, rel=1e-9)
    # The ratio lands within an ULP of exactly 1.0, so the strict >= boolean can
    # fall either side of the boundary purely on rounding. Feasibility just
    # inside the boundary is asserted in the next test instead.
    assert abs(scaled.minimum_life_ratio - 1.0) < 1.0e-12


def test_just_inside_and_outside_the_allowable_scale() -> None:
    result = allowable_cycle_scale_for_fatigue(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, SELECTED_GEOMETRY, ADHESIVE, CURVES, REQUIREMENT
    )

    def feasible_at(scale: float) -> bool:
        return assess_thermal_cycle_fatigue(
            MEMBER_AL, MEMBER_TI, scale_environment(ENVIRONMENT, scale), SELECTED_GEOMETRY,
            ADHESIVE, CURVES, REQUIREMENT, YIELD_BASIS, ADHESIVE_BASIS,
        ).fatigue_feasible

    assert feasible_at(result.allowable_scale * 0.99)
    assert not feasible_at(result.allowable_scale * 1.01)


def test_allowable_scale_exceeds_one_for_a_generous_bondline() -> None:
    """A wide, thick bond survives the full canonical excursion with margin."""
    from thermal_joint.illustrative import radiator_joint_overlap

    result = allowable_cycle_scale_for_fatigue(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, radiator_joint_overlap(40.0, 100.0, 2.0),
        ADHESIVE, CURVES, REQUIREMENT,
    )
    assert result.allowable_scale > 1.0
    assert result.canonical_cycle_passes


# ---------------------------------------------------- fatigue design map

def test_map_shape_ordering_and_determinism() -> None:
    candidates = canonical_map()
    assert len(candidates) == len(MAP_WIDTHS_MM) * len(MAP_THICKNESSES_MM) == 30
    expected = [(w * MM, t * MM) for w in MAP_WIDTHS_MM for t in MAP_THICKNESSES_MM]
    for (want_w, want_t), candidate in zip(expected, candidates):
        assert candidate.bond_width == pytest.approx(want_w, rel=1e-12)
        assert candidate.adhesive_thickness == pytest.approx(want_t, rel=1e-12)
    assert canonical_map() == canonical_map()


def test_map_requires_all_three_screens() -> None:
    for candidate in canonical_map():
        assessment = candidate.assessment
        assert candidate.overall_feasible == (
            assessment.static_yield_feasible
            and assessment.static_adhesive_feasible
            and assessment.fatigue_feasible
        )


def test_map_feasible_count_and_membership() -> None:
    """12 of 30 grid points clear static yield, static shear and fatigue."""
    feasible = [c for c in canonical_map() if c.overall_feasible]
    assert len(feasible) == 12
    combinations = {
        (round(c.bond_width / MM), round(c.adhesive_thickness / MM, 2)) for c in feasible
    }
    assert (30, 2.0) in combinations
    assert (20, 2.0) not in combinations  # fails static shear
    assert (30, 1.5) not in combinations  # passes static, fails fatigue


def test_selection_picks_the_smallest_feasible_bond_area() -> None:
    selected = select_preliminary_fatigue_design(canonical_map())
    assert selected is not None
    assert selected.bond_width == pytest.approx(30.0 * MM, rel=1e-12)
    assert selected.adhesive_thickness == pytest.approx(2.0 * MM, rel=1e-12)
    assert selected.bond_area == min(
        c.bond_area for c in canonical_map() if c.overall_feasible
    )
    assessment = selected.assessment
    assert assessment.overall_feasible
    assert assessment.minimum_life_ratio == pytest.approx(1.0545, rel=1e-3)
    assert assessment.minimum_adhesive_margin == pytest.approx(0.1465, rel=1e-3)


def test_selection_is_order_independent_and_deterministic() -> None:
    candidates = list(canonical_map())
    assert select_preliminary_fatigue_design(candidates) == select_preliminary_fatigue_design(
        list(reversed(candidates))
    )
    rotated = candidates[13:] + candidates[:13]
    assert select_preliminary_fatigue_design(rotated) == select_preliminary_fatigue_design(
        candidates
    )


def test_selection_returns_none_when_nothing_is_feasible() -> None:
    narrow = fatigue_design_map(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASELINE_GEOMETRY, ADHESIVE, CURVES,
        REQUIREMENT, [20.0 * MM], [0.5 * MM, 0.75 * MM], YIELD_BASIS, ADHESIVE_BASIS,
    )
    assert all(not c.overall_feasible for c in narrow)
    assert select_preliminary_fatigue_design(narrow) is None
    assert select_preliminary_fatigue_design([]) is None


def test_fatigue_selection_differs_from_the_milestone_4_static_selection() -> None:
    """The static-selected point is not the fatigue-selected point.

    The Milestone 4 grid stopped at 1.0 mm bondlines; allowing thicker ones lets
    fatigue be met on a NARROWER bond, so the two policies land in different
    places even though the ranking rule is identical.
    """
    static_selected = select_preliminary_joint_design(
        width_thickness_design_map(
            MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASELINE_GEOMETRY, ADHESIVE,
            [w * MM for w in (10, 20, 30, 40, 50, 75, 100)],
            [t * MM for t in (0.1, 0.2, 0.3, 0.5, 0.75, 1.0)],
            YIELD_BASIS, ADHESIVE_BASIS,
        )
    )
    fatigue_selected = select_preliminary_fatigue_design(canonical_map())

    assert static_selected.bond_width == pytest.approx(50.0 * MM, rel=1e-12)
    assert fatigue_selected.bond_width == pytest.approx(30.0 * MM, rel=1e-12)
    assert fatigue_selected.bond_area < static_selected.bond_area

    # The M4 point itself fails the fatigue screen.
    at_static_point = assess_thermal_cycle_fatigue(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, SELECTED_GEOMETRY, ADHESIVE, CURVES,
        REQUIREMENT, YIELD_BASIS, ADHESIVE_BASIS,
    )
    assert at_static_point.static_adhesive_feasible
    assert not at_static_point.fatigue_feasible


# ------------------------------------------------------------ regression

def test_milestone_1_outputs_unchanged() -> None:
    hot = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, 100.0)
    cold = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, -140.0)
    assert hot.member_1_stress == pytest.approx(-62.0278e6, rel=1e-5)
    assert cold.member_1_stress == pytest.approx(+86.8389e6, rel=1e-5)
    assert assess_temperature_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, YIELD_BASIS
    ).minimum_yield_margin == pytest.approx(1.4874, rel=1e-3)


def test_milestone_2_outputs_unchanged() -> None:
    from thermal_joint import RestraintLimitStatus, maximum_allowable_restraint_stiffness

    restrained = assess_restrained_temperature_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT,
        AxialRestraint.from_stiffness_ratio(1.0, MEMBER_AL, MEMBER_TI), YIELD_BASIS,
    )
    assert restrained.cold_result.member_1_stress == pytest.approx(156.1194e6, rel=1e-5)
    assert restrained.minimum_yield_margin == pytest.approx(0.3836, rel=1e-3)
    limit = maximum_allowable_restraint_stiffness(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, YIELD_BASIS
    )
    assert limit.status is RestraintLimitStatus.FINITE_LIMIT
    assert limit.maximum_stiffness_ratio == pytest.approx(13.740543, abs=1e-5)


def test_milestone_3_outputs_unchanged() -> None:
    from thermal_joint import OverlapLimitStatus, required_overlap_length

    cold = solve_shear_lag(MEMBER_AL, MEMBER_TI, BASELINE_GEOMETRY, ADHESIVE, -140.0)
    assert cold.beta == pytest.approx(152.894157, rel=1e-8)
    assert cold.peak_shear_stress == pytest.approx(66.3864e6, rel=1e-5)
    assert assess_shear_lag_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASELINE_GEOMETRY, ADHESIVE, ADHESIVE_BASIS
    ).minimum_margin == pytest.approx(-0.6987, rel=1e-3)
    assert required_overlap_length(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASELINE_GEOMETRY, ADHESIVE, ADHESIVE_BASIS
    ).status is OverlapLimitStatus.NO_FINITE_LENGTH_WITHIN_MODEL


def test_milestone_4_outputs_unchanged() -> None:
    from thermal_joint import (
        maximum_allowable_adhesive_shear_modulus,
        required_adhesive_thickness,
        required_bond_width,
    )

    assert required_bond_width(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASELINE_GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        tolerance=1.0e-9,
    ).required_bond_width == pytest.approx(0.220354, rel=1e-4)
    assert required_adhesive_thickness(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASELINE_GEOMETRY, ADHESIVE, ADHESIVE_BASIS
    ).required_adhesive_thickness == pytest.approx(2.4987e-3, rel=1e-3)
    assert maximum_allowable_adhesive_shear_modulus(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASELINE_GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        tolerance=1.0,
    ).maximum_shear_modulus == pytest.approx(0.080043e9, rel=1e-3)


def test_milestone_5_does_not_perturb_the_earlier_models() -> None:
    before = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, -140.0)
    assess_thermal_cycle_fatigue(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, SELECTED_GEOMETRY, ADHESIVE, CURVES,
        REQUIREMENT, YIELD_BASIS, ADHESIVE_BASIS,
    )
    canonical_map()
    after = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, -140.0)
    assert before == after


def test_requirement_choice_is_reported_not_hidden() -> None:
    """A different requirement changes the map outcome - so it must stay visible."""
    lenient = fatigue_design_map(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASELINE_GEOMETRY, ADHESIVE, CURVES,
        ThermalCycleRequirement(1.0e3, "lenient sweep"),
        [w * MM for w in MAP_WIDTHS_MM], [t * MM for t in MAP_THICKNESSES_MM],
        YIELD_BASIS, ADHESIVE_BASIS,
    )
    assert sum(c.overall_feasible for c in lenient) > sum(
        c.overall_feasible for c in canonical_map()
    )

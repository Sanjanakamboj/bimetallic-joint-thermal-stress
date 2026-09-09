"""Milestone 5 sensitivity sweeps: restraint, temperature, geometry, area ratio."""

from __future__ import annotations

import math

import pytest
from fatigue_cases import (
    ADHESIVE,
    ADHESIVE_BASIS,
    CURVES,
    ENVIRONMENT,
    FREE_ALTERNATING,
    MEMBER_AL,
    MEMBER_TI,
    MM,
    REQUIREMENT,
    SELECTED_GEOMETRY,
    YIELD_BASIS,
)

from thermal_joint import (
    AxialRestraint,
    adhesive_modulus_fatigue_sensitivity,
    adhesive_thickness_fatigue_sensitivity,
    area_ratio_fatigue_sensitivity,
    assess_restrained_temperature_extremes,
    bond_width_fatigue_sensitivity,
    joint_axial_stiffness,
    member_stress_cycles,
    restraint_fatigue_sensitivity,
    scale_environment,
    solve_bimetallic_joint,
    temperature_scale_fatigue_sensitivity,
    zero_stress_restraint_stiffness,
)
from thermal_joint.illustrative import ALUMINIUM_LIKE, ADHESIVE_LIKE, radiator_joint_overlap

ETA_RATIOS = (0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 13.7405)
TEMPERATURE_SCALES = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5)
WIDTHS_MM = (20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 220.0)
THICKNESSES_MM = (0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0)
MODULI_GPA = (0.05, 0.1, 0.2, 0.5, 1.0, 2.0)
AREA_RATIOS = (0.25, 0.5, 1.0, 2.0, 4.0)


def restraint_points():
    return restraint_fatigue_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, CURVES, REQUIREMENT, ETA_RATIOS
    )


# ------------------------------------------------------------- restraint

def test_restraint_zero_reproduces_the_free_joint_cycles() -> None:
    """eta_r = 0 is exactly the Milestone 1 free-joint cycle."""
    point = restraint_points()[0]
    assert point.stiffness_ratio == 0.0
    assert point.member_1_cycle.alternating_stress == pytest.approx(FREE_ALTERNATING, rel=1e-5)
    assert point.member_2_cycle.alternating_stress == pytest.approx(FREE_ALTERNATING, rel=1e-5)
    hot = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, 100.0)
    assert point.member_1_cycle.hot_value == pytest.approx(hot.member_1_stress, rel=1e-13)


def test_restraint_uses_the_verified_milestone_2_stresses() -> None:
    """The eta_r = 1 endpoints are the Milestone 2 restrained stresses."""
    restraint = AxialRestraint.from_stiffness_ratio(1.0, MEMBER_AL, MEMBER_TI)
    expected = assess_restrained_temperature_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, restraint, YIELD_BASIS
    )
    cycle_1, _ = member_stress_cycles(MEMBER_AL, MEMBER_TI, ENVIRONMENT, restraint)
    assert cycle_1.hot_value == pytest.approx(expected.hot_result.member_1_stress, rel=1e-14)
    assert cycle_1.cold_value == pytest.approx(expected.cold_result.member_1_stress, rel=1e-14)
    assert cycle_1.alternating_stress == pytest.approx(133.8167e6, rel=1e-5)
    assert cycle_1.mean_stress == pytest.approx(22.3028e6, rel=1e-4)


def test_restraint_raises_member_1_amplitude_monotonically() -> None:
    """The high-CTE member is loaded harder and harder as restraint rises."""
    amplitudes = [p.member_1_cycle.alternating_stress for p in restraint_points()]
    assert amplitudes == sorted(amplitudes)
    assert amplitudes[0] == pytest.approx(FREE_ALTERNATING, rel=1e-5)
    assert amplitudes[-1] == pytest.approx(185.14e6, rel=1e-3)

    lives = [p.member_1_fatigue.predicted_cycles_to_failure for p in restraint_points()]
    assert lives == sorted(lives, reverse=True)


def test_restraint_does_not_worsen_both_members(  ) -> None:
    """Member 2's amplitude is NOT monotonic - it falls, then rises again.

    Restraint therefore does not degrade both members identically; at
    eta_r = 1 member 1's amplitude is up 80% while member 2's is down 75%.
    """
    points = restraint_points()
    amplitudes = [p.member_2_cycle.alternating_stress for p in points]
    assert amplitudes != sorted(amplitudes)
    assert amplitudes != sorted(amplitudes, reverse=True)
    assert min(amplitudes) < amplitudes[0]
    assert amplitudes[-1] > min(amplitudes)

    at_unity = [p for p in points if p.stiffness_ratio == 1.0][0]
    assert at_unity.member_1_cycle.alternating_stress > FREE_ALTERNATING
    assert at_unity.member_2_cycle.alternating_stress < FREE_ALTERNATING
    assert at_unity.member_2_cycle.alternating_stress == pytest.approx(18.8833e6, rel=1e-4)


def test_member_2_cycle_vanishes_at_the_milestone_2_zero_stress_crossing() -> None:
    """The restraint that nulls member 2's static stress also nulls its cycle.

    The Milestone 2 crossing is dT-independent, so BOTH endpoints vanish there
    and the titanium-like member experiences no thermal cycle at all.
    """
    crossing = zero_stress_restraint_stiffness(MEMBER_AL, MEMBER_TI, 100.0, 2)
    ratio = crossing / joint_axial_stiffness(MEMBER_AL, MEMBER_TI)
    assert ratio == pytest.approx(0.663399, rel=1e-5)

    _, cycle_2 = member_stress_cycles(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, AxialRestraint(crossing)
    )
    assert abs(cycle_2.hot_value) < 1.0
    assert abs(cycle_2.cold_value) < 1.0
    assert cycle_2.alternating_stress < 1.0
    life = CURVES.member_2.life_at_alternating_stress(cycle_2.alternating_stress)
    assert life > 1.0e50


def test_all_restraint_levels_pass_metal_fatigue_for_this_joint() -> None:
    """Reported as computed: metal fatigue never governs across the whole sweep."""
    assert all(point.passes for point in restraint_points())
    assert all(point.minimum_life_ratio > 1.0 for point in restraint_points())


def test_restraint_sweep_is_deterministic() -> None:
    assert restraint_points() == restraint_points()


# ---------------------------------------------------------- temperature

def test_temperature_scale_is_linear_in_amplitude() -> None:
    """Scaling both excursions scales every alternating stress exactly."""
    points = temperature_scale_fatigue_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, SELECTED_GEOMETRY, ADHESIVE, CURVES,
        REQUIREMENT, TEMPERATURE_SCALES, YIELD_BASIS, ADHESIVE_BASIS,
    )
    unit = [p for p in points if p.value == 1.0][0]
    for point in points:
        for component in ("member_1_fatigue", "member_2_fatigue", "adhesive_fatigue"):
            scaled = getattr(point.assessment, component).alternating_stress
            base = getattr(unit.assessment, component).alternating_stress
            assert scaled == pytest.approx(point.value * base, rel=1e-12)


def test_temperature_scale_follows_the_basquin_power_law() -> None:
    """Life ratio scales as scale^(1/b) - verified against the closed form."""
    points = temperature_scale_fatigue_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, SELECTED_GEOMETRY, ADHESIVE, CURVES,
        REQUIREMENT, TEMPERATURE_SCALES, YIELD_BASIS, ADHESIVE_BASIS,
    )
    unit = [p for p in points if p.value == 1.0][0]
    exponent = 1.0 / CURVES.adhesive_shear.exponent_b
    for point in points:
        expected = unit.assessment.adhesive_fatigue.predicted_cycles_to_failure * (
            point.value**exponent
        )
        assert point.assessment.adhesive_fatigue.predicted_cycles_to_failure == pytest.approx(
            expected, rel=1e-10
        )
    lives = [p.assessment.adhesive_fatigue.predicted_cycles_to_failure for p in points]
    assert lives == sorted(lives, reverse=True)


def test_scale_environment_preserves_the_asymmetry_ratio() -> None:
    for scale in TEMPERATURE_SCALES:
        scaled = scale_environment(ENVIRONMENT, scale)
        assert scaled.delta_temperature_hot == pytest.approx(100.0 * scale, rel=1e-12)
        assert scaled.delta_temperature_cold == pytest.approx(-140.0 * scale, rel=1e-12)
        assert scaled.reference_temperature == ENVIRONMENT.reference_temperature
    with pytest.raises(ValueError):
        scale_environment(ENVIRONMENT, -0.5)


# -------------------------------------------------------------- geometry

def test_bond_width_raises_adhesive_life_monotonically() -> None:
    points = bond_width_fatigue_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, radiator_joint_overlap(40.0, 20.0, 1.0),
        ADHESIVE, CURVES, REQUIREMENT, [w * MM for w in WIDTHS_MM],
        YIELD_BASIS, ADHESIVE_BASIS,
    )
    amplitudes = [p.assessment.adhesive_fatigue.alternating_stress for p in points]
    lives = [p.assessment.adhesive_fatigue.predicted_cycles_to_failure for p in points]
    assert amplitudes == sorted(amplitudes, reverse=True)
    assert lives == sorted(lives)


def test_static_and_fatigue_width_boundaries_differ() -> None:
    """The fatigue-critical width is NOT the static-critical width."""
    points = bond_width_fatigue_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, radiator_joint_overlap(40.0, 20.0, 1.0),
        ADHESIVE, CURVES, REQUIREMENT, [w * MM for w in WIDTHS_MM],
        YIELD_BASIS, ADHESIVE_BASIS,
    )
    static_ok = [p.value for p in points if p.assessment.static_adhesive_feasible]
    fatigue_ok = [p.value for p in points if p.assessment.fatigue_feasible]
    assert static_ok and fatigue_ok
    assert min(fatigue_ok) > min(static_ok)  # fatigue needs the wider bond


def test_thickness_raises_adhesive_life_monotonically() -> None:
    points = adhesive_thickness_fatigue_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, radiator_joint_overlap(40.0, 50.0, 1.0),
        ADHESIVE, CURVES, REQUIREMENT, [t * MM for t in THICKNESSES_MM],
        YIELD_BASIS, ADHESIVE_BASIS,
    )
    lives = [p.assessment.adhesive_fatigue.predicted_cycles_to_failure for p in points]
    assert lives == sorted(lives)
    assert points[0].value == pytest.approx(0.2 * MM, rel=1e-12)


def test_lower_modulus_raises_adhesive_life() -> None:
    points = adhesive_modulus_fatigue_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, SELECTED_GEOMETRY, ADHESIVE, CURVES,
        REQUIREMENT, [g * 1.0e9 for g in MODULI_GPA], YIELD_BASIS, ADHESIVE_BASIS,
    )
    lives = [p.assessment.adhesive_fatigue.predicted_cycles_to_failure for p in points]
    assert lives == sorted(lives, reverse=True)
    assert ADHESIVE_LIKE.shear_modulus == 1.0e9  # canonical record untouched


def test_area_ratio_sweep_recomputes_everything_and_reports_governing() -> None:
    """The static-best area ratio need not be the fatigue-best one - reported."""
    points = area_ratio_fatigue_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, SELECTED_GEOMETRY, ADHESIVE, CURVES,
        REQUIREMENT, AREA_RATIOS, YIELD_BASIS, ADHESIVE_BASIS,
    )
    adhesive_amplitudes = [p.assessment.adhesive_fatigue.alternating_stress for p in points]
    assert adhesive_amplitudes == sorted(adhesive_amplitudes)  # rises with area ratio
    for point in points:
        assert point.assessment.governing_fatigue_component == "adhesive shear"
    metal_amplitudes = [p.assessment.member_1_fatigue.alternating_stress for p in points]
    assert len(set(metal_amplitudes)) == len(metal_amplitudes)  # M1 recomputed each time
    assert ALUMINIUM_LIKE.thermal_expansion_coefficient == 23.0e-6


def test_sweeps_do_not_mutate_canonical_inputs() -> None:
    geometry = radiator_joint_overlap(40.0, 50.0, 1.0)
    bond_width_fatigue_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, geometry, ADHESIVE, CURVES, REQUIREMENT,
        [w * MM for w in WIDTHS_MM], YIELD_BASIS, ADHESIVE_BASIS,
    )
    adhesive_thickness_fatigue_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, geometry, ADHESIVE, CURVES, REQUIREMENT,
        [t * MM for t in THICKNESSES_MM], YIELD_BASIS, ADHESIVE_BASIS,
    )
    assert geometry == radiator_joint_overlap(40.0, 50.0, 1.0)
    assert ENVIRONMENT.delta_temperature_cold == -140.0
    assert ADHESIVE_LIKE.shear_strength == 25.0e6
    assert MEMBER_AL.area == pytest.approx(100.0e-6, rel=1e-12)

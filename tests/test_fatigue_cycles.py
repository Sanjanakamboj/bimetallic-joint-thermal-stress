"""Stress-cycle construction: signs, alternating/mean stress and the magnitude trap."""

from __future__ import annotations

import math

import pytest
from fatigue_cases import (
    ADHESIVE,
    BASELINE_GEOMETRY,
    BASELINE_TAU_A,
    BASELINE_TAU_M,
    ENVIRONMENT,
    FREE_ALTERNATING,
    FREE_MEAN_1,
    FREE_MEAN_2,
    MEMBER_AL,
    MEMBER_TI,
    SELECTED_GEOMETRY,
    SELECTED_TAU_A,
    SELECTED_TAU_M,
)

from thermal_joint import (
    adhesive_shear_cycle,
    member_stress_cycles,
    solve_bimetallic_joint,
    solve_shear_lag,
    stress_cycle,
)


def test_cycle_hand_calculation_for_the_free_aluminium_like_member() -> None:
    """hot -62.0278, cold +86.8389 -> sig_a 74.4333, sig_m +12.4056, R -0.7143."""
    cycle = stress_cycle(-62.0278e6, 86.8389e6)
    assert cycle.maximum_stress == pytest.approx(86.8389e6, rel=1e-12)
    assert cycle.minimum_stress == pytest.approx(-62.0278e6, rel=1e-12)
    assert cycle.alternating_stress == pytest.approx(74.43335e6, rel=1e-9)
    assert cycle.mean_stress == pytest.approx(12.40555e6, rel=1e-9)
    assert cycle.stress_ratio == pytest.approx(-62.0278 / 86.8389, rel=1e-12)
    assert cycle.crosses_zero


def test_alternating_stress_is_never_negative() -> None:
    """Whatever the endpoint order or sign, the amplitude is non-negative."""
    for hot, cold in ((-5.0, 3.0), (3.0, -5.0), (-9.0, -2.0), (7.0, 7.0), (0.0, 0.0)):
        assert stress_cycle(hot, cold).alternating_stress >= 0.0


def test_endpoint_order_does_not_change_the_cycle() -> None:
    """max/min are algebraic, so swapping hot and cold leaves sig_a and sig_m alone."""
    forward = stress_cycle(-62.0e6, 87.0e6)
    swapped = stress_cycle(87.0e6, -62.0e6)
    assert forward.alternating_stress == swapped.alternating_stress
    assert forward.mean_stress == swapped.mean_stress
    assert forward.stress_ratio == swapped.stress_ratio


def test_mean_stress_keeps_its_sign() -> None:
    """The two free members have equal amplitude but opposite mean stress."""
    cycle_1, cycle_2 = member_stress_cycles(MEMBER_AL, MEMBER_TI, ENVIRONMENT)
    assert cycle_1.alternating_stress == pytest.approx(cycle_2.alternating_stress, rel=1e-12)
    assert cycle_1.mean_stress == pytest.approx(-cycle_2.mean_stress, rel=1e-12)
    assert cycle_1.mean_stress > 0.0 > cycle_2.mean_stress
    assert cycle_1.alternating_stress == pytest.approx(FREE_ALTERNATING, rel=1e-5)
    assert cycle_1.mean_stress == pytest.approx(FREE_MEAN_1, rel=1e-4)
    assert cycle_2.mean_stress == pytest.approx(FREE_MEAN_2, rel=1e-4)


def test_member_cycles_come_from_the_verified_milestone_1_solver() -> None:
    """The endpoints are exactly the Milestone 1 free-joint stresses."""
    hot = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, 100.0)
    cold = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, -140.0)
    cycle_1, cycle_2 = member_stress_cycles(MEMBER_AL, MEMBER_TI, ENVIRONMENT)
    assert cycle_1.hot_value == pytest.approx(hot.member_1_stress, rel=1e-14)
    assert cycle_1.cold_value == pytest.approx(cold.member_1_stress, rel=1e-14)
    assert cycle_2.hot_value == pytest.approx(hot.member_2_stress, rel=1e-14)
    assert cycle_2.cold_value == pytest.approx(cold.member_2_stress, rel=1e-14)


def test_stress_ratio_is_none_only_when_the_maximum_is_exactly_zero() -> None:
    assert stress_cycle(-5.0, -3.0).stress_ratio is not None
    assert stress_cycle(0.0, -3.0).stress_ratio is None
    assert stress_cycle(0.0, 0.0).stress_ratio is None


def test_zero_amplitude_cycle_is_well_formed() -> None:
    cycle = stress_cycle(50.0e6, 50.0e6)
    assert cycle.alternating_stress == 0.0
    assert cycle.mean_stress == 50.0e6
    assert cycle.stress_ratio == pytest.approx(1.0, rel=1e-15)
    assert not cycle.crosses_zero


# ------------------------------------------------------- the adhesive sign trap

def test_adhesive_cycle_endpoints_have_opposite_signs() -> None:
    """The linear model reverses shear with dT, so the cycle crosses zero."""
    cycle = adhesive_shear_cycle(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASELINE_GEOMETRY, ADHESIVE
    )
    assert cycle.hot_value > 0.0
    assert cycle.cold_value < 0.0
    assert cycle.crosses_zero
    assert abs(cycle.hot_value) != pytest.approx(abs(cycle.cold_value), rel=1e-3)


def test_adhesive_alternating_shear_is_the_half_sum_of_magnitudes() -> None:
    """For opposite-sign endpoints, tau_a = (|tau_hot| + |tau_cold|)/2."""
    for geometry, expected_a, expected_m in (
        (BASELINE_GEOMETRY, BASELINE_TAU_A, BASELINE_TAU_M),
        (SELECTED_GEOMETRY, SELECTED_TAU_A, SELECTED_TAU_M),
    ):
        cycle = adhesive_shear_cycle(
            MEMBER_AL, MEMBER_TI, ENVIRONMENT, geometry, ADHESIVE
        )
        assert cycle.alternating_stress == pytest.approx(
            0.5 * (abs(cycle.hot_value) + abs(cycle.cold_value)), rel=1e-14
        )
        assert cycle.alternating_stress == pytest.approx(expected_a, rel=1e-5)
        assert cycle.mean_stress == pytest.approx(expected_m, rel=1e-4)
        assert cycle.mean_stress < 0.0  # asymmetric excursion -> nonzero mean


def test_using_magnitudes_instead_of_signed_endpoints_destroys_the_cycle() -> None:
    """Guard against the specific mistake of feeding two peak magnitudes in.

    Peak magnitudes are both positive, so the cycle collapses to a small
    amplitude about a large mean - here a 7x under-prediction of the amplitude
    and a badly wrong life.
    """
    cycle = adhesive_shear_cycle(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASELINE_GEOMETRY, ADHESIVE
    )
    wrong = stress_cycle(abs(cycle.hot_value), abs(cycle.cold_value))

    assert wrong.alternating_stress < 0.2 * cycle.alternating_stress
    assert wrong.alternating_stress == pytest.approx(9.4838e6, rel=1e-3)
    assert not wrong.crosses_zero
    assert cycle.crosses_zero
    # The production helper must produce the signed version, not the magnitude one.
    assert cycle.alternating_stress == pytest.approx(BASELINE_TAU_A, rel=1e-5)


def test_adhesive_cycle_reads_both_endpoints_at_the_same_station() -> None:
    """Both endpoints are the free-edge peak, x = L_b, so the pair is comparable."""
    station = BASELINE_GEOMETRY.overlap_length
    hot = solve_shear_lag(
        MEMBER_AL, MEMBER_TI, BASELINE_GEOMETRY, ADHESIVE, 100.0
    ).shear_stress(station)
    cold = solve_shear_lag(
        MEMBER_AL, MEMBER_TI, BASELINE_GEOMETRY, ADHESIVE, -140.0
    ).shear_stress(station)
    cycle = adhesive_shear_cycle(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASELINE_GEOMETRY, ADHESIVE
    )
    assert cycle.hot_value == hot
    assert cycle.cold_value == cold
    assert abs(hot) == pytest.approx(
        solve_shear_lag(MEMBER_AL, MEMBER_TI, BASELINE_GEOMETRY, ADHESIVE, 100.0).peak_shear_stress,
        rel=1e-12,
    )


def test_symmetric_excursion_gives_a_zero_mean_adhesive_cycle() -> None:
    """The nonzero mean is a consequence of the asymmetry, not of the model."""
    from thermal_joint import ThermalEnvironment

    symmetric = ThermalEnvironment(20.0, -80.0, 120.0)
    cycle = adhesive_shear_cycle(
        MEMBER_AL, MEMBER_TI, symmetric, BASELINE_GEOMETRY, ADHESIVE
    )
    assert cycle.mean_stress == pytest.approx(0.0, abs=1e-3)
    assert cycle.alternating_stress > 0.0


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_non_finite_endpoints_rejected(bad: float) -> None:
    with pytest.raises(ValueError):
        stress_cycle(bad, 1.0)
    with pytest.raises(ValueError):
        stress_cycle(1.0, bad)

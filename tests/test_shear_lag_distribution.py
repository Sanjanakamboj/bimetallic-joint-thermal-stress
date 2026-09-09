"""Tests R-Y: the closed-form distributions and the force-transfer identities."""

from __future__ import annotations

import math

import pytest
from shear_lag_cases import (
    HAND_ADHESIVE,
    HAND_AVERAGE_SHEAR,
    HAND_GEOMETRY,
    HAND_MEMBER_1,
    HAND_MEMBER_2,
    HAND_SIGNED_FORCE,
    HAND_TRANSFERRED_FORCE,
)

from thermal_joint import solve_shear_lag
from thermal_joint.illustrative import (
    ADHESIVE_LIKE,
    aluminium_like_member,
    radiator_joint_overlap,
    titanium_like_member,
)

DEMAND = solve_shear_lag(HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, HAND_ADHESIVE, 50.0)
LENGTH = HAND_GEOMETRY.overlap_length
WIDTH = HAND_GEOMETRY.bond_width


def simpson(function, lower: float, upper: float, intervals: int = 2000) -> float:
    """Composite Simpson's rule, implemented here so the check is independent."""
    if intervals % 2:
        intervals += 1
    step = (upper - lower) / intervals
    total = function(lower) + function(upper)
    for index in range(1, intervals):
        weight = 4.0 if index % 2 else 2.0
        total += weight * function(lower + index * step)
    return total * step / 3.0


def test_r_shear_stress_at_hand_calculated_positions() -> None:
    """R. tau(0) = 25e6/sinh(2) and tau(L_b) = 25e6*coth(2), from the hand case."""
    scale = HAND_TRANSFERRED_FORCE * 100.0 / WIDTH  # |N_t| beta / b = 25 MPa
    assert scale == pytest.approx(25.0e6, rel=1e-13)

    assert abs(DEMAND.shear_stress(0.0)) == pytest.approx(
        scale / math.sinh(2.0), rel=1e-12
    )
    assert abs(DEMAND.shear_stress(LENGTH)) == pytest.approx(
        scale * math.cosh(2.0) / math.sinh(2.0), rel=1e-12
    )
    assert abs(DEMAND.shear_stress(0.010)) == pytest.approx(
        scale * math.cosh(1.0) / math.sinh(2.0), rel=1e-12
    )


def test_s_shear_stress_keeps_one_sign_and_grows_toward_the_free_edge() -> None:
    """S. tau(x) ~ cosh(beta x): single-signed and strictly increasing in |tau|."""
    positions = [LENGTH * i / 20 for i in range(21)]
    magnitudes = [abs(DEMAND.shear_stress(x)) for x in positions]
    signs = {DEMAND.shear_stress(x) > 0.0 for x in positions}
    assert len(signs) == 1
    assert magnitudes == sorted(magnitudes)
    assert magnitudes[0] < magnitudes[-1]


def test_s_cooling_reverses_the_whole_distribution() -> None:
    """S. Reversing dT mirrors tau(x) and N_1(x) exactly."""
    cooling = solve_shear_lag(
        HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, HAND_ADHESIVE, -50.0
    )
    for x in (0.0, 0.005, 0.010, LENGTH):
        assert cooling.shear_stress(x) == pytest.approx(-DEMAND.shear_stress(x), rel=1e-13)
        assert cooling.member_1_force(x) == pytest.approx(
            -DEMAND.member_1_force(x), rel=1e-13, abs=1e-9
        )


def test_t_peak_is_at_the_free_edge_and_matches_the_reported_peak() -> None:
    """T. The analytic peak location is x = L_b; verified by scanning the domain."""
    assert DEMAND.peak_location == pytest.approx(LENGTH, rel=1e-15)
    scanned = max(
        abs(DEMAND.shear_stress(LENGTH * i / 5000)) for i in range(5001)
    )
    assert DEMAND.peak_shear_stress == pytest.approx(scanned, rel=1e-12)
    assert DEMAND.peak_shear_stress == pytest.approx(
        abs(DEMAND.shear_stress(DEMAND.peak_location)), rel=1e-13
    )


def test_t_peak_shear_hand_calculation() -> None:
    """T. tau_peak = 25 MPa * coth(2) = 25.9329 MPa."""
    assert DEMAND.peak_shear_stress == pytest.approx(
        25.0e6 / math.tanh(2.0), rel=1e-12
    )
    assert DEMAND.peak_shear_stress == pytest.approx(25.932868e6, rel=1e-6)


def test_u_member_1_force_distribution_hand_calculation() -> None:
    """U. N_1(x) = -5000 [1 - sinh(100x)/sinh(2)]; N_1(0) = -5000, N_1(L_b) = 0."""
    assert DEMAND.member_1_force(0.0) == pytest.approx(HAND_SIGNED_FORCE, rel=1e-13)
    assert DEMAND.member_1_force(LENGTH) == pytest.approx(0.0, abs=1e-9)
    assert DEMAND.member_1_force(0.010) == pytest.approx(
        HAND_SIGNED_FORCE * (1.0 - math.sinh(1.0) / math.sinh(2.0)), rel=1e-12
    )
    assert DEMAND.member_1_force(0.010) == pytest.approx(-3379.9, rel=1e-4)


def test_v_member_2_force_distribution_hand_calculation() -> None:
    """V. N_2(x) = -N_1(x) everywhere."""
    assert DEMAND.member_2_force(0.0) == pytest.approx(-HAND_SIGNED_FORCE, rel=1e-13)
    assert DEMAND.member_2_force(LENGTH) == pytest.approx(0.0, abs=1e-9)
    assert DEMAND.member_2_force(0.010) == pytest.approx(
        -HAND_SIGNED_FORCE * (1.0 - math.sinh(1.0) / math.sinh(2.0)), rel=1e-12
    )


def test_w_member_forces_sum_to_zero_everywhere() -> None:
    """W. The pair is self-equilibrating at every station in the overlap."""
    for index in range(51):
        x = LENGTH * index / 50
        total = DEMAND.member_1_force(x) + DEMAND.member_2_force(x)
        assert total == pytest.approx(0.0, abs=1e-12 * HAND_TRANSFERRED_FORCE)


def test_x_independent_numerical_integration_of_the_shear() -> None:
    """X. Simpson integration of tau(x) reproduces the adherend force change."""
    integral = simpson(DEMAND.shear_stress, 0.0, LENGTH)
    transferred = WIDTH * integral
    force_change = DEMAND.member_1_force(LENGTH) - DEMAND.member_1_force(0.0)
    assert transferred == pytest.approx(force_change, rel=1e-9)
    assert transferred == pytest.approx(-HAND_SIGNED_FORCE, rel=1e-9)


def test_x_independent_integration_over_partial_spans() -> None:
    """X. The identity holds on sub-intervals too, not just end to end."""
    for upper in (0.004, 0.010, 0.016, LENGTH):
        integral = WIDTH * simpson(DEMAND.shear_stress, 0.0, upper)
        change = DEMAND.member_1_force(upper) - DEMAND.member_1_force(0.0)
        assert integral == pytest.approx(change, rel=1e-8, abs=1e-9)


def test_x_independent_integration_for_the_illustrative_joint() -> None:
    """X. Same check on the canonical illustrative overlap at both extremes."""
    member_1, member_2 = aluminium_like_member(), titanium_like_member()
    geometry = radiator_joint_overlap()
    for delta_t in (100.0, -140.0):
        demand = solve_shear_lag(member_1, member_2, geometry, ADHESIVE_LIKE, delta_t)
        integral = geometry.bond_width * simpson(
            demand.shear_stress, 0.0, geometry.overlap_length, 4000
        )
        change = demand.member_1_force(geometry.overlap_length) - demand.member_1_force(0.0)
        assert integral == pytest.approx(change, rel=1e-8)
        assert integral == pytest.approx(-demand.signed_transferred_force, rel=1e-8)


def test_y_force_transfer_residual_is_near_machine_precision() -> None:
    """Y. The model's own force-transfer residual vanishes."""
    for delta_t in (-140.0, -50.0, 12.5, 50.0, 100.0):
        demand = solve_shear_lag(
            HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, HAND_ADHESIVE, delta_t
        )
        scale = max(demand.transferred_force, 1.0)
        assert abs(demand.force_transfer_residual) <= 1.0e-10 * scale


def test_average_shear_and_ratio_hand_calculation() -> None:
    """tau_avg = 5000/(0.02*0.02) = 12.5 MPa; peak/avg = 2 coth(2) = 2.0746."""
    assert DEMAND.average_transfer_shear == pytest.approx(HAND_AVERAGE_SHEAR, rel=1e-13)
    assert DEMAND.peak_to_average_ratio == pytest.approx(
        2.0 / math.tanh(2.0), rel=1e-12
    )
    assert DEMAND.peak_to_average_ratio == pytest.approx(2.0746294, rel=1e-6)


def test_average_shear_uses_the_whole_bond_area() -> None:
    """No factor-of-two ambiguity: tau_avg = |N_t| / (b L_b) over the full overlap."""
    assert DEMAND.average_transfer_shear == pytest.approx(
        DEMAND.transferred_force / HAND_GEOMETRY.bond_area, rel=1e-14
    )


def test_distribution_rejects_coordinates_outside_the_overlap() -> None:
    """Positions must lie inside [0, L_b]."""
    for bad in (-1.0e-6, LENGTH * 1.0001, 1.0, math.nan, math.inf):
        with pytest.raises(ValueError):
            DEMAND.shear_stress(bad)
        with pytest.raises(ValueError):
            DEMAND.member_1_force(bad)
        with pytest.raises(ValueError):
            DEMAND.member_2_force(bad)
    assert DEMAND.shear_stress(0.0) is not None
    assert DEMAND.shear_stress(LENGTH) is not None


def test_distribution_is_numerically_stable_for_a_very_long_overlap() -> None:
    """Overflow-safe hyperbolic ratios: a 1 m overlap (lambda = 100) still evaluates."""
    long_geometry = HAND_GEOMETRY.with_overlap_length(1.0)
    demand = solve_shear_lag(
        HAND_MEMBER_1, HAND_MEMBER_2, long_geometry, HAND_ADHESIVE, 50.0
    )
    assert demand.dimensionless_overlap == pytest.approx(100.0, rel=1e-12)
    assert math.isfinite(demand.peak_shear_stress)
    assert math.isfinite(demand.shear_stress(0.0))
    assert math.isfinite(demand.member_1_force(0.5))
    assert demand.member_1_force(0.0) == pytest.approx(HAND_SIGNED_FORCE, rel=1e-13)
    assert abs(demand.force_transfer_residual) <= 1.0e-9 * demand.transferred_force

"""Tests Z-AH: limiting cases and scaling laws for the shear-lag screen."""

from __future__ import annotations

import math

import pytest
from shear_lag_cases import HAND_ADHESIVE, HAND_GEOMETRY, HAND_MEMBER_1, HAND_MEMBER_2

from thermal_joint import (
    AxialMember,
    ThermoelasticMaterial,
    long_overlap_peak_shear_stress,
    shear_lag_parameter,
    solve_shear_lag,
)
from thermal_joint.illustrative import (
    ADHESIVE_LIKE,
    aluminium_like_member,
    radiator_joint_overlap,
    titanium_like_member,
)

MEMBER_AL = aluminium_like_member(100.0)
MEMBER_TI = titanium_like_member(100.0)
GEOMETRY = radiator_joint_overlap()


def demand_at(delta_t: float, geometry=HAND_GEOMETRY):
    return solve_shear_lag(HAND_MEMBER_1, HAND_MEMBER_2, geometry, HAND_ADHESIVE, delta_t)


def test_z_zero_delta_temperature_gives_a_zero_demand() -> None:
    """Z. dT = 0 means no mismatch force, no shear anywhere."""
    demand = demand_at(0.0)
    assert demand.mismatch_strain == 0.0
    assert demand.transferred_force == 0.0
    assert demand.peak_shear_stress == 0.0
    assert demand.average_transfer_shear == 0.0
    assert demand.shear_stress(0.0) == 0.0
    assert demand.shear_stress(HAND_GEOMETRY.overlap_length) == 0.0
    assert demand.member_1_force(0.0) == 0.0
    assert demand.force_transfer_residual == 0.0


@pytest.mark.parametrize("alpha", [0.0, 8.5e-6, 23.0e-6, -3.0e-6])
def test_aa_equal_cte_gives_a_zero_demand(alpha: float) -> None:
    """AA. Matched CTEs generate no mismatch force, so the bondline is unloaded."""
    member_1 = AxialMember(ThermoelasticMaterial("A", 70.0e9, alpha, 270.0e6), 100.0e-6)
    member_2 = AxialMember(ThermoelasticMaterial("B", 110.0e9, alpha, 830.0e6), 250.0e-6)
    for delta_t in (-140.0, 100.0):
        demand = solve_shear_lag(member_1, member_2, GEOMETRY, ADHESIVE_LIKE, delta_t)
        assert demand.mismatch_strain == 0.0
        assert demand.transferred_force == pytest.approx(0.0, abs=1e-9)
        assert demand.peak_shear_stress == pytest.approx(0.0, abs=1e-3)
        assert demand.shear_stress(GEOMETRY.overlap_length) == pytest.approx(0.0, abs=1e-3)


@pytest.mark.parametrize("delta_t", [30.0, 50.0, 140.0])
def test_ab_hot_and_cold_shear_reverse_sign(delta_t: float) -> None:
    """AB. Heating and cooling give equal and opposite shear distributions."""
    heating = demand_at(delta_t)
    cooling = demand_at(-delta_t)
    for x in (0.0, 0.007, HAND_GEOMETRY.overlap_length):
        assert cooling.shear_stress(x) == pytest.approx(-heating.shear_stress(x), rel=1e-13)
    assert cooling.peak_shear_stress == pytest.approx(heating.peak_shear_stress, rel=1e-13)


@pytest.mark.parametrize("factor", [0.5, 2.0, 3.0])
def test_ac_peak_shear_is_linear_in_delta_temperature(factor: float) -> None:
    """AC. tau_peak scales exactly with |dT| at fixed geometry."""
    base = demand_at(50.0)
    scaled = demand_at(50.0 * factor)
    assert scaled.transferred_force == pytest.approx(factor * base.transferred_force, rel=1e-13)
    assert scaled.peak_shear_stress == pytest.approx(factor * base.peak_shear_stress, rel=1e-13)
    assert scaled.average_transfer_shear == pytest.approx(
        factor * base.average_transfer_shear, rel=1e-13
    )
    assert scaled.peak_to_average_ratio == pytest.approx(base.peak_to_average_ratio, rel=1e-13)


@pytest.mark.parametrize("length", [1.0e-5, 1.0e-4, 5.0e-4])
def test_ad_short_overlap_peak_to_average_approaches_one(length: float) -> None:
    """AD. lambda coth(lambda) -> 1: a very short overlap shears almost uniformly."""
    demand = demand_at(50.0, HAND_GEOMETRY.with_overlap_length(length))
    assert demand.dimensionless_overlap < 0.06
    assert demand.peak_to_average_ratio == pytest.approx(1.0, abs=2.0e-3)
    assert demand.peak_to_average_ratio >= 1.0


def test_ad_short_overlap_shear_is_nearly_flat_across_the_bond() -> None:
    """AD. The distribution itself flattens, not just the ratio."""
    demand = demand_at(50.0, HAND_GEOMETRY.with_overlap_length(1.0e-4))
    values = [abs(demand.shear_stress(1.0e-4 * i / 10)) for i in range(11)]
    assert max(values) / min(values) == pytest.approx(1.0, abs=1.0e-3)


def test_ae_peak_to_average_ratio_increases_with_overlap() -> None:
    """AE. Longer overlaps concentrate shear at the free edge."""
    ratios = [
        demand_at(50.0, HAND_GEOMETRY.with_overlap_length(length)).peak_to_average_ratio
        for length in (1.0e-4, 1.0e-3, 5.0e-3, 0.02, 0.08, 0.32)
    ]
    assert ratios == sorted(ratios)
    assert ratios[0] == pytest.approx(1.0, abs=1e-3)
    assert ratios[-1] > 20.0


def test_af_long_overlap_peak_approaches_a_finite_nonzero_asymptote() -> None:
    """AF. tau_peak -> |N_t| beta / b from above; it never decays to zero."""
    beta = shear_lag_parameter(HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, HAND_ADHESIVE)
    reference = demand_at(50.0)
    asymptote = long_overlap_peak_shear_stress(
        reference.transferred_force, beta, HAND_GEOMETRY.bond_width
    )
    assert asymptote == pytest.approx(25.0e6, rel=1e-12)
    assert asymptote > 0.0

    peaks = [
        demand_at(50.0, HAND_GEOMETRY.with_overlap_length(length)).peak_shear_stress
        for length in (0.02, 0.05, 0.1, 0.5, 1.0)
    ]
    assert peaks == sorted(peaks, reverse=True)
    # Approached strictly from above; at lambda = 100 coth() saturates to 1.0 in
    # double precision, so the longest overlap lands exactly on the asymptote.
    assert all(peak >= asymptote for peak in peaks)
    assert peaks[0] > asymptote
    assert peaks[-1] == pytest.approx(asymptote, rel=1e-12)


def test_af_asymptote_is_analytic_not_a_long_numerical_extrapolation() -> None:
    """AF. The asymptote helper is a closed form, independent of any overlap length."""
    beta = shear_lag_parameter(HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, HAND_ADHESIVE)
    values = {
        round(
            long_overlap_peak_shear_stress(5000.0, beta, HAND_GEOMETRY.bond_width), 6
        )
    }
    assert len(values) == 1
    assert long_overlap_peak_shear_stress(
        5000.0, beta, HAND_GEOMETRY.bond_width
    ) == pytest.approx(5000.0 * beta / HAND_GEOMETRY.bond_width, rel=1e-15)


def test_ag_added_overlap_gives_diminishing_peak_shear_returns() -> None:
    """AG. Each doubling of overlap buys less peak-shear reduction than the last."""
    lengths = [0.005, 0.010, 0.020, 0.040, 0.080, 0.160]
    peaks = [
        demand_at(50.0, HAND_GEOMETRY.with_overlap_length(length)).peak_shear_stress
        for length in lengths
    ]
    reductions = [peaks[i] - peaks[i + 1] for i in range(len(peaks) - 1)]
    assert all(reduction > 0.0 for reduction in reductions)
    assert reductions == sorted(reductions, reverse=True)
    # Beyond a few transfer lengths the benefit is essentially exhausted.
    assert reductions[-1] < 1.0e-4 * reductions[0]


def test_ag_average_shear_keeps_falling_like_one_over_length() -> None:
    """AG. Average shear does decay as 1/L even though the peak does not."""
    for length in (0.02, 0.08, 0.32):
        demand = demand_at(50.0, HAND_GEOMETRY.with_overlap_length(length))
        assert demand.average_transfer_shear == pytest.approx(
            demand.transferred_force / (HAND_GEOMETRY.bond_width * length), rel=1e-13
        )
    short = demand_at(50.0, HAND_GEOMETRY.with_overlap_length(0.02))
    quadruple = demand_at(50.0, HAND_GEOMETRY.with_overlap_length(0.08))
    assert quadruple.average_transfer_shear == pytest.approx(
        short.average_transfer_shear / 4.0, rel=1e-12
    )


def test_ah_cold_to_hot_peak_ratio_is_1_point_4_for_the_canonical_environment() -> None:
    """AH. |dT_cold|/|dT_hot| = 140/100 = 1.4 carries straight into the peak shear."""
    hot = solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, 100.0)
    cold = solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, -140.0)
    assert cold.peak_shear_stress / hot.peak_shear_stress == pytest.approx(1.4, rel=1e-12)
    assert cold.transferred_force / hot.transferred_force == pytest.approx(1.4, rel=1e-12)


def test_canonical_illustrative_peak_shear_values() -> None:
    """Regression lock on the canonical illustrative overlap."""
    hot = solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, 100.0)
    cold = solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, -140.0)

    assert hot.beta == pytest.approx(152.894157, rel=1e-8)
    assert hot.transfer_length == pytest.approx(6.5405e-3, rel=1e-4)
    assert hot.dimensionless_overlap == pytest.approx(6.115766, rel=1e-6)
    assert hot.transferred_force == pytest.approx(6202.7778, rel=1e-6)
    assert cold.transferred_force == pytest.approx(8683.8889, rel=1e-6)
    assert hot.peak_shear_stress == pytest.approx(47.4189e6, rel=1e-5)
    assert cold.peak_shear_stress == pytest.approx(66.3864e6, rel=1e-5)
    assert hot.average_transfer_shear == pytest.approx(7.7535e6, rel=1e-4)
    assert cold.average_transfer_shear == pytest.approx(10.8549e6, rel=1e-4)
    assert cold.peak_to_average_ratio == pytest.approx(6.1158, rel=1e-4)


def test_peak_to_average_ratio_equals_lambda_coth_lambda() -> None:
    """The concentration factor is exactly lambda coth(lambda)."""
    for length in (1.0e-4, 0.005, 0.02, 0.16):
        demand = demand_at(50.0, HAND_GEOMETRY.with_overlap_length(length))
        lam = demand.dimensionless_overlap
        assert demand.peak_to_average_ratio == pytest.approx(
            lam / math.tanh(lam), rel=1e-12
        )

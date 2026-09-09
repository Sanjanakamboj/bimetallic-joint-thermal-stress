"""Tests G-AE: the one-variable design sensitivity sweeps."""

from __future__ import annotations

import math

import pytest
from design_cases import (
    ADHESIVE,
    ADHESIVE_BASIS,
    ALLOWABLE,
    BASELINE_COLD_MARGIN,
    BASELINE_COLD_PEAK,
    BASELINE_YIELD_MARGIN,
    COLD,
    ENVIRONMENT,
    GEOMETRY,
    MEMBER_AL,
    MEMBER_TI,
    MM,
    YIELD_BASIS,
)

from thermal_joint import (
    adhesive_modulus_sensitivity,
    adhesive_thickness_sensitivity,
    area_ratio_sensitivity,
    bond_width_sensitivity,
    cte_mismatch_sensitivity,
    evaluate_joint_design,
    solve_shear_lag,
    thermal_excursion_sensitivity,
)
from thermal_joint.illustrative import ADHESIVE_LIKE, ALUMINIUM_LIKE, TITANIUM_LIKE

WIDTHS_MM = (10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0)
THICKNESSES_MM = (0.05, 0.10, 0.20, 0.30, 0.50, 0.75, 1.00, 1.50)
MODULI_GPA = (0.05, 0.10, 0.20, 0.50, 1.0, 2.0, 5.0)
AREA_RATIOS = (0.25, 0.5, 1.0, 2.0, 4.0)
EXCURSIONS = (20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0)
MISMATCH_SCALES = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5)


def widths():
    return bond_width_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE,
        [w * MM for w in WIDTHS_MM], YIELD_BASIS, ADHESIVE_BASIS,
    )


def thicknesses():
    return adhesive_thickness_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE,
        [t * MM for t in THICKNESSES_MM], YIELD_BASIS, ADHESIVE_BASIS,
    )


def moduli():
    return adhesive_modulus_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE,
        [g * 1.0e9 for g in MODULI_GPA], YIELD_BASIS, ADHESIVE_BASIS,
    )


# --------------------------------------------------------------- width G-K

def test_g_width_increase_lowers_peak_shear() -> None:
    """G. Peak shear falls monotonically with bond width."""
    peaks = [candidate.peak_shear_stress for candidate in widths()]
    assert peaks == sorted(peaks, reverse=True)
    assert peaks[0] > peaks[-1]


def test_h_width_increase_raises_beta_as_sqrt() -> None:
    """H. beta ~ sqrt(b), so width raises beta even as it lowers peak shear."""
    candidates = widths()
    for candidate in candidates:
        expected = candidates[2].beta * math.sqrt(
            candidate.bond_width / candidates[2].bond_width
        )
        assert candidate.beta == pytest.approx(expected, rel=1e-12)
    betas = [candidate.beta for candidate in candidates]
    assert betas == sorted(betas)


def test_i_width_sweep_is_deterministic() -> None:
    """I. Repeating the sweep reproduces identical records."""
    assert widths() == widths()


def test_j_width_sweep_does_not_mutate_canonical_inputs() -> None:
    """J. The shipped geometry, members and adhesive are untouched afterwards."""
    widths()
    assert GEOMETRY.bond_width == pytest.approx(20.0 * MM, rel=1e-12)
    assert GEOMETRY.adhesive_thickness == pytest.approx(0.2 * MM, rel=1e-12)
    assert ADHESIVE_LIKE.shear_modulus == 1.0e9
    assert ADHESIVE_LIKE.shear_strength == 25.0e6
    assert ALUMINIUM_LIKE.thermal_expansion_coefficient == 23.0e-6


@pytest.mark.parametrize("factor", [4.0, 9.0, 16.0])
def test_k_long_overlap_width_power_law_holds_across_scaled_widths(factor: float) -> None:
    """K. At lambda >> 1 the finite-overlap peak follows b^(-1/2) too."""
    base = evaluate_joint_design(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, YIELD_BASIS, ADHESIVE_BASIS
    )
    scaled = bond_width_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE,
        [GEOMETRY.bond_width * factor], YIELD_BASIS, ADHESIVE_BASIS,
    )[0]
    assert scaled.peak_shear_stress == pytest.approx(
        base.peak_shear_stress / math.sqrt(factor), rel=2e-4
    )
    assert scaled.long_overlap_asymptote == pytest.approx(
        base.long_overlap_asymptote / math.sqrt(factor), rel=1e-12
    )


def test_width_sweep_over_the_canonical_range_never_passes() -> None:
    """The 10-100 mm range cannot recover the canonical joint; stated explicitly."""
    assert all(not candidate.adhesive_feasible for candidate in widths())
    assert all(candidate.yield_feasible for candidate in widths())


# ----------------------------------------------------------- thickness L-P

def test_l_thicker_adhesive_lowers_beta() -> None:
    """L. beta ~ t_a^(-1/2)."""
    betas = [candidate.beta for candidate in thicknesses()]
    assert betas == sorted(betas, reverse=True)


def test_m_thicker_adhesive_increases_transfer_length() -> None:
    """M. Transfer length 1/beta grows with bondline thickness."""
    lengths = [candidate.transfer_length for candidate in thicknesses()]
    assert lengths == sorted(lengths)


def test_n_thicker_adhesive_lowers_peak_shear_toward_the_floor() -> None:
    """N. Peak shear falls monotonically but stays above the uniform-shear floor."""
    candidates = thicknesses()
    peaks = [candidate.peak_shear_stress for candidate in candidates]
    assert peaks == sorted(peaks, reverse=True)
    for candidate in candidates:
        assert candidate.peak_shear_stress > candidate.uniform_shear_floor


def test_o_thickness_sweep_is_deterministic() -> None:
    assert thicknesses() == thicknesses()


def test_p_baseline_thickness_reproduces_milestone_3() -> None:
    """P. The 0.2 mm point is exactly the Milestone 3 canonical result."""
    baseline = [c for c in thicknesses() if c.adhesive_thickness == pytest.approx(0.2 * MM)][0]
    assert baseline.cold_peak_shear_stress == pytest.approx(BASELINE_COLD_PEAK, rel=1e-5)
    assert baseline.adhesive_margin == pytest.approx(BASELINE_COLD_MARGIN, rel=1e-3)
    assert baseline.beta == pytest.approx(152.894157, rel=1e-8)


def test_thickness_sweep_over_the_canonical_range_never_passes() -> None:
    """Thickness alone up to 1.5 mm cannot recover the canonical joint at 1 GPa."""
    assert all(not candidate.adhesive_feasible for candidate in thicknesses())


# ------------------------------------------------------------ modulus Q-U

def test_q_lower_modulus_lowers_beta() -> None:
    """Q. beta ~ sqrt(G_a)."""
    betas = [candidate.beta for candidate in moduli()]
    assert betas == sorted(betas)


def test_r_lower_modulus_increases_transfer_length() -> None:
    lengths = [candidate.transfer_length for candidate in moduli()]
    assert lengths == sorted(lengths, reverse=True)


def test_s_lower_modulus_reduces_peak_shear() -> None:
    peaks = [candidate.peak_shear_stress for candidate in moduli()]
    assert peaks == sorted(peaks)


def test_t_baseline_modulus_reproduces_milestone_3() -> None:
    """T. The 1 GPa point is exactly the Milestone 3 canonical result."""
    baseline = [c for c in moduli() if c.adhesive_shear_modulus == pytest.approx(1.0e9)][0]
    assert baseline.cold_peak_shear_stress == pytest.approx(BASELINE_COLD_PEAK, rel=1e-5)
    assert baseline.adhesive_margin == pytest.approx(BASELINE_COLD_MARGIN, rel=1e-3)


def test_u_point_one_gpa_alone_still_fails_at_the_canonical_thickness() -> None:
    """U. Milestone 3's finding is locked: 0.1 GPa at 0.2 mm still fails."""
    soft = [c for c in moduli() if c.adhesive_shear_modulus == pytest.approx(0.1e9)][0]
    assert soft.peak_shear_stress == pytest.approx(21.8893e6, rel=1e-4)
    assert soft.adhesive_margin == pytest.approx(-0.0863, abs=1e-3)
    assert not soft.adhesive_feasible


def test_modulus_sweep_shear_strength_is_never_modified() -> None:
    """The sweep varies stiffness only; the allowable is constant throughout."""
    for candidate in moduli():
        assert candidate.allowable_shear_stress == pytest.approx(ALLOWABLE, rel=1e-15)


def test_modulus_sweep_is_deterministic() -> None:
    assert moduli() == moduli()


# ------------------------------------------------- temperature and CTE V-Z

def test_v_peak_shear_is_linear_in_excursion_magnitude() -> None:
    """V. Both the transferred force and the peak shear scale exactly with |dT|."""
    points = thermal_excursion_sensitivity(
        MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE, EXCURSIONS,
        yield_basis=YIELD_BASIS, adhesive_basis=ADHESIVE_BASIS,
    )
    reference = points[EXCURSIONS.index(20.0)]
    for point, magnitude in zip(points, EXCURSIONS):
        factor = magnitude / 20.0
        assert point.peak_shear_stress == pytest.approx(
            factor * reference.peak_shear_stress, rel=1e-12
        )
        assert point.transferred_force == pytest.approx(
            factor * reference.transferred_force, rel=1e-12
        )


def test_w_allowable_delta_temperature_appears_as_the_sweep_boundary() -> None:
    """W. The 42.18 K allowable shows up as the pass/fail crossing of the sweep."""
    points = thermal_excursion_sensitivity(
        MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE, EXCURSIONS,
        yield_basis=YIELD_BASIS, adhesive_basis=ADHESIVE_BASIS,
    )
    passing = [p.value for p in points if p.adhesive_feasible]
    failing = [p.value for p in points if not p.adhesive_feasible]
    assert max(passing) == 40.0
    assert min(failing) == 60.0

    from thermal_joint import allowable_temperature_change_for_adhesive_shear

    allowable = allowable_temperature_change_for_adhesive_shear(
        MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE, ADHESIVE_BASIS
    )
    assert allowable == pytest.approx(42.1773, rel=1e-5)
    assert max(passing) < allowable < min(failing)


def test_x_zero_cte_mismatch_gives_zero_shear() -> None:
    """X. A scale factor of zero removes the mismatch entirely."""
    points = cte_mismatch_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, MISMATCH_SCALES,
        YIELD_BASIS, ADHESIVE_BASIS,
    )
    zero = points[0]
    assert zero.value == 0.0
    assert zero.transferred_force == pytest.approx(0.0, abs=1e-9)
    assert zero.peak_shear_stress == pytest.approx(0.0, abs=1e-3)
    assert math.isinf(zero.adhesive_margin)
    assert zero.adhesive_feasible


def test_y_peak_shear_is_linear_in_the_cte_mismatch() -> None:
    """Y. tau_peak scales exactly with |alpha_1 - alpha_2| at fixed E, A, geometry, dT."""
    points = cte_mismatch_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, MISMATCH_SCALES,
        YIELD_BASIS, ADHESIVE_BASIS,
    )
    unit = [p for p in points if p.value == 1.0][0]
    for point in points:
        assert point.peak_shear_stress == pytest.approx(
            point.value * unit.peak_shear_stress, rel=1e-11, abs=1e-3
        )
    assert unit.peak_shear_stress == pytest.approx(BASELINE_COLD_PEAK, rel=1e-5)


def test_z_synthetic_materials_do_not_mutate_canonical_records() -> None:
    """Z. The CTE sweep builds copies; the shipped materials are unchanged."""
    cte_mismatch_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, MISMATCH_SCALES,
        YIELD_BASIS, ADHESIVE_BASIS,
    )
    assert ALUMINIUM_LIKE.thermal_expansion_coefficient == 23.0e-6
    assert TITANIUM_LIKE.thermal_expansion_coefficient == 8.5e-6
    assert MEMBER_AL.material.thermal_expansion_coefficient == 23.0e-6
    assert MEMBER_AL.area == pytest.approx(100.0e-6, rel=1e-12)


# ------------------------------------------------------------ area ratio AA-AE

def ratios():
    return area_ratio_sensitivity(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, AREA_RATIOS,
        YIELD_BASIS, ADHESIVE_BASIS,
    )


def test_aa_area_ratio_sweep_is_deterministic() -> None:
    assert ratios() == ratios()


def test_ab_transfer_force_is_recomputed_at_every_ratio() -> None:
    """AB. The Milestone 1 demand is re-derived for each area, never held fixed."""
    from thermal_joint import AxialMember, solve_bimetallic_joint

    points = ratios()
    forces = [point.transferred_force for point in points]
    assert len(set(forces)) == len(forces)
    for point, ratio in zip(points, AREA_RATIOS):
        member = AxialMember(ALUMINIUM_LIKE, ratio * MEMBER_TI.area)
        expected = abs(solve_bimetallic_joint(member, MEMBER_TI, COLD).member_1_force)
        assert point.transferred_force == pytest.approx(expected, rel=1e-12)
    assert forces == sorted(forces)


def test_ac_beta_is_recomputed_at_every_ratio() -> None:
    """AC. A larger aluminium-like area lowers C and therefore beta."""
    betas = [point.beta for point in ratios()]
    assert betas == sorted(betas, reverse=True)
    assert len(set(betas)) == len(betas)


def test_ac_lowest_peak_shear_is_at_the_smallest_area_ratio() -> None:
    """AC. Equal areas are NOT optimal for adhesive shear - determined, not assumed."""
    points = ratios()
    peaks = [point.peak_shear_stress for point in points]
    assert peaks == sorted(peaks)
    best = min(points, key=lambda point: point.peak_shear_stress)
    assert best.value == 0.25
    assert best.value != 1.0


def test_ad_yield_margins_are_reported_separately_at_every_ratio() -> None:
    """AD. The member yield margin is carried alongside, never merged."""
    for point in ratios():
        assert math.isfinite(point.minimum_yield_margin)
        assert point.minimum_yield_margin != point.adhesive_margin
        assert point.yield_feasible
        assert point.overall_feasible == (point.yield_feasible and point.adhesive_feasible)


def test_ae_equal_area_point_reproduces_milestone_3() -> None:
    """AE. The ratio-1 point is the Milestone 3 canonical case."""
    unit = [point for point in ratios() if point.value == 1.0][0]
    assert unit.peak_shear_stress == pytest.approx(BASELINE_COLD_PEAK, rel=1e-5)
    assert unit.adhesive_margin == pytest.approx(BASELINE_COLD_MARGIN, rel=1e-3)
    assert unit.minimum_yield_margin == pytest.approx(BASELINE_YIELD_MARGIN, rel=1e-3)
    assert unit.governing_extreme == "cold"


def test_all_sweeps_report_cold_as_governing_for_the_canonical_environment() -> None:
    """Section 27: computed from the margins at every point, never assumed."""
    for candidate in widths() + thicknesses() + moduli():
        assert candidate.governing_extreme == "cold"
        assert candidate.cold_peak_shear_stress > candidate.hot_peak_shear_stress
    for point in ratios():
        assert point.governing_extreme == "cold"

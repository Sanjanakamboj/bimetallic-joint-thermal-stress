"""Tests A-F: the analytical long-overlap scaling laws."""

from __future__ import annotations

import math

import pytest
from design_cases import (
    ADHESIVE,
    ADHESIVE_BASIS,
    BASELINE_ASYMPTOTE,
    BASELINE_FLOOR,
    COLD,
    ENVIRONMENT,
    GEOMETRY,
    MEMBER_AL,
    MEMBER_TI,
    YIELD_BASIS,
)

from thermal_joint import (
    AdhesiveMaterial,
    BondedOverlapGeometry,
    adherend_compliance_sum,
    long_overlap_peak_shear,
    long_overlap_peak_shear_stress,
    shear_lag_parameter,
    solve_shear_lag,
    thermal_mismatch_strain,
    uniform_shear_floor,
)


def exponent(scaled: float, base: float, factor: float) -> float:
    """Empirical power-law exponent from a single factor change."""
    return math.log(scaled / base) / math.log(factor)


def asymptote_with(width=None, thickness=None, modulus=None, delta_t=COLD) -> float:
    geometry = BondedOverlapGeometry(
        overlap_length=GEOMETRY.overlap_length,
        bond_width=GEOMETRY.bond_width if width is None else width,
        adhesive_thickness=(
            GEOMETRY.adhesive_thickness if thickness is None else thickness
        ),
    )
    adhesive = ADHESIVE if modulus is None else AdhesiveMaterial(
        "variant", modulus, ADHESIVE.shear_strength, "test fixture"
    )
    return long_overlap_peak_shear(MEMBER_AL, MEMBER_TI, geometry, adhesive, delta_t)


@pytest.mark.parametrize("factor", [0.25, 4.0, 9.0])
def test_a_asymptotic_width_scaling_is_inverse_square_root(factor: float) -> None:
    """A. tau_inf ~ b^(-1/2). NOT 1/b - the width also raises beta."""
    base = asymptote_with()
    scaled = asymptote_with(width=GEOMETRY.bond_width * factor)
    assert exponent(scaled, base, factor) == pytest.approx(-0.5, abs=1e-12)
    assert scaled != pytest.approx(base / factor, rel=1e-3)  # explicitly not 1/b


@pytest.mark.parametrize("factor", [0.25, 4.0, 9.0])
def test_b_asymptotic_thickness_scaling_is_inverse_square_root(factor: float) -> None:
    """B. tau_inf ~ t_a^(-1/2)."""
    base = asymptote_with()
    scaled = asymptote_with(thickness=GEOMETRY.adhesive_thickness * factor)
    assert exponent(scaled, base, factor) == pytest.approx(-0.5, abs=1e-12)


@pytest.mark.parametrize("factor", [0.25, 4.0, 9.0])
def test_c_asymptotic_modulus_scaling_is_square_root(factor: float) -> None:
    """C. tau_inf ~ G_a^(+1/2)."""
    base = asymptote_with()
    scaled = asymptote_with(modulus=ADHESIVE.shear_modulus * factor)
    assert exponent(scaled, base, factor) == pytest.approx(+0.5, abs=1e-12)


@pytest.mark.parametrize("factor", [0.5, 2.0, 4.0])
def test_d_transfer_force_and_asymptote_are_linear_in_delta_temperature(factor: float) -> None:
    """D. tau_inf ~ |d_eps|^1, so linear in dT and in the CTE mismatch."""
    base = asymptote_with()
    scaled = asymptote_with(delta_t=COLD * factor)
    assert exponent(scaled, base, factor) == pytest.approx(+1.0, abs=1e-12)

    demand = solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE, COLD)
    doubled = solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE, 2.0 * COLD)
    assert doubled.transferred_force == pytest.approx(2.0 * demand.transferred_force, rel=1e-13)


def test_d_asymptote_matches_the_closed_form_from_raw_scalars() -> None:
    """D. tau_inf = |d_eps| sqrt(G_a / (b t_a C)), evaluated independently."""
    expected = abs(thermal_mismatch_strain(MEMBER_AL, MEMBER_TI, COLD)) * math.sqrt(
        ADHESIVE.shear_modulus
        / (
            GEOMETRY.bond_width
            * GEOMETRY.adhesive_thickness
            * adherend_compliance_sum(MEMBER_AL, MEMBER_TI)
        )
    )
    assert asymptote_with() == pytest.approx(expected, rel=1e-14)


def test_e_helper_matches_the_milestone_3_primitive_and_a_very_long_overlap() -> None:
    """E. The M4 wrapper agrees with the M3 primitive and with numerics at L -> inf."""
    beta = shear_lag_parameter(MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE)
    demand = solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE, COLD)
    primitive = long_overlap_peak_shear_stress(
        demand.transferred_force, beta, GEOMETRY.bond_width
    )
    assert asymptote_with() == pytest.approx(primitive, rel=1e-12)

    for length in (0.2, 1.0, 5.0):
        long_demand = solve_shear_lag(
            MEMBER_AL, MEMBER_TI, GEOMETRY.with_overlap_length(length), ADHESIVE, COLD
        )
        assert long_demand.peak_shear_stress >= asymptote_with()
    assert long_demand.peak_shear_stress == pytest.approx(asymptote_with(), rel=1e-12)


def test_f_canonical_asymptote_and_floor_regression() -> None:
    """F. Canonical cold asymptote 66.386 MPa and uniform-shear floor 10.855 MPa."""
    assert asymptote_with() == pytest.approx(BASELINE_ASYMPTOTE, rel=1e-5)
    floor = uniform_shear_floor(MEMBER_AL, MEMBER_TI, GEOMETRY, COLD)
    assert floor == pytest.approx(BASELINE_FLOOR, rel=1e-5)
    assert floor < asymptote_with()


def test_f_uniform_shear_floor_is_the_beta_to_zero_limit() -> None:
    """F. A vanishingly stiff bondline drives peak shear onto the floor."""
    floor = uniform_shear_floor(MEMBER_AL, MEMBER_TI, GEOMETRY, COLD)
    very_soft = AdhesiveMaterial("very soft", 1.0e3, ADHESIVE.shear_strength, "fixture")
    demand = solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, very_soft, COLD)
    assert demand.peak_shear_stress == pytest.approx(floor, rel=1e-3)
    assert demand.peak_shear_stress >= floor


def test_f_floor_falls_as_one_over_width_while_asymptote_falls_as_sqrt() -> None:
    """F. The two limits scale differently - that is why width has no floor."""
    wide = BondedOverlapGeometry(GEOMETRY.overlap_length, 4.0 * GEOMETRY.bond_width,
                                 GEOMETRY.adhesive_thickness)
    base_floor = uniform_shear_floor(MEMBER_AL, MEMBER_TI, GEOMETRY, COLD)
    wide_floor = uniform_shear_floor(MEMBER_AL, MEMBER_TI, wide, COLD)
    assert wide_floor == pytest.approx(base_floor / 4.0, rel=1e-12)
    assert asymptote_with(width=4.0 * GEOMETRY.bond_width) == pytest.approx(
        asymptote_with() / 2.0, rel=1e-12
    )

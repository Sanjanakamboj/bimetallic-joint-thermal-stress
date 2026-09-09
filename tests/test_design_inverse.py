"""Tests AF-AW: inverse sizing for bond width, adhesive thickness and modulus."""

from __future__ import annotations

import math

import pytest
from design_cases import (
    ADHESIVE,
    ADHESIVE_BASIS,
    ALLOWABLE,
    BASELINE_FLOOR,
    ENVIRONMENT,
    GEOMETRY,
    MEMBER_AL,
    MEMBER_TI,
    MM,
)

from thermal_joint import (
    AdhesiveMaterial,
    AdhesiveModulusStatus,
    AdhesiveShearBasis,
    AdhesiveThicknessStatus,
    BondWidthStatus,
    BondedOverlapGeometry,
    assess_shear_lag_extremes,
    maximum_allowable_adhesive_shear_modulus,
    required_adhesive_thickness,
    required_bond_width,
)


def peak_at(width=None, thickness=None, modulus=None) -> float:
    geometry = BondedOverlapGeometry(
        overlap_length=GEOMETRY.overlap_length,
        bond_width=GEOMETRY.bond_width if width is None else width,
        adhesive_thickness=GEOMETRY.adhesive_thickness if thickness is None else thickness,
    )
    adhesive = ADHESIVE if modulus is None else AdhesiveMaterial(
        "variant", modulus, ADHESIVE.shear_strength, "test fixture"
    )
    return assess_shear_lag_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, geometry, adhesive, ADHESIVE_BASIS
    ).governing_peak_shear_stress


# ------------------------------------------------------------- width AF-AL

def test_af_lower_bound_already_passes() -> None:
    """AF. Starting the search from an already-wide bond reports it immediately."""
    result = required_bond_width(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        minimum_bond_width=0.5, maximum_bond_width=2.0,
    )
    assert result.status is BondWidthStatus.LOWER_BOUND_ALREADY_PASSES
    assert result.required_bond_width == 0.5
    assert result.margin_at_required_width >= 0.0
    assert result.iterations == 0


def test_ag_finite_required_width_for_the_canonical_case() -> None:
    """AG. The canonical joint needs about 220 mm of bond width."""
    result = required_bond_width(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        tolerance=1.0e-9,
    )
    assert result.status is BondWidthStatus.FINITE_REQUIRED_WIDTH
    assert result.has_finite_width
    assert result.governing_extreme == "cold"
    assert result.required_bond_width == pytest.approx(0.220354, rel=1e-4)
    assert result.required_bond_width / GEOMETRY.bond_width == pytest.approx(11.02, rel=1e-3)


def test_ag_required_width_matches_the_analytical_asymptotic_estimate() -> None:
    """AG. b = |d_eps|^2 G_a / (t_a C tau_allow^2) up to the coth correction."""
    from thermal_joint import adherend_compliance_sum, thermal_mismatch_strain

    mismatch = abs(thermal_mismatch_strain(MEMBER_AL, MEMBER_TI, -140.0))
    compliance = adherend_compliance_sum(MEMBER_AL, MEMBER_TI)
    analytic = (
        mismatch**2
        * ADHESIVE.shear_modulus
        / (GEOMETRY.adhesive_thickness * compliance * ALLOWABLE**2)
    )
    result = required_bond_width(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        tolerance=1.0e-9,
    )
    assert analytic == pytest.approx(0.220353, rel=1e-4)
    assert result.required_bond_width == pytest.approx(analytic, rel=1e-3)


def test_ah_returned_width_gives_a_near_zero_margin() -> None:
    """AH. MS is approximately zero and non-negative at the returned width."""
    result = required_bond_width(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        tolerance=1.0e-9,
    )
    assert result.margin_at_required_width == pytest.approx(0.0, abs=1e-6)
    assert result.margin_at_required_width >= 0.0
    assert peak_at(width=result.required_bond_width) == pytest.approx(ALLOWABLE, rel=1e-6)


def test_ai_just_below_the_boundary_fails() -> None:
    """AI. A slightly narrower bond exceeds the allowable."""
    result = required_bond_width(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        tolerance=1.0e-9,
    )
    assert peak_at(width=result.required_bond_width * 0.99) > ALLOWABLE


def test_aj_just_above_the_boundary_passes() -> None:
    """AJ. A slightly wider bond meets the allowable."""
    result = required_bond_width(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        tolerance=1.0e-9,
    )
    assert peak_at(width=result.required_bond_width * 1.01) < ALLOWABLE


def test_ak_insufficient_upper_bound_is_reported_not_expanded() -> None:
    """AK. A too-narrow bracket is reported; bounds are never widened silently."""
    result = required_bond_width(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        minimum_bond_width=1.0 * MM, maximum_bond_width=0.1,
    )
    assert result.status is BondWidthStatus.NO_BOUNDARY_WITHIN_SEARCH_BOUNDS
    assert result.required_bond_width is None
    assert result.maximum_bond_width == 0.1
    assert result.peak_shear_at_upper_bound > result.allowable_shear_stress


def test_al_width_search_is_deterministic() -> None:
    """AL. Repeating the search reproduces the identical result object."""
    kwargs = dict(minimum_bond_width=1.0 * MM, maximum_bond_width=2.0, tolerance=1.0e-8)
    first = required_bond_width(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS, **kwargs
    )
    second = required_bond_width(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS, **kwargs
    )
    assert first == second


def test_width_search_never_mutates_the_supplied_geometry() -> None:
    geometry = BondedOverlapGeometry(0.040, 0.020, 0.0002)
    required_bond_width(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, geometry, ADHESIVE, ADHESIVE_BASIS
    )
    assert geometry == BondedOverlapGeometry(0.040, 0.020, 0.0002)


# --------------------------------------------------------- thickness AM-AQ

def test_am_thickness_lower_bound_already_passes() -> None:
    """AM. A bracket that starts thick enough reports immediately."""
    result = required_adhesive_thickness(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        minimum_adhesive_thickness=0.010, maximum_adhesive_thickness=0.05,
    )
    assert result.status is AdhesiveThicknessStatus.LOWER_BOUND_ALREADY_PASSES
    assert result.required_adhesive_thickness == 0.010
    assert result.iterations == 0


def test_an_finite_required_thickness_for_the_canonical_case() -> None:
    """AN. The canonical joint needs about a 2.5 mm bondline at 1 GPa."""
    result = required_adhesive_thickness(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS
    )
    assert result.status is AdhesiveThicknessStatus.FINITE_REQUIRED_THICKNESS
    assert result.has_finite_thickness
    assert result.required_adhesive_thickness == pytest.approx(2.4987 * MM, rel=1e-3)
    assert result.required_adhesive_thickness / GEOMETRY.adhesive_thickness == pytest.approx(
        12.49, rel=1e-2
    )
    assert result.uniform_shear_floor == pytest.approx(BASELINE_FLOOR, rel=1e-4)
    assert result.uniform_shear_floor < result.allowable_shear_stress


def test_ao_thickness_boundary_round_trip() -> None:
    """AO. MS is approximately zero there; thinner fails, thicker passes."""
    result = required_adhesive_thickness(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS
    )
    thickness = result.required_adhesive_thickness
    assert result.margin_at_required_thickness == pytest.approx(0.0, abs=1e-6)
    assert result.margin_at_required_thickness >= 0.0
    assert peak_at(thickness=thickness * 0.99) > ALLOWABLE
    assert peak_at(thickness=thickness * 1.01) < ALLOWABLE


def test_no_finite_thickness_when_the_uniform_shear_floor_exceeds_the_allowable() -> None:
    """The floor makes thickness useless below a certain bond width - analytically."""
    narrow = BondedOverlapGeometry(GEOMETRY.overlap_length, 5.0 * MM, GEOMETRY.adhesive_thickness)
    result = required_adhesive_thickness(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, narrow, ADHESIVE, ADHESIVE_BASIS
    )
    assert result.status is AdhesiveThicknessStatus.NO_FINITE_THICKNESS_WITHIN_MODEL
    assert result.required_adhesive_thickness is None
    assert result.uniform_shear_floor > result.allowable_shear_stress
    assert result.iterations == 0


def test_ap_thickness_insufficient_upper_bound() -> None:
    """AP. A bracket ending below the boundary is reported, not expanded."""
    result = required_adhesive_thickness(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        minimum_adhesive_thickness=1.0e-5, maximum_adhesive_thickness=1.0 * MM,
    )
    assert result.status is AdhesiveThicknessStatus.NO_BOUNDARY_WITHIN_SEARCH_BOUNDS
    assert result.required_adhesive_thickness is None
    assert result.uniform_shear_floor < result.allowable_shear_stress


def test_aq_thickness_search_is_deterministic() -> None:
    args = (MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS)
    assert required_adhesive_thickness(*args) == required_adhesive_thickness(*args)


# ----------------------------------------------------------- modulus AR-AW

def test_ar_finite_maximum_modulus_for_the_canonical_case() -> None:
    """AR. At 20 mm / 0.2 mm the adhesive must be softer than about 0.080 GPa."""
    result = maximum_allowable_adhesive_shear_modulus(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        tolerance=1.0,
    )
    assert result.status is AdhesiveModulusStatus.FINITE_MAXIMUM_MODULUS
    assert result.has_finite_modulus
    assert result.maximum_shear_modulus == pytest.approx(0.080043e9, rel=1e-3)
    assert result.maximum_shear_modulus < ADHESIVE.shear_modulus


def test_as_modulus_boundary_gives_a_near_zero_margin() -> None:
    """AS. MS is approximately zero and non-negative at the returned modulus."""
    result = maximum_allowable_adhesive_shear_modulus(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        tolerance=1.0,
    )
    assert result.margin_at_maximum_modulus == pytest.approx(0.0, abs=1e-6)
    assert result.margin_at_maximum_modulus >= 0.0


def test_at_modulus_just_above_the_boundary_fails() -> None:
    """AT. Stiffer than the boundary exceeds the allowable - the direction is a maximum."""
    result = maximum_allowable_adhesive_shear_modulus(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        tolerance=1.0,
    )
    assert peak_at(modulus=result.maximum_shear_modulus * 1.01) > ALLOWABLE


def test_au_modulus_just_below_the_boundary_passes() -> None:
    """AU. Softer than the boundary meets the allowable."""
    result = maximum_allowable_adhesive_shear_modulus(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        tolerance=1.0,
    )
    assert peak_at(modulus=result.maximum_shear_modulus * 0.99) < ALLOWABLE


def test_av_modulus_upper_bound_already_passes() -> None:
    """AV. On a wide, thick bond even a stiff adhesive passes; modulus does not size it."""
    generous = BondedOverlapGeometry(GEOMETRY.overlap_length, 0.200, 2.0 * MM)
    result = maximum_allowable_adhesive_shear_modulus(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, generous, ADHESIVE, ADHESIVE_BASIS,
        minimum_shear_modulus=1.0e6, maximum_shear_modulus=2.0e9,
    )
    assert result.status is AdhesiveModulusStatus.UPPER_BOUND_ALREADY_PASSES
    assert result.maximum_shear_modulus == 2.0e9
    assert result.iterations == 0


def test_av_modulus_lower_bound_already_fails() -> None:
    """AV. A bracket that starts too stiff is reported, not expanded downward."""
    result = maximum_allowable_adhesive_shear_modulus(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        minimum_shear_modulus=0.5e9, maximum_shear_modulus=5.0e9,
    )
    assert result.status is AdhesiveModulusStatus.LOWER_BOUND_ALREADY_FAILS
    assert result.maximum_shear_modulus is None
    assert result.uniform_shear_floor < result.allowable_shear_stress


def test_av_no_finite_modulus_when_the_floor_exceeds_the_allowable() -> None:
    """AV. On a narrow bond no modulus works, because the floor is already too high."""
    narrow = BondedOverlapGeometry(GEOMETRY.overlap_length, 5.0 * MM, GEOMETRY.adhesive_thickness)
    result = maximum_allowable_adhesive_shear_modulus(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, narrow, ADHESIVE, ADHESIVE_BASIS
    )
    assert result.status is AdhesiveModulusStatus.NO_FINITE_MODULUS_WITHIN_MODEL
    assert result.maximum_shear_modulus is None
    assert result.uniform_shear_floor > result.allowable_shear_stress


def test_aw_modulus_search_is_deterministic() -> None:
    args = (MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS)
    assert maximum_allowable_adhesive_shear_modulus(
        *args
    ) == maximum_allowable_adhesive_shear_modulus(*args)


# ------------------------------------------------------------- validation

@pytest.mark.parametrize("bad", [0.0, -1.0, math.nan, math.inf])
def test_all_searches_validate_their_bounds_and_tolerance(bad: float) -> None:
    base = (MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS)
    for kwargs in ({"minimum_bond_width": bad}, {"maximum_bond_width": bad},
                   {"tolerance": bad}):
        with pytest.raises(ValueError):
            required_bond_width(*base, **kwargs)
    for kwargs in ({"minimum_adhesive_thickness": bad},
                   {"maximum_adhesive_thickness": bad}, {"tolerance": bad}):
        with pytest.raises(ValueError):
            required_adhesive_thickness(*base, **kwargs)
    for kwargs in ({"minimum_shear_modulus": bad}, {"maximum_shear_modulus": bad},
                   {"tolerance": bad}):
        with pytest.raises(ValueError):
            maximum_allowable_adhesive_shear_modulus(*base, **kwargs)


def test_all_searches_reject_inverted_brackets_and_bad_iteration_caps() -> None:
    base = (MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS)
    with pytest.raises(ValueError):
        required_bond_width(*base, minimum_bond_width=1.0, maximum_bond_width=0.5)
    with pytest.raises(ValueError):
        required_adhesive_thickness(
            *base, minimum_adhesive_thickness=0.01, maximum_adhesive_thickness=0.001
        )
    with pytest.raises(ValueError):
        maximum_allowable_adhesive_shear_modulus(
            *base, minimum_shear_modulus=1.0e9, maximum_shear_modulus=1.0e8
        )
    with pytest.raises(ValueError):
        required_bond_width(*base, max_iterations=0)


def test_iteration_caps_are_honoured() -> None:
    result = required_bond_width(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS,
        tolerance=1.0e-15, max_iterations=6,
    )
    assert result.iterations == 6


def test_design_factor_sensitivity_of_the_required_width() -> None:
    """Section 26: the design basis materially changes the sizing outcome."""
    widths = {}
    for factor in (1.0, 1.25, 1.5, 2.0):
        result = required_bond_width(
            MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE,
            AdhesiveShearBasis(factor), maximum_bond_width=5.0, tolerance=1.0e-9,
        )
        assert result.status is BondWidthStatus.FINITE_REQUIRED_WIDTH
        widths[factor] = result.required_bond_width

    assert widths[1.0] == pytest.approx(0.141027, rel=1e-4)
    assert widths[1.25] == pytest.approx(0.220354, rel=1e-4)
    assert widths[1.5] == pytest.approx(0.317310, rel=1e-4)
    assert widths[2.0] == pytest.approx(0.564106, rel=1e-4)
    # Required width scales as 1/tau_allow^2, i.e. as design_factor^2.
    assert widths[2.0] / widths[1.0] == pytest.approx(4.0, rel=1e-3)

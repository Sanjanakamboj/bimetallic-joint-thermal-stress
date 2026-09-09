"""Tests AP-AW: allowable temperature change and the required-overlap search."""

from __future__ import annotations

import math

import pytest
from shear_lag_cases import HAND_ADHESIVE, HAND_GEOMETRY, HAND_MEMBER_1, HAND_MEMBER_2

from thermal_joint import (
    AdhesiveMaterial,
    AdhesiveShearBasis,
    AxialMember,
    OverlapLimitStatus,
    ThermalEnvironment,
    ThermoelasticMaterial,
    adherend_compliance_sum,
    allowable_temperature_change_for_adhesive_shear,
    assess_shear_lag_extremes,
    required_overlap_length,
    shear_lag_parameter,
    solve_shear_lag,
)
from thermal_joint.illustrative import (
    ADHESIVE_LIKE,
    aluminium_like_member,
    radiator_joint_environment,
    radiator_joint_overlap,
    titanium_like_member,
)

MEMBER_AL = aluminium_like_member(100.0)
MEMBER_TI = titanium_like_member(100.0)
ENVIRONMENT = radiator_joint_environment()
GEOMETRY = radiator_joint_overlap()
BASIS = AdhesiveShearBasis(1.25)

# A synthetic high-strength adhesive, used only to exercise branches the
# illustrative canonical case does not reach.
STRONG_ADHESIVE = AdhesiveMaterial(
    "Synthetic high-strength (test only)", 1.0e9, 200.0e6, "synthetic test fixture"
)


def test_ap_allowable_delta_temperature_hand_calculation() -> None:
    """AP. dT_allow = tau_allow b C / (|d_alpha| beta coth(beta L_b))."""
    allowable = allowable_temperature_change_for_adhesive_shear(
        MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, BASIS
    )
    beta = shear_lag_parameter(MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE)
    expected = (
        BASIS.allowable_shear_stress(ADHESIVE_LIKE)
        * GEOMETRY.bond_width
        * adherend_compliance_sum(MEMBER_AL, MEMBER_TI)
        / (
            abs(
                MEMBER_AL.material.thermal_expansion_coefficient
                - MEMBER_TI.material.thermal_expansion_coefficient
            )
            * beta
            / math.tanh(beta * GEOMETRY.overlap_length)
        )
    )
    assert allowable == pytest.approx(expected, rel=1e-13)
    assert allowable == pytest.approx(42.1773, rel=1e-5)
    assert allowable > 0.0


def test_aq_allowable_delta_temperature_round_trip() -> None:
    """AQ. At the returned dT the peak shear sits exactly on the allowable."""
    allowable = allowable_temperature_change_for_adhesive_shear(
        MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, BASIS
    )
    for sign in (+1.0, -1.0):
        demand = solve_shear_lag(
            MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, sign * allowable
        )
        assert demand.peak_shear_stress == pytest.approx(
            BASIS.allowable_shear_stress(ADHESIVE_LIKE), rel=1e-12
        )
        assert BASIS.margin_of_safety(
            ADHESIVE_LIKE, demand.peak_shear_stress
        ) == pytest.approx(0.0, abs=1e-12)


def test_aq_just_inside_and_outside_the_allowable_delta_temperature() -> None:
    """AQ. Slightly cooler passes, slightly hotter fails."""
    allowable = allowable_temperature_change_for_adhesive_shear(
        MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, BASIS
    )
    inside = solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, 0.99 * allowable)
    outside = solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, 1.01 * allowable)
    assert BASIS.margin_of_safety(ADHESIVE_LIKE, inside.peak_shear_stress) > 0.0
    assert BASIS.margin_of_safety(ADHESIVE_LIKE, outside.peak_shear_stress) < 0.0


def test_ap_allowable_delta_temperature_is_infinite_without_mismatch() -> None:
    """AP. Matched CTEs generate no shear, so no temperature limit applies."""
    alpha = 15.0e-6
    member_1 = AxialMember(ThermoelasticMaterial("A", 70.0e9, alpha, 270.0e6), 100.0e-6)
    member_2 = AxialMember(ThermoelasticMaterial("B", 110.0e9, alpha, 830.0e6), 100.0e-6)
    assert math.isinf(
        allowable_temperature_change_for_adhesive_shear(
            member_1, member_2, GEOMETRY, ADHESIVE_LIKE, BASIS
        )
    )


def test_ar_lower_bound_already_passes() -> None:
    """AR. With a strong adhesive and a long lower bound, nothing more is needed."""
    result = required_overlap_length(
        MEMBER_AL,
        MEMBER_TI,
        ENVIRONMENT,
        GEOMETRY,
        STRONG_ADHESIVE,
        AdhesiveShearBasis(1.0),
        minimum_overlap_length=0.050,
        maximum_overlap_length=1.0,
    )
    assert result.status is OverlapLimitStatus.LOWER_BOUND_ALREADY_PASSES
    assert result.required_overlap_length == 0.050
    assert result.margin_at_required_length >= 0.0
    assert result.iterations == 0


def test_as_finite_required_overlap_length() -> None:
    """AS. Between the short-overlap peak and the asymptote a finite length exists."""
    result = required_overlap_length(
        MEMBER_AL,
        MEMBER_TI,
        ENVIRONMENT,
        GEOMETRY,
        STRONG_ADHESIVE,
        AdhesiveShearBasis(1.0),
        minimum_overlap_length=1.0e-3,
        maximum_overlap_length=1.0,
    )
    assert result.status is OverlapLimitStatus.FINITE_REQUIRED_LENGTH
    assert result.has_finite_length
    assert result.peak_shear_at_lower_bound > result.allowable_shear_stress
    assert result.long_overlap_asymptote < result.allowable_shear_stress
    assert 1.0e-3 < result.required_overlap_length < 1.0


def test_as_finite_length_matches_an_independent_analytical_solution() -> None:
    """AS. tau_inf coth(beta L) = tau_allow solved directly, no bisection."""
    result = required_overlap_length(
        MEMBER_AL,
        MEMBER_TI,
        ENVIRONMENT,
        GEOMETRY,
        STRONG_ADHESIVE,
        AdhesiveShearBasis(1.0),
        minimum_overlap_length=1.0e-3,
        maximum_overlap_length=1.0,
        length_tolerance=1.0e-9,
    )
    beta = shear_lag_parameter(MEMBER_AL, MEMBER_TI, GEOMETRY, STRONG_ADHESIVE)
    cold_force = solve_shear_lag(
        MEMBER_AL, MEMBER_TI, GEOMETRY, STRONG_ADHESIVE, -140.0
    ).transferred_force
    asymptote = cold_force * beta / GEOMETRY.bond_width
    allowable = STRONG_ADHESIVE.shear_strength
    expected = math.atanh(asymptote / allowable) / beta
    assert result.required_overlap_length == pytest.approx(expected, abs=1e-7)


def test_at_returned_overlap_gives_approximately_zero_margin() -> None:
    """AT. The boundary length sits on MS = 0: passing there, failing just below."""
    result = required_overlap_length(
        MEMBER_AL,
        MEMBER_TI,
        ENVIRONMENT,
        GEOMETRY,
        STRONG_ADHESIVE,
        AdhesiveShearBasis(1.0),
        minimum_overlap_length=1.0e-3,
        maximum_overlap_length=1.0,
        length_tolerance=1.0e-9,
    )
    length = result.required_overlap_length
    assert result.margin_at_required_length == pytest.approx(0.0, abs=1e-6)
    assert result.margin_at_required_length >= 0.0

    def margin_at(candidate: float) -> float:
        return assess_shear_lag_extremes(
            MEMBER_AL,
            MEMBER_TI,
            ENVIRONMENT,
            GEOMETRY.with_overlap_length(candidate),
            STRONG_ADHESIVE,
            AdhesiveShearBasis(1.0),
        ).minimum_margin

    assert margin_at(length * 1.05) > 0.0
    assert margin_at(length * 0.95) < 0.0


def test_au_no_finite_length_within_model_for_the_canonical_case() -> None:
    """AU. The illustrative allowable is below the asymptote: no overlap can pass."""
    result = required_overlap_length(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE, BASIS
    )
    assert result.status is OverlapLimitStatus.NO_FINITE_LENGTH_WITHIN_MODEL
    assert result.required_overlap_length is None
    assert result.margin_at_required_length is None
    assert result.allowable_shear_stress < result.long_overlap_asymptote
    assert result.long_overlap_asymptote == pytest.approx(66.3858e6, rel=1e-5)
    assert result.iterations == 0
    assert not result.has_finite_length


def test_av_insufficient_upper_bound_is_reported_not_expanded() -> None:
    """AV. A boundary above the bracket is reported; bounds are never widened."""
    marginal = AdhesiveMaterial(
        "Synthetic marginal (test only)", 1.0e9, 67.0e6, "synthetic test fixture"
    )
    result = required_overlap_length(
        MEMBER_AL,
        MEMBER_TI,
        ENVIRONMENT,
        GEOMETRY,
        marginal,
        AdhesiveShearBasis(1.0),
        minimum_overlap_length=1.0e-3,
        maximum_overlap_length=5.0e-3,
    )
    assert result.status is OverlapLimitStatus.NO_BOUNDARY_WITHIN_SEARCH_BOUNDS
    assert result.required_overlap_length is None
    assert result.maximum_overlap_length == 5.0e-3
    assert result.allowable_shear_stress > result.long_overlap_asymptote
    assert result.peak_shear_at_upper_bound > result.allowable_shear_stress


def test_av_widening_the_bound_then_finds_the_boundary() -> None:
    """AV. The same case succeeds once the caller supplies a big enough bracket."""
    marginal = AdhesiveMaterial(
        "Synthetic marginal (test only)", 1.0e9, 67.0e6, "synthetic test fixture"
    )
    result = required_overlap_length(
        MEMBER_AL,
        MEMBER_TI,
        ENVIRONMENT,
        GEOMETRY,
        marginal,
        AdhesiveShearBasis(1.0),
        minimum_overlap_length=1.0e-3,
        maximum_overlap_length=0.5,
    )
    assert result.status is OverlapLimitStatus.FINITE_REQUIRED_LENGTH
    assert result.required_overlap_length > 5.0e-3


def test_aw_inverse_result_is_deterministic() -> None:
    """AW. Repeating the search reproduces the identical result object."""
    kwargs = dict(
        minimum_overlap_length=1.0e-3,
        maximum_overlap_length=1.0,
        length_tolerance=1.0e-8,
    )
    first = required_overlap_length(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, STRONG_ADHESIVE,
        AdhesiveShearBasis(1.0), **kwargs
    )
    second = required_overlap_length(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, STRONG_ADHESIVE,
        AdhesiveShearBasis(1.0), **kwargs
    )
    assert first == second


def test_tighter_tolerance_costs_more_iterations_and_respects_the_cap() -> None:
    """The tolerance and iteration cap are explicit and honoured."""
    coarse = required_overlap_length(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, STRONG_ADHESIVE,
        AdhesiveShearBasis(1.0), length_tolerance=1.0e-3,
    )
    fine = required_overlap_length(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, STRONG_ADHESIVE,
        AdhesiveShearBasis(1.0), length_tolerance=1.0e-9,
    )
    assert coarse.iterations < fine.iterations
    capped = required_overlap_length(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, STRONG_ADHESIVE,
        AdhesiveShearBasis(1.0), length_tolerance=1.0e-15, max_iterations=4,
    )
    assert capped.iterations == 4


@pytest.mark.parametrize("bad", [0.0, -1.0, math.nan, math.inf])
def test_search_bounds_and_tolerance_are_validated(bad: float) -> None:
    for kwargs in (
        {"minimum_overlap_length": bad},
        {"maximum_overlap_length": bad},
        {"length_tolerance": bad},
    ):
        with pytest.raises(ValueError):
            required_overlap_length(
                MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE, BASIS, **kwargs
            )


def test_bracket_ordering_and_iteration_cap_are_validated() -> None:
    with pytest.raises(ValueError):
        required_overlap_length(
            MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE, BASIS,
            minimum_overlap_length=0.1, maximum_overlap_length=0.05,
        )
    with pytest.raises(ValueError):
        required_overlap_length(
            MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE, BASIS,
            max_iterations=0,
        )


def test_search_does_not_mutate_the_supplied_geometry() -> None:
    """The bisection builds copies; the caller's geometry is untouched."""
    geometry = radiator_joint_overlap()
    required_overlap_length(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, geometry, STRONG_ADHESIVE, AdhesiveShearBasis(1.0)
    )
    assert geometry == radiator_joint_overlap()
    assert geometry.overlap_length == 0.040


def test_environment_with_zero_excursion_passes_at_the_lower_bound() -> None:
    """A zero-dT environment has no demand, so the shortest overlap already passes."""
    result = required_overlap_length(
        MEMBER_AL,
        MEMBER_TI,
        ThermalEnvironment(20.0, 20.0, 20.0),
        GEOMETRY,
        ADHESIVE_LIKE,
        BASIS,
    )
    assert result.status is OverlapLimitStatus.LOWER_BOUND_ALREADY_PASSES

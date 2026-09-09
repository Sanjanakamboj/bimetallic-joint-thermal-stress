"""Sections 23-27: the maximum-allowable-restraint inverse design search."""

from __future__ import annotations

import pytest

from thermal_joint import (
    AxialMember,
    AxialRestraint,
    RestraintLimitStatus,
    ThermalEnvironment,
    ThermoelasticMaterial,
    YieldBasis,
    assess_restrained_temperature_extremes,
    joint_axial_stiffness,
    maximum_allowable_restraint_stiffness,
    rigid_restraint_minimum_margin,
    rigid_restraint_stresses,
)
from thermal_joint.illustrative import (
    aluminium_like_member,
    radiator_joint_environment,
    titanium_like_member,
)

MEMBER_AL = aluminium_like_member(100.0)
MEMBER_TI = titanium_like_member(100.0)
ENVIRONMENT = radiator_joint_environment()
BASIS = YieldBasis(1.25)


def min_margin_at(ratio: float) -> float:
    restraint = AxialRestraint.from_stiffness_ratio(ratio, MEMBER_AL, MEMBER_TI)
    return assess_restrained_temperature_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, restraint, BASIS
    ).minimum_yield_margin


def test_rigid_asymptote_helper_matches_the_rigid_stresses() -> None:
    """Section 27. The asymptotic margin is built from E_i (eps_ref - alpha_i dT)."""
    margin = rigid_restraint_minimum_margin(MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS)
    expected = []
    for delta_t in (ENVIRONMENT.delta_temperature_cold, ENVIRONMENT.delta_temperature_hot):
        stress_1, stress_2 = rigid_restraint_stresses(MEMBER_AL, MEMBER_TI, delta_t)
        expected.append(BASIS.allowable_stress(MEMBER_AL.material) / abs(stress_1) - 1.0)
        expected.append(BASIS.allowable_stress(MEMBER_TI.material) / abs(stress_2) - 1.0)
    assert margin == pytest.approx(min(expected), rel=1e-14)


def test_canonical_case_has_a_finite_limit() -> None:
    """Section 26 outcome A. The illustrative joint does bound the restraint."""
    result = maximum_allowable_restraint_stiffness(MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS)
    assert result.status is RestraintLimitStatus.FINITE_LIMIT
    assert result.has_finite_limit
    assert result.free_joint_minimum_margin > 0.0
    assert result.rigid_restraint_minimum_margin < 0.0


def test_canonical_limit_matches_an_independent_analytical_derivation() -> None:
    """Section 25. Independent closed form for the boundary, no bisection involved.

    The cold aluminium-like member governs, so the boundary satisfies
    ``E_1 (S dT / (K_joint + K_r) - alpha_1 dT) = sigma_allow_1``, giving
    ``K_joint + K_r = S dT / (sigma_allow_1 / E_1 + alpha_1 dT)``.
    """
    result = maximum_allowable_restraint_stiffness(MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS)

    k_joint = joint_axial_stiffness(MEMBER_AL, MEMBER_TI)
    s_term = (
        MEMBER_AL.axial_stiffness * MEMBER_AL.material.thermal_expansion_coefficient
        + MEMBER_TI.axial_stiffness * MEMBER_TI.material.thermal_expansion_coefficient
    )
    delta_t = ENVIRONMENT.delta_temperature_cold
    allowable = BASIS.allowable_stress(MEMBER_AL.material)
    modulus = MEMBER_AL.material.elastic_modulus
    alpha = MEMBER_AL.material.thermal_expansion_coefficient

    total = s_term * delta_t / (allowable / modulus + alpha * delta_t)
    expected_ratio = (total - k_joint) / k_joint

    assert expected_ratio == pytest.approx(13.740543735, rel=1e-9)
    assert result.maximum_stiffness_ratio == pytest.approx(expected_ratio, abs=1e-5)
    assert result.maximum_restraint_stiffness == pytest.approx(
        expected_ratio * k_joint, rel=1e-6
    )


def test_round_trip_at_the_returned_boundary() -> None:
    """Section 25. The margin is ~0 at the limit, positive below, negative above."""
    result = maximum_allowable_restraint_stiffness(MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS)
    limit = result.maximum_stiffness_ratio

    assert result.minimum_margin_at_limit == pytest.approx(0.0, abs=1e-6)
    assert result.minimum_margin_at_limit >= 0.0
    assert min_margin_at(limit) >= 0.0
    assert min_margin_at(limit * 0.99) > 0.0
    assert min_margin_at(limit * 1.01) < 0.0


def test_returned_stiffness_and_ratio_are_consistent() -> None:
    """K_r and eta_r in the result describe the same restraint."""
    result = maximum_allowable_restraint_stiffness(MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS)
    assert result.joint_axial_stiffness == pytest.approx(
        joint_axial_stiffness(MEMBER_AL, MEMBER_TI), rel=1e-15
    )
    assert result.maximum_restraint_stiffness == pytest.approx(
        result.maximum_stiffness_ratio * result.joint_axial_stiffness, rel=1e-14
    )


def test_free_joint_already_fails_status() -> None:
    """Section 24 outcome C. A joint that already yields unrestrained is reported."""
    weak = AxialMember(
        ThermoelasticMaterial("very low yield", 70.0e9, 23.0e-6, 50.0e6), 100.0e-6
    )
    result = maximum_allowable_restraint_stiffness(weak, MEMBER_TI, ENVIRONMENT, YieldBasis())
    assert result.status is RestraintLimitStatus.FREE_JOINT_ALREADY_FAILS
    assert result.maximum_stiffness_ratio is None
    assert result.maximum_restraint_stiffness is None
    assert result.minimum_margin_at_limit is None
    assert result.free_joint_minimum_margin < 0.0
    assert not result.has_finite_limit


def test_no_finite_limit_status_for_a_mild_environment() -> None:
    """Section 24 outcome B. If even rigid restraint passes, no number is invented."""
    mild = ThermalEnvironment(20.0, 0.0, 40.0)
    result = maximum_allowable_restraint_stiffness(MEMBER_AL, MEMBER_TI, mild, BASIS)
    assert result.status is RestraintLimitStatus.NO_FINITE_LIMIT_WITHIN_MODEL
    assert result.maximum_stiffness_ratio is None
    assert result.maximum_restraint_stiffness is None
    assert result.rigid_restraint_minimum_margin >= 0.0


def test_no_boundary_within_search_bounds_is_reported_not_expanded() -> None:
    """Section 25. A too-small bracket is reported, never silently widened."""
    result = maximum_allowable_restraint_stiffness(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS, upper_stiffness_ratio=1.0
    )
    assert result.status is RestraintLimitStatus.NO_BOUNDARY_WITHIN_SEARCH_BOUNDS
    assert result.maximum_stiffness_ratio is None
    assert result.search_upper_ratio == 1.0
    assert result.rigid_restraint_minimum_margin < 0.0


def test_synthetic_case_exercises_the_bisection_independently() -> None:
    """Section 26. A second, unrelated joint also produces a verified finite limit."""
    member_1 = AxialMember(ThermoelasticMaterial("A", 100.0e9, 20.0e-6, 200.0e6), 200.0e-6)
    member_2 = AxialMember(ThermoelasticMaterial("B", 200.0e9, 10.0e-6, 900.0e6), 100.0e-6)
    environment = ThermalEnvironment(20.0, -100.0, 150.0)
    basis = YieldBasis(1.0)

    result = maximum_allowable_restraint_stiffness(member_1, member_2, environment, basis)
    assert result.status is RestraintLimitStatus.FINITE_LIMIT

    limit = result.maximum_stiffness_ratio

    def margin(ratio: float) -> float:
        restraint = AxialRestraint.from_stiffness_ratio(ratio, member_1, member_2)
        return assess_restrained_temperature_extremes(
            member_1, member_2, environment, restraint, basis
        ).minimum_yield_margin

    assert margin(limit) >= 0.0
    assert margin(limit * 1.05) < 0.0
    assert result.minimum_margin_at_limit == pytest.approx(0.0, abs=1e-5)


def test_nonzero_reference_strain_is_supported_by_the_search() -> None:
    """Section 22/23. The inverse search honours a nonzero restraint reference strain."""
    result = maximum_allowable_restraint_stiffness(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS, reference_strain=5.0e-4
    )
    assert result.status in {
        RestraintLimitStatus.FINITE_LIMIT,
        RestraintLimitStatus.NO_FINITE_LIMIT_WITHIN_MODEL,
    }
    if result.has_finite_limit:
        restraint = AxialRestraint.from_stiffness_ratio(
            result.maximum_stiffness_ratio, MEMBER_AL, MEMBER_TI, reference_strain=5.0e-4
        )
        assert (
            assess_restrained_temperature_extremes(
                MEMBER_AL, MEMBER_TI, ENVIRONMENT, restraint, BASIS
            ).minimum_yield_margin
            >= 0.0
        )


def test_search_is_deterministic() -> None:
    """Section 23. Repeating the search reproduces the identical result object."""
    first = maximum_allowable_restraint_stiffness(MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS)
    second = maximum_allowable_restraint_stiffness(MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS)
    assert first == second


def test_tighter_tolerance_costs_more_iterations_and_stays_bounded() -> None:
    """Section 25. Tolerance and iteration cap are explicit and respected."""
    coarse = maximum_allowable_restraint_stiffness(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS, tolerance=1.0e-2
    )
    fine = maximum_allowable_restraint_stiffness(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS, tolerance=1.0e-9
    )
    assert coarse.iterations < fine.iterations
    assert fine.iterations <= 200
    assert fine.maximum_stiffness_ratio == pytest.approx(
        coarse.maximum_stiffness_ratio, abs=1.0e-2
    )


def test_iteration_cap_is_honoured() -> None:
    """Section 25. The search never exceeds max_iterations."""
    result = maximum_allowable_restraint_stiffness(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS, tolerance=1.0e-15, max_iterations=5
    )
    assert result.iterations == 5


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf")])
def test_search_bounds_and_tolerance_are_validated(bad: float) -> None:
    with pytest.raises(ValueError):
        maximum_allowable_restraint_stiffness(
            MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS, upper_stiffness_ratio=bad
        )
    with pytest.raises(ValueError):
        maximum_allowable_restraint_stiffness(
            MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS, tolerance=bad
        )


def test_max_iterations_is_validated() -> None:
    with pytest.raises(ValueError):
        maximum_allowable_restraint_stiffness(
            MEMBER_AL, MEMBER_TI, ENVIRONMENT, BASIS, max_iterations=0
        )


def test_status_values_are_stable_strings() -> None:
    """Section 24. The status vocabulary is explicit and stable."""
    assert RestraintLimitStatus.FREE_JOINT_ALREADY_FAILS.value == "free_joint_already_fails"
    assert RestraintLimitStatus.FINITE_LIMIT.value == "finite_limit"
    assert (
        RestraintLimitStatus.NO_FINITE_LIMIT_WITHIN_MODEL.value
        == "no_finite_limit_within_model"
    )

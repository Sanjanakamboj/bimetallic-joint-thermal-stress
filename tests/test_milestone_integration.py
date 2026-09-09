"""Tests AX-BC: Milestone 1/2 regression and the combined screening result."""

from __future__ import annotations

import pytest

from thermal_joint import (
    AdhesiveShearBasis,
    AxialRestraint,
    YieldBasis,
    assess_restrained_temperature_extremes,
    assess_temperature_extremes,
    screen_joint,
    solve_bimetallic_joint,
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
YIELD_BASIS = YieldBasis(1.25)
ADHESIVE_BASIS = AdhesiveShearBasis(1.25)


@pytest.mark.parametrize("delta_t", [-140.0, -30.0, 12.5, 100.0])
def test_ax_transferred_force_is_the_milestone_1_member_force(delta_t: float) -> None:
    """AX. The shear-lag demand is taken from Milestone 1, bit for bit."""
    free_joint = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, delta_t)
    demand = solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, delta_t)
    assert demand.signed_transferred_force == free_joint.member_1_force
    assert demand.transferred_force == abs(free_joint.member_1_force)
    assert demand.member_1_force(0.0) == pytest.approx(free_joint.member_1_force, rel=1e-15)


@pytest.mark.parametrize("delta_t", [-140.0, 100.0])
def test_ax_demand_matches_the_independent_mismatch_identity(delta_t: float) -> None:
    """AX. N_t = -d_eps / C, an identity independent of the Milestone 1 code path."""
    from thermal_joint import adherend_compliance_sum, thermal_mismatch_strain

    demand = solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, delta_t)
    expected = -thermal_mismatch_strain(MEMBER_AL, MEMBER_TI, delta_t) / (
        adherend_compliance_sum(MEMBER_AL, MEMBER_TI)
    )
    assert demand.signed_transferred_force == pytest.approx(expected, rel=1e-12)


def test_ay_milestone_1_outputs_are_unchanged() -> None:
    """AY. The Milestone 1 free-joint results are exactly as before."""
    hot = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, 100.0)
    cold = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, -140.0)
    assert hot.member_1_stress == pytest.approx(-62.0278e6, rel=1e-5)
    assert hot.member_2_stress == pytest.approx(+62.0278e6, rel=1e-5)
    assert cold.member_1_stress == pytest.approx(+86.8389e6, rel=1e-5)
    assert cold.member_2_stress == pytest.approx(-86.8389e6, rel=1e-5)

    assessment = assess_temperature_extremes(MEMBER_AL, MEMBER_TI, ENVIRONMENT, YIELD_BASIS)
    assert assessment.governing_extreme == "cold"
    assert assessment.governing_member == 1
    assert assessment.minimum_yield_margin == pytest.approx(1.4874, rel=1e-3)
    assert assessment.hot_assessment.minimum_margin == pytest.approx(2.4823, rel=1e-3)
    assert assessment.passes


def test_az_milestone_2_outputs_are_unchanged() -> None:
    """AZ. The Milestone 2 restrained results are exactly as before."""
    restraint = AxialRestraint.from_stiffness_ratio(1.0, MEMBER_AL, MEMBER_TI)
    assessment = assess_restrained_temperature_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, restraint, YIELD_BASIS
    )
    assert assessment.cold_result.member_1_stress == pytest.approx(+156.1194e6, rel=1e-5)
    assert assessment.cold_result.member_2_stress == pytest.approx(+22.0306e6, rel=1e-5)
    assert assessment.minimum_yield_margin == pytest.approx(0.3836, rel=1e-3)
    assert assessment.governing_extreme == "cold"


def test_az_milestone_2_inverse_result_is_unchanged() -> None:
    """AZ. The Milestone 2 restraint boundary is untouched by Milestone 3."""
    from thermal_joint import (
        RestraintLimitStatus,
        maximum_allowable_restraint_stiffness,
        zero_stress_restraint_stiffness,
    )

    limit = maximum_allowable_restraint_stiffness(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, YIELD_BASIS
    )
    assert limit.status is RestraintLimitStatus.FINITE_LIMIT
    assert limit.maximum_stiffness_ratio == pytest.approx(13.740543, abs=1e-5)
    assert limit.maximum_restraint_stiffness == pytest.approx(2.4733e8, rel=1e-4)

    crossing = zero_stress_restraint_stiffness(MEMBER_AL, MEMBER_TI, 100.0, 2)
    assert crossing / 1.8e7 == pytest.approx(0.663399, rel=1e-5)


def test_ba_milestone_3_does_not_perturb_the_global_models() -> None:
    """BA. Solving a shear-lag problem leaves the global results identical."""
    before = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, -140.0)
    solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, -140.0)
    required = assess_temperature_extremes(MEMBER_AL, MEMBER_TI, ENVIRONMENT, YIELD_BASIS)
    after = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, -140.0)
    assert before == after
    assert required.minimum_yield_margin == pytest.approx(1.4874, rel=1e-3)


def test_bb_overall_feasibility_is_a_boolean_and() -> None:
    """BB. Overall feasibility is the AND of the two independent screens."""
    result = screen_joint(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE, YIELD_BASIS, ADHESIVE_BASIS
    )
    assert result.overall_feasible == (result.yield_feasible and result.adhesive_feasible)
    assert result.yield_feasible is True
    assert result.adhesive_feasible is False
    assert result.overall_feasible is False
    assert result.governing_mode == "adhesive shear"


def test_bb_both_screens_passing_gives_overall_pass() -> None:
    """BB. A mild environment passes both screens, so the AND is True."""
    from thermal_joint import ThermalEnvironment

    mild = ThermalEnvironment(20.0, 10.0, 30.0)
    result = screen_joint(
        MEMBER_AL, MEMBER_TI, mild, GEOMETRY, ADHESIVE_LIKE, YIELD_BASIS, ADHESIVE_BASIS
    )
    assert result.yield_feasible and result.adhesive_feasible
    assert result.overall_feasible is True
    assert result.governing_mode == "none - both screens pass"


def test_bb_both_screens_failing_is_reported_as_both() -> None:
    """BB. When both fail, both are named; neither is silently dropped."""
    from thermal_joint import AxialMember, ThermoelasticMaterial

    weak = AxialMember(
        ThermoelasticMaterial("very low yield", 70.0e9, 23.0e-6, 20.0e6), 100.0e-6
    )
    result = screen_joint(
        weak, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE, YIELD_BASIS, ADHESIVE_BASIS
    )
    assert not result.yield_feasible
    assert not result.adhesive_feasible
    assert result.overall_feasible is False
    assert "member yield" in result.governing_mode
    assert "adhesive shear" in result.governing_mode


def test_bc_yield_and_adhesive_margins_stay_separate() -> None:
    """BC. The two margins are reported side by side, never combined numerically."""
    result = screen_joint(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE, YIELD_BASIS, ADHESIVE_BASIS
    )
    assert result.minimum_yield_margin == pytest.approx(1.4874, rel=1e-3)
    assert result.minimum_adhesive_margin == pytest.approx(-0.6987, rel=1e-3)
    assert result.minimum_yield_margin != result.minimum_adhesive_margin

    # The two screens are the untouched Milestone 1 and Milestone 3 assessments.
    assert result.yield_assessment == assess_temperature_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, YIELD_BASIS
    )
    assert result.adhesive_assessment.minimum_margin == result.minimum_adhesive_margin

    # No combined scalar margin is exposed anywhere on the result.
    attributes = set(dir(result))
    assert not {"combined_margin", "overall_margin", "minimum_margin"} & attributes


def test_screen_joint_defaults_to_unity_bases() -> None:
    """Omitting either basis falls back to a design factor of 1.0."""
    assert screen_joint(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE_LIKE
    ) == screen_joint(
        MEMBER_AL,
        MEMBER_TI,
        ENVIRONMENT,
        GEOMETRY,
        ADHESIVE_LIKE,
        YieldBasis(1.0),
        AdhesiveShearBasis(1.0),
    )

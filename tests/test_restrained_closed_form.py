"""Tests G-M: the restrained closed form, hand calculations and back substitution."""

from __future__ import annotations

import pytest
from cases import HAND_MATERIAL_1, HAND_MATERIAL_2
from restrained_cases import (
    HAND_DELTA_T,
    HAND_MEMBER_1,
    HAND_MEMBER_2,
    HAND_OFFSET_COMMON_STRAIN,
    HAND_OFFSET_RESTRAINT,
    HAND_OFFSET_RESTRAINT_FORCE,
    HAND_OFFSET_STRESS_1,
    HAND_OFFSET_STRESS_2,
    HAND_RESTRAINED_COMMON_STRAIN,
    HAND_RESTRAINED_FORCE_1,
    HAND_RESTRAINED_FORCE_2,
    HAND_RESTRAINED_STRESS_1,
    HAND_RESTRAINED_STRESS_2,
    HAND_RESTRAINT,
    HAND_RESTRAINT_FORCE,
)

from thermal_joint import AxialMember, AxialRestraint, solve_restrained_joint

RESULT = solve_restrained_joint(HAND_MEMBER_1, HAND_MEMBER_2, HAND_DELTA_T, HAND_RESTRAINT)

RESTRAINTS = [
    AxialRestraint(0.0),
    AxialRestraint(1.0e6),
    AxialRestraint(4.0e7),
    AxialRestraint(4.0e9),
    AxialRestraint(1.0e7, reference_strain=1.0e-3),
    AxialRestraint(2.5e7, reference_strain=-7.5e-4),
]
DELTA_TS = [-140.0, -30.0, 12.5, 100.0]


def test_g_restrained_common_strain_hand_calculation() -> None:
    """G. eps = (600 * 50 + 0) / (4.0e7 + 1.0e7) = 6.0e-4."""
    assert RESULT.common_strain == pytest.approx(HAND_RESTRAINED_COMMON_STRAIN, rel=1e-14)


def test_h_member_1_stress_hand_calculation() -> None:
    """H. sigma_1 = 100 GPa (6.0e-4 - 1.0e-3) = -40 MPa."""
    assert RESULT.member_1_stress == pytest.approx(HAND_RESTRAINED_STRESS_1, rel=1e-12)


def test_i_member_2_stress_hand_calculation() -> None:
    """I. sigma_2 = 200 GPa (6.0e-4 - 5.0e-4) = +20 MPa."""
    assert RESULT.member_2_stress == pytest.approx(HAND_RESTRAINED_STRESS_2, rel=1e-12)


def test_j_member_force_hand_calculation() -> None:
    """J. N_1 = -8000 N and N_2 = +2000 N."""
    assert RESULT.member_1_force == pytest.approx(HAND_RESTRAINED_FORCE_1, rel=1e-12)
    assert RESULT.member_2_force == pytest.approx(HAND_RESTRAINED_FORCE_2, rel=1e-12)


def test_k_restraint_force_hand_calculation() -> None:
    """K. N_r = K_r (eps_ref - eps) = 1.0e7 (0 - 6.0e-4) = -6000 N."""
    assert RESULT.restraint_force == pytest.approx(HAND_RESTRAINT_FORCE, rel=1e-12)


def test_k_restraint_force_sign_convention_is_documented_and_signed() -> None:
    """K. Positive N_r acts on the joint in the tensile sense; it is not a magnitude."""
    heating = solve_restrained_joint(HAND_MEMBER_1, HAND_MEMBER_2, +50.0, HAND_RESTRAINT)
    cooling = solve_restrained_joint(HAND_MEMBER_1, HAND_MEMBER_2, -50.0, HAND_RESTRAINT)
    assert heating.restraint_force < 0.0  # joint expanded, restraint pulls it back
    assert cooling.restraint_force > 0.0  # joint contracted, restraint holds it out
    assert cooling.restraint_force == pytest.approx(-heating.restraint_force, rel=1e-14)


@pytest.mark.parametrize("restraint", RESTRAINTS)
@pytest.mark.parametrize("delta_t", DELTA_TS)
def test_l_full_equilibrium_residual_vanishes(restraint: AxialRestraint, delta_t: float) -> None:
    """L. N_1 + N_2 - N_r = 0 to machine precision, at every restraint level."""
    result = solve_restrained_joint(HAND_MEMBER_1, HAND_MEMBER_2, delta_t, restraint)
    scale = max(
        abs(result.member_1_force), abs(result.member_2_force), abs(result.restraint_force), 1.0
    )
    assert abs(result.force_equilibrium_residual) <= 1.0e-12 * scale
    assert result.force_equilibrium_residual == (
        result.member_1_force + result.member_2_force - result.restraint_force
    )


@pytest.mark.parametrize("restraint", RESTRAINTS)
@pytest.mark.parametrize("delta_t", DELTA_TS)
def test_m_independent_direct_equilibrium_back_substitution(
    restraint: AxialRestraint, delta_t: float
) -> None:
    """M. Back-substitute the computed strain into the governing equation.

    Checks ``E1 A1 (eps - a1 dT) + E2 A2 (eps - a2 dT) + K_r (eps - eps_ref) = 0``
    from the raw scalars, independently of the result's own force fields.
    """
    result = solve_restrained_joint(HAND_MEMBER_1, HAND_MEMBER_2, delta_t, restraint)
    eps = result.common_strain

    e1_a1 = HAND_MATERIAL_1.elastic_modulus * HAND_MEMBER_1.area
    e2_a2 = HAND_MATERIAL_2.elastic_modulus * HAND_MEMBER_2.area
    residual = (
        e1_a1 * (eps - HAND_MATERIAL_1.thermal_expansion_coefficient * delta_t)
        + e2_a2 * (eps - HAND_MATERIAL_2.thermal_expansion_coefficient * delta_t)
        + restraint.stiffness * (eps - restraint.reference_strain)
    )
    scale = max(abs(e1_a1), abs(e2_a2), restraint.stiffness) * max(abs(eps), 1e-12)
    assert abs(residual) <= 1.0e-11 * scale


@pytest.mark.parametrize("restraint", RESTRAINTS)
@pytest.mark.parametrize("delta_t", DELTA_TS)
def test_m_independent_closed_form_for_the_common_strain(
    restraint: AxialRestraint, delta_t: float
) -> None:
    """M. Independent evaluation of the closed form from raw scalars."""
    result = solve_restrained_joint(HAND_MEMBER_1, HAND_MEMBER_2, delta_t, restraint)
    e1_a1 = HAND_MATERIAL_1.elastic_modulus * HAND_MEMBER_1.area
    e2_a2 = HAND_MATERIAL_2.elastic_modulus * HAND_MEMBER_2.area
    expected = (
        e1_a1 * HAND_MATERIAL_1.thermal_expansion_coefficient * delta_t
        + e2_a2 * HAND_MATERIAL_2.thermal_expansion_coefficient * delta_t
        + restraint.stiffness * restraint.reference_strain
    ) / (e1_a1 + e2_a2 + restraint.stiffness)
    assert result.common_strain == pytest.approx(expected, rel=1e-14, abs=1e-24)


def test_result_reports_its_restraint_inputs_and_free_reference() -> None:
    """The result carries the restraint inputs and the free-joint strain."""
    assert RESULT.restraint_stiffness == 1.0e7
    assert RESULT.restraint_reference_strain == 0.0
    assert RESULT.delta_temperature == HAND_DELTA_T
    assert RESULT.free_joint_common_strain == pytest.approx(7.5e-4, rel=1e-14)
    assert abs(RESULT.common_strain) < abs(RESULT.free_joint_common_strain)


def test_nonzero_reference_strain_hand_calculation() -> None:
    """Section 22: the generalised formula with eps_ref = 1.0e-3."""
    result = solve_restrained_joint(
        HAND_MEMBER_1, HAND_MEMBER_2, HAND_DELTA_T, HAND_OFFSET_RESTRAINT
    )
    assert result.common_strain == pytest.approx(HAND_OFFSET_COMMON_STRAIN, rel=1e-14)
    assert result.member_1_stress == pytest.approx(HAND_OFFSET_STRESS_1, rel=1e-12)
    assert result.member_2_stress == pytest.approx(HAND_OFFSET_STRESS_2, rel=1e-12)
    assert result.restraint_force == pytest.approx(HAND_OFFSET_RESTRAINT_FORCE, rel=1e-12)


def test_solver_defaults_to_no_restraint() -> None:
    """Omitting the restraint is the same as passing a zero-stiffness restraint."""
    assert solve_restrained_joint(HAND_MEMBER_1, HAND_MEMBER_2, 50.0) == solve_restrained_joint(
        HAND_MEMBER_1, HAND_MEMBER_2, 50.0, AxialRestraint(0.0)
    )


def test_solver_rejects_wrong_argument_types() -> None:
    """Type checking mirrors the Milestone 1 solver."""
    with pytest.raises(TypeError):
        solve_restrained_joint(HAND_MATERIAL_1, HAND_MEMBER_2, 50.0)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        solve_restrained_joint(HAND_MEMBER_1, HAND_MEMBER_2, 50.0, 1.0e7)  # type: ignore[arg-type]


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_solver_rejects_non_finite_delta_temperature(bad: float) -> None:
    with pytest.raises(ValueError):
        solve_restrained_joint(HAND_MEMBER_1, HAND_MEMBER_2, bad, HAND_RESTRAINT)


def test_area_and_modulus_still_drive_the_restrained_solution() -> None:
    """Sanity: the restrained solution still depends on both members' E*A."""
    stiff = AxialMember(HAND_MATERIAL_1, 400.0e-6)
    baseline = solve_restrained_joint(HAND_MEMBER_1, HAND_MEMBER_2, 50.0, HAND_RESTRAINT)
    heavier = solve_restrained_joint(stiff, HAND_MEMBER_2, 50.0, HAND_RESTRAINT)
    assert heavier.common_strain != baseline.common_strain

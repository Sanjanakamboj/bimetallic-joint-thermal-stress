"""Tests N-V: zero-restraint regression to Milestone 1 and the rigid-restraint limit."""

from __future__ import annotations

import pytest
from cases import HAND_MEMBER_1, HAND_MEMBER_2

from thermal_joint import (
    AxialMember,
    AxialRestraint,
    ThermoelasticMaterial,
    rigid_restraint_stresses,
    solve_bimetallic_joint,
    solve_restrained_joint,
)
from thermal_joint.illustrative import aluminium_like_member, titanium_like_member

DELTA_TS = [-140.0, -30.0, 0.0, 12.5, 100.0]
MEMBER_AL = aluminium_like_member(100.0)
MEMBER_TI = titanium_like_member(100.0)

# Restraint stiffness far above the joint stiffness, standing in for K_r -> inf.
VERY_LARGE_RATIO = 1.0e9


@pytest.mark.parametrize("delta_t", DELTA_TS)
def test_n_zero_restraint_reproduces_milestone_1_common_strain(delta_t: float) -> None:
    """N. K_r = 0 gives the Milestone 1 common strain to machine precision."""
    free = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, delta_t)
    restrained = solve_restrained_joint(
        HAND_MEMBER_1, HAND_MEMBER_2, delta_t, AxialRestraint(0.0)
    )
    assert restrained.common_strain == pytest.approx(free.common_strain, rel=1e-15, abs=1e-30)


@pytest.mark.parametrize("delta_t", DELTA_TS)
def test_o_zero_restraint_reproduces_milestone_1_stresses(delta_t: float) -> None:
    """O. K_r = 0 gives the Milestone 1 member stresses."""
    free = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, delta_t)
    restrained = solve_restrained_joint(
        HAND_MEMBER_1, HAND_MEMBER_2, delta_t, AxialRestraint(0.0)
    )
    assert restrained.member_1_stress == pytest.approx(free.member_1_stress, rel=1e-14, abs=1e-9)
    assert restrained.member_2_stress == pytest.approx(free.member_2_stress, rel=1e-14, abs=1e-9)


@pytest.mark.parametrize("delta_t", DELTA_TS)
def test_p_zero_restraint_reproduces_milestone_1_forces(delta_t: float) -> None:
    """P. K_r = 0 gives the Milestone 1 member forces, still equal and opposite."""
    free = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, delta_t)
    restrained = solve_restrained_joint(
        HAND_MEMBER_1, HAND_MEMBER_2, delta_t, AxialRestraint(0.0)
    )
    assert restrained.member_1_force == pytest.approx(free.member_1_force, rel=1e-14, abs=1e-12)
    assert restrained.member_2_force == pytest.approx(free.member_2_force, rel=1e-14, abs=1e-12)
    assert restrained.member_1_force == pytest.approx(-restrained.member_2_force, rel=1e-12, abs=1e-12)


@pytest.mark.parametrize("delta_t", DELTA_TS)
def test_q_zero_restraint_gives_zero_restraint_force(delta_t: float) -> None:
    """Q. A zero-stiffness restraint carries exactly no force."""
    restrained = solve_restrained_joint(
        HAND_MEMBER_1, HAND_MEMBER_2, delta_t, AxialRestraint(0.0)
    )
    assert restrained.restraint_force == 0.0
    free = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, delta_t)
    assert restrained.force_equilibrium_residual == pytest.approx(
        free.force_equilibrium_residual, abs=1e-9
    )


@pytest.mark.parametrize("delta_t", [-140.0, 100.0])
@pytest.mark.parametrize("reference_strain", [0.0, 5.0e-4, -3.0e-4])
def test_r_very_large_restraint_approaches_the_reference_strain(
    delta_t: float, reference_strain: float
) -> None:
    """R. K_r >> K_joint drives the common strain to the restraint reference strain."""
    restraint = AxialRestraint.from_stiffness_ratio(
        VERY_LARGE_RATIO, HAND_MEMBER_1, HAND_MEMBER_2, reference_strain=reference_strain
    )
    result = solve_restrained_joint(HAND_MEMBER_1, HAND_MEMBER_2, delta_t, restraint)
    assert result.common_strain == pytest.approx(reference_strain, rel=1e-8, abs=1e-11)


def test_r_common_strain_moves_monotonically_toward_the_reference_strain() -> None:
    """R. Increasing restraint monotonically shrinks the deviation from eps_ref."""
    deviations = []
    for ratio in (0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0):
        restraint = AxialRestraint.from_stiffness_ratio(ratio, MEMBER_AL, MEMBER_TI)
        result = solve_restrained_joint(MEMBER_AL, MEMBER_TI, -140.0, restraint)
        deviations.append(abs(result.common_strain))
    assert deviations == sorted(deviations, reverse=True)


@pytest.mark.parametrize("delta_t", [-140.0, 100.0])
@pytest.mark.parametrize("reference_strain", [0.0, 5.0e-4])
def test_s_very_large_restraint_approaches_rigid_stress_member_1(
    delta_t: float, reference_strain: float
) -> None:
    """S. sigma_1 -> E_1 (eps_ref - alpha_1 dT)."""
    restraint = AxialRestraint.from_stiffness_ratio(
        VERY_LARGE_RATIO, MEMBER_AL, MEMBER_TI, reference_strain=reference_strain
    )
    result = solve_restrained_joint(MEMBER_AL, MEMBER_TI, delta_t, restraint)
    expected_1, _ = rigid_restraint_stresses(MEMBER_AL, MEMBER_TI, delta_t, reference_strain)
    assert result.member_1_stress == pytest.approx(expected_1, rel=1e-7)


@pytest.mark.parametrize("delta_t", [-140.0, 100.0])
@pytest.mark.parametrize("reference_strain", [0.0, 5.0e-4])
def test_t_very_large_restraint_approaches_rigid_stress_member_2(
    delta_t: float, reference_strain: float
) -> None:
    """T. sigma_2 -> E_2 (eps_ref - alpha_2 dT)."""
    restraint = AxialRestraint.from_stiffness_ratio(
        VERY_LARGE_RATIO, MEMBER_AL, MEMBER_TI, reference_strain=reference_strain
    )
    result = solve_restrained_joint(MEMBER_AL, MEMBER_TI, delta_t, restraint)
    _, expected_2 = rigid_restraint_stresses(MEMBER_AL, MEMBER_TI, delta_t, reference_strain)
    assert result.member_2_stress == pytest.approx(expected_2, rel=1e-7)


def test_rigid_reference_formula_is_the_classic_fully_restrained_result() -> None:
    """For eps_ref = 0 the rigid stresses reduce to -E_i alpha_i dT."""
    for delta_t in (-140.0, 100.0):
        stress_1, stress_2 = rigid_restraint_stresses(MEMBER_AL, MEMBER_TI, delta_t)
        assert stress_1 == pytest.approx(
            -MEMBER_AL.material.elastic_modulus
            * MEMBER_AL.material.thermal_expansion_coefficient
            * delta_t,
            rel=1e-14,
        )
        assert stress_2 == pytest.approx(
            -MEMBER_TI.material.elastic_modulus
            * MEMBER_TI.material.thermal_expansion_coefficient
            * delta_t,
            rel=1e-14,
        )


def test_rigid_limit_drives_both_members_to_the_same_sign() -> None:
    """The qualitative Milestone 2 transition: self-equilibrating -> same-sign.

    The free joint always has one member in tension and the other in
    compression. Under strong restraint with eps_ref = 0, heating puts both in
    compression and cooling puts both in tension.
    """
    strong = AxialRestraint.from_stiffness_ratio(VERY_LARGE_RATIO, MEMBER_AL, MEMBER_TI)

    free_hot = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, 100.0)
    assert free_hot.member_1_stress * free_hot.member_2_stress < 0.0

    hot = solve_restrained_joint(MEMBER_AL, MEMBER_TI, 100.0, strong)
    assert hot.member_1_stress < 0.0 and hot.member_2_stress < 0.0
    assert hot.members_carry_same_sign_stress

    cold = solve_restrained_joint(MEMBER_AL, MEMBER_TI, -140.0, strong)
    assert cold.member_1_stress > 0.0 and cold.member_2_stress > 0.0
    assert cold.members_carry_same_sign_stress


@pytest.mark.parametrize("delta_t", [-140.0, 100.0])
@pytest.mark.parametrize("alpha", [23.0e-6, 8.5e-6, -2.0e-6])
def test_u_identical_cte_with_finite_restraint_gives_nonzero_stress(
    delta_t: float, alpha: float
) -> None:
    """U. Under restraint, matching CTE no longer implies zero stress.

    Milestone 1's zero-stress result for alpha_1 = alpha_2 depended on the joint
    being externally free. With a finite restraint and eps_ref = 0 both members
    develop same-sign restraint stress.
    """
    member_1 = AxialMember(ThermoelasticMaterial("A", 70.0e9, alpha, 270.0e6), 100.0e-6)
    member_2 = AxialMember(ThermoelasticMaterial("B", 110.0e9, alpha, 830.0e6), 250.0e-6)
    restraint = AxialRestraint.from_stiffness_ratio(1.0, member_1, member_2)
    result = solve_restrained_joint(member_1, member_2, delta_t, restraint)

    assert abs(result.member_1_stress) > 1.0e3
    assert abs(result.member_2_stress) > 1.0e3
    assert result.members_carry_same_sign_stress
    # Heating with a positive CTE compresses both; cooling tensions both.
    expected_sign = -1.0 if alpha * delta_t > 0 else 1.0
    assert result.member_1_stress * expected_sign > 0.0
    assert result.member_2_stress * expected_sign > 0.0


@pytest.mark.parametrize("delta_t", [-140.0, 100.0])
@pytest.mark.parametrize("alpha", [23.0e-6, 8.5e-6, -2.0e-6])
def test_v_identical_cte_with_zero_restraint_still_gives_zero_stress(
    delta_t: float, alpha: float
) -> None:
    """V. The Milestone 1 zero-stress result is untouched at K_r = 0."""
    member_1 = AxialMember(ThermoelasticMaterial("A", 70.0e9, alpha, 270.0e6), 100.0e-6)
    member_2 = AxialMember(ThermoelasticMaterial("B", 110.0e9, alpha, 830.0e6), 250.0e-6)
    result = solve_restrained_joint(member_1, member_2, delta_t, AxialRestraint(0.0))

    assert result.member_1_stress == pytest.approx(0.0, abs=1e-6)
    assert result.member_2_stress == pytest.approx(0.0, abs=1e-6)
    assert result.restraint_force == 0.0
    assert result.common_strain == pytest.approx(alpha * delta_t, rel=1e-15, abs=1e-30)


def test_zero_delta_temperature_with_restraint_and_zero_reference_is_a_zero_state() -> None:
    """dT = 0 with eps_ref = 0 leaves the restrained joint completely unloaded."""
    restraint = AxialRestraint.from_stiffness_ratio(2.0, MEMBER_AL, MEMBER_TI)
    result = solve_restrained_joint(MEMBER_AL, MEMBER_TI, 0.0, restraint)
    assert result.common_strain == 0.0
    assert result.member_1_stress == 0.0
    assert result.member_2_stress == 0.0
    assert result.restraint_force == 0.0
    assert result.force_equilibrium_residual == 0.0


def test_matched_free_strain_and_reference_strain_unloads_everything() -> None:
    """Section 22: if both free strains equal eps_ref, nothing is loaded."""
    alpha, delta_t = 15.0e-6, 80.0
    member_1 = AxialMember(ThermoelasticMaterial("A", 70.0e9, alpha, 270.0e6), 100.0e-6)
    member_2 = AxialMember(ThermoelasticMaterial("B", 110.0e9, alpha, 830.0e6), 250.0e-6)
    restraint = AxialRestraint(5.0e7, reference_strain=alpha * delta_t)
    result = solve_restrained_joint(member_1, member_2, delta_t, restraint)

    assert result.common_strain == pytest.approx(alpha * delta_t, rel=1e-15)
    assert result.member_1_stress == pytest.approx(0.0, abs=1e-6)
    assert result.member_2_stress == pytest.approx(0.0, abs=1e-6)
    assert result.restraint_force == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize("reference_strain", [-1.0e-3, -2.0e-4, 0.0, 4.0e-4, 1.5e-3])
def test_changing_the_reference_strain_shifts_the_common_strain_toward_it(
    reference_strain: float,
) -> None:
    """Section 22: the common strain tracks the reference strain monotonically."""
    restraint = AxialRestraint.from_stiffness_ratio(
        1.0, MEMBER_AL, MEMBER_TI, reference_strain=reference_strain
    )
    result = solve_restrained_joint(MEMBER_AL, MEMBER_TI, 100.0, restraint)
    baseline = solve_restrained_joint(
        MEMBER_AL,
        MEMBER_TI,
        100.0,
        AxialRestraint.from_stiffness_ratio(1.0, MEMBER_AL, MEMBER_TI),
    )
    # eps_common is an increasing affine function of eps_ref at fixed stiffness.
    assert (result.common_strain - baseline.common_strain) * reference_strain >= 0.0

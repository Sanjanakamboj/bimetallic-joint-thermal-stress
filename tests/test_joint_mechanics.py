"""Tests M-U: common strain, member stresses, equilibrium and area identities."""

from __future__ import annotations

import pytest
from cases import (
    HAND_COMMON_STRAIN,
    HAND_DELTA_T,
    HAND_FORCE_1,
    HAND_FORCE_2,
    HAND_MATERIAL_1,
    HAND_MATERIAL_2,
    HAND_MEMBER_1,
    HAND_MEMBER_2,
    HAND_STRESS_1,
    HAND_STRESS_2,
)
from reference_solution import reference_stresses

from thermal_joint import AxialMember, common_joint_strain, solve_bimetallic_joint

RESULT = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, HAND_DELTA_T)


def test_m_common_strain_hand_calculation() -> None:
    """M. eps_common = (2e7*20e-6 + 2e7*10e-6)/4e7 * 50 = 7.5e-4."""
    assert RESULT.common_strain == pytest.approx(HAND_COMMON_STRAIN, rel=1e-14)
    assert common_joint_strain(
        HAND_MEMBER_1, HAND_MEMBER_2, HAND_DELTA_T
    ) == pytest.approx(HAND_COMMON_STRAIN, rel=1e-14)


def test_n_stress_1_hand_calculation() -> None:
    """N. sigma_1 = 100 GPa * (7.5e-4 - 1.0e-3) = -25 MPa (compression)."""
    assert RESULT.member_1_stress == pytest.approx(HAND_STRESS_1, rel=1e-12)
    assert RESULT.member_1_stress < 0.0


def test_o_stress_2_hand_calculation() -> None:
    """O. sigma_2 = 200 GPa * (7.5e-4 - 5.0e-4) = +50 MPa (tension)."""
    assert RESULT.member_2_stress == pytest.approx(HAND_STRESS_2, rel=1e-12)
    assert RESULT.member_2_stress > 0.0


def test_free_thermal_strains_are_reported() -> None:
    """The result carries both unrestrained thermal strains."""
    assert RESULT.member_1_free_thermal_strain == pytest.approx(1.0e-3, rel=1e-14)
    assert RESULT.member_2_free_thermal_strain == pytest.approx(5.0e-4, rel=1e-14)


AREA_PAIRS = [
    (200.0e-6, 100.0e-6),
    (100.0e-6, 100.0e-6),
    (25.0e-6, 400.0e-6),
    (1.0e-6, 1000.0e-6),
]
DELTA_TS = [-140.0, -30.0, 12.5, 100.0]


@pytest.mark.parametrize("area_1, area_2", AREA_PAIRS)
@pytest.mark.parametrize("delta_t", DELTA_TS)
def test_p_independent_closed_form_stress_1(area_1: float, area_2: float, delta_t: float) -> None:
    """P. sigma_1 matches [E1 E2 A2/(E1A1+E2A2)] (alpha_2-alpha_1) dT."""
    member_1 = AxialMember(HAND_MATERIAL_1, area_1)
    member_2 = AxialMember(HAND_MATERIAL_2, area_2)
    result = solve_bimetallic_joint(member_1, member_2, delta_t)
    expected_1, _ = reference_stresses(member_1, member_2, delta_t)
    assert result.member_1_stress == pytest.approx(expected_1, rel=1e-12)


@pytest.mark.parametrize("area_1, area_2", AREA_PAIRS)
@pytest.mark.parametrize("delta_t", DELTA_TS)
def test_q_independent_closed_form_stress_2(area_1: float, area_2: float, delta_t: float) -> None:
    """Q. sigma_2 matches [E1 E2 A1/(E1A1+E2A2)] (alpha_1-alpha_2) dT."""
    member_1 = AxialMember(HAND_MATERIAL_1, area_1)
    member_2 = AxialMember(HAND_MATERIAL_2, area_2)
    result = solve_bimetallic_joint(member_1, member_2, delta_t)
    _, expected_2 = reference_stresses(member_1, member_2, delta_t)
    assert result.member_2_stress == pytest.approx(expected_2, rel=1e-12)


def test_r_internal_forces_hand_calculation() -> None:
    """R. N_1 = -5000 N, N_2 = +5000 N for the hand-calculation case."""
    assert RESULT.member_1_force == pytest.approx(HAND_FORCE_1, rel=1e-12)
    assert RESULT.member_2_force == pytest.approx(HAND_FORCE_2, rel=1e-12)


@pytest.mark.parametrize("area_1, area_2", AREA_PAIRS)
@pytest.mark.parametrize("delta_t", DELTA_TS)
def test_r_force_equilibrium_residual_at_machine_precision(
    area_1: float, area_2: float, delta_t: float
) -> None:
    """R. N_1 + N_2 vanishes to machine precision and is reported, not hidden."""
    member_1 = AxialMember(HAND_MATERIAL_1, area_1)
    member_2 = AxialMember(HAND_MATERIAL_2, area_2)
    result = solve_bimetallic_joint(member_1, member_2, delta_t)
    scale = max(abs(result.member_1_force), abs(result.member_2_force))
    assert abs(result.force_equilibrium_residual) <= 1.0e-12 * scale
    assert result.force_equilibrium_residual == (
        result.member_1_force + result.member_2_force
    )


@pytest.mark.parametrize("delta_t", DELTA_TS)
def test_s_equal_areas_give_equal_and_opposite_stresses(delta_t: float) -> None:
    """S. With A_1 = A_2 the member stresses are equal and opposite."""
    member_1 = AxialMember(HAND_MATERIAL_1, 150.0e-6)
    member_2 = AxialMember(HAND_MATERIAL_2, 150.0e-6)
    result = solve_bimetallic_joint(member_1, member_2, delta_t)
    assert result.member_1_stress == pytest.approx(-result.member_2_stress, rel=1e-12)
    assert result.member_1_stress / result.member_2_stress == pytest.approx(-1.0, rel=1e-12)


@pytest.mark.parametrize("area_1, area_2", AREA_PAIRS)
def test_t_unequal_area_stress_ratio_identity(area_1: float, area_2: float) -> None:
    """T. Equilibrium forces sigma_1/sigma_2 = -A_2/A_1, not -1 in general."""
    member_1 = AxialMember(HAND_MATERIAL_1, area_1)
    member_2 = AxialMember(HAND_MATERIAL_2, area_2)
    result = solve_bimetallic_joint(member_1, member_2, 75.0)
    assert result.member_1_stress / result.member_2_stress == pytest.approx(
        -area_2 / area_1, rel=1e-12
    )


@pytest.mark.parametrize("area_1, area_2", AREA_PAIRS)
@pytest.mark.parametrize("delta_t", DELTA_TS)
def test_u_common_strain_lies_between_the_free_strains(
    area_1: float, area_2: float, delta_t: float
) -> None:
    """U. The stiffness-weighted strain is bracketed by the two free strains."""
    member_1 = AxialMember(HAND_MATERIAL_1, area_1)
    member_2 = AxialMember(HAND_MATERIAL_2, area_2)
    result = solve_bimetallic_joint(member_1, member_2, delta_t)
    lower = min(result.member_1_free_thermal_strain, result.member_2_free_thermal_strain)
    upper = max(result.member_1_free_thermal_strain, result.member_2_free_thermal_strain)
    tolerance = 1.0e-15 * max(abs(lower), abs(upper), 1.0e-12)
    assert lower - tolerance <= result.common_strain <= upper + tolerance


def test_solver_rejects_non_member_arguments() -> None:
    """The solver requires AxialMember inputs."""
    with pytest.raises(TypeError):
        solve_bimetallic_joint(HAND_MATERIAL_1, HAND_MEMBER_2, 50.0)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        solve_bimetallic_joint(HAND_MEMBER_1, HAND_MATERIAL_2, 50.0)  # type: ignore[arg-type]

"""Tests V-AC: limiting cases, symmetry, scaling laws and stiffness weighting."""

from __future__ import annotations

import pytest
from cases import HAND_MATERIAL_1, HAND_MATERIAL_2, HAND_MEMBER_1, HAND_MEMBER_2

from thermal_joint import AxialMember, ThermoelasticMaterial, solve_bimetallic_joint

DELTA_TS = [-140.0, -30.0, 12.5, 100.0]
AREA_RATIOS = [0.25, 0.5, 1.0, 2.0, 4.0]


@pytest.mark.parametrize("delta_t", DELTA_TS)
@pytest.mark.parametrize("alpha", [0.0, 8.5e-6, 23.0e-6, -3.0e-6])
def test_v_identical_cte_gives_zero_stress(delta_t: float, alpha: float) -> None:
    """V. Equal CTE means no mismatch, so no thermal stress at any dT."""
    member_1 = AxialMember(ThermoelasticMaterial("A", 70.0e9, alpha, 270.0e6), 100.0e-6)
    member_2 = AxialMember(ThermoelasticMaterial("B", 110.0e9, alpha, 830.0e6), 250.0e-6)
    result = solve_bimetallic_joint(member_1, member_2, delta_t)

    assert result.common_strain == pytest.approx(alpha * delta_t, rel=1e-15, abs=1e-30)
    assert result.member_1_stress == pytest.approx(0.0, abs=1e-6)
    assert result.member_2_stress == pytest.approx(0.0, abs=1e-6)
    assert result.member_1_force == pytest.approx(0.0, abs=1e-9)
    assert result.member_2_force == pytest.approx(0.0, abs=1e-9)


def test_w_zero_delta_temperature_gives_a_zero_state() -> None:
    """W. dT = 0 gives exactly zero strain, stress and force."""
    result = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, 0.0)
    assert result.common_strain == 0.0
    assert result.member_1_free_thermal_strain == 0.0
    assert result.member_2_free_thermal_strain == 0.0
    assert result.member_1_stress == 0.0
    assert result.member_2_stress == 0.0
    assert result.member_1_force == 0.0
    assert result.member_2_force == 0.0
    assert result.force_equilibrium_residual == 0.0


@pytest.mark.parametrize("delta_t", [30.0, 100.0, 140.0])
def test_x_heating_and_cooling_reverse_every_sign(delta_t: float) -> None:
    """X. Reversing dT reverses the sign of strain, stress and force."""
    heating = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, delta_t)
    cooling = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, -delta_t)

    assert cooling.common_strain == pytest.approx(-heating.common_strain, rel=1e-14)
    assert cooling.member_1_stress == pytest.approx(-heating.member_1_stress, rel=1e-14)
    assert cooling.member_2_stress == pytest.approx(-heating.member_2_stress, rel=1e-14)
    assert cooling.member_1_force == pytest.approx(-heating.member_1_force, rel=1e-14)
    assert cooling.member_2_force == pytest.approx(-heating.member_2_force, rel=1e-14)
    assert heating.member_1_stress * cooling.member_1_stress < 0.0


@pytest.mark.parametrize("delta_t", [30.0, 100.0, 140.0])
def test_y_stress_magnitudes_are_symmetric_in_delta_temperature(delta_t: float) -> None:
    """Y. Only the signs differ between +dT and -dT; magnitudes are identical."""
    heating = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, delta_t)
    cooling = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, -delta_t)
    assert abs(cooling.member_1_stress) == pytest.approx(abs(heating.member_1_stress), rel=1e-14)
    assert abs(cooling.member_2_stress) == pytest.approx(abs(heating.member_2_stress), rel=1e-14)


@pytest.mark.parametrize("delta_t", [25.0, -70.0])
def test_z_doubling_delta_temperature_doubles_strain_and_stress(delta_t: float) -> None:
    """Z. Stress and common strain are proportional to dT, signs preserved."""
    single = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, delta_t)
    double = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, 2.0 * delta_t)
    assert double.common_strain == pytest.approx(2.0 * single.common_strain, rel=1e-14)
    assert double.member_1_stress == pytest.approx(2.0 * single.member_1_stress, rel=1e-14)
    assert double.member_2_stress == pytest.approx(2.0 * single.member_2_stress, rel=1e-14)


@pytest.mark.parametrize("delta_t", [25.0, -70.0])
def test_aa_doubling_delta_temperature_doubles_internal_forces(delta_t: float) -> None:
    """AA. Internal forces are proportional to dT, signs preserved."""
    single = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, delta_t)
    double = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, 2.0 * delta_t)
    assert double.member_1_force == pytest.approx(2.0 * single.member_1_force, rel=1e-14)
    assert double.member_2_force == pytest.approx(2.0 * single.member_2_force, rel=1e-14)
    assert double.member_1_force * single.member_1_force > 0.0


@pytest.mark.parametrize("mismatch", [1.0e-6, 5.0e-6, 14.5e-6, 30.0e-6])
def test_ab_stress_magnitude_is_linear_in_cte_mismatch(mismatch: float) -> None:
    """AB. |sigma| scales linearly with |alpha_1 - alpha_2| at fixed E, A, dT."""
    alpha_2 = 8.5e-6
    delta_t = 100.0
    reference_mismatch = 1.0e-6

    def stress_pair(delta_alpha: float) -> tuple[float, float]:
        member_1 = AxialMember(
            ThermoelasticMaterial("A", 70.0e9, alpha_2 + delta_alpha, 270.0e6), 100.0e-6
        )
        member_2 = AxialMember(
            ThermoelasticMaterial("B", 110.0e9, alpha_2, 830.0e6), 250.0e-6
        )
        result = solve_bimetallic_joint(member_1, member_2, delta_t)
        return result.member_1_stress, result.member_2_stress

    base_1, base_2 = stress_pair(reference_mismatch)
    scaled_1, scaled_2 = stress_pair(mismatch)
    factor = mismatch / reference_mismatch
    assert scaled_1 == pytest.approx(factor * base_1, rel=1e-12)
    assert scaled_2 == pytest.approx(factor * base_2, rel=1e-12)


def _strain_weight_toward_member_1(result) -> float:
    """Fraction of the way from member 2's free strain to member 1's free strain."""
    gap = result.member_1_free_thermal_strain - result.member_2_free_thermal_strain
    return (result.common_strain - result.member_2_free_thermal_strain) / gap


@pytest.mark.parametrize("delta_t", [100.0, -140.0])
def test_ac_common_strain_shifts_toward_the_stiffer_ea_member_by_area(delta_t: float) -> None:
    """AC. Growing A_1/A_2 pulls the common strain toward member 1's free strain."""
    weights = []
    for ratio in AREA_RATIOS:
        member_1 = AxialMember(HAND_MATERIAL_1, ratio * 100.0e-6)
        member_2 = AxialMember(HAND_MATERIAL_2, 100.0e-6)
        result = solve_bimetallic_joint(member_1, member_2, delta_t)
        expected = member_1.axial_stiffness / (
            member_1.axial_stiffness + member_2.axial_stiffness
        )
        assert _strain_weight_toward_member_1(result) == pytest.approx(expected, rel=1e-12)
        weights.append(_strain_weight_toward_member_1(result))
    assert weights == sorted(weights)
    assert weights[0] < 0.5 < weights[-1]


@pytest.mark.parametrize("delta_t", [100.0, -140.0])
def test_ac_common_strain_shifts_toward_the_stiffer_ea_member_by_modulus(delta_t: float) -> None:
    """AC. Same behaviour when the modulus rather than the area is varied."""
    weights = []
    for modulus in (25.0e9, 50.0e9, 100.0e9, 200.0e9, 400.0e9):
        material_1 = ThermoelasticMaterial("A", modulus, 20.0e-6, 300.0e6)
        member_1 = AxialMember(material_1, 100.0e-6)
        member_2 = AxialMember(HAND_MATERIAL_2, 100.0e-6)
        result = solve_bimetallic_joint(member_1, member_2, delta_t)
        weights.append(_strain_weight_toward_member_1(result))
    assert weights == sorted(weights)
    assert weights[0] < 0.5 < weights[-1]


@pytest.mark.parametrize("ratio", AREA_RATIOS)
def test_area_ratio_redistributes_stress_but_keeps_forces_balanced(ratio: float) -> None:
    """Stress is not area-independent; the internal forces stay equal and opposite."""
    member_1 = AxialMember(HAND_MATERIAL_1, ratio * 100.0e-6)
    member_2 = AxialMember(HAND_MATERIAL_2, 100.0e-6)
    result = solve_bimetallic_joint(member_1, member_2, 100.0)
    assert result.member_1_force == pytest.approx(-result.member_2_force, rel=1e-12)


def test_area_ratio_changes_stress_magnitudes() -> None:
    """Doubling one area genuinely changes both member stresses."""
    small = solve_bimetallic_joint(
        AxialMember(HAND_MATERIAL_1, 50.0e-6), AxialMember(HAND_MATERIAL_2, 100.0e-6), 100.0
    )
    large = solve_bimetallic_joint(
        AxialMember(HAND_MATERIAL_1, 400.0e-6), AxialMember(HAND_MATERIAL_2, 100.0e-6), 100.0
    )
    assert abs(large.member_1_stress) < abs(small.member_1_stress)
    assert abs(large.member_2_stress) > abs(small.member_2_stress)

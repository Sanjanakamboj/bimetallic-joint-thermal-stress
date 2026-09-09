"""Sections 12, 15, 16 and 17: load sharing, sensitivity and the sign transition."""

from __future__ import annotations

import pytest

from thermal_joint import (
    AxialRestraint,
    joint_axial_stiffness,
    solve_bimetallic_joint,
    solve_restrained_joint,
    zero_stress_restraint_stiffness,
)
from thermal_joint.illustrative import aluminium_like_member, titanium_like_member

MEMBER_AL = aluminium_like_member(100.0)
MEMBER_TI = titanium_like_member(100.0)
K_JOINT = joint_axial_stiffness(MEMBER_AL, MEMBER_TI)
RATIOS = (0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0)
EXTREMES = (100.0, -140.0)


def solve_at(ratio: float, delta_t: float, reference_strain: float = 0.0):
    restraint = AxialRestraint.from_stiffness_ratio(
        ratio, MEMBER_AL, MEMBER_TI, reference_strain=reference_strain
    )
    return solve_restrained_joint(MEMBER_AL, MEMBER_TI, delta_t, restraint)


@pytest.mark.parametrize("ratio", [r for r in RATIOS if r > 0.0])
@pytest.mark.parametrize("delta_t", EXTREMES)
def test_member_forces_are_no_longer_equal_and_opposite(ratio: float, delta_t: float) -> None:
    """Section 12. Under restraint the two member forces stop balancing each other."""
    result = solve_at(ratio, delta_t)
    assert result.member_1_force + result.member_2_force != 0.0
    assert abs(result.member_1_force + result.member_2_force) > 1.0


@pytest.mark.parametrize("ratio", RATIOS)
@pytest.mark.parametrize("delta_t", EXTREMES)
def test_member_forces_balance_the_restraint_reaction(ratio: float, delta_t: float) -> None:
    """Section 12. Instead, N_1 + N_2 balances N_r under the documented convention."""
    result = solve_at(ratio, delta_t)
    total = result.member_1_force + result.member_2_force
    scale = max(abs(total), abs(result.restraint_force), 1.0)
    assert total == pytest.approx(result.restraint_force, rel=1e-11, abs=1e-11 * scale)


@pytest.mark.parametrize("delta_t", EXTREMES)
def test_common_strain_magnitude_decreases_with_restraint(delta_t: float) -> None:
    """Section 15. With eps_ref = 0 the common strain shrinks toward zero."""
    magnitudes = [abs(solve_at(r, delta_t).common_strain) for r in RATIOS]
    assert magnitudes == sorted(magnitudes, reverse=True)
    assert magnitudes[-1] < 0.05 * magnitudes[0]


@pytest.mark.parametrize("delta_t", EXTREMES)
def test_high_cte_member_stress_grows_monotonically_with_restraint(delta_t: float) -> None:
    """Section 15. The aluminium-like member is loaded harder as restraint rises."""
    magnitudes = [abs(solve_at(r, delta_t).member_1_stress) for r in RATIOS]
    assert magnitudes == sorted(magnitudes)


@pytest.mark.parametrize("delta_t", EXTREMES)
def test_low_cte_member_stress_is_not_monotonic_in_restraint(delta_t: float) -> None:
    """Section 15. The titanium-like member is the counter-example, as warned.

    Its stress magnitude falls to zero at the crossing and then grows again with
    the opposite sign, so stress magnitude is *not* monotone in every member.
    """
    magnitudes = [abs(solve_at(r, delta_t).member_2_stress) for r in RATIOS]
    assert magnitudes != sorted(magnitudes)
    assert min(magnitudes) < magnitudes[0]
    assert magnitudes[-1] > min(magnitudes)


@pytest.mark.parametrize("delta_t", EXTREMES)
def test_restraint_force_magnitude_grows_with_restraint(delta_t: float) -> None:
    """Section 15. More restraint means more reaction carried by the structure."""
    magnitudes = [abs(solve_at(r, delta_t).restraint_force) for r in RATIOS]
    assert magnitudes == sorted(magnitudes)
    assert magnitudes[0] == 0.0


def test_sign_transition_of_the_low_cte_member_on_heating() -> None:
    """Section 16. On heating the titanium-like member crosses tension -> compression."""
    free = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, 100.0)
    assert free.member_1_stress < 0.0 and free.member_2_stress > 0.0

    strong = solve_at(100.0, 100.0)
    assert strong.member_1_stress < 0.0 and strong.member_2_stress < 0.0


def test_sign_transition_of_the_low_cte_member_on_cooling() -> None:
    """Section 16. On cooling it crosses compression -> tension."""
    free = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, -140.0)
    assert free.member_1_stress > 0.0 and free.member_2_stress < 0.0

    strong = solve_at(100.0, -140.0)
    assert strong.member_1_stress > 0.0 and strong.member_2_stress > 0.0


def test_high_cte_member_never_changes_sign() -> None:
    """Section 16. The aluminium-like member keeps its sign at every restraint level."""
    for delta_t in EXTREMES:
        signs = {solve_at(r, delta_t).member_1_stress > 0.0 for r in RATIOS}
        assert len(signs) == 1


@pytest.mark.parametrize("delta_t", [100.0, -140.0, 37.0, -5.0])
def test_analytical_zero_stress_restraint_for_the_low_cte_member(delta_t: float) -> None:
    """Section 17. K_r = S/alpha_2 - K_joint drives sigma_2 exactly to zero."""
    crossing = zero_stress_restraint_stiffness(MEMBER_AL, MEMBER_TI, delta_t, 2)
    assert crossing is not None

    stiffness_sum = (
        MEMBER_AL.axial_stiffness * MEMBER_AL.material.thermal_expansion_coefficient
        + MEMBER_TI.axial_stiffness * MEMBER_TI.material.thermal_expansion_coefficient
    )
    expected = stiffness_sum / MEMBER_TI.material.thermal_expansion_coefficient - K_JOINT
    assert crossing == pytest.approx(expected, rel=1e-12)

    result = solve_restrained_joint(MEMBER_AL, MEMBER_TI, delta_t, AxialRestraint(crossing))
    assert result.member_2_stress == pytest.approx(0.0, abs=1e-3)
    assert result.member_1_stress != 0.0


def test_zero_stress_crossing_is_independent_of_delta_temperature() -> None:
    """Section 17. For eps_ref = 0 the crossing stiffness does not depend on dT."""
    crossings = {
        round(zero_stress_restraint_stiffness(MEMBER_AL, MEMBER_TI, delta_t, 2), 6)
        for delta_t in (100.0, -140.0, 1.0, -0.5, 250.0)
    }
    assert len(crossings) == 1
    ratio = crossings.pop() / K_JOINT
    assert ratio == pytest.approx(0.663399, rel=1e-5)


def test_no_crossing_exists_for_the_high_cte_member() -> None:
    """Section 17. The aluminium-like member has no admissible zero-stress restraint."""
    for delta_t in (100.0, -140.0):
        assert zero_stress_restraint_stiffness(MEMBER_AL, MEMBER_TI, delta_t, 1) is None


def test_crossing_brackets_the_observed_sign_change() -> None:
    """Section 17. The analytic crossing agrees with the numerically observed one."""
    crossing = zero_stress_restraint_stiffness(MEMBER_AL, MEMBER_TI, 100.0, 2)
    below = solve_restrained_joint(
        MEMBER_AL, MEMBER_TI, 100.0, AxialRestraint(0.99 * crossing)
    )
    above = solve_restrained_joint(
        MEMBER_AL, MEMBER_TI, 100.0, AxialRestraint(1.01 * crossing)
    )
    assert below.member_2_stress > 0.0
    assert above.member_2_stress < 0.0


def test_zero_stress_crossing_with_nonzero_reference_strain() -> None:
    """Section 17. The generalised crossing formula also holds for eps_ref != 0."""
    reference_strain = 4.0e-4
    crossing = zero_stress_restraint_stiffness(
        MEMBER_AL, MEMBER_TI, 100.0, 2, reference_strain=reference_strain
    )
    assert crossing is not None
    result = solve_restrained_joint(
        MEMBER_AL,
        MEMBER_TI,
        100.0,
        AxialRestraint(crossing, reference_strain=reference_strain),
    )
    assert result.member_2_stress == pytest.approx(0.0, abs=1e-3)


def test_zero_stress_helper_validates_member_index() -> None:
    with pytest.raises(ValueError):
        zero_stress_restraint_stiffness(MEMBER_AL, MEMBER_TI, 100.0, 3)


def test_zero_stress_helper_returns_none_for_a_degenerate_target() -> None:
    """A member with zero free thermal strain and eps_ref = 0 has no finite crossing."""
    from thermal_joint import AxialMember, ThermoelasticMaterial

    athermal = AxialMember(ThermoelasticMaterial("athermal", 100.0e9, 0.0, 300.0e6), 100.0e-6)
    assert zero_stress_restraint_stiffness(MEMBER_AL, athermal, 100.0, 2) is None

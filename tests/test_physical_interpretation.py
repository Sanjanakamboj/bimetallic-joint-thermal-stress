"""Section 29: the computed signs must reproduce the expected physics.

The aluminium-like member has the higher CTE. On heating it wants to expand
more than the titanium-like member, so the bond restrains it into compression
while the titanium-like member is pulled into tension. Cooling reverses both.

Every assertion below reads the sign off the computed result; nothing is
hard-coded into the model.
"""

from __future__ import annotations

import pytest

from thermal_joint import YieldBasis, assess_temperature_extremes, solve_bimetallic_joint
from thermal_joint.illustrative import (
    ALUMINIUM_LIKE,
    TITANIUM_LIKE,
    aluminium_like_member,
    radiator_joint_environment,
    titanium_like_member,
)

MEMBER_AL = aluminium_like_member(100.0)
MEMBER_TI = titanium_like_member(100.0)


def test_aluminium_like_member_has_the_higher_cte() -> None:
    """The premise of the interpretation, asserted rather than assumed."""
    assert (
        ALUMINIUM_LIKE.thermal_expansion_coefficient
        > TITANIUM_LIKE.thermal_expansion_coefficient
    )


def test_heating_compresses_the_high_cte_member() -> None:
    """On heating the restrained high-CTE member goes into compression."""
    result = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, 100.0)
    assert result.member_1_stress < 0.0
    assert result.member_2_stress > 0.0
    assert result.member_1_free_thermal_strain > result.common_strain
    assert result.member_2_free_thermal_strain < result.common_strain


def test_cooling_tensions_the_high_cte_member() -> None:
    """On cooling the high-CTE member wants to contract more, so it is tensioned."""
    result = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, -140.0)
    assert result.member_1_stress > 0.0
    assert result.member_2_stress < 0.0
    assert result.member_1_free_thermal_strain < result.common_strain
    assert result.member_2_free_thermal_strain > result.common_strain


def test_member_ordering_does_not_change_the_physics() -> None:
    """Swapping the member order swaps the reported stresses, nothing else."""
    forward = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, 100.0)
    reversed_order = solve_bimetallic_joint(MEMBER_TI, MEMBER_AL, 100.0)
    assert reversed_order.member_1_stress == pytest.approx(forward.member_2_stress, rel=1e-12)
    assert reversed_order.member_2_stress == pytest.approx(forward.member_1_stress, rel=1e-12)
    assert reversed_order.common_strain == pytest.approx(forward.common_strain, rel=1e-14)


def test_representative_sanity_case_regression() -> None:
    """Regression lock on the representative equal-area illustrative case.

    Values are the closed-form results for the illustrative inputs:
    70/110 GPa, 23e-6 / 8.5e-6 per K, 100 mm^2 each, dT_hot = +100 K,
    dT_cold = -140 K, design factor 1.25.
    """
    assessment = assess_temperature_extremes(
        MEMBER_AL, MEMBER_TI, radiator_joint_environment(), YieldBasis(1.25)
    )

    assert assessment.hot_result.member_1_stress == pytest.approx(-62.0278e6, rel=1e-5)
    assert assessment.hot_result.member_2_stress == pytest.approx(+62.0278e6, rel=1e-5)
    assert assessment.cold_result.member_1_stress == pytest.approx(+86.8389e6, rel=1e-5)
    assert assessment.cold_result.member_2_stress == pytest.approx(-86.8389e6, rel=1e-5)

    assert assessment.governing_extreme == "cold"
    assert assessment.governing_member == 1
    assert assessment.minimum_yield_margin == pytest.approx(1.4874, rel=1e-3)
    assert assessment.passes


def test_representative_case_forces_balance() -> None:
    """Equal areas mean equal and opposite stresses and forces."""
    assessment = assess_temperature_extremes(
        MEMBER_AL, MEMBER_TI, radiator_joint_environment(), YieldBasis(1.25)
    )
    for result in (assessment.cold_result, assessment.hot_result):
        scale = max(abs(result.member_1_force), abs(result.member_2_force))
        assert abs(result.force_equilibrium_residual) <= 1.0e-12 * scale


def test_illustrative_inputs_are_labelled_as_non_allowable() -> None:
    """Illustrative properties must never read as qualified design allowables."""
    from thermal_joint import illustrative

    assert "NOT DESIGN ALLOWABLE" in illustrative.ILLUSTRATIVE_NOTE
    assert "illustrative" in ALUMINIUM_LIKE.name.lower()
    assert "illustrative" in TITANIUM_LIKE.name.lower()

"""Tests AD-AK: preliminary elastic yield allowables, margins and governing member."""

from __future__ import annotations

import math

import pytest
from cases import HAND_DELTA_T, HAND_MEMBER_1, HAND_MEMBER_2

from thermal_joint import (
    AxialMember,
    ThermalJointResult,
    ThermoelasticMaterial,
    YieldBasis,
    assess_yield,
    solve_bimetallic_joint,
)

MATERIAL = ThermoelasticMaterial("M", 70.0e9, 23.0e-6, 300.0e6)


def test_ad_allowable_hand_calculation() -> None:
    """AD. sigma_allow = 300 MPa / 1.5 = 200 MPa."""
    assert YieldBasis(1.5).allowable_stress(MATERIAL) == pytest.approx(200.0e6, rel=1e-15)


def test_ad_default_design_factor_is_unity() -> None:
    """AD. The default basis leaves the yield strength untouched."""
    basis = YieldBasis()
    assert basis.design_factor == 1.0
    assert basis.allowable_stress(MATERIAL) == MATERIAL.yield_strength


@pytest.mark.parametrize("bad_factor", [0.0, 0.999, -1.0, -2.5])
def test_ae_design_factor_below_one_rejected(bad_factor: float) -> None:
    """AE. A design factor below 1 would raise the allowable above yield."""
    with pytest.raises(ValueError):
        YieldBasis(bad_factor)


@pytest.mark.parametrize("bad_factor", [math.nan, math.inf, -math.inf])
def test_ae_non_finite_design_factor_rejected(bad_factor: float) -> None:
    """AE. The design factor must be finite."""
    with pytest.raises(ValueError):
        YieldBasis(bad_factor)


@pytest.mark.parametrize("good_factor", [1.0, 1.25, 1.5, 2.0])
def test_ae_valid_design_factors_accepted(good_factor: float) -> None:
    """AE. Any finite factor >= 1 is a valid design basis."""
    assert YieldBasis(good_factor).design_factor == good_factor


@pytest.mark.parametrize("stress", [-50.0e6, 50.0e6])
def test_af_margin_hand_calculation(stress: float) -> None:
    """AF. MS = (300/1.5) / 50 - 1 = 3, for tension and compression alike."""
    assert YieldBasis(1.5).margin_of_safety(MATERIAL, stress) == pytest.approx(3.0, rel=1e-14)


def test_ag_exact_boundary_passes() -> None:
    """AG. A stress exactly at the allowable gives MS = 0 and passes."""
    boundary_material = ThermoelasticMaterial("boundary", 70.0e9, 23.0e-6, 120.0e6)
    basis = YieldBasis()
    margin = basis.margin_of_safety(boundary_material, -120.0e6)
    assert margin == 0.0
    assert margin >= 0.0


def test_ag_exact_boundary_passes_through_the_assessment() -> None:
    """AG. The boundary also passes when reached through a solved joint."""
    result = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, HAND_DELTA_T)
    tuned_1 = AxialMember(
        ThermoelasticMaterial(
            "tuned-1",
            HAND_MEMBER_1.material.elastic_modulus,
            HAND_MEMBER_1.material.thermal_expansion_coefficient,
            abs(result.member_1_stress),
        ),
        HAND_MEMBER_1.area,
    )
    assessment = assess_yield(result, tuned_1, HAND_MEMBER_2, YieldBasis())
    assert assessment.member_1_margin.margin_of_safety == 0.0
    assert assessment.member_1_margin.passes
    assert assessment.passes


def test_ah_overstressed_member_fails() -> None:
    """AH. A stress above the allowable gives a negative margin and fails."""
    basis = YieldBasis(1.5)
    margin = basis.margin_of_safety(MATERIAL, -250.0e6)
    assert margin == pytest.approx(200.0 / 250.0 - 1.0, rel=1e-14)
    assert margin < 0.0

    weak = AxialMember(ThermoelasticMaterial("weak", 100.0e9, 20.0e-6, 10.0e6), 200.0e-6)
    result = solve_bimetallic_joint(weak, HAND_MEMBER_2, HAND_DELTA_T)
    assessment = assess_yield(result, weak, HAND_MEMBER_2, YieldBasis())
    assert assessment.minimum_margin < 0.0
    assert not assessment.passes
    assert not assessment.member_1_margin.passes


def test_ai_zero_stress_gives_infinite_non_governing_margin() -> None:
    """AI. Zero stress never divides by zero; it yields an infinite margin."""
    assert YieldBasis(2.0).margin_of_safety(MATERIAL, 0.0) == math.inf

    alpha = 15.0e-6
    member_1 = AxialMember(ThermoelasticMaterial("A", 70.0e9, alpha, 270.0e6), 100.0e-6)
    member_2 = AxialMember(ThermoelasticMaterial("B", 110.0e9, alpha, 830.0e6), 300.0e-6)
    result = solve_bimetallic_joint(member_1, member_2, 100.0)
    assessment = assess_yield(result, member_1, member_2, YieldBasis())
    assert assessment.minimum_margin == math.inf
    assert assessment.passes


def test_aj_governing_member_is_the_one_with_the_smaller_margin() -> None:
    """AJ. The higher-yield member can govern; the lower yield does not decide."""
    strong_but_small = AxialMember(
        ThermoelasticMaterial("low yield, large area", 70.0e9, 23.0e-6, 100.0e6), 1000.0e-6
    )
    weak_but_stressed = AxialMember(
        ThermoelasticMaterial("high yield, small area", 110.0e9, 8.5e-6, 500.0e6), 10.0e-6
    )
    result = solve_bimetallic_joint(strong_but_small, weak_but_stressed, 100.0)
    assessment = assess_yield(result, strong_but_small, weak_but_stressed, YieldBasis())

    assert abs(result.member_2_stress) > abs(result.member_1_stress)
    assert weak_but_stressed.material.yield_strength > strong_but_small.material.yield_strength
    assert assessment.governing_member == 2
    assert assessment.minimum_margin == assessment.member_2_margin.margin_of_safety
    assert assessment.governing_margin is assessment.member_2_margin


def test_aj_governing_member_tie_resolves_to_member_one() -> None:
    """AJ. A bit-exact tie in the two margins resolves to member 1."""
    material = ThermoelasticMaterial("symmetric", 100.0e9, 20.0e-6, 300.0e6)
    member_1 = AxialMember(material, 100.0e-6, label="tie-1")
    member_2 = AxialMember(material, 100.0e-6, label="tie-2")
    tied = ThermalJointResult(
        delta_temperature=100.0,
        common_strain=0.0,
        member_1_free_thermal_strain=0.0,
        member_2_free_thermal_strain=0.0,
        member_1_stress=-100.0e6,
        member_2_stress=+100.0e6,
        member_1_force=-10.0e3,
        member_2_force=+10.0e3,
        force_equilibrium_residual=0.0,
    )
    assessment = assess_yield(tied, member_1, member_2, YieldBasis())
    assert (
        assessment.member_1_margin.margin_of_safety
        == assessment.member_2_margin.margin_of_safety
    )
    assert assessment.governing_member == 1


def test_ak_material_carries_no_hidden_design_factor() -> None:
    """AK. Design policy lives in YieldBasis only; the material stores raw yield."""
    assert YieldBasis(1.0).allowable_stress(MATERIAL) == MATERIAL.yield_strength
    fields = set(vars(MATERIAL))
    assert not any("factor" in name or "safety" in name or "allow" in name for name in fields)
    assert fields == {
        "name",
        "elastic_modulus",
        "thermal_expansion_coefficient",
        "yield_strength",
        "density",
    }


def test_ak_design_factor_scales_the_allowable_only() -> None:
    """AK. Changing the basis changes the allowable, never the stored yield."""
    assert YieldBasis(2.0).allowable_stress(MATERIAL) == pytest.approx(
        0.5 * YieldBasis(1.0).allowable_stress(MATERIAL), rel=1e-15
    )
    assert MATERIAL.yield_strength == 300.0e6


def test_assessment_records_labels_and_allowables() -> None:
    """Per-member records carry the label, signed stress and allowable used."""
    result = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, HAND_DELTA_T)
    assessment = assess_yield(result, HAND_MEMBER_1, HAND_MEMBER_2, YieldBasis(1.5))
    assert assessment.member_1_margin.label == "hand-1"
    assert assessment.member_1_margin.stress == result.member_1_stress
    assert assessment.member_1_margin.allowable_stress == pytest.approx(200.0e6, rel=1e-14)


def test_assess_yield_defaults_to_unity_basis() -> None:
    """Omitting the basis is the same as passing YieldBasis()."""
    result = solve_bimetallic_joint(HAND_MEMBER_1, HAND_MEMBER_2, HAND_DELTA_T)
    assert assess_yield(result, HAND_MEMBER_1, HAND_MEMBER_2) == assess_yield(
        result, HAND_MEMBER_1, HAND_MEMBER_2, YieldBasis(1.0)
    )

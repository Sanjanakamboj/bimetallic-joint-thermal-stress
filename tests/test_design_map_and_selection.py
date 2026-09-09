"""Tests AX-BM: the bounded design map, the selection policy and prior-milestone regression."""

from __future__ import annotations

import pytest
from design_cases import (
    ADHESIVE,
    ADHESIVE_BASIS,
    BASELINE_COLD_MARGIN,
    BASELINE_COLD_PEAK,
    BASELINE_HOT_PEAK,
    BASELINE_YIELD_MARGIN,
    ENVIRONMENT,
    GEOMETRY,
    MEMBER_AL,
    MEMBER_TI,
    MM,
    YIELD_BASIS,
)

from thermal_joint import (
    AdhesiveMaterial,
    AxialRestraint,
    BondedOverlapGeometry,
    DesignSelectionPolicy,
    JointDesignCandidate,
    assess_restrained_temperature_extremes,
    assess_shear_lag_extremes,
    assess_temperature_extremes,
    evaluate_joint_design,
    select_preliminary_joint_design,
    solve_bimetallic_joint,
    solve_shear_lag,
    width_thickness_design_map,
)
from thermal_joint.illustrative import ADHESIVE_LIKE

MAP_WIDTHS_MM = (10.0, 20.0, 30.0, 40.0, 50.0, 75.0, 100.0)
MAP_THICKNESSES_MM = (0.1, 0.2, 0.3, 0.5, 0.75, 1.0)


def canonical_map(adhesive=ADHESIVE):
    return width_thickness_design_map(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, adhesive,
        [w * MM for w in MAP_WIDTHS_MM],
        [t * MM for t in MAP_THICKNESSES_MM],
        YIELD_BASIS, ADHESIVE_BASIS,
    )


def test_ax_every_map_point_matches_a_direct_package_evaluation() -> None:
    """AX. No point is shortcut; each is the full verified evaluation."""
    for candidate in canonical_map():
        geometry = BondedOverlapGeometry(
            GEOMETRY.overlap_length, candidate.bond_width, candidate.adhesive_thickness
        )
        expected = evaluate_joint_design(
            MEMBER_AL, MEMBER_TI, ENVIRONMENT, geometry, ADHESIVE, YIELD_BASIS, ADHESIVE_BASIS
        )
        assert candidate == expected
        direct = assess_shear_lag_extremes(
            MEMBER_AL, MEMBER_TI, ENVIRONMENT, geometry, ADHESIVE, ADHESIVE_BASIS
        )
        assert candidate.peak_shear_stress == direct.governing_peak_shear_stress


def test_ax_map_shape_and_ordering_are_deterministic() -> None:
    """AX. Row-major over width then thickness, one point per grid cell."""
    candidates = canonical_map()
    assert len(candidates) == len(MAP_WIDTHS_MM) * len(MAP_THICKNESSES_MM) == 42
    expected = [
        (w * MM, t * MM) for w in MAP_WIDTHS_MM for t in MAP_THICKNESSES_MM
    ]
    actual = [(c.bond_width, c.adhesive_thickness) for c in candidates]
    for (want_w, want_t), (got_w, got_t) in zip(expected, actual):
        assert got_w == pytest.approx(want_w, rel=1e-12)
        assert got_t == pytest.approx(want_t, rel=1e-12)
    assert canonical_map() == canonical_map()


def test_ay_canonical_baseline_point_is_present_and_unchanged() -> None:
    """AY. The 20 mm / 0.2 mm cell reproduces the Milestone 3 baseline exactly."""
    baseline = [
        c for c in canonical_map()
        if c.bond_width == pytest.approx(20.0 * MM)
        and c.adhesive_thickness == pytest.approx(0.2 * MM)
    ][0]
    assert baseline.cold_peak_shear_stress == pytest.approx(BASELINE_COLD_PEAK, rel=1e-5)
    assert baseline.hot_peak_shear_stress == pytest.approx(BASELINE_HOT_PEAK, rel=1e-5)
    assert baseline.adhesive_margin == pytest.approx(BASELINE_COLD_MARGIN, rel=1e-3)
    assert baseline.minimum_yield_margin == pytest.approx(BASELINE_YIELD_MARGIN, rel=1e-3)
    assert not baseline.overall_feasible


def test_az_map_never_modifies_the_adhesive_strength() -> None:
    """AZ. Only geometry varies; the allowable is identical at every point."""
    allowables = {c.allowable_shear_stress for c in canonical_map()}
    assert allowables == {ADHESIVE_BASIS.allowable_shear_stress(ADHESIVE)}
    assert ADHESIVE_LIKE.shear_strength == 25.0e6
    assert ADHESIVE_LIKE.shear_modulus == 1.0e9


def test_ba_feasibility_is_yield_and_adhesive() -> None:
    """BA. Overall feasibility is the boolean AND of the two screens."""
    for candidate in canonical_map():
        assert candidate.overall_feasible == (
            candidate.yield_feasible and candidate.adhesive_feasible
        )
        assert candidate.yield_feasible is True
    assert any(not c.adhesive_feasible for c in canonical_map())


def test_bb_margins_stay_separate_on_every_candidate() -> None:
    """BB. No combined scalar margin is exposed anywhere on the candidate."""
    candidate = canonical_map()[0]
    assert candidate.minimum_yield_margin != candidate.adhesive_margin
    attributes = set(dir(candidate))
    assert not {"combined_margin", "overall_margin", "total_margin"} & attributes


def test_canonical_map_feasible_point_count_and_membership() -> None:
    """Six of the 42 grid points pass at G_a = 1.0 GPa - computed, not assumed."""
    feasible = [c for c in canonical_map() if c.overall_feasible]
    assert len(feasible) == 6
    combinations = {
        (round(c.bond_width / MM), round(c.adhesive_thickness / MM, 2)) for c in feasible
    }
    assert combinations == {
        (50, 1.0), (75, 0.75), (75, 1.0), (100, 0.5), (100, 0.75), (100, 1.0)
    }


def test_bc_selection_is_deterministic() -> None:
    """BC. Repeated selection returns the same candidate."""
    candidates = canonical_map()
    assert select_preliminary_joint_design(candidates) is not None
    first = select_preliminary_joint_design(candidates)
    second = select_preliminary_joint_design(canonical_map())
    assert first == second


def test_bd_selection_is_independent_of_input_order() -> None:
    """BD. The policy ranks on values, so shuffling the grid changes nothing."""
    candidates = list(canonical_map())
    reversed_selection = select_preliminary_joint_design(list(reversed(candidates)))
    rotated = candidates[17:] + candidates[:17]
    assert select_preliminary_joint_design(candidates) == reversed_selection
    assert select_preliminary_joint_design(rotated) == reversed_selection


def test_be_returns_none_when_nothing_is_feasible() -> None:
    """BE. A grid with no feasible point returns None rather than a best-effort pick."""
    narrow = width_thickness_design_map(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE,
        [10.0 * MM, 20.0 * MM], [0.1 * MM, 0.2 * MM], YIELD_BASIS, ADHESIVE_BASIS,
    )
    assert all(not c.overall_feasible for c in narrow)
    assert select_preliminary_joint_design(narrow) is None
    assert select_preliminary_joint_design([]) is None


def test_bf_selected_point_really_passes_both_screens() -> None:
    """BF. The selected candidate is re-verified from the package mechanics."""
    selected = select_preliminary_joint_design(canonical_map())
    assert selected.overall_feasible
    assert selected.adhesive_margin >= 0.0
    assert selected.minimum_yield_margin >= 0.0

    geometry = BondedOverlapGeometry(
        GEOMETRY.overlap_length, selected.bond_width, selected.adhesive_thickness
    )
    shear = assess_shear_lag_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, geometry, ADHESIVE, ADHESIVE_BASIS
    )
    assert shear.passes
    assert shear.minimum_margin == pytest.approx(selected.adhesive_margin, rel=1e-12)


def test_bg_selected_point_obeys_the_stated_tie_break_policy() -> None:
    """BG. Smallest bond area first; the 50 mm / 1.0 mm cell wins."""
    candidates = canonical_map()
    selected = select_preliminary_joint_design(candidates)
    feasible = [c for c in candidates if c.overall_feasible]

    assert selected.bond_width == pytest.approx(50.0 * MM, rel=1e-12)
    assert selected.adhesive_thickness == pytest.approx(1.0 * MM, rel=1e-12)
    assert selected.bond_area == min(c.bond_area for c in feasible)
    assert selected.adhesive_margin == pytest.approx(0.0648, abs=1e-3)
    assert selected.hot_adhesive_margin == pytest.approx(0.4907, abs=1e-3)
    assert selected.cold_adhesive_margin == pytest.approx(0.0648, abs=1e-3)
    assert selected.governing_extreme == "cold"
    assert selected.minimum_yield_margin == pytest.approx(BASELINE_YIELD_MARGIN, rel=1e-3)


def test_bg_thickness_tie_break_applies_at_equal_bond_area() -> None:
    """BG. With equal areas the thinner bondline wins, per the documented policy."""
    geometry_a = BondedOverlapGeometry(GEOMETRY.overlap_length, 0.100, 0.5 * MM)
    geometry_b = BondedOverlapGeometry(GEOMETRY.overlap_length, 0.100, 1.0 * MM)
    thin = evaluate_joint_design(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, geometry_a, ADHESIVE, YIELD_BASIS, ADHESIVE_BASIS
    )
    thick = evaluate_joint_design(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, geometry_b, ADHESIVE, YIELD_BASIS, ADHESIVE_BASIS
    )
    assert thin.bond_area == thick.bond_area
    assert thin.overall_feasible and thick.overall_feasible
    assert thick.adhesive_margin > thin.adhesive_margin  # the thicker point is safer
    # ...yet the policy still prefers the thinner bondline at equal area.
    assert select_preliminary_joint_design([thick, thin]) == thin
    assert select_preliminary_joint_design([thin, thick]) == thin


def test_unsupported_policy_is_rejected() -> None:
    with pytest.raises(ValueError):
        select_preliminary_joint_design(canonical_map(), "not-a-policy")  # type: ignore[arg-type]
    assert DesignSelectionPolicy.MINIMUM_BOND_AREA.value == "minimum_bond_area"


def test_soft_adhesive_map_demonstrates_co_design() -> None:
    """A second illustrative map at 0.1 GPa opens up far more of the grid."""
    soft = AdhesiveMaterial(
        "Structural-adhesive-like, compliant variant (illustrative)",
        0.1e9, ADHESIVE.shear_strength,
        "ILLUSTRATIVE ADHESIVE-LIKE INPUT - NOT DESIGN ALLOWABLE",
    )
    candidates = canonical_map(soft)
    feasible = [c for c in candidates if c.overall_feasible]
    assert len(feasible) == 32
    selected = select_preliminary_joint_design(candidates)
    assert selected.bond_width == pytest.approx(20.0 * MM, rel=1e-12)
    assert selected.adhesive_thickness == pytest.approx(0.3 * MM, rel=1e-12)
    assert selected.adhesive_margin == pytest.approx(0.0717, abs=1e-3)
    # The canonical stiff adhesive is untouched by building the variant.
    assert ADHESIVE_LIKE.shear_modulus == 1.0e9


# ------------------------------------------------------- regression BH-BM

def test_bh_milestone_1_outputs_unchanged() -> None:
    """BH. Free-joint stresses and yield margins are exactly as before."""
    hot = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, 100.0)
    cold = solve_bimetallic_joint(MEMBER_AL, MEMBER_TI, -140.0)
    assert hot.member_1_stress == pytest.approx(-62.0278e6, rel=1e-5)
    assert cold.member_1_stress == pytest.approx(+86.8389e6, rel=1e-5)

    assessment = assess_temperature_extremes(MEMBER_AL, MEMBER_TI, ENVIRONMENT, YIELD_BASIS)
    assert assessment.governing_extreme == "cold"
    assert assessment.minimum_yield_margin == pytest.approx(1.4874, rel=1e-3)


def test_bi_milestone_2_outputs_unchanged() -> None:
    """BI. Restrained results and the restraint boundary are untouched."""
    from thermal_joint import (
        RestraintLimitStatus,
        maximum_allowable_restraint_stiffness,
        zero_stress_restraint_stiffness,
    )

    restraint = AxialRestraint.from_stiffness_ratio(1.0, MEMBER_AL, MEMBER_TI)
    restrained = assess_restrained_temperature_extremes(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, restraint, YIELD_BASIS
    )
    assert restrained.cold_result.member_1_stress == pytest.approx(156.1194e6, rel=1e-5)
    assert restrained.cold_result.member_2_stress == pytest.approx(22.0306e6, rel=1e-5)
    assert restrained.minimum_yield_margin == pytest.approx(0.3836, rel=1e-3)

    limit = maximum_allowable_restraint_stiffness(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, YIELD_BASIS
    )
    assert limit.status is RestraintLimitStatus.FINITE_LIMIT
    assert limit.maximum_stiffness_ratio == pytest.approx(13.740543, abs=1e-5)
    assert zero_stress_restraint_stiffness(
        MEMBER_AL, MEMBER_TI, 100.0, 2
    ) / 1.8e7 == pytest.approx(0.663399, rel=1e-5)


def test_bj_milestone_3_baseline_shear_lag_unchanged() -> None:
    """BJ. The canonical shear-lag numbers are bit-for-bit the Milestone 3 ones."""
    cold = solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE, -140.0)
    assert cold.beta == pytest.approx(152.894157, rel=1e-8)
    assert cold.transfer_length == pytest.approx(6.5405e-3, rel=1e-4)
    assert cold.dimensionless_overlap == pytest.approx(6.115766, rel=1e-6)
    assert cold.transferred_force == pytest.approx(8683.8889, rel=1e-6)
    assert cold.peak_shear_stress == pytest.approx(BASELINE_COLD_PEAK, rel=1e-5)
    assert cold.average_transfer_shear == pytest.approx(10.8549e6, rel=1e-4)

    from thermal_joint import OverlapLimitStatus, required_overlap_length

    overlap = required_overlap_length(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, ADHESIVE_BASIS
    )
    assert overlap.status is OverlapLimitStatus.NO_FINITE_LENGTH_WITHIN_MODEL


def test_bm_combined_compliant_case_reproduces_the_milestone_3_pass() -> None:
    """BM. 20 mm / 1.0 mm / 0.1 GPa still gives 13.435 MPa and MS +0.489."""
    soft = AdhesiveMaterial(
        "compliant variant", 0.1e9, ADHESIVE.shear_strength, "illustrative fixture"
    )
    geometry = GEOMETRY.with_adhesive_thickness(1.0 * MM)
    candidate = evaluate_joint_design(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, geometry, soft, YIELD_BASIS, ADHESIVE_BASIS
    )
    assert candidate.cold_peak_shear_stress == pytest.approx(13.4355e6, rel=1e-4)
    assert candidate.adhesive_margin == pytest.approx(0.4886, abs=1e-3)
    assert candidate.overall_feasible


def test_evaluate_joint_design_defaults_to_unity_bases() -> None:
    """Omitting either basis falls back to a design factor of 1.0."""
    from thermal_joint import AdhesiveShearBasis, YieldBasis

    assert evaluate_joint_design(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE
    ) == evaluate_joint_design(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE,
        YieldBasis(1.0), AdhesiveShearBasis(1.0),
    )


def test_candidate_reports_bond_area_and_adhesive_volume_geometrically() -> None:
    """Section 24/25: geometry diagnostics only - no adhesive mass is invented."""
    candidate = evaluate_joint_design(
        MEMBER_AL, MEMBER_TI, ENVIRONMENT, GEOMETRY, ADHESIVE, YIELD_BASIS, ADHESIVE_BASIS
    )
    assert candidate.bond_area == pytest.approx(GEOMETRY.bond_area, rel=1e-15)
    assert candidate.adhesive_volume == pytest.approx(
        GEOMETRY.bond_area * GEOMETRY.adhesive_thickness, rel=1e-15
    )
    attributes = set(dir(candidate))
    assert not {"adhesive_mass", "mass", "bond_mass"} & attributes
    assert not hasattr(JointDesignCandidate, "adhesive_density")

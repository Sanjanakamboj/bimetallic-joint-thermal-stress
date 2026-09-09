"""Tests BD-BI: overlap, thickness, modulus, area-ratio and dT sensitivity sweeps."""

from __future__ import annotations

import pytest

from thermal_joint import (
    AdhesiveMaterial,
    AdhesiveShearBasis,
    AxialMember,
    shear_lag_parameter,
    solve_shear_lag,
    transfer_length,
)
from thermal_joint.illustrative import (
    ADHESIVE_LIKE,
    ALUMINIUM_LIKE,
    aluminium_like_member,
    radiator_joint_overlap,
    titanium_like_member,
)

MEMBER_AL = aluminium_like_member(100.0)
MEMBER_TI = titanium_like_member(100.0)
GEOMETRY = radiator_joint_overlap()
BASIS = AdhesiveShearBasis(1.25)
COLD = -140.0

OVERLAP_LENGTHS_MM = (5.0, 10.0, 20.0, 40.0, 80.0, 160.0)
THICKNESSES_MM = (0.05, 0.10, 0.20, 0.50, 1.00)
MODULI_GPA = (0.1, 0.25, 0.5, 1.0, 2.0, 5.0)
AREA_RATIOS = (0.25, 0.5, 1.0, 2.0, 4.0)


def test_bd_overlap_length_trend_is_monotonic() -> None:
    """BD. Peak shear falls monotonically with overlap; average falls faster."""
    peaks, averages, ratios = [], [], []
    for length_mm in OVERLAP_LENGTHS_MM:
        demand = solve_shear_lag(
            MEMBER_AL,
            MEMBER_TI,
            GEOMETRY.with_overlap_length(length_mm * 1.0e-3),
            ADHESIVE_LIKE,
            COLD,
        )
        peaks.append(demand.peak_shear_stress)
        averages.append(demand.average_transfer_shear)
        ratios.append(demand.peak_to_average_ratio)

    assert peaks == sorted(peaks, reverse=True)
    assert averages == sorted(averages, reverse=True)
    assert ratios == sorted(ratios)
    # Average shear collapses by 32x over the sweep; the peak barely moves.
    assert averages[0] / averages[-1] == pytest.approx(32.0, rel=1e-9)
    assert peaks[0] / peaks[-1] < 1.6


def test_bd_transferred_force_is_independent_of_overlap_length() -> None:
    """BD. The Milestone 1 demand does not depend on how it is transferred."""
    forces = {
        round(
            solve_shear_lag(
                MEMBER_AL,
                MEMBER_TI,
                GEOMETRY.with_overlap_length(length_mm * 1.0e-3),
                ADHESIVE_LIKE,
                COLD,
            ).transferred_force,
            9,
        )
        for length_mm in OVERLAP_LENGTHS_MM
    }
    assert len(forces) == 1


def test_be_adhesive_thickness_trend() -> None:
    """BE. Thicker bondline -> smaller beta -> longer transfer length -> lower peak."""
    betas, lengths, peaks = [], [], []
    for thickness_mm in THICKNESSES_MM:
        geometry = GEOMETRY.with_adhesive_thickness(thickness_mm * 1.0e-3)
        betas.append(shear_lag_parameter(MEMBER_AL, MEMBER_TI, geometry, ADHESIVE_LIKE))
        lengths.append(transfer_length(MEMBER_AL, MEMBER_TI, geometry, ADHESIVE_LIKE))
        peaks.append(
            solve_shear_lag(MEMBER_AL, MEMBER_TI, geometry, ADHESIVE_LIKE, COLD).peak_shear_stress
        )

    assert betas == sorted(betas, reverse=True)
    assert lengths == sorted(lengths)
    assert peaks == sorted(peaks, reverse=True)
    assert peaks[-1] < peaks[0]


def test_bf_adhesive_modulus_trend() -> None:
    """BF. Stiffer adhesive -> larger beta -> shorter transfer length -> higher peak."""
    betas, lengths, peaks = [], [], []
    for modulus_gpa in MODULI_GPA:
        adhesive = AdhesiveMaterial(
            f"G={modulus_gpa} GPa",
            modulus_gpa * 1.0e9,
            ADHESIVE_LIKE.shear_strength,
            "illustrative sweep fixture",
        )
        betas.append(shear_lag_parameter(MEMBER_AL, MEMBER_TI, GEOMETRY, adhesive))
        lengths.append(transfer_length(MEMBER_AL, MEMBER_TI, GEOMETRY, adhesive))
        peaks.append(
            solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, adhesive, COLD).peak_shear_stress
        )

    assert betas == sorted(betas)
    assert lengths == sorted(lengths, reverse=True)
    assert peaks == sorted(peaks)


def test_bf_bondline_compliance_is_the_effective_lever_but_not_a_single_fix() -> None:
    """BF. Softening the bondline helps a lot, yet no single sweep point passes.

    At the canonical 40 mm overlap every point of the modulus sweep and every
    point of the thickness sweep still fails the illustrative allowable. Only
    combining a soft *and* thick bondline (0.1 GPa over 1.0 mm) gets the peak
    below it. Reported as computed - the levers are real but not individually
    sufficient here.
    """
    stiff_peak = solve_shear_lag(
        MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, COLD
    ).peak_shear_stress
    soft = AdhesiveMaterial(
        "G=0.1 GPa", 0.1e9, ADHESIVE_LIKE.shear_strength, "illustrative sweep fixture"
    )
    soft_peak = solve_shear_lag(MEMBER_AL, MEMBER_TI, GEOMETRY, soft, COLD).peak_shear_stress

    assert soft_peak < stiff_peak / 3.0
    assert BASIS.margin_of_safety(ADHESIVE_LIKE, stiff_peak) < 0.0
    assert BASIS.margin_of_safety(soft, soft_peak) < 0.0  # still fails on its own

    thick_peak = solve_shear_lag(
        MEMBER_AL, MEMBER_TI, GEOMETRY.with_adhesive_thickness(1.0e-3), ADHESIVE_LIKE, COLD
    ).peak_shear_stress
    assert BASIS.margin_of_safety(ADHESIVE_LIKE, thick_peak) < 0.0  # also fails alone

    combined_peak = solve_shear_lag(
        MEMBER_AL, MEMBER_TI, GEOMETRY.with_adhesive_thickness(1.0e-3), soft, COLD
    ).peak_shear_stress
    assert BASIS.margin_of_safety(soft, combined_peak) > 0.0
    assert combined_peak == pytest.approx(13.435e6, rel=1e-3)


def test_bf_no_single_variable_sweep_point_passes_at_the_canonical_overlap() -> None:
    """BF. Stated explicitly so the failing screen cannot be read as a near miss."""
    for modulus_gpa in MODULI_GPA:
        adhesive = AdhesiveMaterial(
            f"G={modulus_gpa} GPa",
            modulus_gpa * 1.0e9,
            ADHESIVE_LIKE.shear_strength,
            "illustrative sweep fixture",
        )
        peak = solve_shear_lag(
            MEMBER_AL, MEMBER_TI, GEOMETRY, adhesive, COLD
        ).peak_shear_stress
        assert BASIS.margin_of_safety(adhesive, peak) < 0.0

    for thickness_mm in THICKNESSES_MM:
        peak = solve_shear_lag(
            MEMBER_AL,
            MEMBER_TI,
            GEOMETRY.with_adhesive_thickness(thickness_mm * 1.0e-3),
            ADHESIVE_LIKE,
            COLD,
        ).peak_shear_stress
        assert BASIS.margin_of_safety(ADHESIVE_LIKE, peak) < 0.0


def test_bg_area_ratio_sweep_is_deterministic_and_couples_to_milestone_1() -> None:
    """BG. Changing the area ratio changes the M1 demand, beta and the peak together."""
    forces, betas, peaks = [], [], []
    for ratio in AREA_RATIOS:
        member_1 = AxialMember(ALUMINIUM_LIKE, ratio * 100.0e-6)
        demand = solve_shear_lag(member_1, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, COLD)
        forces.append(demand.transferred_force)
        betas.append(demand.beta)
        peaks.append(demand.peak_shear_stress)

    # A stiffer aluminium-like member raises the mismatch force it must carry,
    # while lowering beta; the peak shear rises on balance.
    assert forces == sorted(forces)
    assert betas == sorted(betas, reverse=True)
    assert peaks == sorted(peaks)

    repeat = [
        solve_shear_lag(
            AxialMember(ALUMINIUM_LIKE, ratio * 100.0e-6),
            MEMBER_TI,
            GEOMETRY,
            ADHESIVE_LIKE,
            COLD,
        ).peak_shear_stress
        for ratio in AREA_RATIOS
    ]
    assert repeat == peaks


def test_bg_area_ratio_sweep_does_not_hold_the_demand_fixed() -> None:
    """BG. The sweep genuinely re-runs Milestone 1 rather than reusing one force."""
    from thermal_joint import solve_bimetallic_joint

    for ratio in AREA_RATIOS:
        member_1 = AxialMember(ALUMINIUM_LIKE, ratio * 100.0e-6)
        demand = solve_shear_lag(member_1, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, COLD)
        expected = abs(solve_bimetallic_joint(member_1, MEMBER_TI, COLD).member_1_force)
        assert demand.transferred_force == pytest.approx(expected, rel=1e-15)


@pytest.mark.parametrize("delta_t", [10.0, 50.0, 100.0, 140.0, 200.0])
def test_bh_larger_excursions_give_larger_peak_shear(delta_t: float) -> None:
    """BH. Peak shear grows with |dT| at fixed geometry."""
    smaller = solve_shear_lag(
        MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, delta_t
    ).peak_shear_stress
    larger = solve_shear_lag(
        MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, 1.5 * delta_t
    ).peak_shear_stress
    assert larger > smaller
    assert larger == pytest.approx(1.5 * smaller, rel=1e-13)


def test_bh_peak_shear_is_symmetric_in_the_sign_of_delta_temperature() -> None:
    """BH. Only the magnitude of the excursion sets the peak."""
    for delta_t in (50.0, 140.0):
        assert solve_shear_lag(
            MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, delta_t
        ).peak_shear_stress == pytest.approx(
            solve_shear_lag(
                MEMBER_AL, MEMBER_TI, GEOMETRY, ADHESIVE_LIKE, -delta_t
            ).peak_shear_stress,
            rel=1e-14,
        )


def test_bi_sweeps_do_not_mutate_the_canonical_inputs() -> None:
    """BI. After every sweep the shipped members, geometry and adhesive are unchanged."""
    canonical_geometry = radiator_joint_overlap()
    canonical_member = aluminium_like_member(100.0)
    adhesive_before = (
        ADHESIVE_LIKE.name,
        ADHESIVE_LIKE.shear_modulus,
        ADHESIVE_LIKE.shear_strength,
        ADHESIVE_LIKE.source_note,
    )

    for length_mm in OVERLAP_LENGTHS_MM:
        solve_shear_lag(
            MEMBER_AL,
            MEMBER_TI,
            canonical_geometry.with_overlap_length(length_mm * 1.0e-3),
            ADHESIVE_LIKE,
            COLD,
        )
    for thickness_mm in THICKNESSES_MM:
        solve_shear_lag(
            MEMBER_AL,
            MEMBER_TI,
            canonical_geometry.with_adhesive_thickness(thickness_mm * 1.0e-3),
            ADHESIVE_LIKE,
            COLD,
        )
    for ratio in AREA_RATIOS:
        solve_shear_lag(
            AxialMember(ALUMINIUM_LIKE, ratio * 100.0e-6),
            MEMBER_TI,
            canonical_geometry,
            ADHESIVE_LIKE,
            COLD,
        )

    assert canonical_geometry == radiator_joint_overlap()
    assert canonical_geometry.overlap_length == 0.040
    assert canonical_geometry.adhesive_thickness == 0.0002
    assert canonical_member == aluminium_like_member(100.0)
    assert canonical_member.area == pytest.approx(100.0e-6, rel=1e-12)
    assert ALUMINIUM_LIKE.thermal_expansion_coefficient == 23.0e-6
    assert adhesive_before == (
        ADHESIVE_LIKE.name,
        ADHESIVE_LIKE.shear_modulus,
        ADHESIVE_LIKE.shear_strength,
        ADHESIVE_LIKE.source_note,
    )

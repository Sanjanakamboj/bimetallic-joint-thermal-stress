"""Tests A-I: adhesive material, overlap geometry and shear design basis."""

from __future__ import annotations

import math

import pytest

from thermal_joint import AdhesiveMaterial, AdhesiveShearBasis, BondedOverlapGeometry

NON_FINITE = [math.nan, math.inf, -math.inf]
NON_POSITIVE = [0.0, -1.0, -1.0e9]


def test_a_valid_adhesive_stores_inputs() -> None:
    """A. A valid adhesive keeps its properties and provenance."""
    adhesive = AdhesiveMaterial(
        name="  Adhesive-like  ",
        shear_modulus=1.0e9,
        shear_strength=25.0e6,
        source_note="  ILLUSTRATIVE ADHESIVE INPUT - NOT DESIGN ALLOWABLE  ",
        notes="mid-range film adhesive",
    )
    assert adhesive.name == "Adhesive-like"
    assert adhesive.shear_modulus == 1.0e9
    assert adhesive.shear_strength == 25.0e6
    assert adhesive.source_note == "ILLUSTRATIVE ADHESIVE INPUT - NOT DESIGN ALLOWABLE"
    assert adhesive.notes == "mid-range film adhesive"


def test_a_notes_are_optional() -> None:
    """A. Notes may be omitted; provenance may not."""
    assert AdhesiveMaterial("A", 1.0e9, 25.0e6, "illustrative").notes is None


@pytest.mark.parametrize("bad", NON_POSITIVE + NON_FINITE)
def test_b_invalid_shear_modulus_rejected(bad: float) -> None:
    """B. G_a must be finite and strictly positive."""
    with pytest.raises(ValueError):
        AdhesiveMaterial("A", bad, 25.0e6, "illustrative")


@pytest.mark.parametrize("bad", NON_POSITIVE + NON_FINITE)
def test_c_invalid_shear_strength_rejected(bad: float) -> None:
    """C. Shear strength must be a finite positive magnitude."""
    with pytest.raises(ValueError):
        AdhesiveMaterial("A", 1.0e9, bad, "illustrative")


@pytest.mark.parametrize("bad_note", ["", "   ", "\t\n"])
def test_d_provenance_is_required(bad_note: str) -> None:
    """D. An empty provenance string is rejected outright."""
    with pytest.raises(ValueError):
        AdhesiveMaterial("A", 1.0e9, 25.0e6, bad_note)


def test_d_provenance_must_be_a_string() -> None:
    """D. Provenance must be text, not a placeholder object."""
    with pytest.raises(TypeError):
        AdhesiveMaterial("A", 1.0e9, 25.0e6, None)  # type: ignore[arg-type]


def test_d_illustrative_library_entry_is_labelled() -> None:
    """D. The shipped illustrative adhesive says so, in capitals."""
    from thermal_joint.illustrative import ADHESIVE_LIKE, ILLUSTRATIVE_ADHESIVE_NOTE

    assert "NOT DESIGN ALLOWABLE" in ILLUSTRATIVE_ADHESIVE_NOTE
    assert ADHESIVE_LIKE.source_note == ILLUSTRATIVE_ADHESIVE_NOTE
    assert "illustrative" in ADHESIVE_LIKE.name.lower()


def test_e_valid_overlap_geometry() -> None:
    """E. A valid overlap keeps its three dimensions."""
    geometry = BondedOverlapGeometry(0.040, 0.020, 0.0002)
    assert geometry.overlap_length == 0.040
    assert geometry.bond_width == 0.020
    assert geometry.adhesive_thickness == 0.0002


@pytest.mark.parametrize("bad", NON_POSITIVE + NON_FINITE)
def test_f_invalid_overlap_length_rejected(bad: float) -> None:
    """F. Overlap length must be finite and strictly positive."""
    with pytest.raises(ValueError):
        BondedOverlapGeometry(bad, 0.020, 0.0002)


@pytest.mark.parametrize("bad", NON_POSITIVE + NON_FINITE)
def test_g_invalid_bond_width_rejected(bad: float) -> None:
    """G. Bond width must be finite and strictly positive."""
    with pytest.raises(ValueError):
        BondedOverlapGeometry(0.040, bad, 0.0002)


@pytest.mark.parametrize("bad", NON_POSITIVE + NON_FINITE)
def test_h_invalid_adhesive_thickness_rejected(bad: float) -> None:
    """H. Adhesive thickness must be finite and strictly positive."""
    with pytest.raises(ValueError):
        BondedOverlapGeometry(0.040, 0.020, bad)


def test_i_bond_area_hand_calculation() -> None:
    """I. bond_area = b * L_b = 20 mm * 40 mm = 800 mm^2 = 8.0e-4 m^2."""
    assert BondedOverlapGeometry(0.040, 0.020, 0.0002).bond_area == pytest.approx(
        8.0e-4, rel=1e-15
    )


def test_i_bond_area_is_distinct_from_adherend_area() -> None:
    """I. The bonded interface area is not an adherend cross-section."""
    from cases import HAND_MEMBER_1

    geometry = BondedOverlapGeometry(0.040, 0.020, 0.0002)
    assert geometry.bond_area != HAND_MEMBER_1.area
    assert geometry.bond_area == pytest.approx(4.0 * HAND_MEMBER_1.area, rel=1e-12)


def test_geometry_copy_helpers_do_not_mutate_the_original() -> None:
    """Sweeps build copies; the canonical geometry is never edited in place."""
    geometry = BondedOverlapGeometry(0.040, 0.020, 0.0002)
    longer = geometry.with_overlap_length(0.160)
    thicker = geometry.with_adhesive_thickness(0.001)
    assert geometry.overlap_length == 0.040
    assert geometry.adhesive_thickness == 0.0002
    assert longer.overlap_length == 0.160
    assert longer.bond_width == geometry.bond_width
    assert thicker.adhesive_thickness == 0.001


def test_adhesive_and_geometry_are_immutable() -> None:
    """Frozen dataclasses keep validated state from being mutated."""
    adhesive = AdhesiveMaterial("A", 1.0e9, 25.0e6, "illustrative")
    geometry = BondedOverlapGeometry(0.040, 0.020, 0.0002)
    with pytest.raises(Exception):
        adhesive.shear_modulus = -1.0  # type: ignore[misc]
    with pytest.raises(Exception):
        geometry.overlap_length = -1.0  # type: ignore[misc]


def test_shear_basis_rejects_bad_design_factors() -> None:
    """The adhesive basis validates its design factor like the yield basis."""
    for bad in [0.0, 0.5, 0.999, -2.0] + NON_FINITE:
        with pytest.raises(ValueError):
            AdhesiveShearBasis(bad)
    assert AdhesiveShearBasis().design_factor == 1.0


def test_shear_basis_rejects_non_adhesive_arguments() -> None:
    with pytest.raises(TypeError):
        AdhesiveShearBasis().allowable_shear_stress("not-an-adhesive")  # type: ignore[arg-type]

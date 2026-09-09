"""Tests A-G: material and member construction and validation."""

from __future__ import annotations

import math

import pytest

from thermal_joint import AxialMember, ThermoelasticMaterial


def test_a_valid_material_stores_inputs() -> None:
    """A. A valid material keeps its inputs unchanged and strips its name."""
    material = ThermoelasticMaterial(
        name="  Alloy-like  ",
        elastic_modulus=70.0e9,
        thermal_expansion_coefficient=23.0e-6,
        yield_strength=270.0e6,
        density=2700.0,
    )
    assert material.name == "Alloy-like"
    assert material.elastic_modulus == 70.0e9
    assert material.thermal_expansion_coefficient == 23.0e-6
    assert material.yield_strength == 270.0e6
    assert material.density == 2700.0


def test_a_density_is_optional() -> None:
    """A. Density may be omitted entirely."""
    material = ThermoelasticMaterial("M", 70.0e9, 23.0e-6, 270.0e6)
    assert material.density is None


@pytest.mark.parametrize("bad_name", ["", "   "])
def test_a_empty_name_rejected(bad_name: str) -> None:
    """A. The material name must be non-empty."""
    with pytest.raises(ValueError):
        ThermoelasticMaterial(bad_name, 70.0e9, 23.0e-6, 270.0e6)


@pytest.mark.parametrize("bad_modulus", [0.0, -1.0, -70.0e9])
def test_b_nonpositive_modulus_rejected(bad_modulus: float) -> None:
    """B. Young's modulus must be strictly positive."""
    with pytest.raises(ValueError):
        ThermoelasticMaterial("M", bad_modulus, 23.0e-6, 270.0e6)


@pytest.mark.parametrize("bad_modulus", [math.nan, math.inf, -math.inf])
def test_b_non_finite_modulus_rejected(bad_modulus: float) -> None:
    """B. Young's modulus must be finite."""
    with pytest.raises(ValueError):
        ThermoelasticMaterial("M", bad_modulus, 23.0e-6, 270.0e6)


@pytest.mark.parametrize("bad_alpha", [math.nan, math.inf, -math.inf])
def test_c_non_finite_alpha_rejected(bad_alpha: float) -> None:
    """C. The thermal expansion coefficient must be finite."""
    with pytest.raises(ValueError):
        ThermoelasticMaterial("M", 70.0e9, bad_alpha, 270.0e6)


@pytest.mark.parametrize("alpha", [0.0, -1.0e-6, 23.0e-6])
def test_c_finite_alpha_of_any_sign_accepted(alpha: float) -> None:
    """C. Zero and negative CTE are physically meaningful and are accepted."""
    material = ThermoelasticMaterial("M", 70.0e9, alpha, 270.0e6)
    assert material.thermal_expansion_coefficient == alpha


@pytest.mark.parametrize("bad_yield", [0.0, -1.0, -270.0e6, math.nan, math.inf])
def test_d_invalid_yield_strength_rejected(bad_yield: float) -> None:
    """D. Yield strength must be a finite positive magnitude."""
    with pytest.raises(ValueError):
        ThermoelasticMaterial("M", 70.0e9, 23.0e-6, bad_yield)


def test_e_valid_member_defaults_label_to_material_name() -> None:
    """E. A valid member exposes its area, stiffness and a default label."""
    material = ThermoelasticMaterial("Alloy-like", 70.0e9, 23.0e-6, 270.0e6)
    member = AxialMember(material, 100.0e-6)
    assert member.area == 100.0e-6
    assert member.label == "Alloy-like"
    assert member.axial_stiffness == pytest.approx(7.0e6)


def test_e_member_accepts_explicit_label() -> None:
    """E. An explicit label overrides the material name."""
    material = ThermoelasticMaterial("Alloy-like", 70.0e9, 23.0e-6, 270.0e6)
    assert AxialMember(material, 100.0e-6, label="strap").label == "strap"


@pytest.mark.parametrize("bad_area", [0.0, -1.0e-6, math.nan, math.inf, -math.inf])
def test_f_invalid_area_rejected(bad_area: float) -> None:
    """F. Cross-sectional area must be finite and strictly positive."""
    material = ThermoelasticMaterial("M", 70.0e9, 23.0e-6, 270.0e6)
    with pytest.raises(ValueError):
        AxialMember(material, bad_area)


def test_f_member_requires_a_material_instance() -> None:
    """F. A member cannot be built from something that is not a material."""
    with pytest.raises(TypeError):
        AxialMember("not-a-material", 100.0e-6)  # type: ignore[arg-type]


@pytest.mark.parametrize("bad_density", [0.0, -1.0, math.nan, math.inf])
def test_g_invalid_density_rejected(bad_density: float) -> None:
    """G. When density is supplied it must be finite and strictly positive."""
    with pytest.raises(ValueError):
        ThermoelasticMaterial("M", 70.0e9, 23.0e-6, 270.0e6, density=bad_density)


def test_material_and_member_are_immutable() -> None:
    """Frozen dataclasses keep validated state from being mutated afterwards."""
    material = ThermoelasticMaterial("M", 70.0e9, 23.0e-6, 270.0e6)
    member = AxialMember(material, 100.0e-6)
    with pytest.raises(Exception):
        material.elastic_modulus = -1.0  # type: ignore[misc]
    with pytest.raises(Exception):
        member.area = -1.0  # type: ignore[misc]

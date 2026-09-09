"""Tests J-Q: mismatch strain, the shear-lag parameter beta and transfer length."""

from __future__ import annotations

import math

import pytest
from cases import HAND_MATERIAL_1, HAND_MATERIAL_2
from shear_lag_cases import (
    HAND_ADHESIVE,
    HAND_BETA,
    HAND_GEOMETRY,
    HAND_LAMBDA,
    HAND_MEMBER_1,
    HAND_MEMBER_2,
    HAND_MISMATCH_STRAIN,
    HAND_TRANSFER_LENGTH,
)

from thermal_joint import (
    AdhesiveMaterial,
    AxialMember,
    BondedOverlapGeometry,
    ThermoelasticMaterial,
    adherend_compliance_sum,
    dimensionless_overlap,
    shear_lag_parameter,
    thermal_mismatch_strain,
    transfer_length,
)


def test_j_mismatch_strain_hand_calculation() -> None:
    """J. (20e-6 - 10e-6) * 50 K = 5.0e-4."""
    assert thermal_mismatch_strain(HAND_MEMBER_1, HAND_MEMBER_2, 50.0) == pytest.approx(
        HAND_MISMATCH_STRAIN, rel=1e-14
    )


def test_j_mismatch_strain_ordering_is_member_1_minus_member_2() -> None:
    """J. The documented ordering is alpha_1 - alpha_2."""
    swapped = thermal_mismatch_strain(HAND_MEMBER_2, HAND_MEMBER_1, 50.0)
    assert swapped == pytest.approx(-HAND_MISMATCH_STRAIN, rel=1e-14)


@pytest.mark.parametrize("delta_t", [30.0, 100.0, 140.0])
def test_k_mismatch_strain_reverses_sign_with_delta_temperature(delta_t: float) -> None:
    """K. Heating and cooling give equal and opposite mismatch strain."""
    heating = thermal_mismatch_strain(HAND_MEMBER_1, HAND_MEMBER_2, delta_t)
    cooling = thermal_mismatch_strain(HAND_MEMBER_1, HAND_MEMBER_2, -delta_t)
    assert cooling == pytest.approx(-heating, rel=1e-15)
    assert heating * cooling < 0.0


@pytest.mark.parametrize("alpha", [0.0, 8.5e-6, 23.0e-6, -3.0e-6])
def test_l_mismatch_strain_vanishes_for_equal_cte(alpha: float) -> None:
    """L. Equal CTE means no mismatch at any temperature."""
    member_1 = AxialMember(ThermoelasticMaterial("A", 70.0e9, alpha, 270.0e6), 100.0e-6)
    member_2 = AxialMember(ThermoelasticMaterial("B", 110.0e9, alpha, 830.0e6), 250.0e-6)
    for delta_t in (-140.0, 0.0, 100.0):
        assert thermal_mismatch_strain(member_1, member_2, delta_t) == 0.0


def test_l_mismatch_strain_vanishes_at_zero_delta_temperature() -> None:
    """L. dT = 0 gives exactly zero mismatch."""
    assert thermal_mismatch_strain(HAND_MEMBER_1, HAND_MEMBER_2, 0.0) == 0.0


@pytest.mark.parametrize("factor", [0.5, 2.0, 3.0, -1.0])
def test_mismatch_strain_is_linear_in_delta_temperature_and_in_delta_alpha(
    factor: float,
) -> None:
    """Linear scaling in both dT and the CTE difference."""
    base = thermal_mismatch_strain(HAND_MEMBER_1, HAND_MEMBER_2, 50.0)
    assert thermal_mismatch_strain(
        HAND_MEMBER_1, HAND_MEMBER_2, 50.0 * factor
    ) == pytest.approx(factor * base, rel=1e-14)

    scaled_alpha = ThermoelasticMaterial(
        "scaled",
        HAND_MATERIAL_1.elastic_modulus,
        HAND_MATERIAL_2.thermal_expansion_coefficient
        + factor
        * (
            HAND_MATERIAL_1.thermal_expansion_coefficient
            - HAND_MATERIAL_2.thermal_expansion_coefficient
        ),
        HAND_MATERIAL_1.yield_strength,
    )
    scaled_member = AxialMember(scaled_alpha, HAND_MEMBER_1.area)
    assert thermal_mismatch_strain(
        scaled_member, HAND_MEMBER_2, 50.0
    ) == pytest.approx(factor * base, rel=1e-13)


def test_m_beta_hand_calculation() -> None:
    """M. beta^2 = (1e9 * 0.02 / 2e-4) * (1/2e7 + 1/2e7) = 1e4, so beta = 100 /m."""
    assert adherend_compliance_sum(HAND_MEMBER_1, HAND_MEMBER_2) == pytest.approx(
        1.0e-7, rel=1e-15
    )
    beta = shear_lag_parameter(HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, HAND_ADHESIVE)
    assert beta == pytest.approx(HAND_BETA, rel=1e-13)


def test_m_beta_matches_an_independent_expression() -> None:
    """M. Independent re-evaluation from raw scalars, over several geometries."""
    for length, width, thickness, modulus in (
        (0.04, 0.02, 2.0e-4, 1.0e9),
        (0.01, 0.05, 1.0e-3, 5.0e8),
        (0.16, 0.008, 5.0e-5, 2.0e9),
    ):
        geometry = BondedOverlapGeometry(length, width, thickness)
        adhesive = AdhesiveMaterial("A", modulus, 25.0e6, "illustrative")
        expected = math.sqrt(
            (modulus * width / thickness)
            * (
                1.0 / (HAND_MATERIAL_1.elastic_modulus * HAND_MEMBER_1.area)
                + 1.0 / (HAND_MATERIAL_2.elastic_modulus * HAND_MEMBER_2.area)
            )
        )
        assert shear_lag_parameter(
            HAND_MEMBER_1, HAND_MEMBER_2, geometry, adhesive
        ) == pytest.approx(expected, rel=1e-14)


def test_m_beta_does_not_depend_on_overlap_length() -> None:
    """M. beta is a material/geometry property of the bondline, not of L_b."""
    betas = {
        round(
            shear_lag_parameter(
                HAND_MEMBER_1,
                HAND_MEMBER_2,
                HAND_GEOMETRY.with_overlap_length(length),
                HAND_ADHESIVE,
            ),
            9,
        )
        for length in (0.005, 0.02, 0.08, 0.5)
    }
    assert len(betas) == 1


@pytest.mark.parametrize("factor", [0.25, 4.0, 9.0])
def test_n_beta_scales_with_sqrt_of_shear_modulus(factor: float) -> None:
    """N. beta ~ sqrt(G_a)."""
    base = shear_lag_parameter(HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, HAND_ADHESIVE)
    scaled = AdhesiveMaterial(
        "scaled", factor * HAND_ADHESIVE.shear_modulus, 25.0e6, "illustrative"
    )
    assert shear_lag_parameter(
        HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, scaled
    ) == pytest.approx(math.sqrt(factor) * base, rel=1e-13)


@pytest.mark.parametrize("factor", [0.25, 4.0, 9.0])
def test_o_beta_scales_with_inverse_sqrt_of_adhesive_thickness(factor: float) -> None:
    """O. beta ~ 1/sqrt(t_a): a thicker bondline spreads load further."""
    base = shear_lag_parameter(HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, HAND_ADHESIVE)
    thicker = HAND_GEOMETRY.with_adhesive_thickness(
        factor * HAND_GEOMETRY.adhesive_thickness
    )
    assert shear_lag_parameter(
        HAND_MEMBER_1, HAND_MEMBER_2, thicker, HAND_ADHESIVE
    ) == pytest.approx(base / math.sqrt(factor), rel=1e-13)


def test_p_beta_includes_both_adherend_compliances() -> None:
    """P. Both 1/(EA) terms enter; stiffening either adherend lowers beta."""
    stiffer_1 = AxialMember(HAND_MATERIAL_1, 10.0 * HAND_MEMBER_1.area)
    stiffer_2 = AxialMember(HAND_MATERIAL_2, 10.0 * HAND_MEMBER_2.area)
    base = shear_lag_parameter(HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, HAND_ADHESIVE)
    assert (
        shear_lag_parameter(stiffer_1, HAND_MEMBER_2, HAND_GEOMETRY, HAND_ADHESIVE) < base
    )
    assert (
        shear_lag_parameter(HAND_MEMBER_1, stiffer_2, HAND_GEOMETRY, HAND_ADHESIVE) < base
    )


def test_p_beta_is_symmetric_in_the_two_adherends() -> None:
    """P. The compliance sum is symmetric, so member order cannot matter."""
    assert shear_lag_parameter(
        HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, HAND_ADHESIVE
    ) == pytest.approx(
        shear_lag_parameter(HAND_MEMBER_2, HAND_MEMBER_1, HAND_GEOMETRY, HAND_ADHESIVE),
        rel=1e-15,
    )


def test_p_a_rigid_adherend_limit_reduces_to_the_other_compliance() -> None:
    """P. As one adherend becomes rigid, beta tends to the single-compliance value."""
    nearly_rigid = AxialMember(HAND_MATERIAL_2, 1.0e6 * HAND_MEMBER_2.area)
    beta = shear_lag_parameter(HAND_MEMBER_1, nearly_rigid, HAND_GEOMETRY, HAND_ADHESIVE)
    expected = math.sqrt(
        (HAND_ADHESIVE.shear_modulus * HAND_GEOMETRY.bond_width / HAND_GEOMETRY.adhesive_thickness)
        / HAND_MEMBER_1.axial_stiffness
    )
    assert beta == pytest.approx(expected, rel=1e-5)


def test_q_transfer_length_is_the_reciprocal_of_beta() -> None:
    """Q. transfer length = 1/beta = 10 mm for the hand case."""
    length = transfer_length(HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, HAND_ADHESIVE)
    assert length == pytest.approx(HAND_TRANSFER_LENGTH, rel=1e-13)
    assert length * shear_lag_parameter(
        HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, HAND_ADHESIVE
    ) == pytest.approx(1.0, rel=1e-15)


def test_q_dimensionless_overlap_is_beta_times_length() -> None:
    """Q. lambda = beta L_b = 100 * 0.020 = 2, the single metric used everywhere."""
    assert dimensionless_overlap(
        HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, HAND_ADHESIVE
    ) == pytest.approx(HAND_LAMBDA, rel=1e-13)


def test_parameter_helpers_type_check_their_arguments() -> None:
    with pytest.raises(TypeError):
        shear_lag_parameter(HAND_MEMBER_1, HAND_MEMBER_2, "nope", HAND_ADHESIVE)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        shear_lag_parameter(HAND_MEMBER_1, HAND_MEMBER_2, HAND_GEOMETRY, "nope")  # type: ignore[arg-type]


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_mismatch_strain_rejects_non_finite_delta_temperature(bad: float) -> None:
    with pytest.raises(ValueError):
        thermal_mismatch_strain(HAND_MEMBER_1, HAND_MEMBER_2, bad)

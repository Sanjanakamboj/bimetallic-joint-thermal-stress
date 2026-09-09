"""Tests A-F: the AxialRestraint object and the dimensionless restraint parameter."""

from __future__ import annotations

import math

import pytest
from cases import HAND_MEMBER_1, HAND_MEMBER_2

from thermal_joint import (
    AxialRestraint,
    joint_axial_stiffness,
    restraint_stiffness_ratio,
)


def test_a_valid_restraint_stores_inputs() -> None:
    """A. A valid restraint keeps its stiffness, reference strain and label."""
    restraint = AxialRestraint(stiffness=1.0e7, reference_strain=2.5e-4, label="  strut  ")
    assert restraint.stiffness == 1.0e7
    assert restraint.reference_strain == 2.5e-4
    assert restraint.label == "strut"


def test_a_defaults_are_the_canonical_case() -> None:
    """A. Reference strain defaults to zero with a descriptive default label."""
    restraint = AxialRestraint(1.0e7)
    assert restraint.reference_strain == 0.0
    assert restraint.label == "external axial restraint"


def test_b_zero_stiffness_is_valid() -> None:
    """B. K_r = 0 must be allowed; it is the free-joint limit."""
    restraint = AxialRestraint(0.0)
    assert restraint.stiffness == 0.0
    assert restraint.stiffness_ratio(HAND_MEMBER_1, HAND_MEMBER_2) == 0.0


@pytest.mark.parametrize("bad", [-1.0, -1.0e7, -1.0e-12])
def test_c_negative_stiffness_rejected(bad: float) -> None:
    """C. A negative restraint stiffness is unphysical."""
    with pytest.raises(ValueError):
        AxialRestraint(bad)


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_d_non_finite_stiffness_rejected(bad: float) -> None:
    """D. The stiffness must be finite; rigidity is a limit, not an input."""
    with pytest.raises(ValueError):
        AxialRestraint(bad)


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_e_non_finite_reference_strain_rejected(bad: float) -> None:
    """E. The reference strain must be finite."""
    with pytest.raises(ValueError):
        AxialRestraint(1.0e7, reference_strain=bad)


@pytest.mark.parametrize("reference_strain", [-1.0e-3, 0.0, 5.0e-4])
def test_e_finite_reference_strain_of_any_sign_accepted(reference_strain: float) -> None:
    """E. Positive, zero and negative reference strains are all meaningful."""
    assert AxialRestraint(1.0e7, reference_strain).reference_strain == reference_strain


@pytest.mark.parametrize("bad_label", ["", "   ", "\t"])
def test_empty_label_rejected(bad_label: str) -> None:
    """A supplied label must be non-empty."""
    with pytest.raises(ValueError):
        AxialRestraint(1.0e7, label=bad_label)


def test_f_stiffness_ratio_hand_calculation() -> None:
    """F. eta_r = K_r / (E1 A1 + E2 A2) = 1.0e7 / 4.0e7 = 0.25."""
    assert joint_axial_stiffness(HAND_MEMBER_1, HAND_MEMBER_2) == pytest.approx(4.0e7, rel=1e-15)
    restraint = AxialRestraint(1.0e7)
    assert restraint.stiffness_ratio(HAND_MEMBER_1, HAND_MEMBER_2) == pytest.approx(0.25, rel=1e-15)
    assert restraint_stiffness_ratio(restraint, HAND_MEMBER_1, HAND_MEMBER_2) == pytest.approx(
        0.25, rel=1e-15
    )


def test_f_joint_axial_stiffness_is_the_sum_of_member_stiffnesses() -> None:
    """F. Validate the normalisation formula independently of the members' API."""
    expected = (
        HAND_MEMBER_1.material.elastic_modulus * HAND_MEMBER_1.area
        + HAND_MEMBER_2.material.elastic_modulus * HAND_MEMBER_2.area
    )
    assert joint_axial_stiffness(HAND_MEMBER_1, HAND_MEMBER_2) == pytest.approx(expected, rel=1e-15)


@pytest.mark.parametrize("ratio", [0.0, 0.1, 0.25, 1.0, 10.0, 100.0])
def test_f_from_stiffness_ratio_round_trips(ratio: float) -> None:
    """F. Building from eta_r and reading eta_r back is an identity."""
    restraint = AxialRestraint.from_stiffness_ratio(ratio, HAND_MEMBER_1, HAND_MEMBER_2)
    assert restraint.stiffness == pytest.approx(ratio * 4.0e7, rel=1e-14)
    assert restraint.stiffness_ratio(HAND_MEMBER_1, HAND_MEMBER_2) == pytest.approx(
        ratio, rel=1e-14, abs=1e-18
    )


def test_from_stiffness_ratio_rejects_negative_ratio() -> None:
    """A negative eta_r is as unphysical as a negative stiffness."""
    with pytest.raises(ValueError):
        AxialRestraint.from_stiffness_ratio(-0.5, HAND_MEMBER_1, HAND_MEMBER_2)


def test_restraint_is_immutable() -> None:
    """Frozen dataclass keeps validated state from being mutated."""
    restraint = AxialRestraint(1.0e7)
    with pytest.raises(Exception):
        restraint.stiffness = -1.0  # type: ignore[misc]

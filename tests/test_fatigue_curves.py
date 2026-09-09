"""Basquin curve, cycle requirement, life inversion and the mean-stress policy."""

from __future__ import annotations

import math

import pytest
from fatigue_cases import (
    ADHESIVE_FATIGUE,
    AL_FATIGUE,
    CURVES,
    REQUIREMENT,
    TI_FATIGUE,
)

from thermal_joint import (
    BasquinFatigueCurve,
    ThermalCycleRequirement,
    assess_fatigue_life,
    stress_cycle,
)

NON_FINITE = [math.nan, math.inf, -math.inf]


def test_valid_curve_stores_inputs_and_provenance() -> None:
    curve = BasquinFatigueCurve(
        "  Test curve  ", 900.0e6, -0.12, "  ILLUSTRATIVE FATIGUE INPUT  ", "note"
    )
    assert curve.name == "Test curve"
    assert curve.coefficient_A == 900.0e6
    assert curve.exponent_b == -0.12
    assert curve.source_note == "ILLUSTRATIVE FATIGUE INPUT"
    assert curve.notes == "note"


@pytest.mark.parametrize("bad", [0.0, -1.0, -900.0e6] + NON_FINITE)
def test_invalid_coefficient_rejected(bad: float) -> None:
    with pytest.raises(ValueError):
        BasquinFatigueCurve("c", bad, -0.12, "illustrative")


@pytest.mark.parametrize("bad", [0.0, 0.12, 1.0] + NON_FINITE)
def test_non_negative_exponent_rejected(bad: float) -> None:
    """b must be strictly negative so that life falls as amplitude rises."""
    with pytest.raises(ValueError):
        BasquinFatigueCurve("c", 900.0e6, bad, "illustrative")


@pytest.mark.parametrize("bad_note", ["", "   ", "\t"])
def test_provenance_is_required(bad_note: str) -> None:
    with pytest.raises(ValueError):
        BasquinFatigueCurve("c", 900.0e6, -0.12, bad_note)


def test_shipped_curves_are_labelled_illustrative() -> None:
    """Every shipped curve says NOT DESIGN ALLOWABLE, in capitals."""
    from thermal_joint.illustrative import ILLUSTRATIVE_FATIGUE_NOTE

    assert "NOT DESIGN ALLOWABLE" in ILLUSTRATIVE_FATIGUE_NOTE
    for curve in (AL_FATIGUE, TI_FATIGUE, ADHESIVE_FATIGUE):
        assert curve.source_note == ILLUSTRATIVE_FATIGUE_NOTE
        assert "illustrative" in curve.name.lower()
        assert "endurance limit" in curve.notes


def test_adhesive_curve_is_declared_in_shear() -> None:
    """No shear-to-von-Mises conversion is applied anywhere in Milestone 5."""
    assert "shear" in ADHESIVE_FATIGUE.name.lower()
    assert "von-Mises" in ADHESIVE_FATIGUE.notes or "von Mises" in ADHESIVE_FATIGUE.notes
    assert ADHESIVE_FATIGUE is not AL_FATIGUE
    assert CURVES.adhesive_shear is ADHESIVE_FATIGUE


def test_curve_hand_calculations() -> None:
    """A N^b at round decade lives, for each shipped curve."""
    assert AL_FATIGUE.alternating_stress_at_life(1.0e4) == pytest.approx(298.018e6, rel=1e-5)
    assert AL_FATIGUE.alternating_stress_at_life(1.0e7) == pytest.approx(130.09e6, rel=1e-4)
    assert TI_FATIGUE.alternating_stress_at_life(1.0e4) == pytest.approx(796.214e6, rel=1e-5)
    assert ADHESIVE_FATIGUE.alternating_stress_at_life(1.0e4) == pytest.approx(15.0713e6, rel=1e-5)
    assert ADHESIVE_FATIGUE.alternating_stress_at_life(1.0e6) == pytest.approx(7.55355e6, rel=1e-5)
    assert AL_FATIGUE.alternating_stress_at_life(1.0) == pytest.approx(900.0e6, rel=1e-14)


def test_life_inversion_hand_calculation() -> None:
    """N_f = (sigma_a / A)^(1/b)."""
    assert ADHESIVE_FATIGUE.life_at_alternating_stress(15.0713e6) == pytest.approx(
        1.0e4, rel=1e-4
    )
    assert AL_FATIGUE.life_at_alternating_stress(74.4333e6) == pytest.approx(1.0486e9, rel=1e-3)


@pytest.mark.parametrize("life", [1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e8])
@pytest.mark.parametrize("curve_name", ["al", "ti", "adhesive"])
def test_round_trip_between_stress_and_life(life: float, curve_name: str) -> None:
    """stress_at_life(life_at_stress(sigma_a)) == sigma_a, and the reverse."""
    curve = {"al": AL_FATIGUE, "ti": TI_FATIGUE, "adhesive": ADHESIVE_FATIGUE}[curve_name]
    stress = curve.alternating_stress_at_life(life)
    assert curve.life_at_alternating_stress(stress) == pytest.approx(life, rel=1e-9)
    assert curve.alternating_stress_at_life(
        curve.life_at_alternating_stress(stress)
    ) == pytest.approx(stress, rel=1e-9)


def test_higher_amplitude_always_lowers_life() -> None:
    for curve in (AL_FATIGUE, TI_FATIGUE, ADHESIVE_FATIGUE):
        lives = [
            curve.life_at_alternating_stress(stress)
            for stress in (5.0e6, 10.0e6, 20.0e6, 50.0e6, 100.0e6)
        ]
        assert lives == sorted(lives, reverse=True)


def test_life_is_strongly_nonlinear_in_amplitude() -> None:
    """A 2x amplitude costs 2^(-1/b) in life - a factor of ~330 for the adhesive."""
    base = ADHESIVE_FATIGUE.life_at_alternating_stress(10.0e6)
    doubled = ADHESIVE_FATIGUE.life_at_alternating_stress(20.0e6)
    assert base / doubled == pytest.approx(2.0 ** (-1.0 / ADHESIVE_FATIGUE.exponent_b), rel=1e-12)
    assert base / doubled == pytest.approx(101.594, rel=1e-4)


def test_zero_amplitude_gives_infinite_life() -> None:
    """A component that never cycles is explicitly non-governing."""
    for curve in (AL_FATIGUE, TI_FATIGUE, ADHESIVE_FATIGUE):
        assert math.isinf(curve.life_at_alternating_stress(0.0))


def test_life_uses_the_magnitude_of_the_amplitude() -> None:
    """A sign on the amplitude cannot change the answer."""
    assert AL_FATIGUE.life_at_alternating_stress(-50.0e6) == (
        AL_FATIGUE.life_at_alternating_stress(50.0e6)
    )


@pytest.mark.parametrize("bad", [0.0, -1.0] + NON_FINITE)
def test_stress_at_life_validates_its_argument(bad: float) -> None:
    with pytest.raises(ValueError):
        AL_FATIGUE.alternating_stress_at_life(bad)


# ------------------------------------------------------------ requirement

def test_requirement_validation() -> None:
    requirement = ThermalCycleRequirement(1.0e4, "illustrative")
    assert requirement.required_cycles == 1.0e4
    assert requirement.label == "illustrative"
    for bad in [0.0, -1.0] + NON_FINITE:
        with pytest.raises(ValueError):
            ThermalCycleRequirement(bad)
    with pytest.raises(ValueError):
        ThermalCycleRequirement(1.0e4, "   ")


def test_canonical_requirement_is_labelled_illustrative() -> None:
    assert REQUIREMENT.required_cycles == 1.0e4
    assert "not a qualification requirement" in REQUIREMENT.label


# --------------------------------------------------- the mean-stress policy

def test_mean_stress_does_not_affect_life() -> None:
    """Milestone 5 applies NO mean-stress correction; life uses amplitude only.

    Three cycles with identical amplitude but very different mean stress must
    give identical predicted life.
    """
    amplitude = 40.0e6
    cycles = [
        stress_cycle(-amplitude, amplitude),               # mean 0
        stress_cycle(0.0, 2.0 * amplitude),                # mean +40 MPa
        stress_cycle(-2.0 * amplitude, 0.0),               # mean -40 MPa
        stress_cycle(100.0e6 - amplitude, 100.0e6 + amplitude),  # mean +100 MPa
    ]
    means = {cycle.mean_stress for cycle in cycles}
    assert len(means) == 4  # the means really do differ

    lives = {
        assess_fatigue_life("c", cycle, AL_FATIGUE, REQUIREMENT).predicted_cycles_to_failure
        for cycle in cycles
    }
    assert len(lives) == 1


def test_mean_stress_is_still_reported() -> None:
    """The omission is visible: mean stress is carried on every result."""
    cycle = stress_cycle(-20.0e6, 100.0e6)
    result = assess_fatigue_life("c", cycle, AL_FATIGUE, REQUIREMENT)
    assert result.mean_stress == pytest.approx(40.0e6, rel=1e-12)
    assert result.mean_stress == cycle.mean_stress
    assert result.alternating_stress == pytest.approx(60.0e6, rel=1e-12)


# ---------------------------------------------------------- life results

def test_life_result_ratio_and_margin() -> None:
    """life_ratio = N_f / N_required and MS = ratio - 1."""
    cycle = stress_cycle(-15.0713e6, 15.0713e6)
    result = assess_fatigue_life("adhesive shear", cycle, ADHESIVE_FATIGUE, REQUIREMENT)
    assert result.predicted_cycles_to_failure == pytest.approx(1.0e4, rel=1e-4)
    assert result.life_ratio == pytest.approx(1.0, rel=1e-4)
    assert result.margin == pytest.approx(0.0, abs=1e-4)
    assert result.passes


def test_exact_boundary_passes() -> None:
    """N_f == N_required is a PASS."""
    stress = ADHESIVE_FATIGUE.alternating_stress_at_life(1.0e4)
    cycle = stress_cycle(-stress, stress)
    result = assess_fatigue_life("adhesive shear", cycle, ADHESIVE_FATIGUE, REQUIREMENT)
    assert result.predicted_cycles_to_failure == pytest.approx(1.0e4, rel=1e-12)
    assert result.passes
    assert result.margin == pytest.approx(0.0, abs=1e-9)


def test_failing_case_reports_a_negative_margin() -> None:
    cycle = stress_cycle(-30.0e6, 30.0e6)
    result = assess_fatigue_life("adhesive shear", cycle, ADHESIVE_FATIGUE, REQUIREMENT)
    assert result.predicted_cycles_to_failure < 1.0e4
    assert result.life_ratio < 1.0
    assert result.margin < 0.0
    assert not result.passes


def test_zero_amplitude_result_is_infinite_and_passing() -> None:
    result = assess_fatigue_life("c", stress_cycle(5.0e6, 5.0e6), AL_FATIGUE, REQUIREMENT)
    assert math.isinf(result.predicted_cycles_to_failure)
    assert math.isinf(result.life_ratio)
    assert result.passes


def test_life_result_type_checking() -> None:
    cycle = stress_cycle(-1.0e6, 1.0e6)
    with pytest.raises(TypeError):
        assess_fatigue_life("c", "not-a-cycle", AL_FATIGUE, REQUIREMENT)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        assess_fatigue_life("c", cycle, "not-a-curve", REQUIREMENT)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        assess_fatigue_life("c", cycle, AL_FATIGUE, 1.0e4)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        assess_fatigue_life("  ", cycle, AL_FATIGUE, REQUIREMENT)

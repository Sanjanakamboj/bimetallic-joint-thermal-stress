r"""Milestone 5: first-order thermal-cycle fatigue screening primitives.

Milestones 1-4 give static states: global thermal mismatch stress, restrained
stress, adhesive shear transfer and static feasibility. Milestone 5 consumes
the already-verified hot and cold endpoints and asks whether repeated cycling
between them is life-limiting.

Cycle definition
----------------
One complete excursion **cold -> hot -> cold** counts as **one thermal cycle**.
Only the two endpoint states are used: no transient path, no dwell time, no
rate effect, and no variable-amplitude spectrum. ``N`` always means *cycles*,
never reversals, everywhere in this package.

For a signed scalar state with endpoints ``sigma_hot`` and ``sigma_cold``::

    sigma_max = max(sigma_hot, sigma_cold)
    sigma_min = min(sigma_hot, sigma_cold)
    sigma_a   = (sigma_max - sigma_min) / 2       alternating
    sigma_m   = (sigma_max + sigma_min) / 2       mean
    R         = sigma_min / sigma_max             (only when sigma_max != 0)

The signs are kept all the way through: for the adhesive the linear model
reverses shear with ``dT``, so the two endpoints have opposite signs and the
cycle crosses zero. Feeding in two magnitudes instead of two signed values
would collapse that cycle and badly under-predict the amplitude.

Fatigue curve
-------------
A Basquin power law in cycles::

    sigma_a = A * N^b          =>       N_f = (sigma_a / A)^(1/b)

with ``A > 0`` and ``b < 0``.

Mean-stress policy
------------------
**Mean stress is reported but no mean-stress correction is applied in
Milestone 5.** Life depends on the alternating stress alone. Goodman, Gerber
and Soderberg corrections are deliberately excluded: with unsourced,
illustrative curves an unsourced correction would only add false authority.
The mean stress is computed and carried so the omission is visible.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

from ._validation import require_finite, require_non_empty_name, require_positive

__all__ = [
    "StressCycle",
    "stress_cycle",
    "BasquinFatigueCurve",
    "ThermalCycleRequirement",
    "FatigueLifeResult",
    "assess_fatigue_life",
]


@dataclass(frozen=True)
class StressCycle:
    """A constant-amplitude cycle between two signed endpoint states.

    Units follow the inputs: [Pa] for both normal stress and shear stress.

    Attributes
    ----------
    hot_value, cold_value:
        The signed endpoint states at the hot and cold extremes.
    maximum_stress, minimum_stress:
        Algebraic max and min of the two endpoints - not magnitudes.
    alternating_stress:
        ``(max - min) / 2``. Always non-negative.
    mean_stress:
        ``(max + min) / 2``. Signed; reported but not used for life.
    stress_ratio:
        ``R = min / max``, or ``None`` when ``max`` is exactly zero.
    """

    hot_value: float
    cold_value: float
    maximum_stress: float
    minimum_stress: float
    alternating_stress: float
    mean_stress: float
    stress_ratio: Optional[float]

    @property
    def crosses_zero(self) -> bool:
        """True when the two endpoints straddle zero (a fully reversing cycle)."""
        return self.maximum_stress > 0.0 > self.minimum_stress


def stress_cycle(hot_value: float, cold_value: float) -> StressCycle:
    """Build a :class:`StressCycle` from two **signed** endpoint states.

    The same helper serves metal normal stress and adhesive shear stress; the
    algebra is identical and keeping one implementation keeps the two
    consistent.

    Parameters
    ----------
    hot_value, cold_value:
        Signed states at the hot and cold extremes [Pa]. Pass signed values,
        never magnitudes - see the module docstring.
    """
    hot = require_finite(hot_value, "hot_value")
    cold = require_finite(cold_value, "cold_value")
    maximum = max(hot, cold)
    minimum = min(hot, cold)
    return StressCycle(
        hot_value=hot,
        cold_value=cold,
        maximum_stress=maximum,
        minimum_stress=minimum,
        alternating_stress=0.5 * (maximum - minimum),
        mean_stress=0.5 * (maximum + minimum),
        stress_ratio=None if maximum == 0.0 else minimum / maximum,
    )


@dataclass(frozen=True)
class BasquinFatigueCurve:
    """An illustrative Basquin S-N curve ``sigma_a = A * N^b``, ``N`` in cycles.

    Attributes
    ----------
    name:
        Non-empty identifier.
    coefficient_A:
        Basquin coefficient ``A`` [Pa]; the alternating stress the curve gives
        at ``N = 1`` cycle. Finite and strictly positive.
    exponent_b:
        Basquin exponent ``b``; finite and strictly **negative**, so that life
        falls as amplitude rises.
    source_note:
        Required non-empty provenance. Unsourced curves must say
        ``ILLUSTRATIVE FATIGUE INPUT - NOT DESIGN ALLOWABLE``.
    notes:
        Optional remarks, e.g. the range of ``N`` over which the fit is meant
        to be read.

    Notes
    -----
    A Basquin power law has no endurance limit and no low-cycle cut-off, so
    extrapolating it to very small or very large ``N`` is not meaningful. The
    curves shipped in :mod:`thermal_joint.illustrative` state an intended range
    in their notes; nothing is clamped, because silently clamping would hide
    the extrapolation rather than flag it.
    """

    name: str
    coefficient_A: float
    exponent_b: float
    source_note: str
    notes: Optional[str] = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", require_non_empty_name(self.name))
        object.__setattr__(
            self, "coefficient_A", require_positive(self.coefficient_A, "coefficient_A")
        )
        exponent = require_finite(self.exponent_b, "exponent_b")
        if exponent >= 0.0:
            raise ValueError(f"exponent_b must be strictly negative, got {exponent!r}")
        object.__setattr__(self, "exponent_b", exponent)
        object.__setattr__(
            self, "source_note", require_non_empty_name(self.source_note, "source_note")
        )
        if self.notes is not None:
            object.__setattr__(self, "notes", require_non_empty_name(self.notes, "notes"))

    def alternating_stress_at_life(self, cycles: float) -> float:
        """Alternating stress the curve allows at ``N`` cycles [Pa].

        ``sigma_a = A * N^b``.
        """
        life = require_positive(cycles, "cycles")
        return self.coefficient_A * life**self.exponent_b

    def life_at_alternating_stress(self, alternating_stress: float) -> float:
        """Predicted cycles to failure at an alternating stress [-].

        ``N_f = (sigma_a / A)^(1/b)``. A zero amplitude returns ``math.inf``:
        an explicitly non-governing, infinite life rather than a division by
        zero. The magnitude of the input is used, so the caller cannot get a
        different answer by flipping the sign.
        """
        amplitude = abs(require_finite(alternating_stress, "alternating_stress"))
        if amplitude == 0.0:
            return math.inf
        return (amplitude / self.coefficient_A) ** (1.0 / self.exponent_b)


@dataclass(frozen=True)
class ThermalCycleRequirement:
    """The number of thermal cycles the joint must survive.

    Attributes
    ----------
    required_cycles:
        Finite, strictly positive cycle count [-].
    label:
        Non-empty description of where the number comes from.
    """

    required_cycles: float
    label: str = "illustrative thermal-cycle requirement"

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "required_cycles", require_positive(self.required_cycles, "required_cycles")
        )
        object.__setattr__(self, "label", require_non_empty_name(self.label, "label"))


@dataclass(frozen=True)
class FatigueLifeResult:
    """A **preliminary fatigue life margin** for one component. Not certification life.

    Attributes
    ----------
    component:
        Label of the thing being assessed.
    cycle:
        The :class:`StressCycle` used.
    alternating_stress, mean_stress:
        Copied out of the cycle for convenience [Pa]. Mean stress is reported
        but plays no part in the life calculation.
    fatigue_curve:
        The :class:`BasquinFatigueCurve` applied.
    predicted_cycles_to_failure:
        ``N_f`` [-]; ``math.inf`` for a zero-amplitude cycle.
    required_cycles:
        ``N_required`` [-].
    life_ratio:
        ``N_f / N_required`` [-].
    margin:
        ``life_ratio - 1`` [-]. A margin of exactly zero passes.
    """

    component: str
    cycle: StressCycle
    alternating_stress: float
    mean_stress: float
    fatigue_curve: BasquinFatigueCurve
    predicted_cycles_to_failure: float
    required_cycles: float
    life_ratio: float
    margin: float

    @property
    def passes(self) -> bool:
        """True when ``N_f >= N_required``; the boundary passes."""
        return self.predicted_cycles_to_failure >= self.required_cycles


def assess_fatigue_life(
    component: str,
    cycle: StressCycle,
    fatigue_curve: BasquinFatigueCurve,
    requirement: ThermalCycleRequirement,
) -> FatigueLifeResult:
    """Score one stress cycle against one curve and one cycle requirement.

    Life uses the alternating stress only; the mean stress is carried through
    unused, per the Milestone 5 mean-stress policy.
    """
    if not isinstance(cycle, StressCycle):
        raise TypeError(f"cycle must be a StressCycle, got {type(cycle).__name__}")
    if not isinstance(fatigue_curve, BasquinFatigueCurve):
        raise TypeError(
            f"fatigue_curve must be a BasquinFatigueCurve, got {type(fatigue_curve).__name__}"
        )
    if not isinstance(requirement, ThermalCycleRequirement):
        raise TypeError(
            f"requirement must be a ThermalCycleRequirement, got {type(requirement).__name__}"
        )

    predicted = fatigue_curve.life_at_alternating_stress(cycle.alternating_stress)
    ratio = predicted / requirement.required_cycles
    return FatigueLifeResult(
        component=require_non_empty_name(component, "component"),
        cycle=cycle,
        alternating_stress=cycle.alternating_stress,
        mean_stress=cycle.mean_stress,
        fatigue_curve=fatigue_curve,
        predicted_cycles_to_failure=predicted,
        required_cycles=requirement.required_cycles,
        life_ratio=ratio,
        margin=ratio - 1.0,
    )

r"""Milestone 5: integrated thermal-cycle fatigue assessment and sensitivity studies.

Consumes the verified Milestone 1 and Milestone 3 hot/cold states and scores
them against illustrative Basquin curves. Static and fatigue margins are kept
**separately visible and are never numerically blended**; only the pass/fail
booleans are combined.

Scoping note on external restraint
----------------------------------
The Milestone 3 shear-lag model takes its demand from the *free* joint, and
Milestone 3 deliberately declined to push Milestone 2 restrained member forces
through the same overlap model, because the external load path is different and
is not defined there. Milestone 5 keeps that boundary: the integrated
assessment is a free-joint assessment, and external restraint is studied for
**metal fatigue only**, through :func:`restraint_fatigue_sensitivity`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Sequence

from ._validation import require_finite, require_positive
from .adhesive import AdhesiveMaterial, AdhesiveShearBasis, BondedOverlapGeometry
from .design_trade import _modulus_variant, evaluate_joint_design
from .environment import ThermalEnvironment
from .extremes import assess_temperature_extremes
from .fatigue import (
    BasquinFatigueCurve,
    FatigueLifeResult,
    StressCycle,
    ThermalCycleRequirement,
    assess_fatigue_life,
    stress_cycle,
)
from .margins import YieldBasis
from .materials import AxialMember, ThermoelasticMaterial
from .restraint import AxialRestraint, solve_restrained_joint
from .shear_lag import solve_shear_lag
from .shear_lag_extremes import assess_shear_lag_extremes

__all__ = [
    "FatigueCurveSet",
    "ThermalCycleFatigueAssessment",
    "RestraintFatiguePoint",
    "FatigueSweepPoint",
    "AllowableCycleScaleResult",
    "member_stress_cycles",
    "adhesive_shear_cycle",
    "assess_thermal_cycle_fatigue",
    "restraint_fatigue_sensitivity",
    "allowable_cycle_scale_for_fatigue",
    "temperature_scale_fatigue_sensitivity",
    "bond_width_fatigue_sensitivity",
    "adhesive_thickness_fatigue_sensitivity",
    "adhesive_modulus_fatigue_sensitivity",
    "area_ratio_fatigue_sensitivity",
    "fatigue_design_map",
]


@dataclass(frozen=True)
class FatigueCurveSet:
    """The three illustrative curves used by an assessment.

    Separate records: two metal S-N curves in normal stress and one adhesive
    curve in **shear**. No shear-to-von-Mises conversion is performed - that is
    out of scope for Milestone 5.
    """

    member_1: BasquinFatigueCurve
    member_2: BasquinFatigueCurve
    adhesive_shear: BasquinFatigueCurve

    def __post_init__(self) -> None:
        for field, value in (
            ("member_1", self.member_1),
            ("member_2", self.member_2),
            ("adhesive_shear", self.adhesive_shear),
        ):
            if not isinstance(value, BasquinFatigueCurve):
                raise TypeError(
                    f"{field} must be a BasquinFatigueCurve, got {type(value).__name__}"
                )


def member_stress_cycles(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    restraint: AxialRestraint | None = None,
) -> tuple[StressCycle, StressCycle]:
    """Hot/cold member stress cycles, free by default or under a given restraint.

    Uses the verified Milestone 1 / Milestone 2 solvers unchanged.
    """
    external = AxialRestraint(0.0) if restraint is None else restraint
    hot = solve_restrained_joint(
        member_1, member_2, environment.delta_temperature_hot, external
    )
    cold = solve_restrained_joint(
        member_1, member_2, environment.delta_temperature_cold, external
    )
    return (
        stress_cycle(hot.member_1_stress, cold.member_1_stress),
        stress_cycle(hot.member_2_stress, cold.member_2_stress),
    )


def adhesive_shear_cycle(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
) -> StressCycle:
    """Signed adhesive shear cycle at the peak location.

    Both endpoints are read at the **same physical station** - the free edge
    ``x = L_b``, where Milestone 3 proved the peak sits - so the signed
    reversal is preserved. The linear model flips the sign with ``dT``, so for
    the canonical asymmetric excursion the two endpoints have opposite signs
    and unequal magnitudes, giving a zero-crossing cycle with a nonzero mean.
    """
    station = geometry.overlap_length
    hot = solve_shear_lag(
        member_1, member_2, geometry, adhesive, environment.delta_temperature_hot
    ).shear_stress(station)
    cold = solve_shear_lag(
        member_1, member_2, geometry, adhesive, environment.delta_temperature_cold
    ).shear_stress(station)
    return stress_cycle(hot, cold)


@dataclass(frozen=True)
class ThermalCycleFatigueAssessment:
    """Static and fatigue screens for one configuration, side by side.

    Attributes
    ----------
    member_1_fatigue, member_2_fatigue, adhesive_fatigue:
        Per-component :class:`~thermal_joint.fatigue.FatigueLifeResult`.
    static_yield_feasible, static_adhesive_feasible:
        The Milestone 1 and Milestone 3 static screens, unchanged.
    minimum_yield_margin, minimum_adhesive_margin:
        Static margins [-], reported separately and never blended with life.
    governing_fatigue_component:
        ``"member 1"``, ``"member 2"`` or ``"adhesive shear"`` - whichever has
        the smallest life ratio. Computed, never assumed. An exact tie resolves
        in that stated order.
    minimum_life_ratio:
        Smallest ``N_f / N_required`` over the three components [-].
    """

    requirement: ThermalCycleRequirement
    member_1_fatigue: FatigueLifeResult
    member_2_fatigue: FatigueLifeResult
    adhesive_fatigue: FatigueLifeResult
    static_yield_feasible: bool
    static_adhesive_feasible: bool
    minimum_yield_margin: float
    minimum_adhesive_margin: float
    governing_fatigue_component: str
    minimum_life_ratio: float

    @property
    def fatigue_results(self) -> tuple[FatigueLifeResult, ...]:
        """The three per-component results, in a stable order."""
        return (self.member_1_fatigue, self.member_2_fatigue, self.adhesive_fatigue)

    @property
    def governing_fatigue_result(self) -> FatigueLifeResult:
        """The result carrying the smallest life ratio."""
        index = {"member 1": 0, "member 2": 1, "adhesive shear": 2}[
            self.governing_fatigue_component
        ]
        return self.fatigue_results[index]

    @property
    def fatigue_feasible(self) -> bool:
        """True when every fatigue component meets the cycle requirement."""
        return all(result.passes for result in self.fatigue_results)

    @property
    def overall_feasible(self) -> bool:
        """Static yield AND static adhesive shear AND all fatigue checks.

        A boolean conjunction of three independent screens - not a blended
        margin.
        """
        return (
            self.static_yield_feasible
            and self.static_adhesive_feasible
            and self.fatigue_feasible
        )


def assess_thermal_cycle_fatigue(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    curves: FatigueCurveSet,
    requirement: ThermalCycleRequirement,
    yield_basis: YieldBasis | None = None,
    adhesive_basis: AdhesiveShearBasis | None = None,
) -> ThermalCycleFatigueAssessment:
    """Score the free-joint configuration on static and fatigue screens together."""
    if not isinstance(curves, FatigueCurveSet):
        raise TypeError(f"curves must be a FatigueCurveSet, got {type(curves).__name__}")

    cycle_1, cycle_2 = member_stress_cycles(member_1, member_2, environment)
    adhesive_cycle = adhesive_shear_cycle(
        member_1, member_2, environment, geometry, adhesive
    )

    result_1 = assess_fatigue_life("member 1", cycle_1, curves.member_1, requirement)
    result_2 = assess_fatigue_life("member 2", cycle_2, curves.member_2, requirement)
    result_a = assess_fatigue_life(
        "adhesive shear", adhesive_cycle, curves.adhesive_shear, requirement
    )

    yield_assessment = assess_temperature_extremes(
        member_1, member_2, environment, yield_basis
    )
    shear_assessment = assess_shear_lag_extremes(
        member_1, member_2, environment, geometry, adhesive, adhesive_basis
    )

    ordered = (("member 1", result_1), ("member 2", result_2), ("adhesive shear", result_a))
    governing_name, governing = min(ordered, key=lambda item: item[1].life_ratio)

    return ThermalCycleFatigueAssessment(
        requirement=requirement,
        member_1_fatigue=result_1,
        member_2_fatigue=result_2,
        adhesive_fatigue=result_a,
        static_yield_feasible=yield_assessment.passes,
        static_adhesive_feasible=shear_assessment.passes,
        minimum_yield_margin=yield_assessment.minimum_yield_margin,
        minimum_adhesive_margin=shear_assessment.minimum_margin,
        governing_fatigue_component=governing_name,
        minimum_life_ratio=governing.life_ratio,
    )


@dataclass(frozen=True)
class RestraintFatiguePoint:
    """Metal fatigue at one external-restraint level. Metal only - see module note."""

    stiffness_ratio: float
    member_1_cycle: StressCycle
    member_2_cycle: StressCycle
    member_1_fatigue: FatigueLifeResult
    member_2_fatigue: FatigueLifeResult

    @property
    def minimum_life_ratio(self) -> float:
        return min(self.member_1_fatigue.life_ratio, self.member_2_fatigue.life_ratio)

    @property
    def passes(self) -> bool:
        return self.member_1_fatigue.passes and self.member_2_fatigue.passes


def restraint_fatigue_sensitivity(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    curves: FatigueCurveSet,
    requirement: ThermalCycleRequirement,
    stiffness_ratios: Sequence[float],
    restraint_reference_strain: float = 0.0,
) -> tuple[RestraintFatiguePoint, ...]:
    """Metal fatigue across external restraint levels ``eta_r``.

    The adhesive is deliberately excluded: the Milestone 3 overlap model is
    posed for the free joint, and reusing it under external restraint would
    extend it past what has been verified.
    """
    points = []
    for ratio in stiffness_ratios:
        restraint = AxialRestraint.from_stiffness_ratio(
            require_finite(ratio, "stiffness_ratio"),
            member_1,
            member_2,
            reference_strain=restraint_reference_strain,
        )
        cycle_1, cycle_2 = member_stress_cycles(member_1, member_2, environment, restraint)
        points.append(
            RestraintFatiguePoint(
                stiffness_ratio=ratio,
                member_1_cycle=cycle_1,
                member_2_cycle=cycle_2,
                member_1_fatigue=assess_fatigue_life(
                    "member 1", cycle_1, curves.member_1, requirement
                ),
                member_2_fatigue=assess_fatigue_life(
                    "member 2", cycle_2, curves.member_2, requirement
                ),
            )
        )
    return tuple(points)


@dataclass(frozen=True)
class AllowableCycleScaleResult:
    """Largest scale on the canonical thermal excursion that still meets the life.

    Under the linear model every alternating stress is exactly proportional to
    the excursion scale, so this inversion is **exact**, not iterative::

        sigma_a,allow = A * N_required^b
        scale_i       = sigma_a,allow,i / sigma_a,i(canonical)

    The governing scale is the smallest of the three. At that scale the minimum
    life ratio is 1 to within floating point.
    """

    requirement: ThermalCycleRequirement
    member_1_scale: float
    member_2_scale: float
    adhesive_scale: float
    governing_component: str
    allowable_scale: float

    @property
    def canonical_cycle_passes(self) -> bool:
        """True when the unscaled canonical excursion already meets the life."""
        return self.allowable_scale >= 1.0


def _component_scale(result: FatigueLifeResult, requirement: ThermalCycleRequirement) -> float:
    allowable = result.fatigue_curve.alternating_stress_at_life(requirement.required_cycles)
    if result.alternating_stress == 0.0:
        return math.inf
    return allowable / result.alternating_stress


def allowable_cycle_scale_for_fatigue(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    curves: FatigueCurveSet,
    requirement: ThermalCycleRequirement,
) -> AllowableCycleScaleResult:
    """Exact allowable scale factor on the canonical thermal excursion."""
    assessment = assess_thermal_cycle_fatigue(
        member_1, member_2, environment, geometry, adhesive, curves, requirement
    )
    scales = {
        "member 1": _component_scale(assessment.member_1_fatigue, requirement),
        "member 2": _component_scale(assessment.member_2_fatigue, requirement),
        "adhesive shear": _component_scale(assessment.adhesive_fatigue, requirement),
    }
    governing = min(scales, key=lambda name: scales[name])
    return AllowableCycleScaleResult(
        requirement=requirement,
        member_1_scale=scales["member 1"],
        member_2_scale=scales["member 2"],
        adhesive_scale=scales["adhesive shear"],
        governing_component=governing,
        allowable_scale=scales[governing],
    )


def scale_environment(environment: ThermalEnvironment, scale: float) -> ThermalEnvironment:
    """Scale both excursions about ``T_ref``, preserving their asymmetry ratio."""
    factor = require_finite(scale, "scale")
    if factor < 0.0:
        raise ValueError(f"scale must be >= 0, got {factor!r}")
    reference = environment.reference_temperature
    return ThermalEnvironment(
        reference_temperature=reference,
        cold_temperature=reference + factor * (environment.cold_temperature - reference),
        hot_temperature=reference + factor * (environment.hot_temperature - reference),
    )


@dataclass(frozen=True)
class FatigueSweepPoint:
    """One point of a one-variable fatigue sweep."""

    variable: str
    value: float
    assessment: ThermalCycleFatigueAssessment

    @property
    def minimum_life_ratio(self) -> float:
        return self.assessment.minimum_life_ratio

    @property
    def overall_feasible(self) -> bool:
        return self.assessment.overall_feasible


def _sweep_point(variable, value, member_1, member_2, environment, geometry, adhesive,
                 curves, requirement, yield_basis, adhesive_basis) -> FatigueSweepPoint:
    return FatigueSweepPoint(
        variable=variable,
        value=value,
        assessment=assess_thermal_cycle_fatigue(
            member_1, member_2, environment, geometry, adhesive, curves, requirement,
            yield_basis, adhesive_basis,
        ),
    )


def temperature_scale_fatigue_sensitivity(
    member_1, member_2, environment, geometry, adhesive, curves, requirement,
    scales: Sequence[float], yield_basis=None, adhesive_basis=None,
) -> tuple[FatigueSweepPoint, ...]:
    """Scale both excursions by each factor, preserving the hot/cold asymmetry."""
    return tuple(
        _sweep_point(
            "temperature_scale", scale, member_1, member_2,
            scale_environment(environment, scale), geometry, adhesive,
            curves, requirement, yield_basis, adhesive_basis,
        )
        for scale in scales
    )


def bond_width_fatigue_sensitivity(
    member_1, member_2, environment, geometry, adhesive, curves, requirement,
    bond_widths: Sequence[float], yield_basis=None, adhesive_basis=None,
) -> tuple[FatigueSweepPoint, ...]:
    """Sweep bond width [m]; copies of the geometry are built for each point."""
    return tuple(
        _sweep_point(
            "bond_width", width, member_1, member_2, environment,
            BondedOverlapGeometry(
                overlap_length=geometry.overlap_length,
                bond_width=require_positive(width, "bond_width"),
                adhesive_thickness=geometry.adhesive_thickness,
            ),
            adhesive, curves, requirement, yield_basis, adhesive_basis,
        )
        for width in bond_widths
    )


def adhesive_thickness_fatigue_sensitivity(
    member_1, member_2, environment, geometry, adhesive, curves, requirement,
    adhesive_thicknesses: Sequence[float], yield_basis=None, adhesive_basis=None,
) -> tuple[FatigueSweepPoint, ...]:
    """Sweep bondline thickness [m]."""
    return tuple(
        _sweep_point(
            "adhesive_thickness", thickness, member_1, member_2, environment,
            geometry.with_adhesive_thickness(thickness), adhesive,
            curves, requirement, yield_basis, adhesive_basis,
        )
        for thickness in adhesive_thicknesses
    )


def adhesive_modulus_fatigue_sensitivity(
    member_1, member_2, environment, geometry, adhesive, curves, requirement,
    shear_moduli: Sequence[float], yield_basis=None, adhesive_basis=None,
) -> tuple[FatigueSweepPoint, ...]:
    """Sweep adhesive shear modulus [Pa]; strength and curve are held fixed."""
    return tuple(
        _sweep_point(
            "adhesive_shear_modulus", modulus, member_1, member_2, environment,
            geometry, _modulus_variant(adhesive, modulus),
            curves, requirement, yield_basis, adhesive_basis,
        )
        for modulus in shear_moduli
    )


def area_ratio_fatigue_sensitivity(
    member_1, member_2, environment, geometry, adhesive, curves, requirement,
    area_ratios: Sequence[float], yield_basis=None, adhesive_basis=None,
) -> tuple[FatigueSweepPoint, ...]:
    """Sweep ``A_1 / A_2`` end to end, recomputing every upstream quantity."""
    points = []
    for ratio in area_ratios:
        scaled = AxialMember(
            member_1.material,
            require_positive(ratio, "area_ratio") * member_2.area,
            label=member_1.label,
        )
        points.append(
            _sweep_point(
                "area_ratio", ratio, scaled, member_2, environment, geometry,
                adhesive, curves, requirement, yield_basis, adhesive_basis,
            )
        )
    return tuple(points)


@dataclass(frozen=True)
class FatigueDesignCandidate:
    """One grid point of the bounded fatigue design map."""

    bond_width: float
    adhesive_thickness: float
    bond_area: float
    assessment: ThermalCycleFatigueAssessment

    @property
    def overall_feasible(self) -> bool:
        return self.assessment.overall_feasible


def fatigue_design_map(
    member_1, member_2, environment, geometry, adhesive, curves, requirement,
    bond_widths: Sequence[float], adhesive_thicknesses: Sequence[float],
    yield_basis=None, adhesive_basis=None,
) -> tuple[FatigueDesignCandidate, ...]:
    """A **bounded preliminary design map** requiring static *and* fatigue feasibility.

    Not an optimization: a deterministic row-major evaluation of an explicitly
    bounded grid, using only the verified mechanics plus the Milestone 5
    screens.
    """
    candidates = []
    for width in bond_widths:
        for thickness in adhesive_thicknesses:
            local = BondedOverlapGeometry(
                overlap_length=geometry.overlap_length,
                bond_width=require_positive(width, "bond_width"),
                adhesive_thickness=require_positive(thickness, "adhesive_thickness"),
            )
            candidates.append(
                FatigueDesignCandidate(
                    bond_width=local.bond_width,
                    adhesive_thickness=local.adhesive_thickness,
                    bond_area=local.bond_area,
                    assessment=assess_thermal_cycle_fatigue(
                        member_1, member_2, environment, local, adhesive, curves,
                        requirement, yield_basis, adhesive_basis,
                    ),
                )
            )
    return tuple(candidates)


def select_preliminary_fatigue_design(
    candidates: Sequence[FatigueDesignCandidate],
) -> Optional[FatigueDesignCandidate]:
    """Pick one feasible point using the Milestone 4 geometry policy.

    Same deterministic rule as
    :func:`~thermal_joint.design_trade.select_preliminary_joint_design`:
    feasible first, then smallest bond area, then smaller adhesive thickness,
    then larger minimum life ratio, then input order.

    Returns ``None`` when nothing in the grid passes all three screens. The
    result is a preliminary feasible point within the bounded study grid, not
    an optimized flight-joint design.
    """
    feasible = [
        (index, candidate)
        for index, candidate in enumerate(candidates)
        if candidate.overall_feasible
    ]
    if not feasible:
        return None

    def key(item):
        index, candidate = item
        return (
            candidate.bond_area,
            candidate.adhesive_thickness,
            -candidate.assessment.minimum_life_ratio,
            index,
        )

    return min(feasible, key=key)[1]


__all__.extend(["FatigueDesignCandidate", "select_preliminary_fatigue_design",
                "scale_environment"])

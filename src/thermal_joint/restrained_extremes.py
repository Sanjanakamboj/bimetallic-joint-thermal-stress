"""Hot/cold extreme assessment of the bonded joint under external restraint.

This is the Milestone 2 counterpart of
:mod:`thermal_joint.extremes`. It reuses the Milestone 1 yield machinery
(:class:`~thermal_joint.margins.YieldBasis` and
:func:`~thermal_joint.margins.assess_yield`) unchanged - the margin definition,
the zero-stress handling and the governing-member rule are all identical.
"""

from __future__ import annotations

from dataclasses import dataclass

from .environment import ThermalEnvironment
from .margins import JointYieldAssessment, YieldBasis, assess_yield
from .materials import AxialMember
from .restraint import (
    AxialRestraint,
    RestrainedThermalJointResult,
    solve_restrained_joint,
)

__all__ = [
    "RestrainedTemperatureExtremeAssessment",
    "assess_restrained_temperature_extremes",
]


@dataclass(frozen=True)
class RestrainedTemperatureExtremeAssessment:
    """Restrained joint state and yield margins at both temperature extremes.

    Attributes
    ----------
    environment:
        The :class:`~thermal_joint.environment.ThermalEnvironment` evaluated.
    restraint:
        The :class:`~thermal_joint.restraint.AxialRestraint` applied.
    basis:
        The :class:`~thermal_joint.margins.YieldBasis` applied.
    cold_result, hot_result:
        Restrained joint results at each extreme.
    cold_assessment, hot_assessment:
        Per-extreme yield assessments.
    governing_extreme:
        ``"cold"`` or ``"hot"``, computed from the actual margins. An exact tie
        resolves to ``"cold"``, matching Milestone 1.
    governing_member:
        ``1`` or ``2`` at the governing extreme.
    minimum_yield_margin:
        Smallest preliminary elastic yield margin over both extremes and both
        members [-].
    """

    environment: ThermalEnvironment
    restraint: AxialRestraint
    basis: YieldBasis
    cold_result: RestrainedThermalJointResult
    hot_result: RestrainedThermalJointResult
    cold_assessment: JointYieldAssessment
    hot_assessment: JointYieldAssessment
    governing_extreme: str
    governing_member: int
    minimum_yield_margin: float

    @property
    def governing_result(self) -> RestrainedThermalJointResult:
        """Restrained joint result at the governing extreme."""
        return self.cold_result if self.governing_extreme == "cold" else self.hot_result

    @property
    def governing_assessment(self) -> JointYieldAssessment:
        """Yield assessment at the governing extreme."""
        return (
            self.cold_assessment
            if self.governing_extreme == "cold"
            else self.hot_assessment
        )

    @property
    def passes(self) -> bool:
        """True when the minimum margin over both extremes is ``>= 0``."""
        return self.minimum_yield_margin >= 0.0


def assess_restrained_temperature_extremes(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    restraint: AxialRestraint | None = None,
    basis: YieldBasis | None = None,
) -> RestrainedTemperatureExtremeAssessment:
    """Solve and assess the restrained joint at the cold and hot extremes.

    Parameters
    ----------
    member_1, member_2:
        The two bonded members.
    environment:
        Reference, cold and hot temperatures.
    restraint:
        External restraint; defaults to ``AxialRestraint(0.0)`` (free joint).
    basis:
        Design basis; defaults to ``YieldBasis()`` (design factor 1.0).

    Returns
    -------
    RestrainedTemperatureExtremeAssessment
        Both results, both assessments, and the governing extreme, governing
        member and minimum margin computed from the actual margins.
    """
    if not isinstance(environment, ThermalEnvironment):
        raise TypeError(
            f"environment must be a ThermalEnvironment, got {type(environment).__name__}"
        )
    external = AxialRestraint(0.0) if restraint is None else restraint
    design_basis = YieldBasis() if basis is None else basis

    cold_result = solve_restrained_joint(
        member_1, member_2, environment.delta_temperature_cold, external
    )
    hot_result = solve_restrained_joint(
        member_1, member_2, environment.delta_temperature_hot, external
    )

    cold_assessment = assess_yield(cold_result, member_1, member_2, design_basis)
    hot_assessment = assess_yield(hot_result, member_1, member_2, design_basis)

    if hot_assessment.minimum_margin < cold_assessment.minimum_margin:
        governing_extreme = "hot"
        governing_assessment = hot_assessment
    else:
        governing_extreme = "cold"
        governing_assessment = cold_assessment

    return RestrainedTemperatureExtremeAssessment(
        environment=environment,
        restraint=external,
        basis=design_basis,
        cold_result=cold_result,
        hot_result=hot_result,
        cold_assessment=cold_assessment,
        hot_assessment=hot_assessment,
        governing_extreme=governing_extreme,
        governing_member=governing_assessment.governing_member,
        minimum_yield_margin=governing_assessment.minimum_margin,
    )

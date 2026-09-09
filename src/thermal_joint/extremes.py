"""Evaluation of a bonded joint at both cold and hot temperature extremes."""

from __future__ import annotations

from dataclasses import dataclass

from .environment import ThermalEnvironment
from .joint import ThermalJointResult, solve_bimetallic_joint
from .margins import JointYieldAssessment, YieldBasis, assess_yield
from .materials import AxialMember

__all__ = ["TemperatureExtremeAssessment", "assess_temperature_extremes"]


@dataclass(frozen=True)
class TemperatureExtremeAssessment:
    """Joint state and preliminary yield margins at both temperature extremes.

    Attributes
    ----------
    environment:
        The :class:`~thermal_joint.environment.ThermalEnvironment` evaluated.
    basis:
        The :class:`~thermal_joint.margins.YieldBasis` applied.
    cold_result, hot_result:
        Solved :class:`~thermal_joint.joint.ThermalJointResult` at each extreme.
    cold_assessment, hot_assessment:
        Per-extreme :class:`~thermal_joint.margins.JointYieldAssessment`.
    governing_extreme:
        ``"cold"`` or ``"hot"`` - whichever yields the smaller minimum margin.
        Computed from the margins, never assumed. An exact tie (for instance
        symmetric excursions with temperature-independent properties) resolves
        to ``"cold"``.
    governing_member:
        ``1`` or ``2`` - the governing member at the governing extreme.
    minimum_yield_margin:
        The smallest preliminary elastic yield margin over both extremes and
        both members [-].
    """

    environment: ThermalEnvironment
    basis: YieldBasis
    cold_result: ThermalJointResult
    hot_result: ThermalJointResult
    cold_assessment: JointYieldAssessment
    hot_assessment: JointYieldAssessment
    governing_extreme: str
    governing_member: int
    minimum_yield_margin: float

    @property
    def governing_result(self) -> ThermalJointResult:
        """Joint result at the governing extreme."""
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


def assess_temperature_extremes(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    basis: YieldBasis | None = None,
) -> TemperatureExtremeAssessment:
    """Solve and assess the joint at the cold and hot temperature extremes.

    Parameters
    ----------
    member_1, member_2:
        The two bonded members.
    environment:
        Reference, cold and hot temperatures.
    basis:
        Design basis; defaults to ``YieldBasis()`` (design factor 1.0).

    Returns
    -------
    TemperatureExtremeAssessment
        Both joint results, both yield assessments, and the governing extreme,
        governing member and minimum margin computed from the actual margins.
    """
    if not isinstance(environment, ThermalEnvironment):
        raise TypeError(
            "environment must be a ThermalEnvironment, got "
            f"{type(environment).__name__}"
        )
    design_basis = YieldBasis() if basis is None else basis

    cold_result = solve_bimetallic_joint(
        member_1, member_2, environment.delta_temperature_cold
    )
    hot_result = solve_bimetallic_joint(
        member_1, member_2, environment.delta_temperature_hot
    )

    cold_assessment = assess_yield(cold_result, member_1, member_2, design_basis)
    hot_assessment = assess_yield(hot_result, member_1, member_2, design_basis)

    if hot_assessment.minimum_margin < cold_assessment.minimum_margin:
        governing_extreme = "hot"
        governing_assessment = hot_assessment
    else:
        governing_extreme = "cold"
        governing_assessment = cold_assessment

    return TemperatureExtremeAssessment(
        environment=environment,
        basis=design_basis,
        cold_result=cold_result,
        hot_result=hot_result,
        cold_assessment=cold_assessment,
        hot_assessment=hot_assessment,
        governing_extreme=governing_extreme,
        governing_member=governing_assessment.governing_member,
        minimum_yield_margin=governing_assessment.minimum_margin,
    )

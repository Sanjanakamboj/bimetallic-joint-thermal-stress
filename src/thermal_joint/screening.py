"""Combined global + local preliminary screening.

Pairs the Milestone 1 member yield assessment with the Milestone 3 adhesive
shear assessment. The two margins are reported **separately and are never
numerically combined** - they measure different failure modes against different
allowables, and averaging or minimising across them would be meaningless. Only
the pass/fail booleans are combined, with a plain boolean AND.
"""

from __future__ import annotations

from dataclasses import dataclass

from .adhesive import AdhesiveMaterial, AdhesiveShearBasis, BondedOverlapGeometry
from .environment import ThermalEnvironment
from .extremes import TemperatureExtremeAssessment, assess_temperature_extremes
from .margins import YieldBasis
from .materials import AxialMember
from .shear_lag_extremes import ShearLagExtremeAssessment, assess_shear_lag_extremes

__all__ = ["PreliminaryScreeningResult", "screen_joint"]


@dataclass(frozen=True)
class PreliminaryScreeningResult:
    """Global yield and local adhesive shear screens, reported side by side.

    Attributes
    ----------
    yield_assessment:
        Milestone 1 free-joint member yield assessment over both extremes.
    adhesive_assessment:
        Milestone 3 adhesive shear assessment over both extremes.
    """

    yield_assessment: TemperatureExtremeAssessment
    adhesive_assessment: ShearLagExtremeAssessment

    @property
    def minimum_yield_margin(self) -> float:
        """Minimum preliminary elastic yield margin [-]. Kept separate."""
        return self.yield_assessment.minimum_yield_margin

    @property
    def minimum_adhesive_margin(self) -> float:
        """Minimum preliminary adhesive shear margin [-]. Kept separate."""
        return self.adhesive_assessment.minimum_margin

    @property
    def yield_feasible(self) -> bool:
        """True when member yield passes at both extremes."""
        return self.yield_assessment.passes

    @property
    def adhesive_feasible(self) -> bool:
        """True when adhesive shear passes at both extremes."""
        return self.adhesive_assessment.passes

    @property
    def overall_feasible(self) -> bool:
        """Boolean AND of the two screens - not a combined margin."""
        return self.yield_feasible and self.adhesive_feasible

    @property
    def governing_mode(self) -> str:
        """Which screen fails, as text: ``"yield"``, ``"adhesive shear"``, both or none.

        A label for reporting only; it does not rank the two margins numerically.
        """
        if self.overall_feasible:
            return "none - both screens pass"
        failing = []
        if not self.yield_feasible:
            failing.append("member yield")
        if not self.adhesive_feasible:
            failing.append("adhesive shear")
        return " and ".join(failing)


def screen_joint(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    yield_basis: YieldBasis | None = None,
    adhesive_basis: AdhesiveShearBasis | None = None,
) -> PreliminaryScreeningResult:
    """Run the global yield screen and the local adhesive shear screen together.

    Both use the free-joint (Milestone 1) thermal state, so the adhesive screen
    sees exactly the mismatch force the yield screen is based on.
    """
    return PreliminaryScreeningResult(
        yield_assessment=assess_temperature_extremes(
            member_1, member_2, environment, yield_basis
        ),
        adhesive_assessment=assess_shear_lag_extremes(
            member_1, member_2, environment, geometry, adhesive, adhesive_basis
        ),
    )

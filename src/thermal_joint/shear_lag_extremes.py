"""Hot/cold extreme assessment of the adhesive shear-lag screen.

Mirrors the Milestone 1/2 extreme-assessment pattern: both extremes are solved,
both margins are computed, and the governing extreme is derived from the actual
margins rather than assumed. An exact tie resolves to ``"cold"``, matching the
earlier milestones.
"""

from __future__ import annotations

from dataclasses import dataclass

from .adhesive import AdhesiveMaterial, AdhesiveShearBasis, BondedOverlapGeometry
from .environment import ThermalEnvironment
from .materials import AxialMember
from .shear_lag import ShearLagDemand, solve_shear_lag

__all__ = ["ShearLagExtremeAssessment", "assess_shear_lag_extremes"]


@dataclass(frozen=True)
class ShearLagExtremeAssessment:
    """Adhesive shear demand and margins at both temperature extremes.

    Attributes
    ----------
    environment, geometry, adhesive, basis:
        The inputs used.
    cold_demand, hot_demand:
        :class:`~thermal_joint.shear_lag.ShearLagDemand` at each extreme.
    cold_margin, hot_margin:
        Preliminary adhesive shear margins at each extreme [-].
    governing_extreme:
        ``"cold"`` or ``"hot"``, computed from the margins.
    minimum_margin:
        The smaller of the two margins [-].
    """

    environment: ThermalEnvironment
    geometry: BondedOverlapGeometry
    adhesive: AdhesiveMaterial
    basis: AdhesiveShearBasis
    cold_demand: ShearLagDemand
    hot_demand: ShearLagDemand
    cold_margin: float
    hot_margin: float
    governing_extreme: str
    minimum_margin: float

    @property
    def governing_demand(self) -> ShearLagDemand:
        """The :class:`~thermal_joint.shear_lag.ShearLagDemand` that governs."""
        return self.cold_demand if self.governing_extreme == "cold" else self.hot_demand

    @property
    def governing_peak_shear_stress(self) -> float:
        """Peak adhesive shear stress at the governing extreme [Pa]."""
        return self.governing_demand.peak_shear_stress

    @property
    def allowable_shear_stress(self) -> float:
        """Adhesive shear allowable used [Pa]."""
        return self.basis.allowable_shear_stress(self.adhesive)

    @property
    def passes(self) -> bool:
        """True when the minimum preliminary adhesive shear margin is ``>= 0``."""
        return self.minimum_margin >= 0.0


def assess_shear_lag_extremes(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    basis: AdhesiveShearBasis | None = None,
) -> ShearLagExtremeAssessment:
    """Solve and assess the bonded overlap at the cold and hot extremes.

    Parameters
    ----------
    member_1, member_2:
        The two bonded adherends.
    environment:
        Reference, cold and hot temperatures.
    geometry:
        Bonded overlap geometry.
    adhesive:
        Adhesive material.
    basis:
        Shear design basis; defaults to ``AdhesiveShearBasis()`` (factor 1.0).
    """
    if not isinstance(environment, ThermalEnvironment):
        raise TypeError(
            f"environment must be a ThermalEnvironment, got {type(environment).__name__}"
        )
    design_basis = AdhesiveShearBasis() if basis is None else basis
    if not isinstance(design_basis, AdhesiveShearBasis):
        raise TypeError(
            f"basis must be an AdhesiveShearBasis, got {type(design_basis).__name__}"
        )

    cold_demand = solve_shear_lag(
        member_1, member_2, geometry, adhesive, environment.delta_temperature_cold
    )
    hot_demand = solve_shear_lag(
        member_1, member_2, geometry, adhesive, environment.delta_temperature_hot
    )

    cold_margin = design_basis.margin_of_safety(adhesive, cold_demand.peak_shear_stress)
    hot_margin = design_basis.margin_of_safety(adhesive, hot_demand.peak_shear_stress)

    if hot_margin < cold_margin:
        governing_extreme = "hot"
        minimum_margin = hot_margin
    else:
        governing_extreme = "cold"
        minimum_margin = cold_margin

    return ShearLagExtremeAssessment(
        environment=environment,
        geometry=geometry,
        adhesive=adhesive,
        basis=design_basis,
        cold_demand=cold_demand,
        hot_demand=hot_demand,
        cold_margin=cold_margin,
        hot_margin=hot_margin,
        governing_extreme=governing_extreme,
        minimum_margin=minimum_margin,
    )

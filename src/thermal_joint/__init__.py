"""Axial thermal-stress model for a perfectly bonded bimetallic joint.

Milestone 1 scope: two dissimilar isotropic members bonded so that they share a
common axial strain under a uniform temperature change, with no external axial
load. Provides the closed-form common strain, member thermal stresses, internal
forces, an explicit force-equilibrium residual, preliminary elastic yield
margins and hot/cold extreme evaluation.

Milestone 2 adds finite external axial restraint from the surrounding
structure: :class:`AxialRestraint`, :func:`solve_restrained_joint`,
:func:`assess_restrained_temperature_extremes` and a preliminary inverse-design
search for the maximum tolerable restraint stiffness. The Milestone 1
free-joint solution is preserved exactly as the zero-restraint limit.

Milestone 3 adds a *local* load-transfer layer: a one-dimensional linear-elastic
adhesive shear-lag screen for a finite bonded overlap
(:class:`AdhesiveMaterial`, :class:`BondedOverlapGeometry`,
:func:`solve_shear_lag`, :func:`assess_shear_lag_extremes`), plus allowable-dT
and minimum-overlap inverse design. The shear-lag demand is imported from the
Milestone 1 free-joint solution, never re-derived.

Milestone 4 adds a bounded preliminary *design trade* on top of those verified
mechanics: analytical asymptote scaling, one-variable sensitivity sweeps,
inverse sizing for bond width, adhesive thickness and adhesive shear modulus, a
bounded width x thickness design map, and a deterministic selection policy. It
does not optimize a flight joint.

SI units throughout: E [Pa], alpha [1/K], area [m^2], stress [Pa], force [N],
strain [-]. Temperature differences may be K or degC increments.

Sign convention: tension positive, compression negative, positive dT = heating.
"""

from __future__ import annotations

from .adhesive import AdhesiveMaterial, AdhesiveShearBasis, BondedOverlapGeometry
from .design_inverse import (
    AdhesiveModulusStatus,
    AdhesiveThicknessStatus,
    BondWidthStatus,
    MaximumAdhesiveModulusResult,
    RequiredAdhesiveThicknessResult,
    RequiredBondWidthResult,
    maximum_allowable_adhesive_shear_modulus,
    required_adhesive_thickness,
    required_bond_width,
)
from .design_trade import (
    DesignSelectionPolicy,
    JointDesignCandidate,
    SensitivityPoint,
    adhesive_modulus_sensitivity,
    adhesive_thickness_sensitivity,
    area_ratio_sensitivity,
    bond_width_sensitivity,
    cte_mismatch_sensitivity,
    evaluate_joint_design,
    long_overlap_peak_shear,
    select_preliminary_joint_design,
    thermal_excursion_sensitivity,
    uniform_shear_floor,
    width_thickness_design_map,
)
from .environment import ThermalEnvironment
from .extremes import TemperatureExtremeAssessment, assess_temperature_extremes
from .inverse import (
    MaximumRestraintResult,
    RestraintLimitStatus,
    maximum_allowable_restraint_stiffness,
    rigid_restraint_minimum_margin,
)
from .joint import ThermalJointResult, common_joint_strain, solve_bimetallic_joint
from .margins import (
    JointYieldAssessment,
    MemberStressState,
    MemberYieldMargin,
    YieldBasis,
    assess_yield,
)
from .materials import AxialMember, ThermoelasticMaterial
from .restrained_extremes import (
    RestrainedTemperatureExtremeAssessment,
    assess_restrained_temperature_extremes,
)
from .screening import PreliminaryScreeningResult, screen_joint
from .shear_lag import (
    ShearLagDemand,
    adherend_compliance_sum,
    dimensionless_overlap,
    long_overlap_peak_shear_stress,
    shear_lag_parameter,
    solve_shear_lag,
    thermal_mismatch_strain,
    transfer_length,
)
from .shear_lag_extremes import ShearLagExtremeAssessment, assess_shear_lag_extremes
from .shear_lag_inverse import (
    OverlapLimitStatus,
    RequiredOverlapResult,
    allowable_temperature_change_for_adhesive_shear,
    required_overlap_length,
)
from .restraint import (
    AxialRestraint,
    RestrainedThermalJointResult,
    joint_axial_stiffness,
    restraint_stiffness_ratio,
    rigid_restraint_stresses,
    solve_restrained_joint,
    zero_stress_restraint_stiffness,
)

__all__ = [
    # Milestone 1 - free bonded joint
    "AxialMember",
    "JointYieldAssessment",
    "MemberStressState",
    "MemberYieldMargin",
    "TemperatureExtremeAssessment",
    "ThermalEnvironment",
    "ThermalJointResult",
    "ThermoelasticMaterial",
    "YieldBasis",
    "assess_temperature_extremes",
    "assess_yield",
    "common_joint_strain",
    "solve_bimetallic_joint",
    # Milestone 2 - external axial restraint
    "AxialRestraint",
    "MaximumRestraintResult",
    "RestrainedTemperatureExtremeAssessment",
    "RestrainedThermalJointResult",
    "RestraintLimitStatus",
    "assess_restrained_temperature_extremes",
    "joint_axial_stiffness",
    "maximum_allowable_restraint_stiffness",
    "restraint_stiffness_ratio",
    "rigid_restraint_minimum_margin",
    "rigid_restraint_stresses",
    "solve_restrained_joint",
    "zero_stress_restraint_stiffness",
    # Milestone 3 - adhesive shear-lag over a finite bonded overlap
    "AdhesiveMaterial",
    "AdhesiveShearBasis",
    "BondedOverlapGeometry",
    "OverlapLimitStatus",
    "PreliminaryScreeningResult",
    "RequiredOverlapResult",
    "ShearLagDemand",
    "ShearLagExtremeAssessment",
    "adherend_compliance_sum",
    "allowable_temperature_change_for_adhesive_shear",
    "assess_shear_lag_extremes",
    "dimensionless_overlap",
    "long_overlap_peak_shear_stress",
    "required_overlap_length",
    "screen_joint",
    "shear_lag_parameter",
    "solve_shear_lag",
    "thermal_mismatch_strain",
    "transfer_length",
    # Milestone 4 - bounded preliminary design trade
    "AdhesiveModulusStatus",
    "AdhesiveThicknessStatus",
    "BondWidthStatus",
    "DesignSelectionPolicy",
    "JointDesignCandidate",
    "MaximumAdhesiveModulusResult",
    "RequiredAdhesiveThicknessResult",
    "RequiredBondWidthResult",
    "SensitivityPoint",
    "adhesive_modulus_sensitivity",
    "adhesive_thickness_sensitivity",
    "area_ratio_sensitivity",
    "bond_width_sensitivity",
    "cte_mismatch_sensitivity",
    "evaluate_joint_design",
    "long_overlap_peak_shear",
    "maximum_allowable_adhesive_shear_modulus",
    "required_adhesive_thickness",
    "required_bond_width",
    "select_preliminary_joint_design",
    "thermal_excursion_sensitivity",
    "uniform_shear_floor",
    "width_thickness_design_map",
]

__version__ = "0.4.0"

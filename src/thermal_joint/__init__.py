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

SI units throughout: E [Pa], alpha [1/K], area [m^2], stress [Pa], force [N],
strain [-]. Temperature differences may be K or degC increments.

Sign convention: tension positive, compression negative, positive dT = heating.
"""

from __future__ import annotations

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
]

__version__ = "0.2.0"

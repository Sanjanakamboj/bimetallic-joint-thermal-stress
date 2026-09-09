"""Axial thermal-stress model for a perfectly bonded bimetallic joint.

Milestone 1 scope: two dissimilar isotropic members bonded so that they share a
common axial strain under a uniform temperature change, with no external axial
load. Provides the closed-form common strain, member thermal stresses, internal
forces, an explicit force-equilibrium residual, preliminary elastic yield
margins and hot/cold extreme evaluation.

SI units throughout: E [Pa], alpha [1/K], area [m^2], stress [Pa], force [N],
strain [-]. Temperature differences may be K or degC increments.

Sign convention: tension positive, compression negative, positive dT = heating.
"""

from __future__ import annotations

from .environment import ThermalEnvironment
from .extremes import TemperatureExtremeAssessment, assess_temperature_extremes
from .joint import ThermalJointResult, common_joint_strain, solve_bimetallic_joint
from .margins import (
    JointYieldAssessment,
    MemberYieldMargin,
    YieldBasis,
    assess_yield,
)
from .materials import AxialMember, ThermoelasticMaterial

__all__ = [
    "AxialMember",
    "JointYieldAssessment",
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
]

__version__ = "0.1.0"

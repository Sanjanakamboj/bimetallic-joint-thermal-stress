"""Shared canonical fixtures for the Milestone 5 thermal-cycle fatigue tests.

Canonical configuration is the untouched Milestone 1-4 one, plus the
illustrative fatigue curves and a 1e4-cycle illustrative requirement.

Locked canonical cycle values (all verified independently in the tests):

    free member 1   hot -62.0278  cold +86.8389  -> sig_a 74.4333, sig_m +12.4056
    free member 2   hot +62.0278  cold -86.8389  -> sig_a 74.4333, sig_m -12.4056
    M3 adhesive     hot +47.4189  cold -66.3864  -> tau_a 56.9027, tau_m  -9.4838
    M4 adhesive     hot +13.4167  cold -18.7833  -> tau_a 16.1000, tau_m  -2.6833

Predicted lives at those amplitudes:

    member 1  1.049e9    member 2  1.962e14    M3 adhesive  1.42    M4 adhesive  6439
"""

from __future__ import annotations

from thermal_joint import AdhesiveShearBasis, YieldBasis
from thermal_joint.illustrative import (
    ADHESIVE_LIKE,
    ADHESIVE_SHEAR_FATIGUE,
    ALUMINIUM_LIKE_FATIGUE,
    TITANIUM_LIKE_FATIGUE,
    aluminium_like_member,
    illustrative_cycle_requirement,
    illustrative_fatigue_curves,
    radiator_joint_environment,
    radiator_joint_overlap,
    titanium_like_member,
)

__all__ = [
    "MEMBER_AL", "MEMBER_TI", "ENVIRONMENT", "ADHESIVE", "CURVES", "REQUIREMENT",
    "YIELD_BASIS", "ADHESIVE_BASIS", "BASELINE_GEOMETRY", "SELECTED_GEOMETRY", "MM",
    "FREE_ALTERNATING", "FREE_MEAN_1", "FREE_MEAN_2",
    "BASELINE_TAU_A", "BASELINE_TAU_M", "SELECTED_TAU_A", "SELECTED_TAU_M",
    "LIFE_MEMBER_1", "LIFE_MEMBER_2", "LIFE_BASELINE_ADHESIVE", "LIFE_SELECTED_ADHESIVE",
    "AL_FATIGUE", "TI_FATIGUE", "ADHESIVE_FATIGUE",
]

MM = 1.0e-3

MEMBER_AL = aluminium_like_member(100.0)
MEMBER_TI = titanium_like_member(100.0)
ENVIRONMENT = radiator_joint_environment()
ADHESIVE = ADHESIVE_LIKE
CURVES = illustrative_fatigue_curves()
REQUIREMENT = illustrative_cycle_requirement()
YIELD_BASIS = YieldBasis(1.25)
ADHESIVE_BASIS = AdhesiveShearBasis(1.25)

AL_FATIGUE = ALUMINIUM_LIKE_FATIGUE
TI_FATIGUE = TITANIUM_LIKE_FATIGUE
ADHESIVE_FATIGUE = ADHESIVE_SHEAR_FATIGUE

BASELINE_GEOMETRY = radiator_joint_overlap()                    # 40 x 20 mm, 0.2 mm
SELECTED_GEOMETRY = radiator_joint_overlap(40.0, 50.0, 1.0)     # M4 selected point

FREE_ALTERNATING = 74.4333e6
FREE_MEAN_1 = +12.4056e6
FREE_MEAN_2 = -12.4056e6

BASELINE_TAU_A = 56.9027e6
BASELINE_TAU_M = -9.4838e6
SELECTED_TAU_A = 16.1000e6
SELECTED_TAU_M = -2.6833e6

LIFE_MEMBER_1 = 1.0486e9
LIFE_MEMBER_2 = 1.9620e14
LIFE_BASELINE_ADHESIVE = 1.4238
LIFE_SELECTED_ADHESIVE = 6439.0

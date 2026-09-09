"""Shared canonical fixtures for the Milestone 4 design-trade tests.

Everything here is the untouched Milestone 3 canonical configuration:

    members    aluminium-like / titanium-like, 100 mm^2 each
    environment T_ref = +20 degC, T_cold = -120 degC, T_hot = +120 degC
    overlap    L_b = 40 mm, b = 20 mm, t_a = 0.2 mm
    adhesive   G_a = 1.0 GPa, shear strength 25 MPa
    bases      yield DF 1.25, adhesive DF 1.25 -> tau_allow = 20 MPa

Baseline results locked by Milestone 3:

    beta = 152.894 /m, transfer length 6.5405 mm, lambda = 6.1158
    cold tau_peak = 66.3864 MPa, adhesive margin -0.6987  (FAIL)
    hot  tau_peak = 47.4189 MPa, adhesive margin -0.5782  (FAIL)
    minimum member yield margin +1.4874                   (PASS)
    long-overlap asymptote 66.3858 MPa, uniform-shear floor 10.8549 MPa
"""

from __future__ import annotations

from thermal_joint import AdhesiveShearBasis, YieldBasis
from thermal_joint.illustrative import (
    ADHESIVE_LIKE,
    aluminium_like_member,
    radiator_joint_environment,
    radiator_joint_overlap,
    titanium_like_member,
)

__all__ = [
    "MEMBER_AL",
    "MEMBER_TI",
    "ENVIRONMENT",
    "GEOMETRY",
    "ADHESIVE",
    "YIELD_BASIS",
    "ADHESIVE_BASIS",
    "COLD",
    "HOT",
    "ALLOWABLE",
    "BASELINE_COLD_PEAK",
    "BASELINE_HOT_PEAK",
    "BASELINE_COLD_MARGIN",
    "BASELINE_ASYMPTOTE",
    "BASELINE_FLOOR",
    "BASELINE_YIELD_MARGIN",
    "MM",
]

MM = 1.0e-3

MEMBER_AL = aluminium_like_member(100.0)
MEMBER_TI = titanium_like_member(100.0)
ENVIRONMENT = radiator_joint_environment()
GEOMETRY = radiator_joint_overlap()
ADHESIVE = ADHESIVE_LIKE
YIELD_BASIS = YieldBasis(1.25)
ADHESIVE_BASIS = AdhesiveShearBasis(1.25)

COLD = -140.0
HOT = +100.0
ALLOWABLE = 20.0e6

BASELINE_COLD_PEAK = 66.3864e6
BASELINE_HOT_PEAK = 47.4189e6
BASELINE_COLD_MARGIN = -0.6987
BASELINE_ASYMPTOTE = 66.3858e6
BASELINE_FLOOR = 10.8549e6
BASELINE_YIELD_MARGIN = 1.4874

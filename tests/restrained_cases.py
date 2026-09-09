"""Shared hand-checkable inputs for the Milestone 2 restrained-joint tests.

Reuses the Milestone 1 hand-calculation members:

    E1 = 100 GPa, A1 = 200 mm^2, alpha_1 = 20e-6 /K   -> E1 A1 = 2.0e7 N
    E2 = 200 GPa, A2 = 100 mm^2, alpha_2 = 10e-6 /K   -> E2 A2 = 2.0e7 N
    K_joint = 4.0e7 N,  S = E1 A1 a1 + E2 A2 a2 = 600

With dT = +50 K and eps_ref = 0, the restrained closed form
``eps = (S dT + K eps_ref) / (K_joint + K)`` gives round numbers:

  K_r = 1.0e7 N (eta_r = 0.25):
      eps      = 30000 / 5.0e7 = 6.0e-4
      sigma_1  = 100e9 (6.0e-4 - 1.0e-3) = -40 MPa
      sigma_2  = 200e9 (6.0e-4 - 5.0e-4) = +20 MPa
      N_1      = -8000 N,  N_2 = +2000 N,  N_1 + N_2 = -6000 N
      N_r      = 1.0e7 (0 - 6.0e-4) = -6000 N
      residual = N_1 + N_2 - N_r = 0

  K_r = 2.0e7 N (eta_r = 0.5) is exactly the member-2 zero-stress crossing:
      eps      = 30000 / 6.0e7 = 5.0e-4 = alpha_2 dT  ->  sigma_2 = 0 exactly

  K_r = 1.0e7 N with eps_ref = 1.0e-3:
      eps      = (30000 + 10000) / 5.0e7 = 8.0e-4
      sigma_1  = -20 MPa,  sigma_2 = +60 MPa
      N_1      = -4000 N,  N_2 = +6000 N,  N_1 + N_2 = +2000 N
      N_r      = 1.0e7 (1.0e-3 - 8.0e-4) = +2000 N
"""

from __future__ import annotations

from cases import HAND_DELTA_T, HAND_MEMBER_1, HAND_MEMBER_2

from thermal_joint import AxialRestraint

__all__ = [
    "HAND_MEMBER_1",
    "HAND_MEMBER_2",
    "HAND_DELTA_T",
    "HAND_JOINT_STIFFNESS",
    "HAND_RESTRAINT",
    "HAND_RESTRAINED_COMMON_STRAIN",
    "HAND_RESTRAINED_STRESS_1",
    "HAND_RESTRAINED_STRESS_2",
    "HAND_RESTRAINED_FORCE_1",
    "HAND_RESTRAINED_FORCE_2",
    "HAND_RESTRAINT_FORCE",
    "HAND_CROSSING_RESTRAINT",
    "HAND_OFFSET_RESTRAINT",
    "HAND_OFFSET_COMMON_STRAIN",
    "HAND_OFFSET_STRESS_1",
    "HAND_OFFSET_STRESS_2",
    "HAND_OFFSET_RESTRAINT_FORCE",
]

HAND_JOINT_STIFFNESS = 4.0e7

HAND_RESTRAINT = AxialRestraint(stiffness=1.0e7, reference_strain=0.0, label="hand restraint")
HAND_RESTRAINED_COMMON_STRAIN = 6.0e-4
HAND_RESTRAINED_STRESS_1 = -40.0e6
HAND_RESTRAINED_STRESS_2 = +20.0e6
HAND_RESTRAINED_FORCE_1 = -8000.0
HAND_RESTRAINED_FORCE_2 = +2000.0
HAND_RESTRAINT_FORCE = -6000.0

HAND_CROSSING_RESTRAINT = AxialRestraint(stiffness=2.0e7)

HAND_OFFSET_RESTRAINT = AxialRestraint(stiffness=1.0e7, reference_strain=1.0e-3)
HAND_OFFSET_COMMON_STRAIN = 8.0e-4
HAND_OFFSET_STRESS_1 = -20.0e6
HAND_OFFSET_STRESS_2 = +60.0e6
HAND_OFFSET_RESTRAINT_FORCE = +2000.0

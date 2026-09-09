"""Shared hand-checkable inputs for the Milestone 3 shear-lag tests.

Reuses the Milestone 1 hand-calculation members, chosen so that beta comes out
exactly round:

    E1 A1 = E2 A2 = 2.0e7 N   ->  C = 1/2e7 + 1/2e7 = 1.0e-7 1/N
    G_a = 1.0 GPa, b = 20 mm, t_a = 0.2 mm  ->  G_a b / t_a = 1.0e11 N/m^2
    beta^2 = 1.0e11 * 1.0e-7 = 1.0e4      ->  beta = 100 /m exactly
    transfer length = 1/beta = 10 mm

With an overlap L_b = 20 mm the dimensionless overlap is lambda = beta L_b = 2.

At dT = +50 K (the Milestone 1 hand case):

    mismatch strain d_eps = (20e-6 - 10e-6) * 50 = 5.0e-4
    N_t = -d_eps / C = -5000 N        (equals the Milestone 1 member-1 force)
    tau_avg  = 5000 / (0.02 * 0.02)          = 12.5 MPa exactly
    tau_peak = 5000 * 100 * coth(2) / 0.02   = 25e6 * coth(2) = 25.9329 MPa
    tau(0)   = 25e6 / sinh(2)                = 6.8930 MPa
    peak/avg = 2 coth(2)                     = 2.0746
    N_1(10 mm) = -5000 (1 - sinh(1)/sinh(2)) = -3379.9 N
"""

from __future__ import annotations

from cases import HAND_DELTA_T, HAND_MEMBER_1, HAND_MEMBER_2

from thermal_joint import AdhesiveMaterial, AdhesiveShearBasis, BondedOverlapGeometry

__all__ = [
    "HAND_MEMBER_1",
    "HAND_MEMBER_2",
    "HAND_DELTA_T",
    "HAND_ADHESIVE",
    "HAND_GEOMETRY",
    "HAND_BETA",
    "HAND_TRANSFER_LENGTH",
    "HAND_LAMBDA",
    "HAND_MISMATCH_STRAIN",
    "HAND_TRANSFERRED_FORCE",
    "HAND_SIGNED_FORCE",
    "HAND_AVERAGE_SHEAR",
]

HAND_ADHESIVE = AdhesiveMaterial(
    name="Hand-calc adhesive",
    shear_modulus=1.0e9,
    shear_strength=25.0e6,
    source_note="ILLUSTRATIVE ADHESIVE INPUT - NOT DESIGN ALLOWABLE (test fixture)",
)

HAND_GEOMETRY = BondedOverlapGeometry(
    overlap_length=0.020,
    bond_width=0.020,
    adhesive_thickness=0.0002,
)

HAND_BASIS = AdhesiveShearBasis(1.0)

HAND_BETA = 100.0
HAND_TRANSFER_LENGTH = 0.010
HAND_LAMBDA = 2.0
HAND_MISMATCH_STRAIN = 5.0e-4
HAND_SIGNED_FORCE = -5000.0
HAND_TRANSFERRED_FORCE = 5000.0
HAND_AVERAGE_SHEAR = 12.5e6

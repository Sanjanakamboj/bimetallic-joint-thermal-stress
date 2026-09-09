"""Shared, deliberately hand-checkable inputs for the verification tests.

The primary hand-calculation case uses round numbers so the expected values can
be written down exactly:

    E1 = 100 GPa, A1 = 200 mm^2, alpha_1 = 20e-6 /K
    E2 = 200 GPa, A2 = 100 mm^2, alpha_2 = 10e-6 /K
    dT = +50 K

    E1 A1 = E2 A2 = 2.0e7 N,  sum = 4.0e7 N
    eps_common = (2e7*20e-6 + 2e7*10e-6) / 4e7 * 50 = 7.5e-4
    sigma_1 = 100e9 * (7.5e-4 - 1.0e-3) = -25 MPa
    sigma_2 = 200e9 * (7.5e-4 - 5.0e-4) = +50 MPa
    N_1 = -25e6 * 200e-6 = -5000 N,  N_2 = +50e6 * 100e-6 = +5000 N
"""

from __future__ import annotations

from thermal_joint import AxialMember, ThermoelasticMaterial

__all__ = [
    "HAND_MATERIAL_1",
    "HAND_MATERIAL_2",
    "HAND_MEMBER_1",
    "HAND_MEMBER_2",
    "HAND_DELTA_T",
    "HAND_COMMON_STRAIN",
    "HAND_STRESS_1",
    "HAND_STRESS_2",
    "HAND_FORCE_1",
    "HAND_FORCE_2",
]

HAND_MATERIAL_1 = ThermoelasticMaterial(
    name="Hand-calc material 1",
    elastic_modulus=100.0e9,
    thermal_expansion_coefficient=20.0e-6,
    yield_strength=300.0e6,
)
HAND_MATERIAL_2 = ThermoelasticMaterial(
    name="Hand-calc material 2",
    elastic_modulus=200.0e9,
    thermal_expansion_coefficient=10.0e-6,
    yield_strength=600.0e6,
)

HAND_MEMBER_1 = AxialMember(HAND_MATERIAL_1, 200.0e-6, label="hand-1")
HAND_MEMBER_2 = AxialMember(HAND_MATERIAL_2, 100.0e-6, label="hand-2")

HAND_DELTA_T = 50.0
HAND_COMMON_STRAIN = 7.5e-4
HAND_STRESS_1 = -25.0e6
HAND_STRESS_2 = +50.0e6
HAND_FORCE_1 = -5000.0
HAND_FORCE_2 = +5000.0

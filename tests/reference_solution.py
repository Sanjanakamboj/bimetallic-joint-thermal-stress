"""Independent closed-form reference solution used to verify the production code.

These formulas are derived separately from the production implementation: the
production code computes the stiffness-weighted common strain first and then
applies the constitutive relation, whereas these expressions give the member
stresses directly in terms of the CTE mismatch. They deliberately do not call
any function from :mod:`thermal_joint`.

    sigma_1 = [E1 E2 A2 / (E1 A1 + E2 A2)] * (alpha_2 - alpha_1) * dT
    sigma_2 = [E1 E2 A1 / (E1 A1 + E2 A2)] * (alpha_1 - alpha_2) * dT
"""

from __future__ import annotations

__all__ = ["reference_stress_1", "reference_stress_2", "reference_stresses"]


def reference_stress_1(
    e_1: float,
    a_1: float,
    alpha_1: float,
    e_2: float,
    a_2: float,
    alpha_2: float,
    delta_temperature: float,
) -> float:
    """Independent closed form for the member 1 axial thermal stress [Pa]."""
    return (
        (e_1 * e_2 * a_2 / (e_1 * a_1 + e_2 * a_2))
        * (alpha_2 - alpha_1)
        * delta_temperature
    )


def reference_stress_2(
    e_1: float,
    a_1: float,
    alpha_1: float,
    e_2: float,
    a_2: float,
    alpha_2: float,
    delta_temperature: float,
) -> float:
    """Independent closed form for the member 2 axial thermal stress [Pa]."""
    return (
        (e_1 * e_2 * a_1 / (e_1 * a_1 + e_2 * a_2))
        * (alpha_1 - alpha_2)
        * delta_temperature
    )


def reference_stresses(member_1, member_2, delta_temperature: float) -> tuple[float, float]:
    """Both reference stresses for a pair of ``AxialMember`` instances [Pa]."""
    args = (
        member_1.material.elastic_modulus,
        member_1.area,
        member_1.material.thermal_expansion_coefficient,
        member_2.material.elastic_modulus,
        member_2.area,
        member_2.material.thermal_expansion_coefficient,
        delta_temperature,
    )
    return reference_stress_1(*args), reference_stress_2(*args)

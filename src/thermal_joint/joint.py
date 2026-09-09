"""Closed-form axial thermal stress in a perfectly bonded two-member joint.

Model assumptions (Milestone 1)
-------------------------------
* Two members, perfectly bonded, forced to the same total axial strain.
* Uniform temperature in both members, common reference state ``T_ref``.
* No externally applied axial force on the joint.
* Linear elasticity, no bending, no interface slip, no adhesive compliance,
  no thermal gradient, temperature-independent properties.

Governing equations
-------------------
Compatibility (perfect bond)::

    eps_1 = eps_2 = eps_common

Constitutive relation for member i::

    sigma_i = E_i * (eps_common - alpha_i * dT)

Axial force equilibrium for a free joint::

    sigma_1 * A_1 + sigma_2 * A_2 = 0

Solving the three together gives the stiffness-weighted common strain::

    eps_common = (E_1 A_1 alpha_1 + E_2 A_2 alpha_2) / (E_1 A_1 + E_2 A_2) * dT

Units: stress [Pa], strain [-], force [N], area [m^2], dT [K or degC increment].
"""

from __future__ import annotations

from dataclasses import dataclass

from ._validation import require_finite
from .materials import AxialMember

__all__ = ["ThermalJointResult", "common_joint_strain", "solve_bimetallic_joint"]


def common_joint_strain(
    member_1: AxialMember, member_2: AxialMember, delta_temperature: float
) -> float:
    """Stiffness-weighted common axial strain of the bonded joint [-].

    ``eps_common = (E1 A1 a1 + E2 A2 a2) / (E1 A1 + E2 A2) * dT``

    This is the exact closed-form solution of compatibility plus axial force
    equilibrium; it is evaluated directly rather than solved numerically.
    """
    d_t = require_finite(delta_temperature, "delta_temperature")
    ea_1 = member_1.axial_stiffness
    ea_2 = member_2.axial_stiffness
    a_1 = member_1.material.thermal_expansion_coefficient
    a_2 = member_2.material.thermal_expansion_coefficient
    return (ea_1 * a_1 + ea_2 * a_2) / (ea_1 + ea_2) * d_t


@dataclass(frozen=True)
class ThermalJointResult:
    """Thermal state of a bonded two-member joint at one temperature change.

    Attributes
    ----------
    delta_temperature:
        ``T - T_ref`` [K or degC increment]; positive is heating.
    common_strain:
        Total axial strain shared by both members [-].
    member_1_free_thermal_strain, member_2_free_thermal_strain:
        Unrestrained thermal strains ``alpha_i * dT`` [-].
    member_1_stress, member_2_stress:
        Axial thermal stresses [Pa]; tension positive, compression negative.
    member_1_force, member_2_force:
        Internal axial forces ``sigma_i * A_i`` [N].
    force_equilibrium_residual:
        ``N_1 + N_2`` [N]. Zero to machine precision for a joint with no
        external axial load; reported rather than discarded so that the
        equilibrium check is visible to the caller.
    """

    delta_temperature: float
    common_strain: float
    member_1_free_thermal_strain: float
    member_2_free_thermal_strain: float
    member_1_stress: float
    member_2_stress: float
    member_1_force: float
    member_2_force: float
    force_equilibrium_residual: float

    @property
    def stresses(self) -> tuple[float, float]:
        """Member stresses as ``(sigma_1, sigma_2)`` [Pa]."""
        return (self.member_1_stress, self.member_2_stress)

    @property
    def forces(self) -> tuple[float, float]:
        """Member internal forces as ``(N_1, N_2)`` [N]."""
        return (self.member_1_force, self.member_2_force)


def solve_bimetallic_joint(
    member_1: AxialMember, member_2: AxialMember, delta_temperature: float
) -> ThermalJointResult:
    """Solve the bonded two-member joint for a uniform temperature change.

    Parameters
    ----------
    member_1, member_2:
        The two bonded :class:`~thermal_joint.materials.AxialMember` instances.
    delta_temperature:
        ``T - T_ref`` [K or degC increment]; positive is heating.

    Returns
    -------
    ThermalJointResult
        Common strain, free thermal strains, member stresses, internal forces
        and the axial force equilibrium residual.
    """
    for member, name in ((member_1, "member_1"), (member_2, "member_2")):
        if not isinstance(member, AxialMember):
            raise TypeError(
                f"{name} must be an AxialMember, got {type(member).__name__}"
            )

    d_t = require_finite(delta_temperature, "delta_temperature")

    eps_free_1 = member_1.free_thermal_strain(d_t)
    eps_free_2 = member_2.free_thermal_strain(d_t)
    eps_common = common_joint_strain(member_1, member_2, d_t)

    sigma_1 = member_1.material.elastic_modulus * (eps_common - eps_free_1)
    sigma_2 = member_2.material.elastic_modulus * (eps_common - eps_free_2)

    force_1 = sigma_1 * member_1.area
    force_2 = sigma_2 * member_2.area

    return ThermalJointResult(
        delta_temperature=d_t,
        common_strain=eps_common,
        member_1_free_thermal_strain=eps_free_1,
        member_2_free_thermal_strain=eps_free_2,
        member_1_stress=sigma_1,
        member_2_stress=sigma_2,
        member_1_force=force_1,
        member_2_force=force_2,
        force_equilibrium_residual=force_1 + force_2,
    )

"""Milestone 2: the bonded two-member joint under finite external axial restraint.

Milestone 1 solved the *free* joint, where the two members equilibrate only
against each other::

    N_1 + N_2 = 0

Milestone 2 generalises this to a joint attached to surrounding spacecraft
structure of finite axial stiffness, while preserving the Milestone 1 solution
exactly as the zero-restraint limit.

Restraint representation
------------------------
The surrounding structure is a single linear axial spring acting on the common
joint strain. Its stiffness ``K_r`` has units of **newtons** [N] - force per
unit strain, i.e. an equivalent ``E*A`` for the surrounding load path (``EA/L``
if a reference length is introduced). It is deliberately *not* a translational
stiffness in N/m, because no reference length is defined in this model.

The restraint is force-free when ``eps_common == reference_strain``.

Sign convention (used identically in the code, the tests and the README)
------------------------------------------------------------------------
Member forces are tension-positive, exactly as in Milestone 1. The restraint is
a parallel load path carrying its own tension-positive force::

    T_restraint_element = K_r * (eps_common - reference_strain)

and the force the restraint exerts **on the joint** is the reaction to that::

    N_r = -T_restraint_element = K_r * (reference_strain - eps_common)

``restraint_force`` in the result is ``N_r``, the force on the joint, positive
when the restraint acts on the joint in the tensile (positive axial) sense.
Joint equilibrium is then::

    N_1 + N_2 - N_r = 0

Closed form
-----------
Substituting ``sigma_i = E_i (eps - alpha_i dT)`` into the equilibrium statement
``E_1 A_1 (eps - a_1 dT) + E_2 A_2 (eps - a_2 dT) + K_r (eps - eps_ref) = 0``
gives::

                 E_1 A_1 alpha_1 dT + E_2 A_2 alpha_2 dT + K_r eps_ref
    eps_common = -----------------------------------------------------
                            E_1 A_1 + E_2 A_2 + K_r

The denominator is a **sum** of non-negative stiffnesses, so it is strictly
positive and the solution is unconditionally stable. A denominator of the form
``E_1 A_1 + E_2 A_2 - K_r`` would indicate an inconsistent force-direction
convention.

Limits (both enforced by tests):

* ``K_r = 0``   -> the exact Milestone 1 free-joint solution
* ``K_r -> inf`` -> ``eps_common -> reference_strain`` (fully restrained), with
  ``sigma_i -> E_i (eps_ref - alpha_i dT)``, which for ``eps_ref = 0`` is
  ``sigma_i -> -E_i alpha_i dT`` - both members driven to the *same* sign,
  unlike the self-equilibrating free joint.

Units: K_r [N], strains [-], stress [Pa], force [N], dT [K or degC increment].
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

from ._validation import require_finite, require_non_empty_name
from .joint import common_joint_strain
from .materials import AxialMember

__all__ = [
    "AxialRestraint",
    "RestrainedThermalJointResult",
    "joint_axial_stiffness",
    "restraint_stiffness_ratio",
    "rigid_restraint_stresses",
    "solve_restrained_joint",
    "zero_stress_restraint_stiffness",
]


def joint_axial_stiffness(member_1: AxialMember, member_2: AxialMember) -> float:
    """Combined axial stiffness of the bonded pair ``E_1 A_1 + E_2 A_2`` [N]."""
    return member_1.axial_stiffness + member_2.axial_stiffness


@dataclass(frozen=True)
class AxialRestraint:
    """A linear axial restraint representing the surrounding structure.

    Attributes
    ----------
    stiffness:
        Axial restraint stiffness ``K_r`` [N], i.e. force per unit strain. Must
        be finite and ``>= 0``; zero is allowed and recovers the free joint.
    reference_strain:
        The joint strain at which the restraint carries no force [-]. The
        canonical Milestone 2 case is ``0.0``, meaning the surrounding structure
        tends to hold the joint at its reference-temperature length. Nonzero
        values are supported for verification and future use.
    label:
        Optional non-empty label; defaults to ``"external axial restraint"``.
    """

    stiffness: float
    reference_strain: float = 0.0
    label: Optional[str] = None

    def __post_init__(self) -> None:
        stiffness = require_finite(self.stiffness, "stiffness")
        if stiffness < 0.0:
            raise ValueError(f"stiffness must be >= 0, got {stiffness!r}")
        object.__setattr__(self, "stiffness", stiffness)
        object.__setattr__(
            self,
            "reference_strain",
            require_finite(self.reference_strain, "reference_strain"),
        )
        if self.label is None:
            object.__setattr__(self, "label", "external axial restraint")
        else:
            object.__setattr__(self, "label", require_non_empty_name(self.label, "label"))

    def stiffness_ratio(self, member_1: AxialMember, member_2: AxialMember) -> float:
        """Dimensionless restraint parameter ``eta_r = K_r / (E_1 A_1 + E_2 A_2)`` [-]."""
        return self.stiffness / joint_axial_stiffness(member_1, member_2)

    @classmethod
    def from_stiffness_ratio(
        cls,
        stiffness_ratio: float,
        member_1: AxialMember,
        member_2: AxialMember,
        reference_strain: float = 0.0,
        label: Optional[str] = None,
    ) -> "AxialRestraint":
        """Build a restraint from ``eta_r`` rather than an absolute stiffness.

        ``K_r = eta_r * (E_1 A_1 + E_2 A_2)``.
        """
        ratio = require_finite(stiffness_ratio, "stiffness_ratio")
        if ratio < 0.0:
            raise ValueError(f"stiffness_ratio must be >= 0, got {ratio!r}")
        return cls(
            stiffness=ratio * joint_axial_stiffness(member_1, member_2),
            reference_strain=reference_strain,
            label=label,
        )


def restraint_stiffness_ratio(
    restraint: AxialRestraint, member_1: AxialMember, member_2: AxialMember
) -> float:
    """Dimensionless restraint parameter ``eta_r = K_r / (E_1 A_1 + E_2 A_2)`` [-].

    ``eta_r = 0`` is a free joint, ``eta_r ~ 1`` means the surrounding structure
    is about as axially stiff as the joint itself, and ``eta_r >> 1`` is a
    strongly restrained joint.
    """
    return restraint.stiffness_ratio(member_1, member_2)


def rigid_restraint_stresses(
    member_1: AxialMember,
    member_2: AxialMember,
    delta_temperature: float,
    reference_strain: float = 0.0,
) -> tuple[float, float]:
    """Fully restrained (``K_r -> inf``) member stresses ``(sigma_1, sigma_2)`` [Pa].

    ``sigma_i = E_i * (eps_ref - alpha_i * dT)``, which for ``eps_ref = 0``
    reduces to the classic fully restrained result ``sigma_i = -E_i alpha_i dT``.
    Both members take the same sign, unlike the free bimetallic joint.
    """
    d_t = require_finite(delta_temperature, "delta_temperature")
    eps_ref = require_finite(reference_strain, "reference_strain")
    return (
        member_1.material.elastic_modulus * (eps_ref - member_1.free_thermal_strain(d_t)),
        member_2.material.elastic_modulus * (eps_ref - member_2.free_thermal_strain(d_t)),
    )


@dataclass(frozen=True)
class RestrainedThermalJointResult:
    """State of the bonded joint under external restraint at one temperature change.

    Attributes
    ----------
    delta_temperature:
        ``T - T_ref`` [K or degC increment]; positive is heating.
    restraint_stiffness:
        ``K_r`` [N] used for this solution.
    restraint_reference_strain:
        The force-free strain of the restraint [-].
    common_strain:
        Total axial strain shared by both members [-].
    member_1_free_thermal_strain, member_2_free_thermal_strain:
        Unrestrained thermal strains ``alpha_i * dT`` [-].
    member_1_stress, member_2_stress:
        Axial thermal stresses [Pa]; tension positive.
    member_1_force, member_2_force:
        Internal axial member forces ``sigma_i A_i`` [N]; tension positive.
    restraint_force:
        ``N_r = K_r (eps_ref - eps_common)`` [N] - the force the restraint
        exerts **on the joint**, positive in the tensile (positive axial) sense.
    force_equilibrium_residual:
        ``N_1 + N_2 - N_r`` [N], assembled from the signed forces actually used
        by the model. Zero to machine precision; reported, never discarded.
    free_joint_common_strain:
        The Milestone 1 free-joint common strain [-] for the same members and
        ``dT``, provided for direct comparison. No ratio against it is exposed,
        because that ratio is undefined whenever the free strain is zero.
    """

    delta_temperature: float
    restraint_stiffness: float
    restraint_reference_strain: float
    common_strain: float
    member_1_free_thermal_strain: float
    member_2_free_thermal_strain: float
    member_1_stress: float
    member_2_stress: float
    member_1_force: float
    member_2_force: float
    restraint_force: float
    force_equilibrium_residual: float
    free_joint_common_strain: float

    @property
    def stresses(self) -> tuple[float, float]:
        """Member stresses as ``(sigma_1, sigma_2)`` [Pa]."""
        return (self.member_1_stress, self.member_2_stress)

    @property
    def forces(self) -> tuple[float, float]:
        """Member internal forces as ``(N_1, N_2)`` [N], excluding the restraint."""
        return (self.member_1_force, self.member_2_force)

    @property
    def members_carry_same_sign_stress(self) -> bool:
        """True when both members are in tension or both in compression.

        The free bimetallic joint is always self-equilibrating (opposite signs);
        a strongly restrained joint drives both members to the same sign.
        """
        return self.member_1_stress * self.member_2_stress > 0.0


def solve_restrained_joint(
    member_1: AxialMember,
    member_2: AxialMember,
    delta_temperature: float,
    restraint: AxialRestraint | None = None,
) -> RestrainedThermalJointResult:
    """Solve the bonded joint under finite external axial restraint.

    Parameters
    ----------
    member_1, member_2:
        The two bonded members.
    delta_temperature:
        ``T - T_ref`` [K or degC increment].
    restraint:
        The external restraint; defaults to ``AxialRestraint(0.0)``, which
        reproduces the Milestone 1 free-joint solution exactly.

    Returns
    -------
    RestrainedThermalJointResult
    """
    for member, name in ((member_1, "member_1"), (member_2, "member_2")):
        if not isinstance(member, AxialMember):
            raise TypeError(f"{name} must be an AxialMember, got {type(member).__name__}")

    external = AxialRestraint(0.0) if restraint is None else restraint
    if not isinstance(external, AxialRestraint):
        raise TypeError(
            f"restraint must be an AxialRestraint, got {type(external).__name__}"
        )

    d_t = require_finite(delta_temperature, "delta_temperature")

    eps_free_1 = member_1.free_thermal_strain(d_t)
    eps_free_2 = member_2.free_thermal_strain(d_t)

    ea_1 = member_1.axial_stiffness
    ea_2 = member_2.axial_stiffness
    k_r = external.stiffness
    eps_ref = external.reference_strain

    eps_common = (ea_1 * eps_free_1 + ea_2 * eps_free_2 + k_r * eps_ref) / (
        ea_1 + ea_2 + k_r
    )

    sigma_1 = member_1.material.elastic_modulus * (eps_common - eps_free_1)
    sigma_2 = member_2.material.elastic_modulus * (eps_common - eps_free_2)
    force_1 = sigma_1 * member_1.area
    force_2 = sigma_2 * member_2.area
    restraint_force = k_r * (eps_ref - eps_common)

    return RestrainedThermalJointResult(
        delta_temperature=d_t,
        restraint_stiffness=k_r,
        restraint_reference_strain=eps_ref,
        common_strain=eps_common,
        member_1_free_thermal_strain=eps_free_1,
        member_2_free_thermal_strain=eps_free_2,
        member_1_stress=sigma_1,
        member_2_stress=sigma_2,
        member_1_force=force_1,
        member_2_force=force_2,
        restraint_force=restraint_force,
        force_equilibrium_residual=force_1 + force_2 - restraint_force,
        free_joint_common_strain=common_joint_strain(member_1, member_2, d_t),
    )


def zero_stress_restraint_stiffness(
    member_1: AxialMember,
    member_2: AxialMember,
    delta_temperature: float,
    member_index: int,
    reference_strain: float = 0.0,
) -> Optional[float]:
    """Restraint stiffness at which one member's stress crosses through zero.

    Member ``i`` is unstressed when ``eps_common == alpha_i * dT``. Substituting
    the restrained closed form and solving for ``K_r``::

        K_r = dT * (alpha_i * K_joint - S) / (eps_ref - alpha_i * dT)

    with ``K_joint = E_1 A_1 + E_2 A_2`` and
    ``S = E_1 A_1 alpha_1 + E_2 A_2 alpha_2``.

    For ``eps_ref = 0`` this collapses to the dT-independent
    ``K_r = S / alpha_i - K_joint``.

    Parameters
    ----------
    member_index:
        ``1`` or ``2``, selecting which member is driven to zero stress.

    Returns
    -------
    float or None
        The crossing stiffness [N] when a strictly positive, finite one exists
        (a physically reachable restraint), otherwise ``None``. ``None`` means
        this member never passes through zero stress for any admissible
        ``K_r >= 0`` - for example the high-CTE member on heating with
        ``eps_ref = 0``, which is compressive both free and fully restrained.
    """
    if member_index not in (1, 2):
        raise ValueError(f"member_index must be 1 or 2, got {member_index!r}")

    d_t = require_finite(delta_temperature, "delta_temperature")
    eps_ref = require_finite(reference_strain, "reference_strain")

    target = member_1 if member_index == 1 else member_2
    alpha = target.material.thermal_expansion_coefficient

    k_joint = joint_axial_stiffness(member_1, member_2)
    s_term = (
        member_1.axial_stiffness * member_1.material.thermal_expansion_coefficient
        + member_2.axial_stiffness * member_2.material.thermal_expansion_coefficient
    )

    denominator = eps_ref - alpha * d_t
    if denominator == 0.0:
        return None

    stiffness = d_t * (alpha * k_joint - s_term) / denominator
    if not math.isfinite(stiffness) or stiffness <= 0.0:
        return None
    return stiffness

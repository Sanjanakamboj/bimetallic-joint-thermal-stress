"""Preliminary elastic yield margins for the bonded joint.

These are *preliminary elastic yield margins*, not certification margins: they
compare a linear-elastic thermal stress against a yield allowable derived from
an explicit design factor. No knockdowns, statistical basis, joint efficiency,
fatigue or plasticity effects are included.

Convention::

    sigma_allow_i = yield_strength_i / design_factor
    MS_yield_i    = sigma_allow_i / abs(sigma_i) - 1

A margin of exactly zero passes (``MS >= 0``). A member carrying zero stress
has an infinite, non-governing margin (``math.inf``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from ._validation import require_finite
from .materials import AxialMember, ThermoelasticMaterial

__all__ = [
    "YieldBasis",
    "MemberYieldMargin",
    "JointYieldAssessment",
    "MemberStressState",
    "assess_yield",
]


@runtime_checkable
class MemberStressState(Protocol):
    """Anything carrying a signed axial stress for each of the two members.

    Both :class:`~thermal_joint.joint.ThermalJointResult` (Milestone 1, free
    joint) and
    :class:`~thermal_joint.restraint.RestrainedThermalJointResult`
    (Milestone 2, externally restrained) satisfy this, so the yield-margin
    machinery is written once and reused unchanged by both.
    """

    member_1_stress: float
    member_2_stress: float


@dataclass(frozen=True)
class YieldBasis:
    """Explicit design basis for yield allowables.

    Attributes
    ----------
    design_factor:
        Factor applied to the material yield strength to obtain the allowable
        stress. Must be finite and ``>= 1``. The default of 1.0 means the
        allowable equals the raw yield strength.

    Notes
    -----
    The design factor lives here, never inside
    :class:`~thermal_joint.materials.ThermoelasticMaterial`, so that material
    data and design policy stay separable.
    """

    design_factor: float = 1.0

    def __post_init__(self) -> None:
        factor = require_finite(self.design_factor, "design_factor")
        if factor < 1.0:
            raise ValueError(f"design_factor must be >= 1, got {factor!r}")
        object.__setattr__(self, "design_factor", factor)

    def allowable_stress(self, material: ThermoelasticMaterial) -> float:
        """Allowable stress magnitude ``yield_strength / design_factor`` [Pa]."""
        if not isinstance(material, ThermoelasticMaterial):
            raise TypeError(
                "material must be a ThermoelasticMaterial, got "
                f"{type(material).__name__}"
            )
        return material.yield_strength / self.design_factor

    def margin_of_safety(self, material: ThermoelasticMaterial, stress: float) -> float:
        """Preliminary elastic yield margin of safety [-].

        ``MS = sigma_allow / abs(sigma) - 1``; returns ``math.inf`` for zero
        stress, which is non-governing by construction.
        """
        sigma = require_finite(stress, "stress")
        magnitude = abs(sigma)
        if magnitude == 0.0:
            return math.inf
        return self.allowable_stress(material) / magnitude - 1.0


@dataclass(frozen=True)
class MemberYieldMargin:
    """Yield margin bookkeeping for a single member.

    Attributes
    ----------
    label:
        Member label.
    stress:
        Signed axial stress [Pa].
    allowable_stress:
        Allowable stress magnitude [Pa].
    margin_of_safety:
        ``sigma_allow / abs(sigma) - 1`` [-]; ``math.inf`` when stress is zero.
    """

    label: str
    stress: float
    allowable_stress: float
    margin_of_safety: float

    @property
    def passes(self) -> bool:
        """True when ``MS >= 0``; the boundary ``MS == 0`` passes."""
        return self.margin_of_safety >= 0.0


@dataclass(frozen=True)
class JointYieldAssessment:
    """Yield assessment of both members at one temperature condition.

    Attributes
    ----------
    member_1_margin, member_2_margin:
        Per-member :class:`MemberYieldMargin` records.
    governing_member:
        ``1`` or ``2`` - the member with the smaller margin. Computed, never
        assumed from which material has the lower yield strength. An exact tie
        resolves to member 1.
    minimum_margin:
        The smaller of the two margins [-].
    """

    member_1_margin: MemberYieldMargin
    member_2_margin: MemberYieldMargin
    governing_member: int
    minimum_margin: float

    @property
    def governing_margin(self) -> MemberYieldMargin:
        """The :class:`MemberYieldMargin` record of the governing member."""
        return self.member_1_margin if self.governing_member == 1 else self.member_2_margin

    @property
    def passes(self) -> bool:
        """True when both members have ``MS >= 0``."""
        return self.minimum_margin >= 0.0


def assess_yield(
    result: MemberStressState,
    member_1: AxialMember,
    member_2: AxialMember,
    basis: YieldBasis | None = None,
) -> JointYieldAssessment:
    """Compute preliminary elastic yield margins for both members.

    Parameters
    ----------
    result:
        Any solved result exposing ``member_1_stress`` and ``member_2_stress``:
        a Milestone 1 :class:`~thermal_joint.joint.ThermalJointResult` or a
        Milestone 2
        :class:`~thermal_joint.restraint.RestrainedThermalJointResult`.
    member_1, member_2:
        The same members used to produce ``result``, in the same order.
    basis:
        Design basis; defaults to ``YieldBasis()`` (design factor 1.0).
    """
    design_basis = YieldBasis() if basis is None else basis
    if not isinstance(design_basis, YieldBasis):
        raise TypeError(
            f"basis must be a YieldBasis, got {type(design_basis).__name__}"
        )

    margins = []
    for member, stress in (
        (member_1, result.member_1_stress),
        (member_2, result.member_2_stress),
    ):
        margins.append(
            MemberYieldMargin(
                label=member.label,
                stress=stress,
                allowable_stress=design_basis.allowable_stress(member.material),
                margin_of_safety=design_basis.margin_of_safety(member.material, stress),
            )
        )

    margin_1, margin_2 = margins
    governing_member = 2 if margin_2.margin_of_safety < margin_1.margin_of_safety else 1
    minimum_margin = min(margin_1.margin_of_safety, margin_2.margin_of_safety)

    return JointYieldAssessment(
        member_1_margin=margin_1,
        member_2_margin=margin_2,
        governing_member=governing_member,
        minimum_margin=minimum_margin,
    )

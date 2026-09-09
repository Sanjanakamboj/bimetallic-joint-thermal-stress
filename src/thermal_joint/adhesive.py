"""Adhesive material, bonded-overlap geometry and adhesive shear design basis.

Milestone 3 adds a *local* load-transfer layer on top of the Milestone 1/2
*global* models. Nothing here changes the global mechanics: the adhesive data,
the overlap geometry and the shear allowable are new, additive objects.

Units are SI throughout: shear modulus and shear strength [Pa], lengths [m],
bond area [m^2].
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

from ._validation import require_finite, require_non_empty_name, require_positive

__all__ = ["AdhesiveMaterial", "BondedOverlapGeometry", "AdhesiveShearBasis"]


@dataclass(frozen=True)
class AdhesiveMaterial:
    """A linear-elastic adhesive characterised by shear stiffness and strength.

    Attributes
    ----------
    name:
        Non-empty identifier. No commercial adhesive is named unless the
        properties are genuinely sourced.
    shear_modulus:
        Adhesive shear modulus ``G_a`` [Pa]. Finite and strictly positive.
    shear_strength:
        Adhesive shear strength magnitude [Pa]. Finite and strictly positive.
        Stored as a positive magnitude with no embedded design factor - the
        design factor lives in :class:`AdhesiveShearBasis`.
    source_note:
        Required non-empty provenance string. For unsourced values this must
        say ``ILLUSTRATIVE ADHESIVE INPUT - NOT DESIGN ALLOWABLE`` or similar;
        provenance is mandatory so a screening number can never be mistaken for
        a qualified allowable.
    notes:
        Optional free-text remarks.
    """

    name: str
    shear_modulus: float
    shear_strength: float
    source_note: str
    notes: Optional[str] = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", require_non_empty_name(self.name))
        object.__setattr__(
            self, "shear_modulus", require_positive(self.shear_modulus, "shear_modulus")
        )
        object.__setattr__(
            self, "shear_strength", require_positive(self.shear_strength, "shear_strength")
        )
        object.__setattr__(
            self, "source_note", require_non_empty_name(self.source_note, "source_note")
        )
        if self.notes is not None:
            object.__setattr__(self, "notes", require_non_empty_name(self.notes, "notes"))


@dataclass(frozen=True)
class BondedOverlapGeometry:
    """Geometry of the bonded overlap.

    Attributes
    ----------
    overlap_length:
        Bonded overlap length ``L_b`` [m], measured along the load path.
    bond_width:
        Bond width ``b`` [m], transverse to the load path.
    adhesive_thickness:
        Adhesive layer thickness ``t_a`` [m].

    Notes
    -----
    ``bond_area = bond_width * overlap_length`` is the *bonded interface* area.
    It is deliberately distinct from the adherends' cross-sectional areas, which
    live on :class:`~thermal_joint.materials.AxialMember`.
    """

    overlap_length: float
    bond_width: float
    adhesive_thickness: float

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "overlap_length", require_positive(self.overlap_length, "overlap_length")
        )
        object.__setattr__(
            self, "bond_width", require_positive(self.bond_width, "bond_width")
        )
        object.__setattr__(
            self,
            "adhesive_thickness",
            require_positive(self.adhesive_thickness, "adhesive_thickness"),
        )

    @property
    def bond_area(self) -> float:
        """Bonded interface area ``b * L_b`` [m^2]."""
        return self.bond_width * self.overlap_length

    def with_overlap_length(self, overlap_length: float) -> "BondedOverlapGeometry":
        """Return a copy with a different overlap length; the original is untouched."""
        return BondedOverlapGeometry(
            overlap_length=overlap_length,
            bond_width=self.bond_width,
            adhesive_thickness=self.adhesive_thickness,
        )

    def with_adhesive_thickness(self, adhesive_thickness: float) -> "BondedOverlapGeometry":
        """Return a copy with a different adhesive thickness."""
        return BondedOverlapGeometry(
            overlap_length=self.overlap_length,
            bond_width=self.bond_width,
            adhesive_thickness=adhesive_thickness,
        )


@dataclass(frozen=True)
class AdhesiveShearBasis:
    """Explicit design basis for the adhesive shear allowable.

    ``tau_allow = shear_strength / design_factor`` and

    ``MS_adh = tau_allow / abs(tau_peak) - 1``

    This is a **preliminary adhesive shear margin**, not a certification margin.
    A margin of exactly zero passes; zero demand returns ``math.inf``, an
    explicitly non-governing value rather than a division by zero.
    """

    design_factor: float = 1.0

    def __post_init__(self) -> None:
        factor = require_finite(self.design_factor, "design_factor")
        if factor < 1.0:
            raise ValueError(f"design_factor must be >= 1, got {factor!r}")
        object.__setattr__(self, "design_factor", factor)

    def allowable_shear_stress(self, adhesive: AdhesiveMaterial) -> float:
        """Allowable adhesive shear stress ``shear_strength / design_factor`` [Pa]."""
        if not isinstance(adhesive, AdhesiveMaterial):
            raise TypeError(
                f"adhesive must be an AdhesiveMaterial, got {type(adhesive).__name__}"
            )
        return adhesive.shear_strength / self.design_factor

    def margin_of_safety(self, adhesive: AdhesiveMaterial, shear_stress: float) -> float:
        """Preliminary adhesive shear margin of safety [-]."""
        demand = abs(require_finite(shear_stress, "shear_stress"))
        if demand == 0.0:
            return math.inf
        return self.allowable_shear_stress(adhesive) / demand - 1.0

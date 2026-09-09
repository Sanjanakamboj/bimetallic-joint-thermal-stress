"""Isotropic thermoelastic material and axial member representations.

Sign conventions used throughout this package
---------------------------------------------
* Tensile stress is positive; compressive stress is negative.
* ``delta_temperature = T - T_ref``; positive means heating.
* A positive coefficient of thermal expansion means the material expands
  on heating, so the free thermal strain ``alpha * delta_T`` is positive.
* ``yield_strength`` is stored as a positive magnitude and is compared
  against ``abs(stress)``.

Temperature *differences* may be supplied in K or in degrees Celsius; the two
increments are numerically identical, so no conversion is performed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from ._validation import require_finite, require_non_empty_name, require_positive

__all__ = ["ThermoelasticMaterial", "AxialMember"]


@dataclass(frozen=True)
class ThermoelasticMaterial:
    """A linear-elastic, isotropic material with a single thermal expansion coefficient.

    Attributes
    ----------
    name:
        Non-empty identifier for the material.
    elastic_modulus:
        Young's modulus E [Pa]. Must be finite and strictly positive.
    thermal_expansion_coefficient:
        Coefficient of thermal expansion alpha [1/K]. Must be finite. Negative
        values are permitted (some materials contract on heating); zero is
        permitted (a notionally athermal member).
    yield_strength:
        Yield strength magnitude [Pa]. Must be finite and strictly positive.
        Stored as a positive magnitude with no embedded design or safety
        factor - design factors belong to :class:`~thermal_joint.margins.YieldBasis`.
    density:
        Optional mass density [kg/m^3]. If given it must be finite and positive.
        Unused by the Milestone 1 mechanics; carried for later mass bookkeeping.
    """

    name: str
    elastic_modulus: float
    thermal_expansion_coefficient: float
    yield_strength: float
    density: Optional[float] = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", require_non_empty_name(self.name))
        object.__setattr__(
            self,
            "elastic_modulus",
            require_positive(self.elastic_modulus, "elastic_modulus"),
        )
        object.__setattr__(
            self,
            "thermal_expansion_coefficient",
            require_finite(
                self.thermal_expansion_coefficient, "thermal_expansion_coefficient"
            ),
        )
        object.__setattr__(
            self,
            "yield_strength",
            require_positive(self.yield_strength, "yield_strength"),
        )
        if self.density is not None:
            object.__setattr__(self, "density", require_positive(self.density, "density"))

    def free_thermal_strain(self, delta_temperature: float) -> float:
        """Unrestrained thermal strain ``alpha * delta_T`` [-].

        Parameters
        ----------
        delta_temperature:
            ``T - T_ref`` [K or degC increment].
        """
        d_t = require_finite(delta_temperature, "delta_temperature")
        return self.thermal_expansion_coefficient * d_t


@dataclass(frozen=True)
class AxialMember:
    """One axially loaded member of the joint.

    Attributes
    ----------
    material:
        The member's :class:`ThermoelasticMaterial`.
    area:
        Cross-sectional area A [m^2]. Must be finite and strictly positive.
    label:
        Optional human-readable label. Defaults to the material name.

    Notes
    -----
    Both members of a Milestone 1 joint share the same reference length, so
    length cancels out of the compatibility equation and is deliberately not
    modelled here.
    """

    material: ThermoelasticMaterial
    area: float
    label: Optional[str] = None

    def __post_init__(self) -> None:
        if not isinstance(self.material, ThermoelasticMaterial):
            raise TypeError(
                "material must be a ThermoelasticMaterial, got "
                f"{type(self.material).__name__}"
            )
        object.__setattr__(self, "area", require_positive(self.area, "area"))
        if self.label is None:
            object.__setattr__(self, "label", self.material.name)
        else:
            object.__setattr__(self, "label", require_non_empty_name(self.label, "label"))

    @property
    def axial_stiffness(self) -> float:
        """Axial stiffness product ``E * A`` [N] (per unit strain)."""
        return self.material.elastic_modulus * self.area

    def free_thermal_strain(self, delta_temperature: float) -> float:
        """Unrestrained thermal strain of this member [-]."""
        return self.material.free_thermal_strain(delta_temperature)

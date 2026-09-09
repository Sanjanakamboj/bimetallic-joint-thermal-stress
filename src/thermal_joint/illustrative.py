"""Illustrative material and joint inputs for the Milestone 1 sanity study.

ILLUSTRATIVE MATERIAL INPUT - NOT DESIGN ALLOWABLE.

Every value below is a clean, round, order-of-magnitude figure chosen to make
the mechanics legible. None of it is traceable to a specific alloy, temper,
product form, statistical basis or qualification programme, so no alloy
designation is claimed. Do not use these numbers for design.
"""

from __future__ import annotations

from .environment import ThermalEnvironment
from .materials import AxialMember, ThermoelasticMaterial

__all__ = [
    "ILLUSTRATIVE_NOTE",
    "ALUMINIUM_LIKE",
    "TITANIUM_LIKE",
    "SQUARE_MILLIMETRE",
    "aluminium_like_member",
    "titanium_like_member",
    "radiator_joint_environment",
]

ILLUSTRATIVE_NOTE = "ILLUSTRATIVE MATERIAL INPUT - NOT DESIGN ALLOWABLE"

SQUARE_MILLIMETRE = 1.0e-6
"""Conversion factor: 1 mm^2 in m^2."""

ALUMINIUM_LIKE = ThermoelasticMaterial(
    name="Aluminium-like (illustrative)",
    elastic_modulus=70.0e9,
    thermal_expansion_coefficient=23.0e-6,
    yield_strength=270.0e6,
    density=2700.0,
)
"""High-CTE, low-modulus member. ILLUSTRATIVE MATERIAL INPUT - NOT DESIGN ALLOWABLE."""

TITANIUM_LIKE = ThermoelasticMaterial(
    name="Titanium-like (illustrative)",
    elastic_modulus=110.0e9,
    thermal_expansion_coefficient=8.5e-6,
    yield_strength=830.0e6,
    density=4430.0,
)
"""Low-CTE, high-modulus member. ILLUSTRATIVE MATERIAL INPUT - NOT DESIGN ALLOWABLE."""


def aluminium_like_member(area_mm2: float = 100.0) -> AxialMember:
    """Aluminium-like member with cross-sectional area given in mm^2."""
    return AxialMember(
        material=ALUMINIUM_LIKE,
        area=area_mm2 * SQUARE_MILLIMETRE,
        label="Member 1 - aluminium-like",
    )


def titanium_like_member(area_mm2: float = 100.0) -> AxialMember:
    """Titanium-like member with cross-sectional area given in mm^2."""
    return AxialMember(
        material=TITANIUM_LIKE,
        area=area_mm2 * SQUARE_MILLIMETRE,
        label="Member 2 - titanium-like",
    )


def radiator_joint_environment() -> ThermalEnvironment:
    """Illustrative radiator-joint excursion in degrees Celsius.

    ``T_ref = +20 degC``, ``T_cold = -120 degC``, ``T_hot = +120 degC``, giving
    ``dT_cold = -140 K`` and ``dT_hot = +100 K``. These are illustrative
    excursions for a portfolio study, not spacecraft qualification temperatures.
    """
    return ThermalEnvironment(
        reference_temperature=20.0,
        cold_temperature=-120.0,
        hot_temperature=120.0,
    )

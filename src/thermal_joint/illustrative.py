"""Illustrative material and joint inputs for the Milestone 1 sanity study.

ILLUSTRATIVE MATERIAL INPUT - NOT DESIGN ALLOWABLE.

Every value below is a clean, round, order-of-magnitude figure chosen to make
the mechanics legible. None of it is traceable to a specific alloy, temper,
product form, statistical basis or qualification programme, so no alloy
designation is claimed. Do not use these numbers for design.
"""

from __future__ import annotations

from .adhesive import AdhesiveMaterial, BondedOverlapGeometry
from .environment import ThermalEnvironment
from .materials import AxialMember, ThermoelasticMaterial

__all__ = [
    "ILLUSTRATIVE_NOTE",
    "ILLUSTRATIVE_ADHESIVE_NOTE",
    "ALUMINIUM_LIKE",
    "TITANIUM_LIKE",
    "ADHESIVE_LIKE",
    "SQUARE_MILLIMETRE",
    "MILLIMETRE",
    "aluminium_like_member",
    "titanium_like_member",
    "radiator_joint_environment",
    "radiator_joint_overlap",
]

ILLUSTRATIVE_NOTE = "ILLUSTRATIVE MATERIAL INPUT - NOT DESIGN ALLOWABLE"

SQUARE_MILLIMETRE = 1.0e-6
"""Conversion factor: 1 mm^2 in m^2."""

MILLIMETRE = 1.0e-3
"""Conversion factor: 1 mm in m."""

ILLUSTRATIVE_ADHESIVE_NOTE = "ILLUSTRATIVE ADHESIVE-LIKE INPUT - NOT DESIGN ALLOWABLE"

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


ADHESIVE_LIKE = AdhesiveMaterial(
    name="Structural-adhesive-like (illustrative)",
    shear_modulus=1.0e9,
    shear_strength=25.0e6,
    source_note=ILLUSTRATIVE_ADHESIVE_NOTE,
    notes=(
        "Clean mid-range values for a stiff structural film adhesive at room "
        "temperature. Not traceable to any product, cure state, bondline "
        "preparation or test programme, so no commercial adhesive is named. "
        "Chosen before the resulting shear stresses were computed; not tuned "
        "to produce a pass or a failure."
    ),
)
"""Illustrative adhesive. ILLUSTRATIVE ADHESIVE-LIKE INPUT - NOT DESIGN ALLOWABLE."""


def radiator_joint_overlap(
    overlap_length_mm: float = 40.0,
    bond_width_mm: float = 20.0,
    adhesive_thickness_mm: float = 0.2,
) -> BondedOverlapGeometry:
    """Illustrative bonded overlap for the radiator joint, dimensions in mm.

    The default 40 x 20 mm overlap with a 0.2 mm bondline spans about 6.1
    transfer lengths for the illustrative adhesive and members, which is a
    meaningful but not pathological overlap: long enough for the interior to be
    nearly shear-free, short enough that the end effect is still visible.
    """
    return BondedOverlapGeometry(
        overlap_length=overlap_length_mm * MILLIMETRE,
        bond_width=bond_width_mm * MILLIMETRE,
        adhesive_thickness=adhesive_thickness_mm * MILLIMETRE,
    )

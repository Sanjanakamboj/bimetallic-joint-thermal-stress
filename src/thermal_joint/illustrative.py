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
from .fatigue import BasquinFatigueCurve, ThermalCycleRequirement
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
    "ILLUSTRATIVE_FATIGUE_NOTE",
    "ALUMINIUM_LIKE_FATIGUE",
    "TITANIUM_LIKE_FATIGUE",
    "ADHESIVE_SHEAR_FATIGUE",
    "illustrative_fatigue_curves",
    "illustrative_cycle_requirement",
]

ILLUSTRATIVE_NOTE = "ILLUSTRATIVE MATERIAL INPUT - NOT DESIGN ALLOWABLE"

SQUARE_MILLIMETRE = 1.0e-6
"""Conversion factor: 1 mm^2 in m^2."""

MILLIMETRE = 1.0e-3
"""Conversion factor: 1 mm in m."""

ILLUSTRATIVE_ADHESIVE_NOTE = "ILLUSTRATIVE ADHESIVE-LIKE INPUT - NOT DESIGN ALLOWABLE"

ILLUSTRATIVE_FATIGUE_NOTE = "ILLUSTRATIVE FATIGUE INPUT - NOT DESIGN ALLOWABLE"

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


# ---------------------------------------------------------------------------
# Milestone 5: illustrative Basquin fatigue curves, sigma_a = A * N^b, N in cycles.
#
# ILLUSTRATIVE FATIGUE INPUT - NOT DESIGN ALLOWABLE.
#
# Round order-of-magnitude figures chosen to sit in a plausible range for the
# respective material classes, selected after inspecting the canonical
# alternating stresses and before the resulting lives were used for any
# conclusion. None is traceable to an alloy, temper, adhesive product, surface
# condition, environment or test programme, so no such designation is claimed
# and no curve is tuned to force a pass or a failure.
#
# A Basquin power law has no endurance limit and no low-cycle cut-off, so each
# curve states the cycle range over which it is meant to be read. Nothing is
# clamped: silently clamping would hide an extrapolation rather than flag it.
# ---------------------------------------------------------------------------

_FATIGUE_RANGE_NOTE = (
    "Intended to be read over roughly 1e3 to 1e8 cycles. A Basquin fit has no "
    "endurance limit and no low-cycle cut-off, so lives extrapolated outside "
    "that range are not meaningful."
)

ALUMINIUM_LIKE_FATIGUE = BasquinFatigueCurve(
    name="Aluminium-like S-N (illustrative)",
    coefficient_A=900.0e6,
    exponent_b=-0.12,
    source_note=ILLUSTRATIVE_FATIGUE_NOTE,
    notes=(
        "Axial constant-amplitude alternating stress. Gives about 298 MPa at "
        "1e4 cycles and 130 MPa at 1e7. " + _FATIGUE_RANGE_NOTE
    ),
)
"""Illustrative Al-like S-N curve. ILLUSTRATIVE FATIGUE INPUT - NOT DESIGN ALLOWABLE."""

TITANIUM_LIKE_FATIGUE = BasquinFatigueCurve(
    name="Titanium-like S-N (illustrative)",
    coefficient_A=2000.0e6,
    exponent_b=-0.10,
    source_note=ILLUSTRATIVE_FATIGUE_NOTE,
    notes=(
        "Axial constant-amplitude alternating stress; a shallower exponent than "
        "the aluminium-like curve, as titanium alloys generally are. About "
        "796 MPa at 1e4 cycles and 399 MPa at 1e7. " + _FATIGUE_RANGE_NOTE
    ),
)
"""Illustrative Ti-like S-N curve. ILLUSTRATIVE FATIGUE INPUT - NOT DESIGN ALLOWABLE."""

ADHESIVE_SHEAR_FATIGUE = BasquinFatigueCurve(
    name="Structural-adhesive-like shear S-N (illustrative)",
    coefficient_A=60.0e6,
    exponent_b=-0.15,
    source_note=ILLUSTRATIVE_FATIGUE_NOTE,
    notes=(
        "Alternating ADHESIVE SHEAR stress, used directly - no shear-to-von-Mises "
        "conversion is applied anywhere in Milestone 5. The steeper exponent "
        "reflects the greater cyclic sensitivity typical of polymeric adhesives. "
        "About 15.1 MPa at 1e4 cycles and 7.6 MPa at 1e6. " + _FATIGUE_RANGE_NOTE
    ),
)
"""Illustrative adhesive shear S-N curve. ILLUSTRATIVE FATIGUE INPUT - NOT DESIGN ALLOWABLE."""


def illustrative_fatigue_curves():
    """The three illustrative curves as a
    :class:`~thermal_joint.fatigue_assessment.FatigueCurveSet`."""
    from .fatigue_assessment import FatigueCurveSet

    return FatigueCurveSet(
        member_1=ALUMINIUM_LIKE_FATIGUE,
        member_2=TITANIUM_LIKE_FATIGUE,
        adhesive_shear=ADHESIVE_SHEAR_FATIGUE,
    )


def illustrative_cycle_requirement(required_cycles: float = 1.0e4):
    """Illustrative thermal-cycle requirement, default 1e4 cycles.

    1e4 is a round figure of the order of a couple of years of low-Earth-orbit
    thermal cycling. It is an **illustrative study input, not a spacecraft
    qualification requirement**, and it was picked from a short list of round
    values (1e3, 1e4, 1e5, 1e6) rather than tuned: the accompanying study
    reports the full requirement sensitivity across all four so that no single
    choice drives the conclusion.
    """
    return ThermalCycleRequirement(
        required_cycles=required_cycles,
        label="illustrative thermal-cycle requirement - not a qualification requirement",
    )

r"""Milestone 4: a bounded preliminary design trade over the verified mechanics.

Milestone 3 ended on a negative result: the canonical 40 mm overlap fails
adhesive shear, and lengthening the overlap cannot fix it because peak shear
converges on a finite asymptote that already exceeds the allowable. This module
asks which *other* first-order levers can recover feasibility, using nothing but
the Milestone 1 and Milestone 3 equations.

**Milestone 4 does not optimize a flight joint.** It uses the verified
one-dimensional shear-lag equations to identify feasible first-order
combinations of bond width and adhesive compliance within explicitly bounded
design spaces.

Analytical scaling
------------------
Substituting the Milestone 1 demand ``N_t = d_eps / C`` (with
``C = 1/(E_1 A_1) + 1/(E_2 A_2)`` and ``d_eps = (alpha_1 - alpha_2) dT``) into
the Milestone 3 long-overlap asymptote ``tau_inf = |N_t| beta / b``, with
``beta = sqrt(G_a b C / t_a)``::

    tau_inf = |d_eps| * sqrt( G_a / (b * t_a * C) )

so the asymptotic power laws are

===================  ==========================================
lever                exponent on ``tau_inf``
===================  ==========================================
bond width ``b``     ``-1/2``  (**not** ``1/b`` - the width also raises ``beta``)
thickness ``t_a``    ``-1/2``
modulus ``G_a``      ``+1/2``
``|d_eps|``          ``+1``    (so linear in both ``|dT|`` and ``|alpha_1-alpha_2|``)
compliance ``C``     ``-1/2``  (stiffer adherends raise the asymptote)
===================  ==========================================

The practical consequence is that halving the peak shear costs a **factor of
four** in width, thickness or modulus - the square-root makes every single
lever expensive, which is why combinations matter.

Monotonicity (established before any inverse search)
----------------------------------------------------
Write ``tau_peak = |N_t| f(beta) / b`` with ``f(u) = u coth(u L_b)``. With
``y = u L_b``, ``df/dy = (coth y - y csch^2 y)/L_b``, and multiplying through by
``sinh^2 y`` gives ``sinh(2y)/2 - y > 0`` for ``y > 0``. So ``f`` is **strictly
increasing in beta**, and therefore

* ``t_a`` up   -> ``beta`` down -> ``tau_peak`` strictly **down**
* ``G_a`` up   -> ``beta`` up   -> ``tau_peak`` strictly **up**
* ``b`` up     -> ``tau_peak = |N_t| k coth(k L_b sqrt(b)) / sqrt(b)`` with
  ``k = sqrt(G_a C / t_a)``; both factors fall, so strictly **down**, and it
  tends to zero.

The uniform-shear floor
-----------------------
As ``beta -> 0`` (a very thick or very soft bondline) ``f(beta) -> 1/L_b`` and

    tau_peak -> |N_t| / (b L_b) = tau_avg

``tau_avg`` is a hard **floor** on what the thickness and modulus levers can
achieve at fixed width and overlap: if ``tau_avg`` already exceeds the
allowable, no thickness and no modulus can pass. Bond width has no such floor,
because ``tau_avg`` itself falls as ``1/b``. This is what makes width the only
unconditionally effective single lever in this model.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Sequence

from ._validation import require_finite, require_positive
from .adhesive import AdhesiveMaterial, AdhesiveShearBasis, BondedOverlapGeometry
from .environment import ThermalEnvironment
from .extremes import assess_temperature_extremes
from .margins import YieldBasis
from .materials import AxialMember, ThermoelasticMaterial
from .shear_lag import (
    adherend_compliance_sum,
    shear_lag_parameter,
    solve_shear_lag,
    thermal_mismatch_strain,
)
from .shear_lag_extremes import assess_shear_lag_extremes

__all__ = [
    "long_overlap_peak_shear",
    "uniform_shear_floor",
    "JointDesignCandidate",
    "SensitivityPoint",
    "DesignSelectionPolicy",
    "evaluate_joint_design",
    "bond_width_sensitivity",
    "adhesive_thickness_sensitivity",
    "adhesive_modulus_sensitivity",
    "area_ratio_sensitivity",
    "thermal_excursion_sensitivity",
    "cte_mismatch_sensitivity",
    "width_thickness_design_map",
    "select_preliminary_joint_design",
]


def long_overlap_peak_shear(
    member_1: AxialMember,
    member_2: AxialMember,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    delta_temperature: float,
) -> float:
    """Analytic ``L_b -> inf`` peak adhesive shear [Pa] for a whole configuration.

    A convenience wrapper over the Milestone 3 primitives, evaluated as
    ``|d_eps| sqrt(G_a / (b t_a C))``. Independent of the overlap length, and a
    hard lower bound on peak shear at that width, thickness and modulus.
    """
    mismatch = abs(thermal_mismatch_strain(member_1, member_2, delta_temperature))
    compliance = adherend_compliance_sum(member_1, member_2)
    return mismatch * math.sqrt(
        adhesive.shear_modulus
        / (geometry.bond_width * geometry.adhesive_thickness * compliance)
    )


def uniform_shear_floor(
    member_1: AxialMember,
    member_2: AxialMember,
    geometry: BondedOverlapGeometry,
    delta_temperature: float,
) -> float:
    """Uniform-shear floor ``|N_t| / (b L_b)`` [Pa].

    The peak shear that a perfectly compliant bondline would give. Thickness and
    modulus can approach it but never beat it, so it decides whether those two
    levers can work at all at a given width and overlap.
    """
    mismatch = abs(thermal_mismatch_strain(member_1, member_2, delta_temperature))
    transferred = mismatch / adherend_compliance_sum(member_1, member_2)
    return transferred / geometry.bond_area


@dataclass(frozen=True)
class JointDesignCandidate:
    """One fully evaluated bondline configuration.

    Yield and adhesive margins are carried **separately and are never combined
    numerically**; only the feasibility booleans are ANDed.
    """

    bond_width: float
    adhesive_thickness: float
    adhesive_shear_modulus: float
    overlap_length: float
    bond_area: float
    adhesive_volume: float
    beta: float
    transfer_length: float
    dimensionless_overlap: float
    hot_peak_shear_stress: float
    cold_peak_shear_stress: float
    governing_extreme: str
    peak_shear_stress: float
    allowable_shear_stress: float
    hot_adhesive_margin: float
    cold_adhesive_margin: float
    adhesive_margin: float
    minimum_yield_margin: float
    long_overlap_asymptote: float
    uniform_shear_floor: float

    @property
    def adhesive_feasible(self) -> bool:
        """True when the preliminary adhesive shear margin is ``>= 0``."""
        return self.adhesive_margin >= 0.0

    @property
    def yield_feasible(self) -> bool:
        """True when the preliminary member yield margin is ``>= 0``."""
        return self.minimum_yield_margin >= 0.0

    @property
    def overall_feasible(self) -> bool:
        """Boolean AND of the two screens - never a combined margin."""
        return self.adhesive_feasible and self.yield_feasible


def evaluate_joint_design(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    yield_basis: YieldBasis | None = None,
    adhesive_basis: AdhesiveShearBasis | None = None,
) -> JointDesignCandidate:
    """Evaluate one configuration at both temperature extremes.

    Both extremes are always computed and the governing one is taken from the
    actual margins - cold is never assumed.
    """
    shear = assess_shear_lag_extremes(
        member_1, member_2, environment, geometry, adhesive, adhesive_basis
    )
    yield_assessment = assess_temperature_extremes(
        member_1, member_2, environment, yield_basis
    )
    governing_delta_t = shear.governing_demand.delta_temperature

    return JointDesignCandidate(
        bond_width=geometry.bond_width,
        adhesive_thickness=geometry.adhesive_thickness,
        adhesive_shear_modulus=adhesive.shear_modulus,
        overlap_length=geometry.overlap_length,
        bond_area=geometry.bond_area,
        adhesive_volume=geometry.bond_area * geometry.adhesive_thickness,
        beta=shear.cold_demand.beta,
        transfer_length=shear.cold_demand.transfer_length,
        dimensionless_overlap=shear.cold_demand.dimensionless_overlap,
        hot_peak_shear_stress=shear.hot_demand.peak_shear_stress,
        cold_peak_shear_stress=shear.cold_demand.peak_shear_stress,
        governing_extreme=shear.governing_extreme,
        peak_shear_stress=shear.governing_peak_shear_stress,
        allowable_shear_stress=shear.allowable_shear_stress,
        hot_adhesive_margin=shear.hot_margin,
        cold_adhesive_margin=shear.cold_margin,
        adhesive_margin=shear.minimum_margin,
        minimum_yield_margin=yield_assessment.minimum_yield_margin,
        long_overlap_asymptote=long_overlap_peak_shear(
            member_1, member_2, geometry, adhesive, governing_delta_t
        ),
        uniform_shear_floor=uniform_shear_floor(
            member_1, member_2, geometry, governing_delta_t
        ),
    )


@dataclass(frozen=True)
class SensitivityPoint:
    """One point of a one-variable sweep whose variable is not a bondline dimension.

    Used for the area-ratio, thermal-excursion and CTE-mismatch studies, where
    the swept quantity belongs to the members or the environment rather than to
    :class:`JointDesignCandidate`.
    """

    variable: str
    value: float
    beta: float
    transfer_length: float
    transferred_force: float
    peak_shear_stress: float
    governing_extreme: str
    adhesive_margin: float
    minimum_yield_margin: float

    @property
    def adhesive_feasible(self) -> bool:
        return self.adhesive_margin >= 0.0

    @property
    def yield_feasible(self) -> bool:
        return self.minimum_yield_margin >= 0.0

    @property
    def overall_feasible(self) -> bool:
        return self.adhesive_feasible and self.yield_feasible


def _sensitivity_point(
    variable: str,
    value: float,
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    yield_basis: YieldBasis | None,
    adhesive_basis: AdhesiveShearBasis | None,
) -> SensitivityPoint:
    shear = assess_shear_lag_extremes(
        member_1, member_2, environment, geometry, adhesive, adhesive_basis
    )
    yield_assessment = assess_temperature_extremes(
        member_1, member_2, environment, yield_basis
    )
    governing = shear.governing_demand
    return SensitivityPoint(
        variable=variable,
        value=value,
        beta=governing.beta,
        transfer_length=governing.transfer_length,
        transferred_force=governing.transferred_force,
        peak_shear_stress=governing.peak_shear_stress,
        governing_extreme=shear.governing_extreme,
        adhesive_margin=shear.minimum_margin,
        minimum_yield_margin=yield_assessment.minimum_yield_margin,
    )


def bond_width_sensitivity(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    bond_widths: Sequence[float],
    yield_basis: YieldBasis | None = None,
    adhesive_basis: AdhesiveShearBasis | None = None,
) -> tuple[JointDesignCandidate, ...]:
    """Sweep bond width at fixed overlap, thickness, modulus, members and dT.

    Widths are in metres. Copies of the geometry are built for each point; the
    supplied geometry is never mutated.
    """
    return tuple(
        evaluate_joint_design(
            member_1,
            member_2,
            environment,
            BondedOverlapGeometry(
                overlap_length=geometry.overlap_length,
                bond_width=require_positive(width, "bond_width"),
                adhesive_thickness=geometry.adhesive_thickness,
            ),
            adhesive,
            yield_basis,
            adhesive_basis,
        )
        for width in bond_widths
    )


def adhesive_thickness_sensitivity(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    adhesive_thicknesses: Sequence[float],
    yield_basis: YieldBasis | None = None,
    adhesive_basis: AdhesiveShearBasis | None = None,
) -> tuple[JointDesignCandidate, ...]:
    """Sweep adhesive thickness [m] at fixed width, modulus, overlap and members."""
    return tuple(
        evaluate_joint_design(
            member_1,
            member_2,
            environment,
            geometry.with_adhesive_thickness(thickness),
            adhesive,
            yield_basis,
            adhesive_basis,
        )
        for thickness in adhesive_thicknesses
    )


def _modulus_variant(adhesive: AdhesiveMaterial, shear_modulus: float) -> AdhesiveMaterial:
    """A copy of ``adhesive`` with a different shear modulus, same strength and note."""
    return AdhesiveMaterial(
        name=f"{adhesive.name} (G_a sweep variant)",
        shear_modulus=require_positive(shear_modulus, "shear_modulus"),
        shear_strength=adhesive.shear_strength,
        source_note=adhesive.source_note,
        notes=adhesive.notes,
    )


def adhesive_modulus_sensitivity(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    shear_moduli: Sequence[float],
    yield_basis: YieldBasis | None = None,
    adhesive_basis: AdhesiveShearBasis | None = None,
) -> tuple[JointDesignCandidate, ...]:
    """Sweep adhesive shear modulus [Pa]; the shear strength is held fixed.

    Only the stiffness is varied, so the allowable is unchanged across the
    sweep. The canonical adhesive record is never mutated - variants are copies.
    """
    return tuple(
        evaluate_joint_design(
            member_1,
            member_2,
            environment,
            geometry,
            _modulus_variant(adhesive, modulus),
            yield_basis,
            adhesive_basis,
        )
        for modulus in shear_moduli
    )


def area_ratio_sensitivity(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    area_ratios: Sequence[float],
    yield_basis: YieldBasis | None = None,
    adhesive_basis: AdhesiveShearBasis | None = None,
) -> tuple[SensitivityPoint, ...]:
    """Sweep ``A_1 / A_2`` end to end, recomputing the Milestone 1 demand each time.

    Member 2 is held fixed and member 1's area is set to ``ratio * A_2``, so the
    transferred force, ``beta`` and the yield margins are all recomputed.
    """
    points = []
    for ratio in area_ratios:
        scaled = AxialMember(
            member_1.material,
            require_positive(ratio, "area_ratio") * member_2.area,
            label=member_1.label,
        )
        points.append(
            _sensitivity_point(
                "area_ratio", ratio, scaled, member_2, environment, geometry,
                adhesive, yield_basis, adhesive_basis,
            )
        )
    return tuple(points)


def thermal_excursion_sensitivity(
    member_1: AxialMember,
    member_2: AxialMember,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    excursion_magnitudes: Sequence[float],
    reference_temperature: float = 20.0,
    yield_basis: YieldBasis | None = None,
    adhesive_basis: AdhesiveShearBasis | None = None,
) -> tuple[SensitivityPoint, ...]:
    """Sweep the excursion magnitude ``|dT|`` [K] with symmetric hot/cold extremes.

    Each point uses ``T_ref -/+ |dT|``, so both extremes have the same magnitude
    and the governing selection is exercised at every point.
    """
    points = []
    for magnitude in excursion_magnitudes:
        span = require_finite(magnitude, "excursion_magnitude")
        if span < 0.0:
            raise ValueError(f"excursion_magnitude must be >= 0, got {span!r}")
        environment = ThermalEnvironment(
            reference_temperature=reference_temperature,
            cold_temperature=reference_temperature - span,
            hot_temperature=reference_temperature + span,
        )
        points.append(
            _sensitivity_point(
                "excursion_magnitude", span, member_1, member_2, environment,
                geometry, adhesive, yield_basis, adhesive_basis,
            )
        )
    return tuple(points)


def cte_mismatch_sensitivity(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    mismatch_scales: Sequence[float],
    yield_basis: YieldBasis | None = None,
    adhesive_basis: AdhesiveShearBasis | None = None,
) -> tuple[SensitivityPoint, ...]:
    """Scale the CTE mismatch ``alpha_1 - alpha_2`` by each factor.

    Member 2 keeps its CTE and a synthetic copy of member 1's material is built
    with ``alpha_1' = alpha_2 + scale * (alpha_1 - alpha_2)``. Everything else -
    moduli, areas, geometry, temperatures - is held fixed. The canonical
    material records are never mutated.
    """
    alpha_2 = member_2.material.thermal_expansion_coefficient
    base_mismatch = member_1.material.thermal_expansion_coefficient - alpha_2

    points = []
    for scale in mismatch_scales:
        factor = require_finite(scale, "mismatch_scale")
        material = ThermoelasticMaterial(
            name=f"{member_1.material.name} (CTE mismatch x{factor:g})",
            elastic_modulus=member_1.material.elastic_modulus,
            thermal_expansion_coefficient=alpha_2 + factor * base_mismatch,
            yield_strength=member_1.material.yield_strength,
            density=member_1.material.density,
        )
        scaled = AxialMember(material, member_1.area, label=member_1.label)
        points.append(
            _sensitivity_point(
                "mismatch_scale", factor, scaled, member_2, environment, geometry,
                adhesive, yield_basis, adhesive_basis,
            )
        )
    return tuple(points)


def width_thickness_design_map(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    bond_widths: Sequence[float],
    adhesive_thicknesses: Sequence[float],
    yield_basis: YieldBasis | None = None,
    adhesive_basis: AdhesiveShearBasis | None = None,
) -> tuple[JointDesignCandidate, ...]:
    """A **bounded preliminary design map** over bond width x adhesive thickness.

    This is not an optimization. It is a deterministic evaluation of every point
    of an explicitly bounded grid, in row-major order (width outer, thickness
    inner), using only the verified mechanics.
    """
    return tuple(
        evaluate_joint_design(
            member_1,
            member_2,
            environment,
            BondedOverlapGeometry(
                overlap_length=geometry.overlap_length,
                bond_width=require_positive(width, "bond_width"),
                adhesive_thickness=require_positive(thickness, "adhesive_thickness"),
            ),
            adhesive,
            yield_basis,
            adhesive_basis,
        )
        for width in bond_widths
        for thickness in adhesive_thicknesses
    )


class DesignSelectionPolicy(str, Enum):
    """Explicit, deterministic rule for picking one point from a bounded grid."""

    MINIMUM_BOND_AREA = "minimum_bond_area"
    """Feasible on both screens, then smallest bond area.

    Ranking key, applied in order:

    1. must be ``overall_feasible`` (member yield **and** adhesive shear)
    2. smallest ``bond_area`` (= width x overlap; with the overlap fixed in this
       milestone this is simply the narrowest bond)
    3. tie-break: smallest ``adhesive_thickness``
    4. tie-break: largest ``adhesive_margin``
    5. tie-break: earliest position in the input sequence

    This is a screening rule, not an objective function. A thinner adhesive is
    not automatically better for manufacturing or durability, and a wider bond
    is not automatically worse.
    """


def select_preliminary_joint_design(
    candidates: Sequence[JointDesignCandidate],
    policy: DesignSelectionPolicy = DesignSelectionPolicy.MINIMUM_BOND_AREA,
) -> Optional[JointDesignCandidate]:
    """Pick one feasible candidate from a bounded grid under an explicit policy.

    Returns ``None`` when no candidate passes both screens. The result is a
    **preliminary feasible point within the bounded study grid, not an optimized
    flight-joint design** and not a global optimum.
    """
    if policy is not DesignSelectionPolicy.MINIMUM_BOND_AREA:
        raise ValueError(f"unsupported selection policy: {policy!r}")

    feasible = [
        (index, candidate)
        for index, candidate in enumerate(candidates)
        if candidate.overall_feasible
    ]
    if not feasible:
        return None

    def key(item: tuple[int, JointDesignCandidate]):
        index, candidate = item
        return (
            candidate.bond_area,
            candidate.adhesive_thickness,
            -candidate.adhesive_margin,
            index,
        )

    return min(feasible, key=key)[1]

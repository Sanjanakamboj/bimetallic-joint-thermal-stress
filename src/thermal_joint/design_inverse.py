r"""Milestone 4 inverse sizing: bond width, adhesive thickness, adhesive modulus.

Each utility answers a single-lever sizing question against the governing
temperature extreme, using the verified Milestone 3 mechanics. Every search is
deterministic bounded bisection with explicit bounds, tolerance and iteration
cap, and **no bracket is ever expanded silently**.

Feasibility is settled analytically before any search runs, using the two limits
derived in :mod:`thermal_joint.design_trade`:

* ``tau_peak`` is strictly **decreasing** in bond width and in adhesive
  thickness, and strictly **increasing** in adhesive shear modulus.
* Thickness and modulus are both bounded below by the **uniform-shear floor**
  ``tau_avg = |N_t| / (b L_b)``. If that floor exceeds the allowable, no
  thickness and no modulus can pass, and the search is skipped rather than
  returning an invented number.
* Bond width has no floor - ``tau_peak -> 0`` as ``b -> inf`` - so a finite
  required width always exists in principle, and only the search bracket can
  fail to contain it.

These are model-based screening values, not allowables.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Callable, Optional

from ._validation import require_positive
from .adhesive import AdhesiveMaterial, AdhesiveShearBasis, BondedOverlapGeometry
from .design_trade import _modulus_variant, uniform_shear_floor
from .environment import ThermalEnvironment
from .materials import AxialMember
from .shear_lag_extremes import assess_shear_lag_extremes

__all__ = [
    "BondWidthStatus",
    "AdhesiveThicknessStatus",
    "AdhesiveModulusStatus",
    "RequiredBondWidthResult",
    "RequiredAdhesiveThicknessResult",
    "MaximumAdhesiveModulusResult",
    "required_bond_width",
    "required_adhesive_thickness",
    "maximum_allowable_adhesive_shear_modulus",
]


def _validate_bracket(low: float, high: float, tolerance: float, max_iterations: int,
                      low_name: str, high_name: str) -> tuple[float, float, float]:
    lower = require_positive(low, low_name)
    upper = require_positive(high, high_name)
    tol = require_positive(tolerance, "tolerance")
    if lower >= upper:
        raise ValueError(f"{low_name} must be < {high_name}, got {lower!r} and {upper!r}")
    if int(max_iterations) < 1:
        raise ValueError(f"max_iterations must be >= 1, got {max_iterations!r}")
    return lower, upper, tol


def _bisect(
    passes: Callable[[float], bool],
    low: float,
    high: float,
    tolerance: float,
    max_iterations: int,
    keep: str,
) -> tuple[float, int]:
    """Bisect a monotone pass/fail boundary.

    ``keep='high'`` for a decreasing demand (the high end passes, return it);
    ``keep='low'`` for an increasing demand (the low end passes, return it).
    """
    iterations = 0
    while high - low > tolerance and iterations < int(max_iterations):
        middle = 0.5 * (low + high)
        if passes(middle):
            if keep == "high":
                high = middle
            else:
                low = middle
        else:
            if keep == "high":
                low = middle
            else:
                high = middle
        iterations += 1
    return (high if keep == "high" else low), iterations


class BondWidthStatus(str, Enum):
    """Outcome of the minimum-bond-width search."""

    LOWER_BOUND_ALREADY_PASSES = "lower_bound_already_passes"
    """The narrowest width in the bracket already meets the allowable."""

    FINITE_REQUIRED_WIDTH = "finite_required_width"
    """A finite minimum bond width satisfies the adhesive shear allowable."""

    NO_BOUNDARY_WITHIN_SEARCH_BOUNDS = "no_boundary_within_search_bounds"
    """A finite width exists but lies above the requested upper bound.

    Peak shear falls to zero as the bond widens, so a sufficient width always
    exists in this model; only the bracket was too narrow. Re-run with a larger
    ``maximum_bond_width`` - bounds are never expanded silently.
    """


class AdhesiveThicknessStatus(str, Enum):
    """Outcome of the minimum-adhesive-thickness search."""

    LOWER_BOUND_ALREADY_PASSES = "lower_bound_already_passes"
    FINITE_REQUIRED_THICKNESS = "finite_required_thickness"
    NO_FINITE_THICKNESS_WITHIN_MODEL = "no_finite_thickness_within_model"
    """The uniform-shear floor already exceeds the allowable at this width.

    No bondline, however thick or compliant, can pass. The lever has to be bond
    width, overlap, CTE mismatch or temperature excursion instead.
    """
    NO_BOUNDARY_WITHIN_SEARCH_BOUNDS = "no_boundary_within_search_bounds"


class AdhesiveModulusStatus(str, Enum):
    """Outcome of the maximum-adhesive-shear-modulus search."""

    UPPER_BOUND_ALREADY_PASSES = "upper_bound_already_passes"
    """Even the stiffest modulus in the bracket passes; modulus does not size this joint."""

    FINITE_MAXIMUM_MODULUS = "finite_maximum_modulus"
    NO_FINITE_MODULUS_WITHIN_MODEL = "no_finite_modulus_within_model"
    """The uniform-shear floor already exceeds the allowable; even ``G_a -> 0`` fails."""

    LOWER_BOUND_ALREADY_FAILS = "lower_bound_already_fails"
    """The softest modulus in the bracket still fails, though a softer one could pass.

    The boundary lies below the requested lower bound; re-run with a smaller
    ``minimum_shear_modulus``. Bounds are never expanded silently.
    """


@dataclass(frozen=True)
class RequiredBondWidthResult:
    """Result of the minimum-bond-width search. All stresses [Pa], widths [m]."""

    status: BondWidthStatus
    required_bond_width: Optional[float]
    margin_at_required_width: Optional[float]
    governing_extreme: str
    peak_shear_at_lower_bound: float
    peak_shear_at_upper_bound: float
    allowable_shear_stress: float
    minimum_bond_width: float
    maximum_bond_width: float
    tolerance: float
    iterations: int

    @property
    def has_finite_width(self) -> bool:
        return self.status is BondWidthStatus.FINITE_REQUIRED_WIDTH


@dataclass(frozen=True)
class RequiredAdhesiveThicknessResult:
    """Result of the minimum-adhesive-thickness search."""

    status: AdhesiveThicknessStatus
    required_adhesive_thickness: Optional[float]
    margin_at_required_thickness: Optional[float]
    governing_extreme: str
    peak_shear_at_lower_bound: float
    peak_shear_at_upper_bound: float
    uniform_shear_floor: float
    allowable_shear_stress: float
    minimum_adhesive_thickness: float
    maximum_adhesive_thickness: float
    tolerance: float
    iterations: int

    @property
    def has_finite_thickness(self) -> bool:
        return self.status is AdhesiveThicknessStatus.FINITE_REQUIRED_THICKNESS


@dataclass(frozen=True)
class MaximumAdhesiveModulusResult:
    """Result of the maximum-adhesive-shear-modulus search."""

    status: AdhesiveModulusStatus
    maximum_shear_modulus: Optional[float]
    margin_at_maximum_modulus: Optional[float]
    governing_extreme: str
    peak_shear_at_lower_bound: float
    peak_shear_at_upper_bound: float
    uniform_shear_floor: float
    allowable_shear_stress: float
    minimum_shear_modulus: float
    maximum_search_shear_modulus: float
    tolerance: float
    iterations: int

    @property
    def has_finite_modulus(self) -> bool:
        return self.status is AdhesiveModulusStatus.FINITE_MAXIMUM_MODULUS


def _governing(member_1, member_2, environment, geometry, adhesive, basis):
    return assess_shear_lag_extremes(
        member_1, member_2, environment, geometry, adhesive, basis
    )


def required_bond_width(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    basis: AdhesiveShearBasis | None = None,
    minimum_bond_width: float = 1.0e-3,
    maximum_bond_width: float = 2.0,
    tolerance: float = 1.0e-6,
    max_iterations: int = 200,
) -> RequiredBondWidthResult:
    """Narrowest bond width whose peak adhesive shear meets the allowable.

    Overlap length, adhesive thickness, adhesive modulus, members and thermal
    environment are all held fixed. Both extremes are evaluated at every trial
    width and the governing one is used.

    Peak shear is strictly decreasing in width (see
    :mod:`thermal_joint.design_trade`), so the bisection is well posed.
    """
    design_basis = AdhesiveShearBasis() if basis is None else basis
    lower, upper, tol = _validate_bracket(
        minimum_bond_width, maximum_bond_width, tolerance, max_iterations,
        "minimum_bond_width", "maximum_bond_width",
    )
    allowable = design_basis.allowable_shear_stress(adhesive)

    def geometry_at(width: float) -> BondedOverlapGeometry:
        return BondedOverlapGeometry(
            overlap_length=geometry.overlap_length,
            bond_width=width,
            adhesive_thickness=geometry.adhesive_thickness,
        )

    def assess_at(width: float):
        return _governing(
            member_1, member_2, environment, geometry_at(width), adhesive, design_basis
        )

    lower_assessment = assess_at(lower)
    governing_extreme = lower_assessment.governing_extreme
    peak_lower = lower_assessment.governing_peak_shear_stress
    peak_upper = assess_at(upper).governing_peak_shear_stress

    def build(status, width, margin, iterations):
        return RequiredBondWidthResult(
            status=status,
            required_bond_width=width,
            margin_at_required_width=margin,
            governing_extreme=governing_extreme,
            peak_shear_at_lower_bound=peak_lower,
            peak_shear_at_upper_bound=peak_upper,
            allowable_shear_stress=allowable,
            minimum_bond_width=lower,
            maximum_bond_width=upper,
            tolerance=tol,
            iterations=iterations,
        )

    if peak_lower <= allowable:
        return build(
            BondWidthStatus.LOWER_BOUND_ALREADY_PASSES, lower,
            lower_assessment.minimum_margin, 0,
        )
    if peak_upper > allowable:
        return build(BondWidthStatus.NO_BOUNDARY_WITHIN_SEARCH_BOUNDS, None, None, 0)

    width, iterations = _bisect(
        lambda candidate: assess_at(candidate).governing_peak_shear_stress <= allowable,
        lower, upper, tol, max_iterations, keep="high",
    )
    return build(
        BondWidthStatus.FINITE_REQUIRED_WIDTH, width,
        assess_at(width).minimum_margin, iterations,
    )


def required_adhesive_thickness(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    basis: AdhesiveShearBasis | None = None,
    minimum_adhesive_thickness: float = 1.0e-5,
    maximum_adhesive_thickness: float = 0.05,
    tolerance: float = 1.0e-9,
    max_iterations: int = 200,
) -> RequiredAdhesiveThicknessResult:
    """Thinnest adhesive layer whose peak adhesive shear meets the allowable.

    Peak shear is strictly decreasing in thickness but is bounded below by the
    uniform-shear floor, which is checked analytically first.
    """
    design_basis = AdhesiveShearBasis() if basis is None else basis
    lower, upper, tol = _validate_bracket(
        minimum_adhesive_thickness, maximum_adhesive_thickness, tolerance, max_iterations,
        "minimum_adhesive_thickness", "maximum_adhesive_thickness",
    )
    allowable = design_basis.allowable_shear_stress(adhesive)

    def assess_at(thickness: float):
        return _governing(
            member_1, member_2, environment,
            geometry.with_adhesive_thickness(thickness), adhesive, design_basis,
        )

    lower_assessment = assess_at(lower)
    governing_extreme = lower_assessment.governing_extreme
    peak_lower = lower_assessment.governing_peak_shear_stress
    peak_upper = assess_at(upper).governing_peak_shear_stress
    floor = max(
        uniform_shear_floor(member_1, member_2, geometry, delta_t)
        for delta_t in (
            environment.delta_temperature_cold,
            environment.delta_temperature_hot,
        )
    )

    def build(status, thickness, margin, iterations):
        return RequiredAdhesiveThicknessResult(
            status=status,
            required_adhesive_thickness=thickness,
            margin_at_required_thickness=margin,
            governing_extreme=governing_extreme,
            peak_shear_at_lower_bound=peak_lower,
            peak_shear_at_upper_bound=peak_upper,
            uniform_shear_floor=floor,
            allowable_shear_stress=allowable,
            minimum_adhesive_thickness=lower,
            maximum_adhesive_thickness=upper,
            tolerance=tol,
            iterations=iterations,
        )

    if peak_lower <= allowable:
        return build(
            AdhesiveThicknessStatus.LOWER_BOUND_ALREADY_PASSES, lower,
            lower_assessment.minimum_margin, 0,
        )
    if floor > allowable:
        return build(
            AdhesiveThicknessStatus.NO_FINITE_THICKNESS_WITHIN_MODEL, None, None, 0
        )
    if peak_upper > allowable:
        return build(
            AdhesiveThicknessStatus.NO_BOUNDARY_WITHIN_SEARCH_BOUNDS, None, None, 0
        )

    thickness, iterations = _bisect(
        lambda candidate: assess_at(candidate).governing_peak_shear_stress <= allowable,
        lower, upper, tol, max_iterations, keep="high",
    )
    return build(
        AdhesiveThicknessStatus.FINITE_REQUIRED_THICKNESS, thickness,
        assess_at(thickness).minimum_margin, iterations,
    )


def maximum_allowable_adhesive_shear_modulus(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    basis: AdhesiveShearBasis | None = None,
    minimum_shear_modulus: float = 1.0e6,
    maximum_shear_modulus: float = 1.0e10,
    tolerance: float = 1.0e3,
    max_iterations: int = 200,
) -> MaximumAdhesiveModulusResult:
    """Stiffest adhesive whose peak shear still meets the allowable.

    A stiffer adhesive concentrates load transfer, so peak shear is strictly
    **increasing** in ``G_a`` - the boundary is a maximum, not a minimum. The
    shear strength is held fixed across the search, so the allowable does not
    move with the modulus.

    As ``G_a -> 0`` peak shear tends to the uniform-shear floor, which is
    checked analytically before searching.
    """
    design_basis = AdhesiveShearBasis() if basis is None else basis
    lower, upper, tol = _validate_bracket(
        minimum_shear_modulus, maximum_shear_modulus, tolerance, max_iterations,
        "minimum_shear_modulus", "maximum_shear_modulus",
    )
    allowable = design_basis.allowable_shear_stress(adhesive)

    def assess_at(modulus: float):
        return _governing(
            member_1, member_2, environment, geometry,
            _modulus_variant(adhesive, modulus), design_basis,
        )

    lower_assessment = assess_at(lower)
    upper_assessment = assess_at(upper)
    governing_extreme = upper_assessment.governing_extreme
    peak_lower = lower_assessment.governing_peak_shear_stress
    peak_upper = upper_assessment.governing_peak_shear_stress
    floor = max(
        uniform_shear_floor(member_1, member_2, geometry, delta_t)
        for delta_t in (
            environment.delta_temperature_cold,
            environment.delta_temperature_hot,
        )
    )

    def build(status, modulus, margin, iterations):
        return MaximumAdhesiveModulusResult(
            status=status,
            maximum_shear_modulus=modulus,
            margin_at_maximum_modulus=margin,
            governing_extreme=governing_extreme,
            peak_shear_at_lower_bound=peak_lower,
            peak_shear_at_upper_bound=peak_upper,
            uniform_shear_floor=floor,
            allowable_shear_stress=allowable,
            minimum_shear_modulus=lower,
            maximum_search_shear_modulus=upper,
            tolerance=tol,
            iterations=iterations,
        )

    if peak_upper <= allowable:
        return build(
            AdhesiveModulusStatus.UPPER_BOUND_ALREADY_PASSES, upper,
            upper_assessment.minimum_margin, 0,
        )
    if floor > allowable:
        return build(AdhesiveModulusStatus.NO_FINITE_MODULUS_WITHIN_MODEL, None, None, 0)
    if peak_lower > allowable:
        return build(AdhesiveModulusStatus.LOWER_BOUND_ALREADY_FAILS, None, None, 0)

    modulus, iterations = _bisect(
        lambda candidate: assess_at(candidate).governing_peak_shear_stress <= allowable,
        lower, upper, tol, max_iterations, keep="low",
    )
    return build(
        AdhesiveModulusStatus.FINITE_MAXIMUM_MODULUS, modulus,
        assess_at(modulus).minimum_margin, iterations,
    )

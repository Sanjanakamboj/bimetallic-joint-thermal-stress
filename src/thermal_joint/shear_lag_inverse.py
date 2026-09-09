r"""Preliminary inverse design for the adhesive shear-lag screen.

Two questions:

1. **Allowable temperature excursion.** Peak adhesive shear is exactly linear in
   ``dT``, so the allowable excursion has a direct closed form - no search.
2. **Required overlap length.** Peak shear falls monotonically with overlap
   toward the finite asymptote ``tau_inf = |N_t| beta / b``. Whether a finite
   overlap can satisfy the allowable is therefore decided *analytically* against
   that asymptote before any search runs.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from typing import Optional

from ._validation import require_finite, require_positive
from .adhesive import AdhesiveMaterial, AdhesiveShearBasis, BondedOverlapGeometry
from .environment import ThermalEnvironment
from .materials import AxialMember
from .shear_lag import (
    adherend_compliance_sum,
    long_overlap_peak_shear_stress,
    shear_lag_parameter,
    solve_shear_lag,
    thermal_mismatch_strain,
)
from .shear_lag_extremes import assess_shear_lag_extremes

__all__ = [
    "OverlapLimitStatus",
    "RequiredOverlapResult",
    "allowable_temperature_change_for_adhesive_shear",
    "required_overlap_length",
]


def allowable_temperature_change_for_adhesive_shear(
    member_1: AxialMember,
    member_2: AxialMember,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    basis: AdhesiveShearBasis | None = None,
) -> float:
    """Largest ``|dT|`` the adhesive can take at this geometry [K], as a magnitude.

    The transferred force is ``|N_t| = |alpha_1 - alpha_2| |dT| / C`` with
    ``C = 1/(E_1 A_1) + 1/(E_2 A_2)``, and
    ``tau_peak = |N_t| beta coth(beta L_b) / b``. Peak shear is therefore
    strictly proportional to ``|dT|``, and setting ``MS_adh = 0`` gives the
    direct formula::

        |dT|_allow = tau_allow * b * C / ( |alpha_1 - alpha_2| beta coth(beta L_b) )

    No reference ``dT`` and no iteration are needed.

    Returns
    -------
    float
        A positive magnitude [K]. Returns ``math.inf`` when the two members have
        identical CTEs, since no mismatch force is generated at any temperature.
    """
    design_basis = AdhesiveShearBasis() if basis is None else basis
    allowable = design_basis.allowable_shear_stress(adhesive)

    mismatch_per_kelvin = abs(thermal_mismatch_strain(member_1, member_2, 1.0))
    if mismatch_per_kelvin == 0.0:
        return math.inf

    beta = shear_lag_parameter(member_1, member_2, geometry, adhesive)
    compliance = adherend_compliance_sum(member_1, member_2)
    coth = 1.0 / math.tanh(beta * geometry.overlap_length)

    return (
        allowable
        * geometry.bond_width
        * compliance
        / (mismatch_per_kelvin * beta * coth)
    )


class OverlapLimitStatus(str, Enum):
    """Outcome of the required-overlap search."""

    LOWER_BOUND_ALREADY_PASSES = "lower_bound_already_passes"
    """The shortest overlap in the search range already meets the allowable."""

    FINITE_REQUIRED_LENGTH = "finite_required_length"
    """A finite minimum overlap length satisfies the adhesive shear allowable."""

    NO_FINITE_LENGTH_WITHIN_MODEL = "no_finite_length_within_model"
    """The allowable is below the long-overlap asymptote, so no length works.

    Peak shear falls toward ``tau_inf = |N_t| beta / b`` but never below it. If
    ``tau_allow < tau_inf`` the criterion is unreachable at any overlap length -
    the fix has to be a different adhesive, bond width or thermal excursion, not
    more overlap.
    """

    NO_BOUNDARY_WITHIN_SEARCH_BOUNDS = "no_boundary_within_search_bounds"
    """A finite length exists but lies above the requested upper bound.

    Re-run with a larger ``maximum_overlap_length``; bounds are never expanded
    silently.
    """


@dataclass(frozen=True)
class RequiredOverlapResult:
    """Result of the minimum-overlap search.

    Attributes
    ----------
    status:
        An :class:`OverlapLimitStatus`.
    required_overlap_length:
        Minimum overlap length meeting the allowable [m], or ``None`` unless the
        status is ``FINITE_REQUIRED_LENGTH``.
    margin_at_required_length:
        Preliminary adhesive shear margin there [-]; approximately zero and
        non-negative by construction. ``None`` when no length was found.
    governing_extreme:
        The extreme that sets the demand, computed from the margins.
    peak_shear_at_lower_bound, peak_shear_at_upper_bound:
        Governing peak shear at the two search bounds [Pa].
    long_overlap_asymptote:
        ``tau_inf`` at the governing extreme [Pa] - the hard lower bound on peak
        shear.
    allowable_shear_stress:
        The adhesive allowable used [Pa].
    minimum_overlap_length, maximum_overlap_length:
        The explicit search bracket [m].
    length_tolerance:
        Absolute convergence tolerance on length [m].
    iterations:
        Bisection iterations performed.
    """

    status: OverlapLimitStatus
    required_overlap_length: Optional[float]
    margin_at_required_length: Optional[float]
    governing_extreme: str
    peak_shear_at_lower_bound: float
    peak_shear_at_upper_bound: float
    long_overlap_asymptote: float
    allowable_shear_stress: float
    minimum_overlap_length: float
    maximum_overlap_length: float
    length_tolerance: float
    iterations: int

    @property
    def has_finite_length(self) -> bool:
        """True when a finite required overlap length was established."""
        return self.status is OverlapLimitStatus.FINITE_REQUIRED_LENGTH


def required_overlap_length(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    basis: AdhesiveShearBasis | None = None,
    minimum_overlap_length: float = 1.0e-3,
    maximum_overlap_length: float = 1.0,
    length_tolerance: float = 1.0e-6,
    max_iterations: int = 200,
) -> RequiredOverlapResult:
    """Smallest overlap length whose peak adhesive shear meets the allowable.

    Both temperature extremes are evaluated at every trial length and the
    governing one is used, so the answer covers the whole excursion.

    Feasibility is settled analytically before any search: peak shear decreases
    monotonically with overlap length toward ``tau_inf = |N_t| beta / b``, so if
    ``tau_allow < tau_inf`` no finite overlap can ever pass and the search is
    skipped.

    Parameters
    ----------
    minimum_overlap_length, maximum_overlap_length:
        Explicit search bracket [m]; both positive with ``min < max``.
    length_tolerance:
        Absolute convergence tolerance on the length [m].
    max_iterations:
        Hard iteration cap.
    """
    design_basis = AdhesiveShearBasis() if basis is None else basis
    lower = require_positive(minimum_overlap_length, "minimum_overlap_length")
    upper = require_positive(maximum_overlap_length, "maximum_overlap_length")
    tolerance = require_positive(length_tolerance, "length_tolerance")
    if lower >= upper:
        raise ValueError(
            "minimum_overlap_length must be < maximum_overlap_length, got "
            f"{lower!r} and {upper!r}"
        )
    if int(max_iterations) < 1:
        raise ValueError(f"max_iterations must be >= 1, got {max_iterations!r}")

    allowable = design_basis.allowable_shear_stress(adhesive)

    def assess_at(length: float):
        return assess_shear_lag_extremes(
            member_1,
            member_2,
            environment,
            geometry.with_overlap_length(length),
            adhesive,
            design_basis,
        )

    def peak_at(length: float) -> float:
        return assess_at(length).governing_peak_shear_stress

    lower_assessment = assess_at(lower)
    governing_extreme = lower_assessment.governing_extreme
    peak_lower = lower_assessment.governing_peak_shear_stress
    peak_upper = peak_at(upper)

    beta = shear_lag_parameter(member_1, member_2, geometry, adhesive)
    asymptote = max(
        long_overlap_peak_shear_stress(
            solve_shear_lag(
                member_1, member_2, geometry, adhesive, delta_t
            ).transferred_force,
            beta,
            geometry.bond_width,
        )
        for delta_t in (
            environment.delta_temperature_cold,
            environment.delta_temperature_hot,
        )
    )

    def build(
        status: OverlapLimitStatus,
        length: Optional[float],
        margin: Optional[float],
        iterations: int,
    ) -> RequiredOverlapResult:
        return RequiredOverlapResult(
            status=status,
            required_overlap_length=length,
            margin_at_required_length=margin,
            governing_extreme=governing_extreme,
            peak_shear_at_lower_bound=peak_lower,
            peak_shear_at_upper_bound=peak_upper,
            long_overlap_asymptote=asymptote,
            allowable_shear_stress=allowable,
            minimum_overlap_length=lower,
            maximum_overlap_length=upper,
            length_tolerance=tolerance,
            iterations=iterations,
        )

    if peak_lower <= allowable:
        return build(
            OverlapLimitStatus.LOWER_BOUND_ALREADY_PASSES,
            lower,
            lower_assessment.minimum_margin,
            0,
        )

    if allowable < asymptote:
        return build(OverlapLimitStatus.NO_FINITE_LENGTH_WITHIN_MODEL, None, None, 0)

    if peak_upper > allowable:
        return build(OverlapLimitStatus.NO_BOUNDARY_WITHIN_SEARCH_BOUNDS, None, None, 0)

    low, high = lower, upper  # low always fails, high always passes
    iterations = 0
    while high - low > tolerance and iterations < int(max_iterations):
        middle = 0.5 * (low + high)
        if peak_at(middle) <= allowable:
            high = middle
        else:
            low = middle
        iterations += 1

    return build(
        OverlapLimitStatus.FINITE_REQUIRED_LENGTH,
        high,
        assess_at(high).minimum_margin,
        iterations,
    )

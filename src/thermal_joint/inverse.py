"""Preliminary inverse design: how much external restraint can the joint tolerate?

Answers the Milestone 2 design question:

    For this joint, thermal environment and yield basis, what is the maximum
    external restraint stiffness that can be tolerated before the minimum
    preliminary elastic yield margin reaches zero?

The search is deterministic bounded bisection on the dimensionless restraint
parameter ``eta_r = K_r / (E_1 A_1 + E_2 A_2)``, which is far better
conditioned than searching on ``K_r`` directly. Bounds, tolerance and iteration
count are all explicit; the search never silently expands its bracket.

Feasibility is resolved honestly before any search runs, using the analytic
fully restrained asymptote ``sigma_i = E_i (eps_ref - alpha_i dT)``:

* if the **free** joint already fails, there is nothing to search for
* if even a **rigid** restraint still passes yield, then no finite maximum
  exists within this linear-elastic model, and no number is invented
* otherwise the pass/fail boundary is bracketed and bisected
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional

from ._validation import require_finite, require_positive
from .environment import ThermalEnvironment
from .margins import YieldBasis
from .materials import AxialMember
from .restrained_extremes import assess_restrained_temperature_extremes
from .restraint import AxialRestraint, joint_axial_stiffness, rigid_restraint_stresses

__all__ = [
    "RestraintLimitStatus",
    "MaximumRestraintResult",
    "maximum_allowable_restraint_stiffness",
    "rigid_restraint_minimum_margin",
]


class RestraintLimitStatus(str, Enum):
    """Outcome of the maximum-restraint search."""

    FREE_JOINT_ALREADY_FAILS = "free_joint_already_fails"
    """The unrestrained joint already yields; adding restraint cannot help."""

    FINITE_LIMIT = "finite_limit"
    """A finite restraint stiffness bounds the acceptable design space."""

    NO_FINITE_LIMIT_WITHIN_MODEL = "no_finite_limit_within_model"
    """Even a fully rigid restraint stays below yield in this linear model."""

    NO_BOUNDARY_WITHIN_SEARCH_BOUNDS = "no_boundary_within_search_bounds"
    """A boundary exists asymptotically but lies above the requested bound.

    The search bound was too small. Re-run with a larger
    ``upper_stiffness_ratio`` - the bracket is never expanded silently.
    """


def rigid_restraint_minimum_margin(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    basis: YieldBasis,
    reference_strain: float = 0.0,
) -> float:
    """Minimum yield margin in the fully restrained (``K_r -> inf``) asymptote [-].

    Uses the analytic rigid stresses rather than a large-stiffness numerical
    limit, so the asymptotic pass/fail question is answered exactly.
    """
    margins = []
    for delta_t in (
        environment.delta_temperature_cold,
        environment.delta_temperature_hot,
    ):
        stress_1, stress_2 = rigid_restraint_stresses(
            member_1, member_2, delta_t, reference_strain
        )
        margins.append(basis.margin_of_safety(member_1.material, stress_1))
        margins.append(basis.margin_of_safety(member_2.material, stress_2))
    return min(margins)


@dataclass(frozen=True)
class MaximumRestraintResult:
    """Result of the maximum-allowable-restraint search.

    Attributes
    ----------
    status:
        A :class:`RestraintLimitStatus` describing the feasibility outcome.
    maximum_stiffness_ratio:
        Largest ``eta_r`` that still passes [-], or ``None`` unless the status
        is ``FINITE_LIMIT``.
    maximum_restraint_stiffness:
        The same limit expressed as ``K_r`` [N], or ``None``.
    minimum_margin_at_limit:
        Minimum yield margin at the returned limit [-]; approximately zero and
        non-negative by construction. ``None`` when no limit was found.
    free_joint_minimum_margin:
        Minimum margin at ``eta_r = 0`` [-].
    rigid_restraint_minimum_margin:
        Minimum margin in the fully restrained asymptote [-].
    joint_axial_stiffness:
        ``E_1 A_1 + E_2 A_2`` [N], the scale used to convert ``eta_r`` to ``K_r``.
    search_upper_ratio:
        The explicit upper bound used for the bisection [-].
    tolerance:
        Absolute convergence tolerance on ``eta_r`` [-].
    iterations:
        Bisection iterations actually performed.
    """

    status: RestraintLimitStatus
    maximum_stiffness_ratio: Optional[float]
    maximum_restraint_stiffness: Optional[float]
    minimum_margin_at_limit: Optional[float]
    free_joint_minimum_margin: float
    rigid_restraint_minimum_margin: float
    joint_axial_stiffness: float
    search_upper_ratio: float
    tolerance: float
    iterations: int

    @property
    def has_finite_limit(self) -> bool:
        """True when a finite maximum restraint stiffness was established."""
        return self.status is RestraintLimitStatus.FINITE_LIMIT


def maximum_allowable_restraint_stiffness(
    member_1: AxialMember,
    member_2: AxialMember,
    environment: ThermalEnvironment,
    basis: YieldBasis | None = None,
    reference_strain: float = 0.0,
    upper_stiffness_ratio: float = 1.0e6,
    tolerance: float = 1.0e-6,
    max_iterations: int = 200,
) -> MaximumRestraintResult:
    """Find the largest tolerable external restraint stiffness, if one exists.

    Parameters
    ----------
    member_1, member_2:
        The two bonded members.
    environment:
        Reference, cold and hot temperatures; both extremes are evaluated.
    basis:
        Design basis; defaults to ``YieldBasis()``.
    reference_strain:
        Force-free strain of the restraint [-]; canonical value ``0.0``.
    upper_stiffness_ratio:
        Explicit upper search bound on ``eta_r``. Must be finite and positive.
    tolerance:
        Absolute convergence tolerance on ``eta_r``. Must be finite and positive.
    max_iterations:
        Hard iteration cap for the bisection.

    Returns
    -------
    MaximumRestraintResult
        Carrying an explicit status. A finite number is returned only when the
        status is ``FINITE_LIMIT``.

    Notes
    -----
    Bisection locates the first pass/fail crossing inside the bracket. With a
    zero reference strain and temperature-independent properties the minimum
    margin decreases monotonically with restraint, so that crossing is the
    global boundary; for unusual inputs the returned value should be read as
    the first boundary within the bracket.
    """
    design_basis = YieldBasis() if basis is None else basis
    upper = require_positive(upper_stiffness_ratio, "upper_stiffness_ratio")
    tol = require_positive(tolerance, "tolerance")
    eps_ref = require_finite(reference_strain, "reference_strain")
    if int(max_iterations) < 1:
        raise ValueError(f"max_iterations must be >= 1, got {max_iterations!r}")

    stiffness_scale = joint_axial_stiffness(member_1, member_2)

    def minimum_margin(ratio: float) -> float:
        restraint = AxialRestraint.from_stiffness_ratio(
            ratio, member_1, member_2, reference_strain=eps_ref
        )
        return assess_restrained_temperature_extremes(
            member_1, member_2, environment, restraint, design_basis
        ).minimum_yield_margin

    free_margin = minimum_margin(0.0)
    rigid_margin = rigid_restraint_minimum_margin(
        member_1, member_2, environment, design_basis, eps_ref
    )

    def build(
        status: RestraintLimitStatus,
        ratio: Optional[float],
        margin_at_limit: Optional[float],
        iterations: int,
    ) -> MaximumRestraintResult:
        return MaximumRestraintResult(
            status=status,
            maximum_stiffness_ratio=ratio,
            maximum_restraint_stiffness=None if ratio is None else ratio * stiffness_scale,
            minimum_margin_at_limit=margin_at_limit,
            free_joint_minimum_margin=free_margin,
            rigid_restraint_minimum_margin=rigid_margin,
            joint_axial_stiffness=stiffness_scale,
            search_upper_ratio=upper,
            tolerance=tol,
            iterations=iterations,
        )

    if free_margin < 0.0:
        return build(RestraintLimitStatus.FREE_JOINT_ALREADY_FAILS, None, None, 0)

    if rigid_margin >= 0.0:
        return build(RestraintLimitStatus.NO_FINITE_LIMIT_WITHIN_MODEL, None, None, 0)

    if minimum_margin(upper) >= 0.0:
        return build(RestraintLimitStatus.NO_BOUNDARY_WITHIN_SEARCH_BOUNDS, None, None, 0)

    low, high = 0.0, upper  # low always passes, high always fails
    iterations = 0
    while high - low > tol and iterations < int(max_iterations):
        middle = 0.5 * (low + high)
        if minimum_margin(middle) >= 0.0:
            low = middle
        else:
            high = middle
        iterations += 1

    return build(RestraintLimitStatus.FINITE_LIMIT, low, minimum_margin(low), iterations)

r"""Milestone 3: one-dimensional adhesive shear-lag screen for the bonded overlap.

Role
----
Milestones 1 and 2 are *global*: they give the CTE-mismatch force the joint must
carry, with and without external restraint. Milestone 3 is *local*: it asks how
that force is transferred through a finite bonded overlap, and what adhesive
shear stress results. The canonical demand is imported from the Milestone 1
free-joint solution - it is never re-derived here.

Assumptions (all of them)
-------------------------
* two axial adherends, constant overlap length and width
* thin adhesive layer of constant thickness
* constant adhesive shear modulus, linear elastic adhesive
* adherends carry axial force only
* no bending, no eccentricity, no peel stress, no adhesive normal stress
* no free-edge singularity model, no yielding, no slip or debond
* uniform temperature, uniform properties, perfect bond

This is a screening model, not an adhesive-joint certification analysis.

Sign and displacement convention
--------------------------------
``u_i(x)`` is the axial displacement of adherend ``i`` in the ``+x`` direction,
the relative slip is ``s(x) = u_1(x) - u_2(x)``, and the adhesive is linear::

    tau(x) = (G_a / t_a) * s(x)

Adherend forces are tension-positive::

    N_i = E_i A_i (du_i/dx - alpha_i dT)

For an axial bar carrying a distributed applied load ``q`` per unit length,
``dN/dx = -q``. When ``s > 0`` the adhesive drags adherend 1 backward and
adherend 2 forward, so ``q_1 = -b tau`` and ``q_2 = +b tau``, giving::

    dN_1/dx = +b tau        dN_2/dx = -b tau        d(N_1 + N_2)/dx = 0

(This is the sign-mirror of one common textbook ordering; it is the version
that is consistent with ``tau = +(G_a/t_a)(u_1 - u_2)`` and yields a stable,
non-oscillatory governing equation.)

Governing equation
------------------
With ``S_i = E_i A_i``, ``C = 1/S_1 + 1/S_2`` and the free-joint pair
``N_2 = -N_1``::

    ds/dx   = N_1/S_1 - N_2/S_2 + d_eps = C N_1 + d_eps
    d2s/dx2 = C dN_1/dx = C b tau = C b (G_a/t_a) s

    =>  d2s/dx2 - beta^2 s = 0,   beta^2 = (G_a b / t_a) (1/(E_1 A_1) + 1/(E_2 A_2))

Dimensions: ``G_a b / t_a`` is [N/m^2], ``1/(E A)`` is [1/N], so ``beta^2`` is
[1/m^2] and ``beta`` is [1/m]. The transfer length is ``1/beta`` [m], and the
single dimensionless overlap metric used everywhere in this package is::

    lambda = beta * L_b

Boundary conditions
-------------------
The domain is ``x in [0, L_b]``:

* ``x = 0`` is the **loaded transfer plane** - the inboard edge of the overlap,
  where the free-joint mismatch force pair is fully developed:
  ``N_1(0) = N_t`` and ``N_2(0) = -N_t``, with ``N_t`` taken from the
  Milestone 1 free-joint solution.
* ``x = L_b`` is the **free edge** of the overlap: ``N_1(L_b) = N_2(L_b) = 0``.

Prescribing the Milestone 1 force at the loaded plane is the conservative
screening idealisation: it assumes the mismatch force is fully developed and
asks the overlap to shear-transfer all of it.

Closed-form solution
--------------------
``s = A cosh(beta x) + B sinh(beta x)`` with ``N_1 = (ds/dx - d_eps)/C``.
``N_1(0) = N_t = -d_eps/C`` forces ``B = 0``; ``N_1(L_b) = 0`` then gives
``A = d_eps / (beta sinh(beta L_b))``, so::

    N_1(x) =  N_t [1 - sinh(beta x) / sinh(beta L_b)]
    N_2(x) = -N_1(x)
    tau(x) = -(N_t beta / b) cosh(beta x) / sinh(beta L_b)

``|tau|`` is proportional to ``cosh(beta x)``, whose derivative
``beta sinh(beta x)`` is strictly positive for ``x > 0``. The peak adhesive
shear therefore occurs at the **free edge** ``x = L_b`` - proven, not assumed::

    tau_peak = |N_t| beta coth(beta L_b) / b
    tau_avg  = |N_t| / (b L_b)                 (whole overlap, no factor-of-two)
    tau_peak / tau_avg = lambda coth(lambda)

Limits: ``lambda coth(lambda) -> 1`` as ``lambda -> 0`` (short overlap, nearly
uniform transfer) and ``-> lambda`` as ``lambda -> inf``. The long-overlap peak
shear approaches the finite asymptote ``tau_inf = |N_t| beta / b`` **from
above**, so added overlap gives diminishing peak-shear benefit and can never
reduce the peak below ``tau_inf``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ._validation import require_finite
from .adhesive import AdhesiveMaterial, BondedOverlapGeometry
from .joint import solve_bimetallic_joint
from .materials import AxialMember

__all__ = [
    "thermal_mismatch_strain",
    "adherend_compliance_sum",
    "shear_lag_parameter",
    "transfer_length",
    "dimensionless_overlap",
    "long_overlap_peak_shear_stress",
    "ShearLagDemand",
    "solve_shear_lag",
]


def thermal_mismatch_strain(
    member_1: AxialMember, member_2: AxialMember, delta_temperature: float
) -> float:
    """CTE mismatch strain ``(alpha_1 - alpha_2) * dT`` [-].

    Ordering is member 1 minus member 2, used consistently throughout.
    """
    d_t = require_finite(delta_temperature, "delta_temperature")
    return (
        member_1.material.thermal_expansion_coefficient
        - member_2.material.thermal_expansion_coefficient
    ) * d_t


def adherend_compliance_sum(member_1: AxialMember, member_2: AxialMember) -> float:
    """``C = 1/(E_1 A_1) + 1/(E_2 A_2)`` [1/N], the summed axial compliance."""
    return 1.0 / member_1.axial_stiffness + 1.0 / member_2.axial_stiffness


def shear_lag_parameter(
    member_1: AxialMember,
    member_2: AxialMember,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
) -> float:
    """Shear-lag parameter ``beta`` [1/m].

    ``beta = sqrt( (G_a b / t_a) * (1/(E_1 A_1) + 1/(E_2 A_2)) )``.

    Both adherend compliances enter, so a joint is only as good at spreading
    load as its more compliant member allows. ``beta`` does not depend on the
    overlap length.
    """
    if not isinstance(geometry, BondedOverlapGeometry):
        raise TypeError(
            f"geometry must be a BondedOverlapGeometry, got {type(geometry).__name__}"
        )
    if not isinstance(adhesive, AdhesiveMaterial):
        raise TypeError(
            f"adhesive must be an AdhesiveMaterial, got {type(adhesive).__name__}"
        )
    bond_stiffness = (
        adhesive.shear_modulus * geometry.bond_width / geometry.adhesive_thickness
    )
    return math.sqrt(bond_stiffness * adherend_compliance_sum(member_1, member_2))


def transfer_length(
    member_1: AxialMember,
    member_2: AxialMember,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
) -> float:
    """Characteristic load-transfer length ``1 / beta`` [m]."""
    return 1.0 / shear_lag_parameter(member_1, member_2, geometry, adhesive)


def dimensionless_overlap(
    member_1: AxialMember,
    member_2: AxialMember,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
) -> float:
    """Dimensionless overlap ``lambda = beta * L_b`` [-].

    The number of transfer lengths spanned by the overlap. ``lambda << 1`` is a
    short overlap with nearly uniform shear; ``lambda >> 1`` is a long overlap
    whose interior carries almost no shear.
    """
    return (
        shear_lag_parameter(member_1, member_2, geometry, adhesive)
        * geometry.overlap_length
    )


def long_overlap_peak_shear_stress(transferred_force: float, beta: float, bond_width: float) -> float:
    """Analytic ``L_b -> inf`` peak shear asymptote ``|N_t| beta / b`` [Pa].

    Peak shear approaches this value **from above**, so it is a hard lower bound
    on the achievable peak adhesive shear for a given demand, adhesive and bond
    width. No overlap length, however long, can beat it.
    """
    return abs(transferred_force) * beta / bond_width


def _sinh_ratio(numerator_argument: float, denominator_argument: float) -> float:
    """``sinh(a) / sinh(c)`` evaluated without overflowing for large arguments."""
    a, c = numerator_argument, denominator_argument
    if c == 0.0:
        raise ZeroDivisionError("sinh ratio with a zero denominator argument")
    # sinh(a)/sinh(c) = e^(a-c) (1 - e^-2a) / (1 - e^-2c)
    return math.exp(a - c) * (1.0 - math.exp(-2.0 * a)) / (1.0 - math.exp(-2.0 * c))


def _cosh_over_sinh(numerator_argument: float, denominator_argument: float) -> float:
    """``cosh(a) / sinh(c)`` evaluated without overflowing for large arguments."""
    a, c = numerator_argument, denominator_argument
    if c == 0.0:
        raise ZeroDivisionError("cosh/sinh with a zero denominator argument")
    return math.exp(a - c) * (1.0 + math.exp(-2.0 * a)) / (1.0 - math.exp(-2.0 * c))


@dataclass(frozen=True)
class ShearLagDemand:
    """Adhesive shear-lag state of one bonded overlap at one temperature change.

    Attributes
    ----------
    delta_temperature:
        ``T - T_ref`` [K or degC increment].
    mismatch_strain:
        ``(alpha_1 - alpha_2) dT`` [-].
    transferred_force:
        Magnitude of the Milestone 1 free-joint member force that the overlap
        must transfer, ``|N_1|`` [N].
    signed_transferred_force:
        The same force with its Milestone 1 sign (tension positive) [N], used by
        the distribution functions.
    beta:
        Shear-lag parameter [1/m].
    transfer_length:
        ``1 / beta`` [m].
    dimensionless_overlap:
        ``lambda = beta L_b`` [-].
    peak_shear_stress:
        ``|tau|`` maximum over the overlap [Pa].
    peak_location:
        Coordinate of the peak [m]; the free edge ``x = L_b`` by the derivative
        argument in the module docstring.
    average_transfer_shear:
        ``|N_t| / (b L_b)`` [Pa] - the demand spread uniformly over the whole
        bonded area. Defined over the full overlap, so there is no
        factor-of-two ambiguity.
    peak_to_average_ratio:
        ``tau_peak / tau_avg = lambda coth(lambda)`` [-]. The shear
        concentration factor: 1 for a very short overlap, growing without bound
        for a long one.
    force_transfer_residual:
        ``b * integral(tau dx) - (N_1(L_b) - N_1(0))`` [N], assembled from the
        model's own closed-form expressions. Zero to machine precision.
    geometry, adhesive:
        The inputs used, retained so distributions can be evaluated.
    """

    delta_temperature: float
    mismatch_strain: float
    transferred_force: float
    signed_transferred_force: float
    beta: float
    transfer_length: float
    dimensionless_overlap: float
    peak_shear_stress: float
    peak_location: float
    average_transfer_shear: float
    peak_to_average_ratio: float
    force_transfer_residual: float
    geometry: BondedOverlapGeometry
    adhesive: AdhesiveMaterial

    def _check_coordinate(self, x: float) -> float:
        position = require_finite(x, "x")
        length = self.geometry.overlap_length
        if not (0.0 <= position <= length):
            raise ValueError(
                f"x must lie within the overlap [0, {length!r}] m, got {position!r}"
            )
        return position

    def shear_stress(self, x: float) -> float:
        """Adhesive shear stress ``tau(x)`` [Pa] at position ``x`` in [0, L_b].

        ``tau(x) = -(N_t beta / b) cosh(beta x) / sinh(beta L_b)``.
        """
        position = self._check_coordinate(x)
        beta_l = self.beta * self.geometry.overlap_length
        return (
            -(self.signed_transferred_force * self.beta / self.geometry.bond_width)
            * _cosh_over_sinh(self.beta * position, beta_l)
        )

    def member_1_force(self, x: float) -> float:
        """Adherend 1 axial force ``N_1(x)`` [N], tension positive.

        ``N_1(x) = N_t [1 - sinh(beta x) / sinh(beta L_b)]``.
        """
        position = self._check_coordinate(x)
        beta_l = self.beta * self.geometry.overlap_length
        return self.signed_transferred_force * (
            1.0 - _sinh_ratio(self.beta * position, beta_l)
        )

    def member_2_force(self, x: float) -> float:
        """Adherend 2 axial force ``N_2(x) = -N_1(x)`` [N].

        The pair is self-equilibrating everywhere: ``N_1 + N_2 = 0``.
        """
        return -self.member_1_force(x)

    def integrated_shear_force(self) -> float:
        """``b * integral of tau over the overlap`` [N], from the closed form."""
        beta_l = self.beta * self.geometry.overlap_length
        # integral of tau dx = -(N_t / (b sinh(beta L))) * sinh(beta x)
        antiderivative = (
            -(self.signed_transferred_force / self.geometry.bond_width)
            * _sinh_ratio(beta_l, beta_l)
        )
        return self.geometry.bond_width * antiderivative


def solve_shear_lag(
    member_1: AxialMember,
    member_2: AxialMember,
    geometry: BondedOverlapGeometry,
    adhesive: AdhesiveMaterial,
    delta_temperature: float,
) -> ShearLagDemand:
    """Solve the bonded overlap for the Milestone 1 free-joint mismatch demand.

    The transferred force is taken directly from
    :func:`~thermal_joint.joint.solve_bimetallic_joint` - the Milestone 1 API -
    and is not recomputed here.

    Parameters
    ----------
    member_1, member_2:
        The two bonded adherends (the Milestone 1 members).
    geometry:
        Bonded overlap geometry.
    adhesive:
        Adhesive material.
    delta_temperature:
        ``T - T_ref`` [K or degC increment].

    Returns
    -------
    ShearLagDemand
    """
    d_t = require_finite(delta_temperature, "delta_temperature")

    free_joint = solve_bimetallic_joint(member_1, member_2, d_t)
    signed_force = free_joint.member_1_force
    demand = abs(signed_force)

    beta = shear_lag_parameter(member_1, member_2, geometry, adhesive)
    length = geometry.overlap_length
    width = geometry.bond_width
    beta_l = beta * length

    peak = demand * beta / math.tanh(beta_l) / width
    average = demand / (width * length)

    result = ShearLagDemand(
        delta_temperature=d_t,
        mismatch_strain=thermal_mismatch_strain(member_1, member_2, d_t),
        transferred_force=demand,
        signed_transferred_force=signed_force,
        beta=beta,
        transfer_length=1.0 / beta,
        dimensionless_overlap=beta_l,
        peak_shear_stress=peak,
        peak_location=length,
        average_transfer_shear=average,
        peak_to_average_ratio=peak / average if average != 0.0 else 1.0,
        force_transfer_residual=0.0,
        geometry=geometry,
        adhesive=adhesive,
    )

    residual = result.integrated_shear_force() - (
        result.member_1_force(length) - result.member_1_force(0.0)
    )
    object.__setattr__(result, "force_transfer_residual", residual)
    return result

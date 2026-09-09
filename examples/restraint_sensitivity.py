"""Milestone 2 study: the bimetallic radiator joint under external axial restraint.

Sweeps the dimensionless restraint parameter ``eta_r = K_r / (E_1 A_1 + E_2 A_2)``
for the illustrative aluminium-like / titanium-like joint, reports the canonical
restrained case, locates the titanium-like member's zero-stress crossing, shows
the fully restrained asymptote and runs the inverse-design search for the
maximum tolerable restraint.

Run with::

    python examples/restraint_sensitivity.py

All material properties are ILLUSTRATIVE MATERIAL INPUT - NOT DESIGN ALLOWABLE.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from thermal_joint import (  # noqa: E402
    AxialMember,
    AxialRestraint,
    RestraintLimitStatus,
    YieldBasis,
    assess_restrained_temperature_extremes,
    assess_temperature_extremes,
    assess_yield,
    joint_axial_stiffness,
    maximum_allowable_restraint_stiffness,
    rigid_restraint_stresses,
    solve_restrained_joint,
    zero_stress_restraint_stiffness,
)
from thermal_joint.illustrative import (  # noqa: E402
    ILLUSTRATIVE_NOTE,
    SQUARE_MILLIMETRE,
    aluminium_like_member,
    radiator_joint_environment,
    titanium_like_member,
)

MPA = 1.0e6
GPA = 1.0e9
MICROSTRAIN = 1.0e-6
DESIGN_FACTOR = 1.25
CANONICAL_RATIO = 1.0
AREA_MM2 = 100.0
RATIOS = (0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0)
RESTRAINT_LABEL = "illustrative surrounding-structure axial restraint"


def rule(title: str = "") -> None:
    if title:
        print()
        print(title)
    print("-" * 96)


def format_margin(margin: float) -> str:
    return "inf" if math.isinf(margin) else f"{margin:+.4f}"


def print_study_basis(member_1, member_2, environment, basis, k_joint) -> None:
    rule("STUDY BASIS")
    print(f"  {ILLUSTRATIVE_NOTE}")
    print(f"  Restraint is an {RESTRAINT_LABEL}, not a measured spacecraft stiffness.")
    print()
    print(f"  member 1: {member_1.label}")
    print(f"  member 2: {member_2.label}")
    print()
    print(f"  {'Property':<34}{'Member 1':>20}{'Member 2':>20}")
    for caption, value_1, value_2 in (
        ("Young's modulus E [GPa]",
         member_1.material.elastic_modulus / GPA, member_2.material.elastic_modulus / GPA),
        ("CTE alpha [1e-6 / K]",
         member_1.material.thermal_expansion_coefficient / MICROSTRAIN,
         member_2.material.thermal_expansion_coefficient / MICROSTRAIN),
        ("yield strength [MPa]",
         member_1.material.yield_strength / MPA, member_2.material.yield_strength / MPA),
        ("cross-sectional area [mm^2]",
         member_1.area / SQUARE_MILLIMETRE, member_2.area / SQUARE_MILLIMETRE),
        ("axial stiffness E*A [MN]",
         member_1.axial_stiffness / 1.0e6, member_2.axial_stiffness / 1.0e6),
    ):
        print(f"  {caption:<34}{value_1:>20.4g}{value_2:>20.4g}")

    print()
    print(f"  reference temperature T_ref       {environment.reference_temperature:>14.1f} degC")
    print(f"  cold temperature      T_cold      {environment.cold_temperature:>14.1f} degC")
    print(f"  hot temperature       T_hot       {environment.hot_temperature:>14.1f} degC")
    print(f"  cold excursion        dT_cold     {environment.delta_temperature_cold:>14.1f} K")
    print(f"  hot excursion         dT_hot      {environment.delta_temperature_hot:>14.1f} K")
    print(f"  yield design factor                {basis.design_factor:>14.2f}")
    print()
    print(f"  joint axial stiffness K_joint      {k_joint:>14.4e} N   (= E1 A1 + E2 A2)")
    print(f"  canonical restraint    eta_r       {CANONICAL_RATIO:>14.2f}")
    print(f"  canonical restraint    K_r         {CANONICAL_RATIO * k_joint:>14.4e} N")
    print(f"  restraint reference strain eps_ref {0.0:>14.4e}")
    print()
    print("  Sign convention: tension positive. The restraint force N_r is the force")
    print("  exerted ON the joint, N_r = K_r (eps_ref - eps_common), and joint")
    print("  equilibrium is N_1 + N_2 - N_r = 0.")


def print_free_reference(member_1, member_2, environment, basis) -> None:
    rule("FREE-JOINT REFERENCE (Milestone 1, eta_r = 0)")
    free = assess_temperature_extremes(member_1, member_2, environment, basis)
    for tag, result in (("hot ", free.hot_result), ("cold", free.cold_result)):
        print(
            f"  {tag}  dT = {result.delta_temperature:+7.1f} K"
            f"   sigma_1 = {result.member_1_stress / MPA:>9.4f} MPa"
            f"   sigma_2 = {result.member_2_stress / MPA:>9.4f} MPa"
            f"   (self-equilibrating: N_1 + N_2 = 0)"
        )
    print()
    print(f"  governing extreme                  {free.governing_extreme:>14}")
    print(f"  governing member                   {free.governing_member:>14d}")
    print(f"  minimum preliminary yield margin   {format_margin(free.minimum_yield_margin):>14}")
    print(f"  verdict                            {'PASS' if free.passes else 'FAIL':>14}")


def print_case(title: str, result, assessment) -> None:
    rule(title)
    print(f"  temperature change dT                {result.delta_temperature:>16.4f} K")
    print(f"  restraint stiffness K_r              {result.restraint_stiffness:>16.4e} N")
    print(f"  restraint reference strain           {result.restraint_reference_strain:>16.4e}")
    print(
        "  free-joint common strain (eta_r=0)   "
        f"{result.free_joint_common_strain / MICROSTRAIN:>16.4f} microstrain"
    )
    print(
        "  restrained common strain             "
        f"{result.common_strain / MICROSTRAIN:>16.4f} microstrain"
    )
    for index, (stress, force) in enumerate(zip(result.stresses, result.forces), start=1):
        state = "tension" if stress > 0 else "compression" if stress < 0 else "unstressed"
        print(f"  stress, member {index}                     {stress / MPA:>16.4f} MPa  ({state})")
        print(f"  internal force, member {index}             {force:>16.4f} N")
    print(f"  restraint force on joint N_r         {result.restraint_force:>16.4f} N")
    print(
        "  member force sum N_1 + N_2           "
        f"{result.member_1_force + result.member_2_force:>16.4f} N"
    )
    print(
        "  equilibrium residual N_1+N_2-N_r     "
        f"{result.force_equilibrium_residual:>16.3e} N"
    )
    print(
        "  both members same sign?              "
        f"{str(result.members_carry_same_sign_stress):>16}"
    )
    print()
    for index, margin in ((1, assessment.member_1_margin), (2, assessment.member_2_margin)):
        print(
            f"  member {index} yield margin               "
            f"{format_margin(margin.margin_of_safety):>16}"
            f"   {'PASS' if margin.passes else 'FAIL'}"
        )


def print_sensitivity(member_1, member_2, environment, basis, k_joint) -> None:
    rule("RESTRAINT SENSITIVITY")
    print(
        f"  {'eta_r':>7}{'K_r [N]':>13}{'hot s1':>11}{'hot s2':>11}"
        f"{'cold s1':>11}{'cold s2':>11}{'min MS':>11}{'governs':>16}{'':>7}"
    )
    print(f"  {'':>7}{'':>13}{'[MPa]':>11}{'[MPa]':>11}{'[MPa]':>11}{'[MPa]':>11}")
    for ratio in RATIOS:
        restraint = AxialRestraint.from_stiffness_ratio(ratio, member_1, member_2)
        assessment = assess_restrained_temperature_extremes(
            member_1, member_2, environment, restraint, basis
        )
        governs = f"{assessment.governing_extreme}/m{assessment.governing_member}"
        print(
            f"  {ratio:>7.2f}{restraint.stiffness:>13.3e}"
            f"{assessment.hot_result.member_1_stress / MPA:>11.3f}"
            f"{assessment.hot_result.member_2_stress / MPA:>11.3f}"
            f"{assessment.cold_result.member_1_stress / MPA:>11.3f}"
            f"{assessment.cold_result.member_2_stress / MPA:>11.3f}"
            f"{format_margin(assessment.minimum_yield_margin):>11}"
            f"{governs:>16}"
            f"{'  PASS' if assessment.passes else '  FAIL':>7}"
        )
    print()
    print("  Common strain magnitude falls monotonically toward the reference strain as")
    print("  restraint rises. Member 1 (high CTE) is loaded monotonically harder, but")
    print("  member 2 is NOT monotonic: its stress falls through zero and grows again")
    print("  with the opposite sign.")


def print_sign_transition(member_1, member_2, k_joint) -> None:
    rule("SIGN TRANSITION")
    for index, name in ((1, "member 1 (aluminium-like)"), (2, "member 2 (titanium-like)")):
        crossing = zero_stress_restraint_stiffness(member_1, member_2, 100.0, index)
        if crossing is None:
            print(f"  {name}: no admissible K_r >= 0 gives zero stress - sign never changes.")
        else:
            print(
                f"  {name}: crosses zero stress at K_r = {crossing:.4e} N"
                f"  (eta_r = {crossing / k_joint:.6f})"
            )
    print()
    crossing = zero_stress_restraint_stiffness(member_1, member_2, 100.0, 2)
    print("  For eps_ref = 0 the crossing is independent of dT: K_r = S/alpha_2 - K_joint")
    print("  with S = E1 A1 alpha_1 + E2 A2 alpha_2. Verified either side of it:")
    for factor, tag in ((0.9, "below"), (1.1, "above")):
        for delta_t, label in ((100.0, "hot "), (-140.0, "cold")):
            result = solve_restrained_joint(
                member_1, member_2, delta_t, AxialRestraint(factor * crossing)
            )
            print(
                f"    {label} eta_r = {factor * crossing / k_joint:.4f} ({tag:<5})"
                f"  sigma_2 = {result.member_2_stress / MPA:>9.4f} MPa"
            )


def print_rigid_limit(member_1, member_2, environment, basis) -> None:
    rule("RIGID-RESTRAINT LIMIT (K_r -> infinity, eps_ref = 0)")
    print("  sigma_i = E_i (eps_ref - alpha_i dT) = -E_i alpha_i dT")
    print()
    worst = math.inf
    for delta_t, label in (
        (environment.delta_temperature_hot, "hot "),
        (environment.delta_temperature_cold, "cold"),
    ):
        stress_1, stress_2 = rigid_restraint_stresses(member_1, member_2, delta_t)
        margin_1 = basis.margin_of_safety(member_1.material, stress_1)
        margin_2 = basis.margin_of_safety(member_2.material, stress_2)
        worst = min(worst, margin_1, margin_2)
        print(
            f"  {label}  dT = {delta_t:+7.1f} K"
            f"   sigma_1 = {stress_1 / MPA:>9.4f} MPa (MS {format_margin(margin_1)})"
            f"   sigma_2 = {stress_2 / MPA:>9.4f} MPa (MS {format_margin(margin_2)})"
        )
    print()
    print(f"  fully restrained minimum margin      {format_margin(worst):>16}"
          f"   {'PASS' if worst >= 0 else 'FAIL'}")
    print("  Both members take the SAME sign here - compression on heating, tension on")
    print("  cooling - unlike the self-equilibrating free joint.")


def print_inverse_result(member_1, member_2, environment, basis, k_joint) -> None:
    rule("INVERSE RESTRAINT RESULT")
    result = maximum_allowable_restraint_stiffness(member_1, member_2, environment, basis)
    print(f"  status                               {result.status.value:>16}")
    print(f"  free-joint minimum margin            {format_margin(result.free_joint_minimum_margin):>16}")
    print(f"  fully restrained minimum margin      {format_margin(result.rigid_restraint_minimum_margin):>16}")
    print(f"  search upper bound on eta_r          {result.search_upper_ratio:>16.4e}")
    print(f"  tolerance / iterations               {result.tolerance:>10.1e} / {result.iterations:<4d}")
    print()
    if result.status is RestraintLimitStatus.FINITE_LIMIT:
        print(f"  maximum tolerable eta_r              {result.maximum_stiffness_ratio:>16.6f}")
        print(f"  maximum tolerable K_r                {result.maximum_restraint_stiffness:>16.4e} N")
        print(f"  minimum margin at that limit         {result.minimum_margin_at_limit:>16.3e}")
        print()
        print("  Above this restraint stiffness the joint no longer meets the preliminary")
        print("  elastic yield criterion at its governing temperature extreme.")
    elif result.status is RestraintLimitStatus.NO_FINITE_LIMIT_WITHIN_MODEL:
        print("  Even a fully rigid restraint stays below yield, so no finite maximum")
        print("  restraint stiffness exists within this linear-elastic model.")
    elif result.status is RestraintLimitStatus.FREE_JOINT_ALREADY_FAILS:
        print("  The unrestrained joint already fails yield; no restraint limit applies.")
    else:
        print("  A boundary exists asymptotically but lies above the requested search")
        print("  bound. Re-run with a larger upper_stiffness_ratio; bounds are never")
        print("  expanded silently.")


def main() -> None:
    member_1: AxialMember = aluminium_like_member(AREA_MM2)
    member_2: AxialMember = titanium_like_member(AREA_MM2)
    environment = radiator_joint_environment()
    basis = YieldBasis(DESIGN_FACTOR)
    k_joint = joint_axial_stiffness(member_1, member_2)

    print()
    print("=" * 96)
    print("  BIMETALLIC JOINT UNDER EXTERNAL AXIAL RESTRAINT - MILESTONE 2 STUDY")
    print("=" * 96)

    print_study_basis(member_1, member_2, environment, basis, k_joint)
    print_free_reference(member_1, member_2, environment, basis)

    canonical = AxialRestraint.from_stiffness_ratio(
        CANONICAL_RATIO, member_1, member_2, label=RESTRAINT_LABEL
    )
    assessment = assess_restrained_temperature_extremes(
        member_1, member_2, environment, canonical, basis
    )
    print_case(
        f"CANONICAL RESTRAINED CASE - HOT (eta_r = {CANONICAL_RATIO:g})",
        assessment.hot_result,
        assessment.hot_assessment,
    )
    print_case(
        f"CANONICAL RESTRAINED CASE - COLD (eta_r = {CANONICAL_RATIO:g})",
        assessment.cold_result,
        assessment.cold_assessment,
    )

    rule("CANONICAL GOVERNING RESULT")
    print(f"  governing extreme                    {assessment.governing_extreme:>16}")
    print(f"  governing member                     {assessment.governing_member:>16d}")
    print(f"  minimum preliminary yield margin     {format_margin(assessment.minimum_yield_margin):>16}")
    print(f"  verdict                              {'PASS' if assessment.passes else 'FAIL':>16}")

    print_sensitivity(member_1, member_2, environment, basis, k_joint)
    print_sign_transition(member_1, member_2, k_joint)
    print_rigid_limit(member_1, member_2, environment, basis)
    print_inverse_result(member_1, member_2, environment, basis, k_joint)

    rule()
    print(
        "  External restraint can qualitatively change the bimetallic stress state: the\n"
        "  free joint carries self-equilibrating tension/compression, while strong\n"
        "  restraint drives both members toward the same-sign fully restrained thermal\n"
        "  stress state."
    )
    print()
    print(
        "  Perfectly bonded, uniform-temperature axial model with a single linear\n"
        "  restraint; interface shear, finite joint length, thermal gradients,\n"
        "  plasticity and fatigue are not included."
    )
    print(f"  {ILLUSTRATIVE_NOTE}.")
    print()


if __name__ == "__main__":
    main()

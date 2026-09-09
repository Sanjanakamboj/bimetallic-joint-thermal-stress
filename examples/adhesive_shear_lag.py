"""Milestone 3 study: adhesive shear-lag over the bimetallic radiator joint overlap.

Takes the Milestone 1 free-joint CTE-mismatch force and asks what adhesive shear
develops when that force has to be transferred through a finite bonded overlap.
Reports the distributions, the hot/cold governing case, four sensitivity sweeps,
the allowable temperature excursion and the minimum-overlap inverse design.

Run with::

    python examples/adhesive_shear_lag.py

All material properties are ILLUSTRATIVE - NOT DESIGN ALLOWABLE.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from thermal_joint import (  # noqa: E402
    AdhesiveMaterial,
    AdhesiveShearBasis,
    AxialMember,
    OverlapLimitStatus,
    YieldBasis,
    allowable_temperature_change_for_adhesive_shear,
    assess_shear_lag_extremes,
    assess_temperature_extremes,
    long_overlap_peak_shear_stress,
    required_overlap_length,
    screen_joint,
    shear_lag_parameter,
    solve_bimetallic_joint,
    solve_shear_lag,
    transfer_length,
)
from thermal_joint.illustrative import (  # noqa: E402
    ADHESIVE_LIKE,
    ALUMINIUM_LIKE,
    ILLUSTRATIVE_ADHESIVE_NOTE,
    ILLUSTRATIVE_NOTE,
    MILLIMETRE,
    SQUARE_MILLIMETRE,
    aluminium_like_member,
    radiator_joint_environment,
    radiator_joint_overlap,
    titanium_like_member,
)

MPA = 1.0e6
GPA = 1.0e9
MICROSTRAIN = 1.0e-6
YIELD_DESIGN_FACTOR = 1.25
ADHESIVE_DESIGN_FACTOR = 1.25
AREA_MM2 = 100.0

OVERLAP_LENGTHS_MM = (5.0, 10.0, 20.0, 40.0, 80.0, 160.0)
THICKNESSES_MM = (0.05, 0.10, 0.20, 0.50, 1.00)
MODULI_GPA = (0.1, 0.25, 0.5, 1.0, 2.0, 5.0)
AREA_RATIOS = (0.25, 0.5, 1.0, 2.0, 4.0)


def rule(title: str = "") -> None:
    if title:
        print()
        print(title)
    print("-" * 96)


def format_margin(margin: float) -> str:
    return "inf" if math.isinf(margin) else f"{margin:+.4f}"


def verdict(margin: float) -> str:
    return "PASS" if margin >= 0.0 else "FAIL"


def print_study_basis(member_1, member_2, environment, geometry, adhesive, basis) -> None:
    rule("STUDY BASIS")
    print(f"  {ILLUSTRATIVE_NOTE}")
    print(f"  {ILLUSTRATIVE_ADHESIVE_NOTE}")
    print(f"  adhesive provenance: {adhesive.source_note}")
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
    print(f"  reference temperature T_ref        {environment.reference_temperature:>14.1f} degC")
    print(f"  cold temperature      T_cold       {environment.cold_temperature:>14.1f} degC")
    print(f"  hot temperature       T_hot        {environment.hot_temperature:>14.1f} degC")
    print(f"  cold excursion        dT_cold      {environment.delta_temperature_cold:>14.1f} K")
    print(f"  hot excursion         dT_hot       {environment.delta_temperature_hot:>14.1f} K")
    print()
    print(f"  overlap length        L_b          {geometry.overlap_length / MILLIMETRE:>14.2f} mm")
    print(f"  bond width            b            {geometry.bond_width / MILLIMETRE:>14.2f} mm")
    print(f"  adhesive thickness    t_a          {geometry.adhesive_thickness / MILLIMETRE:>14.3f} mm")
    print(f"  bonded area           b*L_b        {geometry.bond_area / SQUARE_MILLIMETRE:>14.1f} mm^2")
    print(f"  adhesive shear modulus G_a         {adhesive.shear_modulus / GPA:>14.3f} GPa")
    print(f"  adhesive shear strength            {adhesive.shear_strength / MPA:>14.2f} MPa")
    print(f"  adhesive design factor             {basis.design_factor:>14.2f}")
    print(f"  adhesive allowable tau_allow       {basis.allowable_shear_stress(adhesive) / MPA:>14.2f} MPa")
    print(f"  member yield design factor         {YIELD_DESIGN_FACTOR:>14.2f}")


def print_free_joint_reference(member_1, member_2, environment) -> None:
    rule("M1 FREE-JOINT REFERENCE (the demand this overlap must transfer)")
    for label, delta_t in (
        ("hot ", environment.delta_temperature_hot),
        ("cold", environment.delta_temperature_cold),
    ):
        result = solve_bimetallic_joint(member_1, member_2, delta_t)
        print(
            f"  {label}  dT = {delta_t:+7.1f} K"
            f"   sigma_1 = {result.member_1_stress / MPA:>9.4f} MPa"
            f"   sigma_2 = {result.member_2_stress / MPA:>9.4f} MPa"
            f"   N_transfer = |N_1| = {abs(result.member_1_force):>9.4f} N"
        )
    print()
    print("  The pair is self-equilibrating (N_1 + N_2 = 0), so the overlap has to")
    print("  transfer the magnitude |N_1| from one adherend into the other.")


def print_parameters(member_1, member_2, geometry, adhesive) -> None:
    rule("SHEAR-LAG PARAMETERS")
    beta = shear_lag_parameter(member_1, member_2, geometry, adhesive)
    length = transfer_length(member_1, member_2, geometry, adhesive)
    print("  beta^2 = (G_a b / t_a) (1/(E1 A1) + 1/(E2 A2))")
    print(f"  beta                               {beta:>14.4f} 1/m")
    print(f"  transfer length 1/beta             {length / MILLIMETRE:>14.4f} mm")
    print(f"  dimensionless overlap lambda=beta*L_b {beta * geometry.overlap_length:>11.4f}")
    print(f"  overlap spans                      {beta * geometry.overlap_length:>14.2f} transfer lengths")
    print()
    print("  Domain x in [0, L_b]:  x=0 is the loaded transfer plane where the")
    print("  mismatch force pair is fully developed; x=L_b is the free edge where")
    print("  both adherend forces vanish. |tau| ~ cosh(beta x) is strictly")
    print("  increasing, so the peak shear is at the FREE EDGE x = L_b.")


def print_case(title: str, demand, basis) -> None:
    rule(title)
    margin = basis.margin_of_safety(demand.adhesive, demand.peak_shear_stress)
    print(f"  temperature change dT                {demand.delta_temperature:>16.4f} K")
    print(f"  mismatch strain (a1-a2)dT            {demand.mismatch_strain / MICROSTRAIN:>16.4f} microstrain")
    print(f"  transferred force N_transfer         {demand.transferred_force:>16.4f} N")
    print(f"  peak adhesive shear tau_peak         {demand.peak_shear_stress / MPA:>16.4f} MPa")
    print(f"  peak location (free edge)            {demand.peak_location / MILLIMETRE:>16.4f} mm")
    print(f"  average transfer shear tau_avg       {demand.average_transfer_shear / MPA:>16.4f} MPa")
    print(f"  peak / average                       {demand.peak_to_average_ratio:>16.4f}")
    print(f"  force-transfer residual              {demand.force_transfer_residual:>16.3e} N")
    print(f"  adhesive allowable tau_allow         {basis.allowable_shear_stress(demand.adhesive) / MPA:>16.4f} MPa")
    print(f"  preliminary adhesive shear margin    {format_margin(margin):>16}   {verdict(margin)}")
    print()
    print(f"  {'x [mm]':>9}{'tau [MPa]':>13}{'N_1 [N]':>13}{'N_2 [N]':>13}{'N_1+N_2 [N]':>14}")
    for index in range(6):
        x = demand.geometry.overlap_length * index / 5
        force_1 = demand.member_1_force(x)
        force_2 = demand.member_2_force(x)
        print(
            f"  {x / MILLIMETRE:>9.2f}{demand.shear_stress(x) / MPA:>13.4f}"
            f"{force_1:>13.3f}{force_2:>13.3f}{force_1 + force_2:>14.2e}"
        )


def print_governing(assessment) -> None:
    rule("GOVERNING EXTREME")
    print(f"  governing extreme                    {assessment.governing_extreme:>16}")
    print(f"  peak adhesive shear                  {assessment.governing_peak_shear_stress / MPA:>16.4f} MPa")
    print(f"  adhesive allowable                   {assessment.allowable_shear_stress / MPA:>16.4f} MPa")
    print(f"  minimum adhesive shear margin        {format_margin(assessment.minimum_margin):>16}")
    print(f"  verdict                              {verdict(assessment.minimum_margin):>16}")
    print()
    print("  Cold governs on magnitude: |dT_cold| = 140 K against |dT_hot| = 100 K,")
    print("  and peak shear is exactly linear in dT, so the ratio is 140/100 = 1.4.")


def print_overlap_sensitivity(member_1, member_2, geometry, adhesive, basis, delta_t) -> None:
    rule(f"OVERLAP SENSITIVITY (governing extreme, dT = {delta_t:+.1f} K)")
    print(
        f"  {'L_b [mm]':>10}{'beta*L_b':>11}{'tau_peak':>12}{'tau_avg':>11}"
        f"{'peak/avg':>11}{'margin':>11}{'':>7}"
    )
    print(f"  {'':>10}{'':>11}{'[MPa]':>12}{'[MPa]':>11}")
    for length_mm in OVERLAP_LENGTHS_MM:
        demand = solve_shear_lag(
            member_1, member_2, geometry.with_overlap_length(length_mm * MILLIMETRE),
            adhesive, delta_t,
        )
        margin = basis.margin_of_safety(adhesive, demand.peak_shear_stress)
        print(
            f"  {length_mm:>10.1f}{demand.dimensionless_overlap:>11.4f}"
            f"{demand.peak_shear_stress / MPA:>12.4f}{demand.average_transfer_shear / MPA:>11.4f}"
            f"{demand.peak_to_average_ratio:>11.4f}{format_margin(margin):>11}"
            f"{'  ' + verdict(margin):>7}"
        )
    asymptote = long_overlap_peak_shear_stress(
        solve_shear_lag(member_1, member_2, geometry, adhesive, delta_t).transferred_force,
        shear_lag_parameter(member_1, member_2, geometry, adhesive),
        geometry.bond_width,
    )
    print()
    print(f"  L_b -> infinity asymptote            {asymptote / MPA:>16.4f} MPa")
    print("  Average shear falls as 1/L_b, but the peak converges to that finite")
    print("  asymptote FROM ABOVE. Beyond a few transfer lengths, more overlap buys")
    print("  essentially nothing in peak shear - it cannot go below the asymptote.")


def print_thickness_sensitivity(member_1, member_2, geometry, adhesive, basis, delta_t) -> None:
    rule(f"ADHESIVE THICKNESS SENSITIVITY (G_a fixed, dT = {delta_t:+.1f} K)")
    print(f"  {'t_a [mm]':>10}{'beta [1/m]':>13}{'1/beta [mm]':>14}{'tau_peak [MPa]':>17}{'margin':>11}{'':>7}")
    for thickness_mm in THICKNESSES_MM:
        local = geometry.with_adhesive_thickness(thickness_mm * MILLIMETRE)
        demand = solve_shear_lag(member_1, member_2, local, adhesive, delta_t)
        margin = basis.margin_of_safety(adhesive, demand.peak_shear_stress)
        print(
            f"  {thickness_mm:>10.2f}{demand.beta:>13.3f}{demand.transfer_length / MILLIMETRE:>14.4f}"
            f"{demand.peak_shear_stress / MPA:>17.4f}{format_margin(margin):>11}"
            f"{'  ' + verdict(margin):>7}"
        )
    print()
    print("  A thicker bondline lowers beta, lengthens the transfer zone and reduces")
    print("  peak shear - computed, not assumed.")


def print_modulus_sensitivity(member_1, member_2, geometry, adhesive, basis, delta_t) -> None:
    rule(f"ADHESIVE MODULUS SENSITIVITY (t_a fixed, dT = {delta_t:+.1f} K)")
    print(f"  {'G_a [GPa]':>11}{'beta [1/m]':>13}{'1/beta [mm]':>14}{'tau_peak [MPa]':>17}{'margin':>11}{'':>7}")
    for modulus_gpa in MODULI_GPA:
        local = AdhesiveMaterial(
            f"G={modulus_gpa} GPa", modulus_gpa * GPA, adhesive.shear_strength,
            "illustrative sweep variant of the canonical adhesive",
        )
        demand = solve_shear_lag(member_1, member_2, geometry, local, delta_t)
        margin = basis.margin_of_safety(local, demand.peak_shear_stress)
        print(
            f"  {modulus_gpa:>11.2f}{demand.beta:>13.3f}{demand.transfer_length / MILLIMETRE:>14.4f}"
            f"{demand.peak_shear_stress / MPA:>17.4f}{format_margin(margin):>11}"
            f"{'  ' + verdict(margin):>7}"
        )
    print()
    print("  A stiffer adhesive raises beta, shortens the transfer zone and raises")
    print("  peak shear. Note that a softer adhesive is not automatically a design")
    print("  improvement: peel, creep, durability and joint deformation are all")
    print("  excluded from this screen.")


def print_area_ratio_sensitivity(member_2, geometry, adhesive, basis, delta_t) -> None:
    rule(f"AREA-RATIO SENSITIVITY (recomputes M1 end to end, dT = {delta_t:+.1f} K)")
    print(
        f"  {'A_Al/A_Ti':>11}{'A_Al [mm^2]':>13}{'M1 |N_1| [N]':>14}{'beta [1/m]':>13}"
        f"{'tau_peak [MPa]':>17}{'margin':>11}{'':>7}"
    )
    for ratio in AREA_RATIOS:
        member_1 = AxialMember(
            ALUMINIUM_LIKE, ratio * AREA_MM2 * SQUARE_MILLIMETRE, label="aluminium-like"
        )
        demand = solve_shear_lag(member_1, member_2, geometry, adhesive, delta_t)
        margin = basis.margin_of_safety(adhesive, demand.peak_shear_stress)
        print(
            f"  {ratio:>11.2f}{member_1.area / SQUARE_MILLIMETRE:>13.1f}"
            f"{demand.transferred_force:>14.3f}{demand.beta:>13.3f}"
            f"{demand.peak_shear_stress / MPA:>17.4f}{format_margin(margin):>11}"
            f"{'  ' + verdict(margin):>7}"
        )
    print()
    print("  The mismatch force is NOT held fixed here: a larger aluminium-like area")
    print("  raises the Milestone 1 demand while lowering beta, and the peak shear")
    print("  rises on balance.")


def print_allowable_delta_t(member_1, member_2, geometry, adhesive, basis) -> None:
    rule("ALLOWABLE DELTA-T")
    allowable = allowable_temperature_change_for_adhesive_shear(
        member_1, member_2, geometry, adhesive, basis
    )
    print("  tau_peak is exactly linear in dT, so this is a direct closed form:")
    print("  |dT|_allow = tau_allow b C / (|a1-a2| beta coth(beta L_b))")
    print()
    print(f"  allowable |dT| at this geometry      {allowable:>16.4f} K")
    check = solve_shear_lag(member_1, member_2, geometry, adhesive, allowable)
    print(f"  peak shear at that dT                {check.peak_shear_stress / MPA:>16.4f} MPa")
    print(f"  margin there                         {format_margin(basis.margin_of_safety(adhesive, check.peak_shear_stress)):>16}")
    print()
    print(f"  For reference the study excursions are -140 K and +100 K, so this")
    print(f"  bondline is limited to roughly +/-{allowable:.0f} K rather than the full")
    print("  excursion.")


def print_required_overlap(member_1, member_2, environment, geometry, adhesive, basis) -> None:
    rule("REQUIRED OVERLAP")
    result = required_overlap_length(
        member_1, member_2, environment, geometry, adhesive, basis
    )
    print(f"  status                               {result.status.value:>16}")
    print(f"  governing extreme                    {result.governing_extreme:>16}")
    print(f"  adhesive allowable                   {result.allowable_shear_stress / MPA:>16.4f} MPa")
    print(f"  peak shear at L_min = {result.minimum_overlap_length / MILLIMETRE:>5.1f} mm        {result.peak_shear_at_lower_bound / MPA:>16.4f} MPa")
    print(f"  peak shear at L_max = {result.maximum_overlap_length / MILLIMETRE:>5.1f} mm        {result.peak_shear_at_upper_bound / MPA:>16.4f} MPa")
    print(f"  long-overlap asymptote               {result.long_overlap_asymptote / MPA:>16.4f} MPa")
    print()
    if result.status is OverlapLimitStatus.FINITE_REQUIRED_LENGTH:
        print(f"  minimum required overlap             {result.required_overlap_length / MILLIMETRE:>16.4f} mm")
        print(f"  margin there                         {format_margin(result.margin_at_required_length):>16}")
    elif result.status is OverlapLimitStatus.LOWER_BOUND_ALREADY_PASSES:
        print("  The shortest overlap in the search range already meets the allowable;")
        print("  adhesive shear does not size this overlap.")
    elif result.status is OverlapLimitStatus.NO_FINITE_LENGTH_WITHIN_MODEL:
        print("  No finite overlap length can satisfy the allowable. Peak shear falls")
        print("  toward the long-overlap asymptote but never below it, and that")
        print("  asymptote already exceeds tau_allow. Adding bonded length is the")
        print("  wrong lever: the fix has to be a more compliant or thicker bondline,")
        print("  a wider bond, a lower CTE mismatch, or a smaller temperature")
        print("  excursion.")
    else:
        print("  A finite length exists but lies above the requested search bound.")
        print("  Re-run with a larger maximum_overlap_length; bounds are never")
        print("  expanded silently.")


def print_combined_screen(member_1, member_2, environment, geometry, adhesive, basis) -> None:
    rule("COMBINED PRELIMINARY SCREEN (global yield AND local adhesive shear)")
    screen = screen_joint(
        member_1, member_2, environment, geometry, adhesive,
        YieldBasis(YIELD_DESIGN_FACTOR), basis,
    )
    print(f"  minimum member yield margin          {format_margin(screen.minimum_yield_margin):>16}"
          f"   {verdict(screen.minimum_yield_margin)}")
    print(f"  minimum adhesive shear margin        {format_margin(screen.minimum_adhesive_margin):>16}"
          f"   {verdict(screen.minimum_adhesive_margin)}")
    print(f"  overall preliminary feasibility      {str(screen.overall_feasible):>16}")
    print(f"  failing screen(s)                    {screen.governing_mode:>16}")
    print()
    print("  The two margins are reported side by side and never combined into one")
    print("  number: they measure different failure modes against different")
    print("  allowables. Only the booleans are ANDed.")


def main() -> None:
    member_1 = aluminium_like_member(AREA_MM2)
    member_2 = titanium_like_member(AREA_MM2)
    environment = radiator_joint_environment()
    geometry = radiator_joint_overlap()
    adhesive = ADHESIVE_LIKE
    basis = AdhesiveShearBasis(ADHESIVE_DESIGN_FACTOR)

    print()
    print("=" * 96)
    print("  ADHESIVE SHEAR-LAG OVER A FINITE BONDED OVERLAP - MILESTONE 3 STUDY")
    print("=" * 96)

    print_study_basis(member_1, member_2, environment, geometry, adhesive, basis)
    print_free_joint_reference(member_1, member_2, environment)
    print_parameters(member_1, member_2, geometry, adhesive)

    assessment = assess_shear_lag_extremes(
        member_1, member_2, environment, geometry, adhesive, basis
    )
    print_case("HOT CASE", assessment.hot_demand, basis)
    print_case("COLD CASE", assessment.cold_demand, basis)
    print_governing(assessment)

    governing_delta_t = assessment.governing_demand.delta_temperature
    print_overlap_sensitivity(member_1, member_2, geometry, adhesive, basis, governing_delta_t)
    print_thickness_sensitivity(member_1, member_2, geometry, adhesive, basis, governing_delta_t)
    print_modulus_sensitivity(member_1, member_2, geometry, adhesive, basis, governing_delta_t)
    print_area_ratio_sensitivity(member_2, geometry, adhesive, basis, governing_delta_t)
    print_allowable_delta_t(member_1, member_2, geometry, adhesive, basis)
    print_required_overlap(member_1, member_2, environment, geometry, adhesive, basis)
    print_combined_screen(member_1, member_2, environment, geometry, adhesive, basis)

    rule()
    print(
        "This is a one-dimensional linear-elastic shear-lag screen. It resolves axial "
        "load transfer and adhesive shear only; peel stress, bending, edge "
        "singularities, nonlinear adhesive behavior and fracture are excluded."
    )
    print()


if __name__ == "__main__":
    main()

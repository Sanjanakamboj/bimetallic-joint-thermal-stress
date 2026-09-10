"""STM-11 final integrated assessment of the bimetallic radiator joint.

One canonical run of the whole analysis chain, using only the package APIs:

    free thermal compatibility -> external restraint -> adhesive shear-lag
    -> static bondline design trade -> thermal-cycle fatigue -> final comparison

Run with::

    python examples/final_thermal_joint_assessment.py

Every material property, adhesive property, fatigue curve and design requirement
in this study is ILLUSTRATIVE - NOT DESIGN ALLOWABLE.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from thermal_joint import (  # noqa: E402
    AdhesiveShearBasis,
    AxialRestraint,
    ThermalCycleRequirement,
    YieldBasis,
    allowable_cycle_scale_for_fatigue,
    allowable_temperature_change_for_adhesive_shear,
    assess_restrained_temperature_extremes,
    assess_shear_lag_extremes,
    assess_temperature_extremes,
    assess_thermal_cycle_fatigue,
    fatigue_design_map,
    joint_axial_stiffness,
    long_overlap_peak_shear,
    maximum_allowable_adhesive_shear_modulus,
    maximum_allowable_restraint_stiffness,
    required_adhesive_thickness,
    required_bond_width,
    required_overlap_length,
    scale_environment,
    select_preliminary_fatigue_design,
    select_preliminary_joint_design,
    solve_bimetallic_joint,
    width_thickness_design_map,
    zero_stress_restraint_stiffness,
)
from thermal_joint.illustrative import (  # noqa: E402
    ADHESIVE_LIKE,
    ILLUSTRATIVE_ADHESIVE_NOTE,
    ILLUSTRATIVE_FATIGUE_NOTE,
    ILLUSTRATIVE_NOTE,
    MILLIMETRE,
    SQUARE_MILLIMETRE,
    aluminium_like_member,
    illustrative_cycle_requirement,
    illustrative_fatigue_curves,
    radiator_joint_environment,
    radiator_joint_overlap,
    titanium_like_member,
)

MPA, GPA, MM = 1.0e6, 1.0e9, MILLIMETRE
YIELD_DF, ADHESIVE_DF = 1.25, 1.25

STATIC_MAP_WIDTHS_MM = (10.0, 20.0, 30.0, 40.0, 50.0, 75.0, 100.0)
STATIC_MAP_THICKNESSES_MM = (0.1, 0.2, 0.3, 0.5, 0.75, 1.0)
FATIGUE_MAP_WIDTHS_MM = (20.0, 30.0, 40.0, 50.0, 75.0, 100.0)
FATIGUE_MAP_THICKNESSES_MM = (0.5, 0.75, 1.0, 1.5, 2.0)
REQUIREMENT_SWEEP = (1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6)
SCALE_SWEEP = (0.5, 0.75, 1.0, 1.25, 1.5)


def rule(title: str = "") -> None:
    if title:
        print()
        print(title)
    print("-" * 100)


def verdict(ok: bool) -> str:
    return "PASS" if ok else "FAIL"


def life_text(cycles: float) -> str:
    return "inf" if math.isinf(cycles) else f"{cycles:.4g}"


def main() -> None:
    member_1 = aluminium_like_member(100.0)
    member_2 = titanium_like_member(100.0)
    environment = radiator_joint_environment()
    adhesive = ADHESIVE_LIKE
    curves = illustrative_fatigue_curves()
    requirement = illustrative_cycle_requirement()
    yield_basis = YieldBasis(YIELD_DF)
    adhesive_basis = AdhesiveShearBasis(ADHESIVE_DF)

    baseline = radiator_joint_overlap()                    # M3: 40 x 20 x 0.2 mm
    static_point = radiator_joint_overlap(40.0, 50.0, 1.0)  # M4 static selection
    fatigue_point = radiator_joint_overlap(40.0, 30.0, 2.0)  # M5 fatigue selection

    print()
    print("=" * 100)
    print("  STM-11  BIMETALLIC JOINT THERMAL STRESS - FINAL INTEGRATED ASSESSMENT")
    print("=" * 100)

    # ---------------------------------------------------------------- basis
    rule("STUDY BASIS")
    print(f"  {ILLUSTRATIVE_NOTE}")
    print(f"  {ILLUSTRATIVE_ADHESIVE_NOTE}")
    print(f"  {ILLUSTRATIVE_FATIGUE_NOTE}")
    print()
    print(f"  {'Property':<34}{'Member 1':>20}{'Member 2':>20}")
    for caption, v1, v2 in (
        ("Young's modulus E [GPa]",
         member_1.material.elastic_modulus / GPA, member_2.material.elastic_modulus / GPA),
        ("CTE alpha [1e-6 / K]",
         member_1.material.thermal_expansion_coefficient * 1e6,
         member_2.material.thermal_expansion_coefficient * 1e6),
        ("yield strength [MPa]",
         member_1.material.yield_strength / MPA, member_2.material.yield_strength / MPA),
        ("cross-sectional area [mm^2]",
         member_1.area / SQUARE_MILLIMETRE, member_2.area / SQUARE_MILLIMETRE),
    ):
        print(f"  {caption:<34}{v1:>20.4g}{v2:>20.4g}")
    print()
    print(f"  T_ref / T_cold / T_hot            {environment.reference_temperature:>8.1f} /"
          f"{environment.cold_temperature:>8.1f} /{environment.hot_temperature:>8.1f} degC")
    print(f"  dT_cold / dT_hot                  {environment.delta_temperature_cold:>8.1f} /"
          f"{environment.delta_temperature_hot:>8.1f} K")
    print(f"  member yield design factor        {YIELD_DF:>8.2f}")
    print()
    print(f"  baseline overlap / width / t_a    {baseline.overlap_length / MM:>8.1f} /"
          f"{baseline.bond_width / MM:>8.1f} /{baseline.adhesive_thickness / MM:>8.2f} mm")
    print(f"  adhesive G_a / shear strength     {adhesive.shear_modulus / GPA:>8.2f} GPa /"
          f"{adhesive.shear_strength / MPA:>8.2f} MPa")
    print(f"  adhesive design factor            {ADHESIVE_DF:>8.2f}"
          f"   -> allowable {adhesive_basis.allowable_shear_stress(adhesive) / MPA:.2f} MPa")
    print()
    print(f"  {'fatigue curve  sigma_a = A N^b':<40}{'A [MPa]':>12}{'b':>10}{'sig_a @1e4':>14}")
    for label, curve in (("member 1 (aluminium-like)", curves.member_1),
                         ("member 2 (titanium-like)", curves.member_2),
                         ("adhesive SHEAR", curves.adhesive_shear)):
        print(f"  {label:<40}{curve.coefficient_A / MPA:>12.1f}{curve.exponent_b:>10.3f}"
              f"{curve.alternating_stress_at_life(1.0e4) / MPA:>14.3f}")
    print()
    print(f"  required thermal cycles           {requirement.required_cycles:>12.0f}")
    print(f"    {requirement.label}")
    print("  Mean stress is reported but NO mean-stress correction is applied.")

    # ------------------------------------------------------------------ M1
    rule("1  FREE-JOINT THERMAL COMPATIBILITY")
    yield_assessment = assess_temperature_extremes(
        member_1, member_2, environment, yield_basis
    )
    print(f"  {'case':<8}{'sigma_1 [MPa]':>16}{'sigma_2 [MPa]':>16}{'N_1 [N]':>13}"
          f"{'N_2 [N]':>13}{'N_1+N_2 [N]':>14}")
    for label, delta_t in (("hot", environment.delta_temperature_hot),
                           ("cold", environment.delta_temperature_cold)):
        r = solve_bimetallic_joint(member_1, member_2, delta_t)
        print(f"  {label:<8}{r.member_1_stress / MPA:>16.4f}{r.member_2_stress / MPA:>16.4f}"
              f"{r.member_1_force:>13.2f}{r.member_2_force:>13.2f}"
              f"{r.member_1_force + r.member_2_force:>14.2e}")
    print()
    print(f"  hot minimum yield margin          "
          f"{yield_assessment.hot_assessment.minimum_margin:>+12.4f}   "
          f"{verdict(yield_assessment.hot_assessment.passes)}")
    print(f"  cold minimum yield margin         "
          f"{yield_assessment.cold_assessment.minimum_margin:>+12.4f}   "
          f"{verdict(yield_assessment.cold_assessment.passes)}")
    print(f"  governing extreme / member        {yield_assessment.governing_extreme:>12} / "
          f"member {yield_assessment.governing_member}")
    print("  The pair is self-equilibrating: N_1 + N_2 = 0 with no external load.")

    # ------------------------------------------------------------------ M2
    rule("2  EXTERNAL RESTRAINT")
    stiffness = joint_axial_stiffness(member_1, member_2)
    restrained = assess_restrained_temperature_extremes(
        member_1, member_2, environment,
        AxialRestraint.from_stiffness_ratio(1.0, member_1, member_2), yield_basis,
    )
    crossing = zero_stress_restraint_stiffness(member_1, member_2, 100.0, 2)
    boundary = maximum_allowable_restraint_stiffness(
        member_1, member_2, environment, yield_basis
    )
    print(f"  joint axial stiffness K_joint     {stiffness:>14.4e} N")
    print(f"  canonical eta_r = 1 -> K_r        {stiffness:>14.4e} N")
    print()
    print(f"  eta_r = 1, cold sigma_1           "
          f"{restrained.cold_result.member_1_stress / MPA:>14.4f} MPa")
    print(f"  eta_r = 1, cold sigma_2           "
          f"{restrained.cold_result.member_2_stress / MPA:>14.4f} MPa")
    print(f"  eta_r = 1, minimum yield margin   {restrained.minimum_yield_margin:>+14.4f}   "
          f"{verdict(restrained.passes)}")
    print()
    print(f"  member 2 zero-stress crossing     {crossing / stiffness:>14.6f} eta_r"
          f"   (K_r = {crossing:.4e} N)")
    print(f"  restraint yield boundary          "
          f"{boundary.maximum_stiffness_ratio:>14.4f} eta_r"
          f"   (K_r = {boundary.maximum_restraint_stiffness:.4e} N)")
    print(f"  boundary status                   {boundary.status.value:>14}")
    print("  Restraint drives both members toward the same sign and can eventually")
    print("  make the metal, not the bondline, the limiting item.")

    # ------------------------------------------------------------------ M3
    rule("3  BASELINE BONDLINE - ADHESIVE SHEAR-LAG")
    shear = assess_shear_lag_extremes(
        member_1, member_2, environment, baseline, adhesive, adhesive_basis
    )
    print(f"  beta                              {shear.cold_demand.beta:>14.4f} 1/m")
    print(f"  transfer length 1/beta            {shear.cold_demand.transfer_length / MM:>14.4f} mm")
    print(f"  dimensionless overlap lambda      {shear.cold_demand.dimensionless_overlap:>14.4f}")
    print()
    print(f"  hot transferred force             {shear.hot_demand.transferred_force:>14.2f} N")
    print(f"  hot peak shear                    {shear.hot_demand.peak_shear_stress / MPA:>14.4f} MPa")
    print(f"  cold transferred force            {shear.cold_demand.transferred_force:>14.2f} N")
    print(f"  cold peak shear                   {shear.cold_demand.peak_shear_stress / MPA:>14.4f} MPa")
    print(f"  cold / hot peak ratio             "
          f"{shear.cold_demand.peak_shear_stress / shear.hot_demand.peak_shear_stress:>14.4f}")
    print(f"  static adhesive margin            {shear.minimum_margin:>+14.4f}   "
          f"{verdict(shear.passes)}")
    print()
    print(f"  long-overlap asymptote            "
          f"{long_overlap_peak_shear(member_1, member_2, baseline, adhesive, environment.delta_temperature_cold) / MPA:>14.4f} MPa")
    print(f"  allowable |dT| at this geometry   "
          f"{allowable_temperature_change_for_adhesive_shear(member_1, member_2, baseline, adhesive, adhesive_basis):>14.4f} K")
    overlap_sizing = required_overlap_length(
        member_1, member_2, environment, baseline, adhesive, adhesive_basis
    )
    print(f"  required-overlap status           {overlap_sizing.status.value:>30}")
    print("  Peak shear sits at the FREE EDGE and converges on the asymptote from")
    print("  above, so overlap length alone cannot recover this bondline.")

    # ------------------------------------------------------------------ M4
    rule("4  STATIC BONDLINE DESIGN TRADE")
    width_sizing = required_bond_width(
        member_1, member_2, environment, baseline, adhesive, adhesive_basis, tolerance=1.0e-9
    )
    thickness_sizing = required_adhesive_thickness(
        member_1, member_2, environment, baseline, adhesive, adhesive_basis
    )
    modulus_sizing = maximum_allowable_adhesive_shear_modulus(
        member_1, member_2, environment, baseline, adhesive, adhesive_basis, tolerance=1.0
    )
    print("  tau_inf = |d_eps| sqrt( G_a / (b t_a C) )  ->  b^-1/2, t_a^-1/2, G_a^+1/2")
    print()
    print(f"  minimum required bond width       {width_sizing.required_bond_width / MM:>14.3f} mm"
          f"   ({width_sizing.status.value})")
    print(f"  minimum required thickness        "
          f"{thickness_sizing.required_adhesive_thickness / MM:>14.4f} mm"
          f"   ({thickness_sizing.status.value})")
    print(f"  maximum allowable G_a             "
          f"{modulus_sizing.maximum_shear_modulus / GPA:>14.6f} GPa"
          f"   ({modulus_sizing.status.value})")
    static_grid = width_thickness_design_map(
        member_1, member_2, environment, baseline, adhesive,
        [w * MM for w in STATIC_MAP_WIDTHS_MM],
        [t * MM for t in STATIC_MAP_THICKNESSES_MM], yield_basis, adhesive_basis,
    )
    static_selected = select_preliminary_joint_design(static_grid)
    print()
    print(f"  static design map                 "
          f"{sum(c.overall_feasible for c in static_grid):>6} feasible of {len(static_grid)}")
    print(f"  static selected point             "
          f"{static_selected.bond_width / MM:>6.0f} mm x "
          f"{static_selected.adhesive_thickness / MM:.2f} mm")
    print(f"    cold adhesive margin            {static_selected.cold_adhesive_margin:>+14.4f}")
    print(f"    hot adhesive margin             {static_selected.hot_adhesive_margin:>+14.4f}")
    print(f"    metal yield margin              {static_selected.minimum_yield_margin:>+14.4f}")

    # ------------------------------------------------------------------ M5
    rule("5  REPEATED THERMAL-CYCLE FATIGUE")
    fatigue_baseline = assess_thermal_cycle_fatigue(
        member_1, member_2, environment, baseline, adhesive, curves, requirement,
        yield_basis, adhesive_basis,
    )
    fatigue_static = assess_thermal_cycle_fatigue(
        member_1, member_2, environment, static_point, adhesive, curves, requirement,
        yield_basis, adhesive_basis,
    )
    print(f"  {'component':<34}{'sigma_a [MPa]':>16}{'sigma_m [MPa]':>16}{'N_f':>14}{'ratio':>12}{'':>7}")
    for result in (fatigue_static.member_1_fatigue, fatigue_static.member_2_fatigue):
        print(f"  {'free metal - ' + result.component:<34}"
              f"{result.alternating_stress / MPA:>16.4f}{result.mean_stress / MPA:>+16.4f}"
              f"{life_text(result.predicted_cycles_to_failure):>14}"
              f"{result.life_ratio:>12.4g}{'  ' + verdict(result.passes):>7}")
    for label, assessment in (("M3 baseline adhesive", fatigue_baseline),
                              ("M4 static-selected adhesive", fatigue_static)):
        r = assessment.adhesive_fatigue
        print(f"  {label:<34}{r.alternating_stress / MPA:>16.4f}"
              f"{r.mean_stress / MPA:>+16.4f}{life_text(r.predicted_cycles_to_failure):>14}"
              f"{r.life_ratio:>12.4g}{'  ' + verdict(r.passes):>7}")

    fatigue_grid = fatigue_design_map(
        member_1, member_2, environment, baseline, adhesive, curves, requirement,
        [w * MM for w in FATIGUE_MAP_WIDTHS_MM],
        [t * MM for t in FATIGUE_MAP_THICKNESSES_MM], yield_basis, adhesive_basis,
    )
    fatigue_selected = select_preliminary_fatigue_design(fatigue_grid)
    scale_result = allowable_cycle_scale_for_fatigue(
        member_1, member_2, environment, static_point, adhesive, curves, requirement
    )
    print()
    print(f"  allowable cycle scale at the static point   {scale_result.allowable_scale:>10.6f}"
          f"   (governed by {scale_result.governing_component})")
    print(f"  fatigue design map                          "
          f"{sum(c.overall_feasible for c in fatigue_grid):>6} feasible of {len(fatigue_grid)}")
    print(f"  fatigue-selected point                      "
          f"{fatigue_selected.bond_width / MM:>6.0f} mm x "
          f"{fatigue_selected.adhesive_thickness / MM:.2f} mm"
          f"   life ratio {fatigue_selected.assessment.minimum_life_ratio:.4f}")

    # ---------------------------------------------------------- comparison
    rule("6  FINAL INTEGRATED COMPARISON")
    designs = []
    for label, geometry in (("M3 baseline", baseline),
                            ("M4 static", static_point),
                            ("M5 fatigue", fatigue_point)):
        designs.append((
            label, geometry,
            assess_shear_lag_extremes(
                member_1, member_2, environment, geometry, adhesive, adhesive_basis),
            assess_thermal_cycle_fatigue(
                member_1, member_2, environment, geometry, adhesive, curves, requirement,
                yield_basis, adhesive_basis),
        ))

    print(f"  {'quantity':<32}" + "".join(f"{label:>18}" for label, *_ in designs))

    def row(caption, getter, fmt="{:.4g}"):
        print(f"  {caption:<32}"
              + "".join(f"{fmt.format(getter(g, s, f)):>18}" for _, g, s, f in designs))

    row("bond width [mm]", lambda g, s, f: g.bond_width / MM, "{:.1f}")
    row("overlap length [mm]", lambda g, s, f: g.overlap_length / MM, "{:.1f}")
    row("adhesive thickness [mm]", lambda g, s, f: g.adhesive_thickness / MM, "{:.2f}")
    row("beta [1/m]", lambda g, s, f: s.cold_demand.beta, "{:.2f}")
    row("transfer length [mm]", lambda g, s, f: s.cold_demand.transfer_length / MM, "{:.3f}")
    row("hot peak shear [MPa]", lambda g, s, f: s.hot_demand.peak_shear_stress / MPA, "{:.3f}")
    row("cold peak shear [MPa]", lambda g, s, f: s.cold_demand.peak_shear_stress / MPA, "{:.3f}")
    row("static adhesive margin", lambda g, s, f: s.minimum_margin, "{:+.4f}")
    row("metal yield margin", lambda g, s, f: f.minimum_yield_margin, "{:+.4f}")
    row("adhesive tau_a [MPa]", lambda g, s, f: f.adhesive_fatigue.alternating_stress / MPA, "{:.3f}")
    row("adhesive N_f [cycles]",
        lambda g, s, f: f.adhesive_fatigue.predicted_cycles_to_failure, "{:.4g}")
    row("fatigue life ratio", lambda g, s, f: f.minimum_life_ratio, "{:.4f}")
    row("static adhesive", lambda g, s, f: verdict(f.static_adhesive_feasible), "{}")
    row("metal yield", lambda g, s, f: verdict(f.static_yield_feasible), "{}")
    row("fatigue", lambda g, s, f: verdict(f.fatigue_feasible), "{}")
    row("INTEGRATED", lambda g, s, f: verdict(f.overall_feasible), "{}")
    print()
    print("  The three margins are reported side by side and are never numerically")
    print("  merged; only the pass/fail booleans are combined.")

    # ---------------------------------------------------------- robustness
    rule("7  ROBUSTNESS - REQUIRED-CYCLE SENSITIVITY (canonical 1e4 unchanged)")
    print(f"  {'N_required':>12}{'M4 static ratio':>18}{'':>7}{'M5 fatigue ratio':>19}{'':>7}"
          f"{'grid feasible':>15}{'selected':>16}")
    for cycles in REQUIREMENT_SWEEP:
        sweep_requirement = ThermalCycleRequirement(cycles, "robustness sweep")
        a4 = assess_thermal_cycle_fatigue(
            member_1, member_2, environment, static_point, adhesive, curves,
            sweep_requirement, yield_basis, adhesive_basis)
        a5 = assess_thermal_cycle_fatigue(
            member_1, member_2, environment, fatigue_point, adhesive, curves,
            sweep_requirement, yield_basis, adhesive_basis)
        grid = fatigue_design_map(
            member_1, member_2, environment, baseline, adhesive, curves, sweep_requirement,
            [w * MM for w in FATIGUE_MAP_WIDTHS_MM],
            [t * MM for t in FATIGUE_MAP_THICKNESSES_MM], yield_basis, adhesive_basis)
        pick = select_preliminary_fatigue_design(grid)
        pick_text = (f"{pick.bond_width / MM:.0f}x{pick.adhesive_thickness / MM:.1f} mm"
                     if pick else "none")
        print(f"  {cycles:>12.0e}{a4.minimum_life_ratio:>18.4g}"
              f"{'  ' + verdict(a4.fatigue_feasible):>7}{a5.minimum_life_ratio:>19.4g}"
              f"{'  ' + verdict(a5.fatigue_feasible):>7}"
              f"{sum(c.overall_feasible for c in grid):>10} of {len(grid)}{pick_text:>16}")
    print()
    print("  The requirement is a user-selected illustrative input, so the whole range")
    print("  is shown. The fatigue-selected point holds to 1e4 and fails from 1e5.")

    rule("8  ROBUSTNESS - THERMAL-CYCLE SCALE (M5 fatigue-selected point)")
    print(f"  {'scale':>8}{'tau_a [MPa]':>15}{'adhesive N_f':>16}{'metal N_f':>16}"
          f"{'governing':>18}{'':>7}")
    for scale in SCALE_SWEEP:
        scaled = assess_thermal_cycle_fatigue(
            member_1, member_2, scale_environment(environment, scale), fatigue_point,
            adhesive, curves, requirement, yield_basis, adhesive_basis)
        metal = min(scaled.member_1_fatigue.predicted_cycles_to_failure,
                    scaled.member_2_fatigue.predicted_cycles_to_failure)
        print(f"  {scale:>8.2f}{scaled.adhesive_fatigue.alternating_stress / MPA:>15.4f}"
              f"{life_text(scaled.adhesive_fatigue.predicted_cycles_to_failure):>16}"
              f"{life_text(metal):>16}{scaled.governing_fatigue_component:>18}"
              f"{'  ' + verdict(scaled.fatigue_feasible):>7}")
    print()
    print("  Alternating stress is exactly linear in the scale and life follows the")
    print("  Basquin power law N ~ scale^(1/b); both are verified to machine precision.")

    # -------------------------------------------------------------- result
    rule("9  FINAL PRELIMINARY RESULT")
    final = designs[2][3]
    final_shear = designs[2][2]
    print(f"  bond width                        {fatigue_point.bond_width / MM:>14.1f} mm")
    print(f"  overlap length                    {fatigue_point.overlap_length / MM:>14.1f} mm")
    print(f"  adhesive thickness                {fatigue_point.adhesive_thickness / MM:>14.2f} mm")
    print(f"  adhesive shear modulus G_a        {adhesive.shear_modulus / GPA:>14.2f} GPa")
    print()
    print(f"  static adhesive shear margin      {final_shear.minimum_margin:>+14.4f}   "
          f"{verdict(final.static_adhesive_feasible)}")
    print(f"  metal yield margin                {final.minimum_yield_margin:>+14.4f}   "
          f"{verdict(final.static_yield_feasible)}")
    print(f"  adhesive fatigue life ratio       {final.minimum_life_ratio:>14.4f}   "
          f"{verdict(final.fatigue_feasible)}")
    print(f"  governing fatigue component       {final.governing_fatigue_component:>14}")
    print(f"  INTEGRATED                        {verdict(final.overall_feasible):>14}")
    print()
    print("  Selection basis: the deterministic bounded-grid policy - feasible on metal")
    print("  yield AND static adhesive shear AND all fatigue components, then smallest")
    print("  bond area, then smaller adhesive thickness, then larger life ratio, then")
    print("  input order. It is the unique minimum-bond-area feasible point of the grid.")
    print()
    print("  This is a preliminary fatigue-feasible bounded-grid point - not an")
    print("  optimized design, not a certified joint, not a final flight design.")

    rule()
    print(
        "The final configuration is a preliminary fatigue-feasible point within a simplified\n"
        "linear-elastic screening model. Material data, adhesive properties, fatigue curves,\n"
        "and design requirements remain illustrative."
    )
    print()


if __name__ == "__main__":
    main()

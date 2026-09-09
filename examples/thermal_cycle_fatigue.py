"""Milestone 5 study: does repeated thermal cycling make the joint fatigue-critical?

Milestone 4 found a statically feasible bondline. This study takes the same
verified hot/cold states, forms constant-amplitude thermal-cycle stress ranges,
and screens them against illustrative Basquin S-N curves.

Run with::

    python examples/thermal_cycle_fatigue.py

All material, adhesive and fatigue inputs are ILLUSTRATIVE - NOT DESIGN ALLOWABLE.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from thermal_joint import (  # noqa: E402
    AdhesiveShearBasis,
    ThermalCycleRequirement,
    YieldBasis,
    adhesive_modulus_fatigue_sensitivity,
    adhesive_shear_cycle,
    adhesive_thickness_fatigue_sensitivity,
    allowable_cycle_scale_for_fatigue,
    area_ratio_fatigue_sensitivity,
    assess_fatigue_life,
    assess_thermal_cycle_fatigue,
    bond_width_fatigue_sensitivity,
    fatigue_design_map,
    member_stress_cycles,
    restraint_fatigue_sensitivity,
    scale_environment,
    select_preliminary_fatigue_design,
    temperature_scale_fatigue_sensitivity,
)
from thermal_joint.illustrative import (  # noqa: E402
    ADHESIVE_LIKE,
    ILLUSTRATIVE_ADHESIVE_NOTE,
    ILLUSTRATIVE_FATIGUE_NOTE,
    ILLUSTRATIVE_NOTE,
    MILLIMETRE,
    aluminium_like_member,
    illustrative_cycle_requirement,
    illustrative_fatigue_curves,
    radiator_joint_environment,
    radiator_joint_overlap,
    titanium_like_member,
)

MPA = 1.0e6
GPA = 1.0e9
MM = MILLIMETRE

ETA_RATIOS = (0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 13.7405)
TEMPERATURE_SCALES = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5)
WIDTHS_MM = (20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 220.0)
THICKNESSES_MM = (0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0)
MODULI_GPA = (0.05, 0.1, 0.2, 0.5, 1.0, 2.0)
AREA_RATIOS = (0.25, 0.5, 1.0, 2.0, 4.0)
REQUIREMENT_CANDIDATES = (1.0e3, 1.0e4, 1.0e5, 1.0e6)
MAP_WIDTHS_MM = (20.0, 30.0, 40.0, 50.0, 75.0, 100.0)
MAP_THICKNESSES_MM = (0.5, 0.75, 1.0, 1.5, 2.0)


def rule(title: str = "") -> None:
    if title:
        print()
        print(title)
    print("-" * 104)


def life_text(cycles: float) -> str:
    return "      inf" if math.isinf(cycles) else f"{cycles:9.4g}"


def ratio_text(value: float) -> str:
    return "      inf" if math.isinf(value) else f"{value:9.4g}"


def verdict(ok: bool) -> str:
    return "PASS" if ok else "FAIL"


def print_study_basis(environment, geometry, adhesive, curves, requirement) -> None:
    rule("STUDY BASIS")
    print(f"  {ILLUSTRATIVE_NOTE}")
    print(f"  {ILLUSTRATIVE_ADHESIVE_NOTE}")
    print(f"  {ILLUSTRATIVE_FATIGUE_NOTE}")
    print()
    print(f"  T_ref / T_cold / T_hot            "
          f"{environment.reference_temperature:>8.1f} /{environment.cold_temperature:>8.1f} /"
          f"{environment.hot_temperature:>8.1f} degC")
    print(f"  dT_cold / dT_hot                  "
          f"{environment.delta_temperature_cold:>8.1f} /{environment.delta_temperature_hot:>8.1f} K")
    print(f"  canonical overlap / width / t_a   "
          f"{geometry.overlap_length / MM:>8.1f} /{geometry.bond_width / MM:>8.1f} /"
          f"{geometry.adhesive_thickness / MM:>8.2f} mm")
    print(f"  adhesive G_a / strength           "
          f"{adhesive.shear_modulus / GPA:>8.2f} GPa /{adhesive.shear_strength / MPA:>8.2f} MPa")
    print()
    print(f"  required thermal cycles           {requirement.required_cycles:>12.0f}")
    print(f"    {requirement.label}")
    print()
    print("  One complete cold -> hot -> cold excursion counts as ONE thermal cycle.")
    print("  Only the two endpoint states are used: no transient path, no dwell time,")
    print("  no rate effect, no variable-amplitude spectrum. N always means cycles.")
    print()
    print(f"  {'fatigue curve  sigma_a = A N^b':<42}{'A [MPa]':>12}{'b':>10}"
          f"{'sig_a @1e4':>14}")
    for label, curve in (("member 1 (aluminium-like)", curves.member_1),
                         ("member 2 (titanium-like)", curves.member_2),
                         ("adhesive SHEAR", curves.adhesive_shear)):
        print(f"  {label:<42}{curve.coefficient_A / MPA:>12.1f}{curve.exponent_b:>10.3f}"
              f"{curve.alternating_stress_at_life(1.0e4) / MPA:>14.4f}")
    print()
    print("  MEAN-STRESS POLICY: mean stress is reported but NO mean-stress correction")
    print("  is applied in Milestone 5. Life uses the alternating stress alone.")
    print("  Goodman / Gerber / Soderberg are deliberately excluded: with unsourced")
    print("  illustrative curves an unsourced correction would only add false authority.")


def print_cycle_row(label, cycle, result) -> None:
    print(
        f"  {label:<26}{cycle.hot_value / MPA:>11.4f}{cycle.cold_value / MPA:>11.4f}"
        f"{cycle.alternating_stress / MPA:>11.4f}{cycle.mean_stress / MPA:>+11.4f}"
        f"{life_text(result.predicted_cycles_to_failure):>11}"
        f"{ratio_text(result.life_ratio):>12}{'  ' + verdict(result.passes):>7}"
    )


def cycle_header() -> None:
    print(f"  {'component':<26}{'hot':>11}{'cold':>11}{'ampl.':>11}{'mean':>11}"
          f"{'N_f':>11}{'life ratio':>12}{'':>7}")
    print(f"  {'':<26}{'[MPa]':>11}{'[MPa]':>11}{'[MPa]':>11}{'[MPa]':>11}")


def print_free_metal(member_1, member_2, environment, curves, requirement) -> None:
    rule("M1 FREE-JOINT METAL CYCLES")
    cycle_header()
    cycle_1, cycle_2 = member_stress_cycles(member_1, member_2, environment)
    for label, cycle, curve in (("member 1 aluminium-like", cycle_1, curves.member_1),
                                ("member 2 titanium-like", cycle_2, curves.member_2)):
        print_cycle_row(label, cycle, assess_fatigue_life(label, cycle, curve, requirement))
    print()
    print("  Equal areas make sigma_1 = -sigma_2, so both members see the SAME")
    print("  alternating stress and equal-and-opposite mean stress.")


def print_restrained_metal(member_1, member_2, environment, curves, requirement) -> None:
    rule("M2 RESTRAINED METAL CYCLES  (eta_r = 1)")
    cycle_header()
    point = restraint_fatigue_sensitivity(
        member_1, member_2, environment, curves, requirement, (1.0,)
    )[0]
    print_cycle_row("member 1 aluminium-like", point.member_1_cycle, point.member_1_fatigue)
    print_cycle_row("member 2 titanium-like", point.member_2_cycle, point.member_2_fatigue)
    free_1, free_2 = member_stress_cycles(member_1, member_2, environment)
    print()
    print(f"  member 1 amplitude {free_1.alternating_stress / MPA:.3f} -> "
          f"{point.member_1_cycle.alternating_stress / MPA:.3f} MPa "
          f"({point.member_1_cycle.alternating_stress / free_1.alternating_stress:.2f}x)")
    print(f"  member 2 amplitude {free_2.alternating_stress / MPA:.3f} -> "
          f"{point.member_2_cycle.alternating_stress / MPA:.3f} MPa "
          f"({point.member_2_cycle.alternating_stress / free_2.alternating_stress:.2f}x)")
    print("  Restraint does NOT degrade both members alike: it loads member 1 harder")
    print("  while unloading member 2. Restraint also shifts the mean stress, but life")
    print("  here depends on amplitude alone.")


def print_adhesive_cycle(label, member_1, member_2, environment, geometry, adhesive,
                         curves, requirement, yield_basis, adhesive_basis) -> None:
    rule(label)
    cycle = adhesive_shear_cycle(member_1, member_2, environment, geometry, adhesive)
    assessment = assess_thermal_cycle_fatigue(
        member_1, member_2, environment, geometry, adhesive, curves, requirement,
        yield_basis, adhesive_basis,
    )
    result = assessment.adhesive_fatigue
    print(f"  bond width / thickness / overlap  "
          f"{geometry.bond_width / MM:>8.1f} /{geometry.adhesive_thickness / MM:>8.2f} /"
          f"{geometry.overlap_length / MM:>8.1f} mm")
    print()
    print(f"  signed peak shear, hot            {cycle.hot_value / MPA:>14.4f} MPa")
    print(f"  signed peak shear, cold           {cycle.cold_value / MPA:>14.4f} MPa")
    print(f"  opposite signs (cycle crosses 0)  {str(cycle.crosses_zero):>14}")
    print(f"  alternating shear tau_a           {cycle.alternating_stress / MPA:>14.4f} MPa")
    print(f"  mean shear tau_m                  {cycle.mean_stress / MPA:>14.4f} MPa")
    print(f"  stress ratio R                    {cycle.stress_ratio:>14.4f}")
    print()
    print(f"  static adhesive shear margin      "
          f"{assessment.minimum_adhesive_margin:>+14.4f}   "
          f"{verdict(assessment.static_adhesive_feasible)}")
    print(f"  predicted life N_f                "
          f"{life_text(result.predicted_cycles_to_failure):>14} cycles")
    print(f"  preliminary fatigue life margin   {result.margin:>+14.4f}   "
          f"{verdict(result.passes)}")
    return assessment


def print_governing(assessment) -> None:
    rule("GOVERNING FATIGUE COMPONENT")
    print(f"  {'component':<26}{'amplitude [MPa]':>18}{'N_f':>14}{'life ratio':>14}{'':>7}")
    for result in assessment.fatigue_results:
        print(f"  {result.component:<26}{result.alternating_stress / MPA:>18.4f}"
              f"{life_text(result.predicted_cycles_to_failure):>14}"
              f"{ratio_text(result.life_ratio):>14}"
              f"{'  ' + verdict(result.passes):>7}")
    print()
    print(f"  governing fatigue component       "
          f"{assessment.governing_fatigue_component:>18}")
    print(f"  minimum life ratio                {assessment.minimum_life_ratio:>18.4f}")
    print()
    print(f"  static metal yield margin         {assessment.minimum_yield_margin:>+18.4f}"
          f"   {verdict(assessment.static_yield_feasible)}")
    print(f"  static adhesive shear margin      {assessment.minimum_adhesive_margin:>+18.4f}"
          f"   {verdict(assessment.static_adhesive_feasible)}")
    print(f"  fatigue feasible                  {str(assessment.fatigue_feasible):>18}")
    print(f"  OVERALL preliminary feasible      {str(assessment.overall_feasible):>18}")
    print()
    print("  Static and fatigue margins are reported side by side and are never")
    print("  numerically blended; only the pass/fail booleans are combined.")


def print_requirement_sensitivity(member_1, member_2, environment, geometry, adhesive,
                                  curves, yield_basis, adhesive_basis) -> None:
    rule("REQUIRED-CYCLE SENSITIVITY  (so the 1e4 choice is not doing hidden work)")
    print(f"  {'N_required':>14}{'adhesive ratio':>18}{'member 1 ratio':>18}"
          f"{'governing':>18}{'':>7}")
    for cycles in REQUIREMENT_CANDIDATES:
        assessment = assess_thermal_cycle_fatigue(
            member_1, member_2, environment, geometry, adhesive, curves,
            ThermalCycleRequirement(cycles, "requirement sweep"), yield_basis, adhesive_basis,
        )
        print(f"  {cycles:>14.0e}{assessment.adhesive_fatigue.life_ratio:>18.4f}"
              f"{assessment.member_1_fatigue.life_ratio:>18.4g}"
              f"{assessment.governing_fatigue_component:>18}"
              f"{'  ' + verdict(assessment.fatigue_feasible):>7}")
    print()
    print("  The M4 selected point passes at 1e3 cycles and fails from 1e4 upward.")
    print("  The requirement is an illustrative study input, not a qualification")
    print("  requirement, so the whole range is shown rather than one chosen number.")


def print_restraint_sensitivity(member_1, member_2, environment, curves, requirement) -> None:
    rule("RESTRAINT SENSITIVITY - METAL FATIGUE ONLY")
    print(f"  {'eta_r':>9}{'m1 ampl':>11}{'m1 mean':>11}{'m1 N_f':>12}"
          f"{'m2 ampl':>11}{'m2 mean':>11}{'m2 N_f':>12}{'':>7}")
    print(f"  {'':>9}{'[MPa]':>11}{'[MPa]':>11}{'':>12}{'[MPa]':>11}{'[MPa]':>11}")
    for point in restraint_fatigue_sensitivity(
        member_1, member_2, environment, curves, requirement, ETA_RATIOS
    ):
        print(
            f"  {point.stiffness_ratio:>9.4f}"
            f"{point.member_1_cycle.alternating_stress / MPA:>11.3f}"
            f"{point.member_1_cycle.mean_stress / MPA:>+11.3f}"
            f"{life_text(point.member_1_fatigue.predicted_cycles_to_failure):>12}"
            f"{point.member_2_cycle.alternating_stress / MPA:>11.3f}"
            f"{point.member_2_cycle.mean_stress / MPA:>+11.3f}"
            f"{life_text(point.member_2_fatigue.predicted_cycles_to_failure):>12}"
            f"{'  ' + verdict(point.passes):>7}"
        )
    print()
    print("  Member 1's amplitude rises monotonically. Member 2's does NOT: it falls to")
    print("  a minimum near eta_r ~ 0.66 and rises again. That is the Milestone 2")
    print("  zero-stress crossing (eta_r = 0.663399), which is dT-independent, so BOTH")
    print("  endpoints vanish there and member 2 sees no thermal cycle at all.")
    print()
    print("  The adhesive is deliberately excluded from this sweep: the Milestone 3")
    print("  overlap model is posed for the free joint, and reusing it under external")
    print("  restraint would extend it past what has been verified.")


def print_sweep(title, caption, points, value_format, note=None, value_scale=1.0) -> None:
    rule(title)
    print(f"  {caption:>14}{'tau_a [MPa]':>14}{'adhesive N_f':>15}{'static MS':>12}"
          f"{'fatigue MS':>13}{'':>7}")
    for point in points:
        assessment = point.assessment
        result = assessment.adhesive_fatigue
        print(
            f"  {value_format.format(point.value / value_scale):>14}"
            f"{result.alternating_stress / MPA:>14.4f}"
            f"{life_text(result.predicted_cycles_to_failure):>15}"
            f"{assessment.minimum_adhesive_margin:>+12.4f}"
            f"{result.margin:>+13.4g}"
            f"{'  ' + verdict(assessment.overall_feasible):>7}"
        )
    if note:
        print()
        print(note)


def print_allowable_scale(member_1, member_2, environment, geometry, adhesive,
                          curves, requirement, yield_basis, adhesive_basis) -> None:
    rule("ALLOWABLE THERMAL-CYCLE AMPLITUDE")
    result = allowable_cycle_scale_for_fatigue(
        member_1, member_2, environment, geometry, adhesive, curves, requirement
    )
    print("  Every alternating stress is exactly proportional to the excursion scale,")
    print("  so this inversion is exact rather than iterative:")
    print("      sigma_a,allow = A N_required^b ,  scale_i = sigma_a,allow,i / sigma_a,i")
    print()
    print(f"  {'component':<26}{'allowable scale':>18}")
    for label, scale in (("member 1", result.member_1_scale),
                         ("member 2", result.member_2_scale),
                         ("adhesive shear", result.adhesive_scale)):
        print(f"  {label:<26}{scale:>18.4f}")
    print()
    print(f"  governing component               {result.governing_component:>18}")
    print(f"  allowable cycle scale             {result.allowable_scale:>18.6f}")
    print(f"  canonical excursion passes        {str(result.canonical_cycle_passes):>18}")
    scaled = scale_environment(environment, result.allowable_scale)
    check = assess_thermal_cycle_fatigue(
        member_1, member_2, scaled, geometry, adhesive, curves, requirement,
        yield_basis, adhesive_basis,
    )
    print()
    print(f"  at that scale: dT_hot {scaled.delta_temperature_hot:+.2f} K, "
          f"dT_cold {scaled.delta_temperature_cold:+.2f} K, "
          f"minimum life ratio {check.minimum_life_ratio:.9f}")
    print(f"  So this bondline needs the excursion cut by about "
          f"{(1 - result.allowable_scale) * 100:.1f}% to reach {requirement.required_cycles:.0f} cycles.")


def print_design_map(candidates, widths_mm, thicknesses_mm) -> None:
    rule("BOUNDED PRELIMINARY FATIGUE DESIGN MAP  (G_a = 1.0 GPa)")
    print("  adhesive life ratio N_f / N_required; '*' passes static yield, static")
    print("  shear AND fatigue")
    print()
    header = "".join(f"{t:>12.2f}" for t in thicknesses_mm)
    print(f"  {'b [mm]':>8}  t_a [mm]:{header}")
    index = 0
    for width_mm in widths_mm:
        cells = []
        for _ in thicknesses_mm:
            candidate = candidates[index]
            marker = "*" if candidate.overall_feasible else " "
            cells.append(f"{candidate.assessment.adhesive_fatigue.life_ratio:>11.3f}{marker}")
            index += 1
        print(f"  {width_mm:>8.0f}           " + "".join(cells))
    feasible = [c for c in candidates if c.overall_feasible]
    print()
    print(f"  feasible on all three screens: {len(feasible)} of {len(candidates)}")
    print("  This is a bounded preliminary design map, not an optimization.")


def print_conflict_and_selection(candidates, selected_assessment, static_geometry) -> None:
    rule("STATIC VERSUS FATIGUE, AND THE PRELIMINARY SELECTED POINT")
    print("  The Milestone 4 selected point (b = 50 mm, t_a = 1.00 mm) passes static")
    print(f"  adhesive shear at MS {selected_assessment.minimum_adhesive_margin:+.4f} but reaches only "
          f"{selected_assessment.adhesive_fatigue.predicted_cycles_to_failure:.0f} cycles,")
    print(f"  a life ratio of {selected_assessment.adhesive_fatigue.life_ratio:.4f}. Static feasibility did not imply")
    print("  fatigue feasibility.")
    print()
    selected = select_preliminary_fatigue_design(candidates)
    if selected is None:
        print("  No point in the bounded grid passes all three screens.")
        return
    print("  Re-selecting over a grid that also allows thicker bondlines, using the")
    print("  same deterministic Milestone 4 policy (feasible first, then smallest bond")
    print("  area, then smaller thickness, then larger life ratio, then input order):")
    print()
    assessment = selected.assessment
    print(f"  {'':<34}{'M4 static pick':>18}{'M5 fatigue pick':>18}")
    print(f"  {'bond width [mm]':<34}{static_geometry.bond_width / MM:>18.2f}"
          f"{selected.bond_width / MM:>18.2f}")
    print(f"  {'adhesive thickness [mm]':<34}{static_geometry.adhesive_thickness / MM:>18.2f}"
          f"{selected.adhesive_thickness / MM:>18.2f}")
    print(f"  {'bond area [mm^2]':<34}{static_geometry.bond_area * 1e6:>18.0f}"
          f"{selected.bond_area * 1e6:>18.0f}")
    print(f"  {'static adhesive margin':<34}"
          f"{selected_assessment.minimum_adhesive_margin:>+18.4f}"
          f"{assessment.minimum_adhesive_margin:>+18.4f}")
    print(f"  {'adhesive life N_f':<34}"
          f"{selected_assessment.adhesive_fatigue.predicted_cycles_to_failure:>18.0f}"
          f"{assessment.adhesive_fatigue.predicted_cycles_to_failure:>18.0f}")
    print(f"  {'minimum life ratio':<34}{selected_assessment.minimum_life_ratio:>18.4f}"
          f"{assessment.minimum_life_ratio:>18.4f}")
    print(f"  {'static metal yield margin':<34}"
          f"{selected_assessment.minimum_yield_margin:>+18.4f}"
          f"{assessment.minimum_yield_margin:>+18.4f}")
    print(f"  {'overall feasible':<34}{str(selected_assessment.overall_feasible):>18}"
          f"{str(assessment.overall_feasible):>18}")
    print()
    print("  The fatigue-feasible pick is actually the NARROWER bond, because the wider")
    print("  grid allows a thicker, more compliant bondline. The selected point is a")
    print("  preliminary feasible point within the bounded study grid, not an optimized")
    print("  flight-joint design, and its life ratio of "
          f"{assessment.minimum_life_ratio:.3f} is slim.")


def main() -> None:
    member_1 = aluminium_like_member(100.0)
    member_2 = titanium_like_member(100.0)
    environment = radiator_joint_environment()
    adhesive = ADHESIVE_LIKE
    curves = illustrative_fatigue_curves()
    requirement = illustrative_cycle_requirement()
    yield_basis = YieldBasis(1.25)
    adhesive_basis = AdhesiveShearBasis(1.25)

    baseline_geometry = radiator_joint_overlap()
    selected_geometry = radiator_joint_overlap(40.0, 50.0, 1.0)

    print()
    print("=" * 104)
    print("  THERMAL-CYCLE FATIGUE SCREEN FOR THE BIMETALLIC JOINT - MILESTONE 5 STUDY")
    print("=" * 104)

    print_study_basis(environment, baseline_geometry, adhesive, curves, requirement)
    print_free_metal(member_1, member_2, environment, curves, requirement)
    print_restrained_metal(member_1, member_2, environment, curves, requirement)

    print_adhesive_cycle(
        "M3 BASELINE ADHESIVE CYCLE  (b = 20 mm, t_a = 0.2 mm - already fails statically)",
        member_1, member_2, environment, baseline_geometry, adhesive, curves,
        requirement, yield_basis, adhesive_basis,
    )
    selected_assessment = print_adhesive_cycle(
        "M4 SELECTED ADHESIVE CYCLE  (b = 50 mm, t_a = 1.0 mm - passes statically)",
        member_1, member_2, environment, selected_geometry, adhesive, curves,
        requirement, yield_basis, adhesive_basis,
    )

    print_governing(selected_assessment)
    print_requirement_sensitivity(
        member_1, member_2, environment, selected_geometry, adhesive, curves,
        yield_basis, adhesive_basis,
    )
    print_restraint_sensitivity(member_1, member_2, environment, curves, requirement)

    print_sweep(
        "TEMPERATURE-RANGE SENSITIVITY  (both excursions scaled, asymmetry preserved)",
        "scale",
        temperature_scale_fatigue_sensitivity(
            member_1, member_2, environment, selected_geometry, adhesive, curves,
            requirement, TEMPERATURE_SCALES, yield_basis, adhesive_basis,
        ),
        "{:.2f}",
        "  Amplitude is linear in the scale, so life follows the Basquin power law\n"
        "  N ~ scale^(1/b): a 25% larger excursion costs a factor of ~4.5 in life.",
    )

    print_sweep(
        "BOND-WIDTH FATIGUE SENSITIVITY  (t_a = 1.0 mm, G_a = 1.0 GPa)",
        "b [mm]",
        bond_width_fatigue_sensitivity(
            member_1, member_2, environment, radiator_joint_overlap(40.0, 20.0, 1.0),
            adhesive, curves, requirement, [w * MM for w in WIDTHS_MM],
            yield_basis, adhesive_basis,
        ),
        "{:.0f}",
        "  The static-critical width and the fatigue-critical width are NOT the same:\n"
        "  fatigue needs the wider bond.",
        value_scale=MM,
    )

    print_sweep(
        "BONDLINE-THICKNESS FATIGUE SENSITIVITY  (b = 50 mm, G_a = 1.0 GPa)",
        "t_a [mm]",
        adhesive_thickness_fatigue_sensitivity(
            member_1, member_2, environment, radiator_joint_overlap(40.0, 50.0, 1.0),
            adhesive, curves, requirement, [t * MM for t in THICKNESSES_MM],
            yield_basis, adhesive_basis,
        ),
        "{:.2f}",
        value_scale=MM,
    )

    print_sweep(
        "ADHESIVE-MODULUS FATIGUE SENSITIVITY  (b = 50 mm, t_a = 1.0 mm)",
        "G_a [GPa]",
        adhesive_modulus_fatigue_sensitivity(
            member_1, member_2, environment, selected_geometry, adhesive, curves,
            requirement, [g * GPA for g in MODULI_GPA], yield_basis, adhesive_basis,
        ),
        "{:.2f}",
        "  A softer adhesive lengthens predicted life here, but creep, peel and\n"
        "  durability are all absent from this model, so low modulus is NOT\n"
        "  universally superior.",
        value_scale=GPA,
    )

    print_sweep(
        "AREA-RATIO FATIGUE SENSITIVITY  (M1 recomputed end to end at each ratio)",
        "A_1/A_2",
        area_ratio_fatigue_sensitivity(
            member_1, member_2, environment, selected_geometry, adhesive, curves,
            requirement, AREA_RATIOS, yield_basis, adhesive_basis,
        ),
        "{:.2f}",
        "  The static-best area ratio is not automatically the fatigue-best one; here\n"
        "  both favour the smaller aluminium-like area, and the adhesive governs at\n"
        "  every ratio.",
    )

    print_allowable_scale(
        member_1, member_2, environment, selected_geometry, adhesive, curves,
        requirement, yield_basis, adhesive_basis,
    )

    candidates = fatigue_design_map(
        member_1, member_2, environment, baseline_geometry, adhesive, curves, requirement,
        [w * MM for w in MAP_WIDTHS_MM], [t * MM for t in MAP_THICKNESSES_MM],
        yield_basis, adhesive_basis,
    )
    print_design_map(candidates, MAP_WIDTHS_MM, MAP_THICKNESSES_MM)
    print_conflict_and_selection(candidates, selected_assessment, selected_geometry)

    rule()
    print(
        "Milestone 5 is a first-order constant-amplitude thermal-cycle fatigue screen built\n"
        "on the verified hot/cold endpoint states. It uses illustrative Basquin curves, no\n"
        "mean-stress correction, and no crack growth, fracture mechanics, rainflow counting,\n"
        "Miner summation, adhesive peel fatigue, multiaxial fatigue, plastic-strain fatigue,\n"
        "creep-fatigue interaction, temperature-dependent fatigue properties, environmental\n"
        "degradation or statistical scatter. It is a screening layer, not certification life."
    )
    print()


if __name__ == "__main__":
    main()

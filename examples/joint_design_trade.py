"""Milestone 4 study: which first-order levers recover adhesive-shear feasibility?

Milestone 3 ended on a negative result - the canonical 40 mm overlap fails
adhesive shear, and more overlap cannot fix it. This study works through the
remaining levers (bond width, adhesive thickness, adhesive modulus, area ratio,
thermal excursion, CTE mismatch), sizes each one inversely, and evaluates a
bounded width x thickness design map.

Run with::

    python examples/joint_design_trade.py

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
    BondedOverlapGeometry,
    YieldBasis,
    adhesive_modulus_sensitivity,
    adhesive_thickness_sensitivity,
    allowable_temperature_change_for_adhesive_shear,
    area_ratio_sensitivity,
    bond_width_sensitivity,
    cte_mismatch_sensitivity,
    evaluate_joint_design,
    long_overlap_peak_shear,
    maximum_allowable_adhesive_shear_modulus,
    required_adhesive_thickness,
    required_bond_width,
    select_preliminary_joint_design,
    thermal_excursion_sensitivity,
    uniform_shear_floor,
    width_thickness_design_map,
)
from thermal_joint.illustrative import (  # noqa: E402
    ADHESIVE_LIKE,
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
MM = MILLIMETRE
YIELD_DESIGN_FACTOR = 1.25
ADHESIVE_DESIGN_FACTOR = 1.25

WIDTHS_MM = (10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0)
THICKNESSES_MM = (0.05, 0.10, 0.20, 0.30, 0.50, 0.75, 1.00, 1.50)
MODULI_GPA = (0.05, 0.10, 0.20, 0.50, 1.0, 2.0, 5.0)
AREA_RATIOS = (0.25, 0.5, 1.0, 2.0, 4.0)
EXCURSIONS = (20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0)
MISMATCH_SCALES = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5)
DESIGN_FACTORS = (1.0, 1.25, 1.5, 2.0)
MAP_WIDTHS_MM = (10.0, 20.0, 30.0, 40.0, 50.0, 75.0, 100.0)
MAP_THICKNESSES_MM = (0.1, 0.2, 0.3, 0.5, 0.75, 1.0)


def rule(title: str = "") -> None:
    if title:
        print()
        print(title)
    print("-" * 100)


def margin_text(margin: float) -> str:
    return "     inf" if math.isinf(margin) else f"{margin:+.4f}"


def verdict(feasible: bool) -> str:
    return "PASS" if feasible else "FAIL"


def print_study_basis(member_1, member_2, environment, geometry, adhesive, basis) -> None:
    rule("STUDY BASIS")
    print(f"  {ILLUSTRATIVE_NOTE}")
    print(f"  {ILLUSTRATIVE_ADHESIVE_NOTE}")
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
    print(f"  T_ref / T_cold / T_hot             "
          f"{environment.reference_temperature:>8.1f} /{environment.cold_temperature:>8.1f} /"
          f"{environment.hot_temperature:>8.1f} degC")
    print(f"  dT_cold / dT_hot                   "
          f"{environment.delta_temperature_cold:>8.1f} /{environment.delta_temperature_hot:>8.1f} K")
    print(f"  overlap length L_b (fixed)         {geometry.overlap_length / MM:>14.2f} mm")
    print(f"  canonical bond width b             {geometry.bond_width / MM:>14.2f} mm")
    print(f"  canonical adhesive thickness t_a   {geometry.adhesive_thickness / MM:>14.3f} mm")
    print(f"  canonical adhesive modulus G_a     {adhesive.shear_modulus / GPA:>14.3f} GPa")
    print(f"  adhesive shear strength            {adhesive.shear_strength / MPA:>14.2f} MPa")
    print(f"  adhesive design factor             {basis.design_factor:>14.2f}")
    print(f"  adhesive allowable tau_allow       "
          f"{basis.allowable_shear_stress(adhesive) / MPA:>14.2f} MPa")
    print(f"  member yield design factor         {YIELD_DESIGN_FACTOR:>14.2f}")


def print_baseline(candidate) -> None:
    rule("BASELINE  (Milestone 3 canonical point)")
    print(f"  bond width / thickness / modulus   "
          f"{candidate.bond_width / MM:>6.1f} mm /{candidate.adhesive_thickness / MM:>6.2f} mm /"
          f"{candidate.adhesive_shear_modulus / GPA:>6.2f} GPa")
    print(f"  beta / transfer length / lambda    "
          f"{candidate.beta:>8.3f} /{candidate.transfer_length / MM:>7.3f} mm /"
          f"{candidate.dimensionless_overlap:>7.4f}")
    print()
    print(f"  hot  peak shear                    {candidate.hot_peak_shear_stress / MPA:>14.4f} MPa"
          f"   margin {margin_text(candidate.hot_adhesive_margin)}")
    print(f"  cold peak shear                    {candidate.cold_peak_shear_stress / MPA:>14.4f} MPa"
          f"   margin {margin_text(candidate.cold_adhesive_margin)}")
    print(f"  governing extreme                  {candidate.governing_extreme:>14}")
    print()
    print(f"  minimum member yield margin        {margin_text(candidate.minimum_yield_margin):>14}"
          f"   {verdict(candidate.yield_feasible)}")
    print(f"  minimum adhesive shear margin      {margin_text(candidate.adhesive_margin):>14}"
          f"   {verdict(candidate.adhesive_feasible)}")
    print(f"  overall preliminary feasibility    {str(candidate.overall_feasible):>14}")
    print()
    print("  The metal passes comfortably; the bondline is what fails. The two")
    print("  margins are kept separate throughout - only the booleans are ANDed.")


def print_scaling(member_1, member_2, geometry, adhesive, candidate) -> None:
    rule("ANALYTICAL SCALING  (long-overlap asymptote)")
    print("  Substituting N_t = d_eps / C into tau_inf = |N_t| beta / b gives")
    print("      tau_inf = |d_eps| * sqrt( G_a / (b * t_a * C) )")
    print()
    print(f"  {'lever':<26}{'exponent on tau_inf':>22}{'to halve tau_inf':>22}")
    for lever, exponent, action in (
        ("bond width b", "-1/2", "4x wider"),
        ("adhesive thickness t_a", "-1/2", "4x thicker"),
        ("adhesive modulus G_a", "+1/2", "4x softer"),
        ("|dT| and |alpha_1-alpha_2|", "+1", "2x smaller"),
    ):
        print(f"  {lever:<26}{exponent:>22}{action:>22}")
    print()
    print(f"  long-overlap asymptote (cold)      {candidate.long_overlap_asymptote / MPA:>14.4f} MPa")
    print(f"  uniform-shear floor |N_t|/(b L_b)  {candidate.uniform_shear_floor / MPA:>14.4f} MPa")
    print(f"  allowable                          {candidate.allowable_shear_stress / MPA:>14.4f} MPa")
    print()
    print("  Why overlap failed: past a few transfer lengths tau_peak sits on the")
    print("  asymptote, which is independent of L_b. Why thickness and modulus are")
    print("  bounded: as the bondline softens, tau_peak falls only to the")
    print("  uniform-shear floor. Width is the one lever with no floor, because the")
    print("  floor itself falls as 1/b. And because every geometric lever enters")
    print("  under a square root, halving the peak shear costs a factor of four.")


def print_candidate_sweep(title, caption, candidates, value_of, value_format) -> None:
    rule(title)
    print(f"  {caption:>12}{'beta [1/m]':>13}{'1/beta [mm]':>13}{'tau_peak':>12}"
          f"{'adh margin':>13}{'yield margin':>14}{'':>7}")
    print(f"  {'':>12}{'':>13}{'':>13}{'[MPa]':>12}")
    for candidate in candidates:
        print(
            f"  {value_format.format(value_of(candidate)):>12}{candidate.beta:>13.3f}"
            f"{candidate.transfer_length / MM:>13.4f}"
            f"{candidate.peak_shear_stress / MPA:>12.4f}"
            f"{margin_text(candidate.adhesive_margin):>13}"
            f"{margin_text(candidate.minimum_yield_margin):>14}"
            f"{'  ' + verdict(candidate.overall_feasible):>7}"
        )


def print_point_sweep(title, caption, points, value_format) -> None:
    rule(title)
    print(f"  {caption:>14}{'M1 force [N]':>15}{'beta [1/m]':>13}{'tau_peak':>12}"
          f"{'adh margin':>13}{'yield margin':>14}{'':>7}")
    print(f"  {'':>14}{'':>15}{'':>13}{'[MPa]':>12}")
    for point in points:
        print(
            f"  {value_format.format(point.value):>14}{point.transferred_force:>15.3f}"
            f"{point.beta:>13.3f}{point.peak_shear_stress / MPA:>12.4f}"
            f"{margin_text(point.adhesive_margin):>13}"
            f"{margin_text(point.minimum_yield_margin):>14}"
            f"{'  ' + verdict(point.overall_feasible):>7}"
        )


def print_design_factor_sensitivity(member_1, member_2, environment, geometry, adhesive) -> None:
    rule("DESIGN-FACTOR SENSITIVITY  (the design basis materially changes the sizing)")
    print(f"  {'factor':>8}{'tau_allow [MPa]':>18}{'cold margin':>14}{'required width [mm]':>22}")
    for factor in DESIGN_FACTORS:
        basis = AdhesiveShearBasis(factor)
        candidate = evaluate_joint_design(
            member_1, member_2, environment, geometry, adhesive,
            YieldBasis(YIELD_DESIGN_FACTOR), basis,
        )
        sizing = required_bond_width(
            member_1, member_2, environment, geometry, adhesive, basis,
            maximum_bond_width=5.0, tolerance=1.0e-9,
        )
        width = (
            f"{sizing.required_bond_width / MM:.3f}"
            if sizing.has_finite_width else sizing.status.value
        )
        print(
            f"  {factor:>8.2f}{basis.allowable_shear_stress(adhesive) / MPA:>18.2f}"
            f"{margin_text(candidate.cold_adhesive_margin):>14}{width:>22}"
        )
    print()
    print("  Required width scales as 1/tau_allow^2, so doubling the design factor")
    print("  quadruples the bond width needed. The factor lives in the basis, never")
    print("  inside the adhesive record.")


def print_inverse(member_1, member_2, environment, geometry, adhesive, basis) -> None:
    rule("INVERSE DESIGN  (single-lever recovery from the canonical baseline)")

    width = required_bond_width(
        member_1, member_2, environment, geometry, adhesive, basis, tolerance=1.0e-9
    )
    print(f"  minimum required bond width        {width.status.value}")
    if width.has_finite_width:
        print(f"    {width.required_bond_width / MM:>10.3f} mm  "
              f"({width.required_bond_width / geometry.bond_width:.2f}x the canonical 20 mm)"
              f"   margin {width.margin_at_required_width:+.2e}")

    thickness = required_adhesive_thickness(
        member_1, member_2, environment, geometry, adhesive, basis
    )
    print(f"  minimum required adhesive thickness {thickness.status.value}")
    if thickness.has_finite_thickness:
        print(f"    {thickness.required_adhesive_thickness / MM:>10.4f} mm "
              f"({thickness.required_adhesive_thickness / geometry.adhesive_thickness:.2f}x the "
              f"canonical 0.2 mm)   margin {thickness.margin_at_required_thickness:+.2e}")
    print(f"    uniform-shear floor {thickness.uniform_shear_floor / MPA:.4f} MPa vs allowable "
          f"{thickness.allowable_shear_stress / MPA:.2f} MPa")

    modulus = maximum_allowable_adhesive_shear_modulus(
        member_1, member_2, environment, geometry, adhesive, basis, tolerance=1.0
    )
    print(f"  maximum allowable adhesive modulus {modulus.status.value}")
    if modulus.has_finite_modulus:
        print(f"    {modulus.maximum_shear_modulus / GPA:>10.6f} GPa "
              f"({adhesive.shear_modulus / modulus.maximum_shear_modulus:.2f}x softer than the "
              f"canonical 1.0 GPa)   margin {modulus.margin_at_maximum_modulus:+.2e}")

    print()
    print("  All three demand roughly the same factor on G_a/(b t_a), because the")
    print("  asymptote depends only on that group. Each is a model-based screening")
    print("  value, not an allowable, and each is individually extreme: a 220 mm")
    print("  bond, a 2.5 mm bondline or an 80 MPa-class adhesive modulus.")


def print_design_map(candidates, widths_mm, thicknesses_mm, label) -> None:
    rule(f"BOUNDED PRELIMINARY DESIGN MAP - {label}")
    print("  cold peak shear [MPa], '*' marks a point passing BOTH screens")
    print()
    header = "".join(f"{t:>11.2f}" for t in thicknesses_mm)
    print(f"  {'b [mm]':>8}   t_a [mm]:{header}")
    index = 0
    for width_mm in widths_mm:
        cells = []
        for _ in thicknesses_mm:
            candidate = candidates[index]
            marker = "*" if candidate.overall_feasible else " "
            cells.append(f"{candidate.cold_peak_shear_stress / MPA:>10.3f}{marker}")
            index += 1
        print(f"  {width_mm:>8.0f}            " + "".join(cells))
    feasible = [c for c in candidates if c.overall_feasible]
    print()
    print(f"  feasible points: {len(feasible)} of {len(candidates)}")
    print("  This is a bounded preliminary design map, not an optimization.")


def print_selection(candidates, baseline) -> None:
    rule("PRELIMINARY SELECTED POINT")
    print("  Policy 'minimum_bond_area', applied in this order:")
    print("    1. must pass member yield AND adhesive shear at both extremes")
    print("    2. smallest bond area (= width x overlap; the overlap is fixed here)")
    print("    3. tie-break: smaller adhesive thickness")
    print("    4. tie-break: larger adhesive margin")
    print("    5. tie-break: earlier position in the grid")
    print()
    selected = select_preliminary_joint_design(candidates)
    if selected is None:
        print("  No point in the bounded grid passes both screens.")
        return

    print(f"  {'':<34}{'baseline':>18}{'selected':>18}")
    for caption, base_value, sel_value, fmt in (
        ("bond width [mm]", baseline.bond_width / MM, selected.bond_width / MM, "{:.2f}"),
        ("adhesive thickness [mm]", baseline.adhesive_thickness / MM,
         selected.adhesive_thickness / MM, "{:.2f}"),
        ("adhesive modulus [GPa]", baseline.adhesive_shear_modulus / GPA,
         selected.adhesive_shear_modulus / GPA, "{:.2f}"),
        ("bond area [mm^2]", baseline.bond_area / SQUARE_MILLIMETRE,
         selected.bond_area / SQUARE_MILLIMETRE, "{:.0f}"),
        ("adhesive volume [mm^3]", baseline.adhesive_volume * 1.0e9,
         selected.adhesive_volume * 1.0e9, "{:.1f}"),
        ("transfer length [mm]", baseline.transfer_length / MM,
         selected.transfer_length / MM, "{:.3f}"),
        ("hot peak shear [MPa]", baseline.hot_peak_shear_stress / MPA,
         selected.hot_peak_shear_stress / MPA, "{:.4f}"),
        ("cold peak shear [MPa]", baseline.cold_peak_shear_stress / MPA,
         selected.cold_peak_shear_stress / MPA, "{:.4f}"),
    ):
        print(f"  {caption:<34}{fmt.format(base_value):>18}{fmt.format(sel_value):>18}")
    print(f"  {'hot adhesive margin':<34}{margin_text(baseline.hot_adhesive_margin):>18}"
          f"{margin_text(selected.hot_adhesive_margin):>18}")
    print(f"  {'cold adhesive margin':<34}{margin_text(baseline.cold_adhesive_margin):>18}"
          f"{margin_text(selected.cold_adhesive_margin):>18}")
    print(f"  {'governing extreme':<34}{baseline.governing_extreme:>18}"
          f"{selected.governing_extreme:>18}")
    print(f"  {'minimum member yield margin':<34}{margin_text(baseline.minimum_yield_margin):>18}"
          f"{margin_text(selected.minimum_yield_margin):>18}")
    print(f"  {'overall feasible':<34}{str(baseline.overall_feasible):>18}"
          f"{str(selected.overall_feasible):>18}")
    print()
    print("  Bond area is a geometric diagnostic only. No adhesive density is")
    print("  supplied, so no adhesive mass is computed or implied.")
    print()
    print("  The selected point is a preliminary feasible point within the bounded")
    print("  study grid, not an optimized flight-joint design, and not a global")
    print("  optimum. A thinner adhesive is not automatically better for")
    print("  manufacturing or durability, and a wider bond is not automatically")
    print("  worse; the policy is a deterministic screening rule only.")


def main() -> None:
    member_1 = aluminium_like_member(100.0)
    member_2 = titanium_like_member(100.0)
    environment = radiator_joint_environment()
    geometry = radiator_joint_overlap()
    adhesive = ADHESIVE_LIKE
    yield_basis = YieldBasis(YIELD_DESIGN_FACTOR)
    basis = AdhesiveShearBasis(ADHESIVE_DESIGN_FACTOR)

    print()
    print("=" * 100)
    print("  BONDLINE DESIGN TRADE FOR THE BIMETALLIC JOINT - MILESTONE 4 STUDY")
    print("=" * 100)

    print_study_basis(member_1, member_2, environment, geometry, adhesive, basis)

    baseline = evaluate_joint_design(
        member_1, member_2, environment, geometry, adhesive, yield_basis, basis
    )
    print_baseline(baseline)
    print_scaling(member_1, member_2, geometry, adhesive, baseline)

    print_candidate_sweep(
        "WIDTH SENSITIVITY  (overlap, thickness, modulus, members, dT all fixed)",
        "b [mm]",
        bond_width_sensitivity(
            member_1, member_2, environment, geometry, adhesive,
            [w * MM for w in WIDTHS_MM], yield_basis, basis,
        ),
        lambda c: c.bond_width / MM, "{:.1f}",
    )
    print()
    print("  tau_peak ~ b^(-1/2) once the overlap is long: widening the bond helps,")
    print("  but only as a square root, so 10-100 mm is not enough on its own.")

    print_candidate_sweep(
        "THICKNESS SENSITIVITY  (G_a fixed at the canonical 1.0 GPa)",
        "t_a [mm]",
        adhesive_thickness_sensitivity(
            member_1, member_2, environment, geometry, adhesive,
            [t * MM for t in THICKNESSES_MM], yield_basis, basis,
        ),
        lambda c: c.adhesive_thickness / MM, "{:.2f}",
    )

    print_candidate_sweep(
        "MODULUS SENSITIVITY  (t_a fixed at the canonical 0.2 mm)",
        "G_a [GPa]",
        adhesive_modulus_sensitivity(
            member_1, member_2, environment, geometry, adhesive,
            [g * GPA for g in MODULI_GPA], yield_basis, basis,
        ),
        lambda c: c.adhesive_shear_modulus / GPA, "{:.2f}",
    )

    print_point_sweep(
        "AREA-RATIO SENSITIVITY  (Milestone 1 recomputed end to end at each ratio)",
        "A_1/A_2",
        area_ratio_sensitivity(
            member_1, member_2, environment, geometry, adhesive, AREA_RATIOS,
            yield_basis, basis,
        ),
        "{:.2f}",
    )
    print()
    print("  Equal areas are NOT optimal for adhesive shear: a smaller aluminium-like")
    print("  area lowers the Milestone 1 mismatch force faster than it raises beta.")

    print_point_sweep(
        "THERMAL-EXCURSION SENSITIVITY  (symmetric +/-|dT| about T_ref)",
        "|dT| [K]",
        thermal_excursion_sensitivity(
            member_1, member_2, geometry, adhesive, EXCURSIONS,
            reference_temperature=environment.reference_temperature,
            yield_basis=yield_basis, adhesive_basis=basis,
        ),
        "{:.0f}",
    )
    allowable_excursion = allowable_temperature_change_for_adhesive_shear(
        member_1, member_2, geometry, adhesive, basis
    )
    print()
    print(f"  Closed-form allowable |dT| = {allowable_excursion:.4f} K, which is exactly where")
    print("  the sweep crosses from PASS to FAIL.")

    print_point_sweep(
        "CTE-MISMATCH SENSITIVITY  (alpha_1 - alpha_2 scaled; synthetic materials)",
        "mismatch x",
        cte_mismatch_sensitivity(
            member_1, member_2, environment, geometry, adhesive, MISMATCH_SCALES,
            yield_basis, basis,
        ),
        "{:.2f}",
    )
    print()
    print("  Peak shear is exactly linear in |alpha_1 - alpha_2| - the only lever here")
    print("  that is first order rather than square-root.")

    print_design_factor_sensitivity(member_1, member_2, environment, geometry, adhesive)
    print_inverse(member_1, member_2, environment, geometry, adhesive, basis)

    canonical_map = width_thickness_design_map(
        member_1, member_2, environment, geometry, adhesive,
        [w * MM for w in MAP_WIDTHS_MM], [t * MM for t in MAP_THICKNESSES_MM],
        yield_basis, basis,
    )
    print_design_map(canonical_map, MAP_WIDTHS_MM, MAP_THICKNESSES_MM,
                     "canonical adhesive, G_a = 1.0 GPa")

    compliant = AdhesiveMaterial(
        name="Structural-adhesive-like, compliant variant (illustrative)",
        shear_modulus=0.1 * GPA,
        shear_strength=adhesive.shear_strength,
        source_note=ILLUSTRATIVE_ADHESIVE_NOTE,
        notes="Illustrative compliant variant used only to demonstrate co-design.",
    )
    compliant_map = width_thickness_design_map(
        member_1, member_2, environment, geometry, compliant,
        [w * MM for w in MAP_WIDTHS_MM], [t * MM for t in MAP_THICKNESSES_MM],
        yield_basis, basis,
    )
    print_design_map(compliant_map, MAP_WIDTHS_MM, MAP_THICKNESSES_MM,
                     "illustrative compliant variant, G_a = 0.1 GPa")
    print("  The canonical adhesive is not replaced by this variant; the second map")
    print("  is shown alongside the first purely to demonstrate co-design.")

    print_selection(canonical_map, baseline)

    rule()
    print(
        "Milestone 4 uses the verified one-dimensional shear-lag model to identify which\n"
        "first-order design levers can recover adhesive-shear feasibility. It remains a\n"
        "screening model and does not include peel stress, edge singularities, nonlinear\n"
        "adhesive response, fracture, creep or fatigue."
    )
    print()


if __name__ == "__main__":
    main()

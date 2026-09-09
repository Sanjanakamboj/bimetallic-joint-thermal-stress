"""Representative Milestone 1 sanity study for a bimetallic radiator joint.

Runs the perfectly bonded, uniform-temperature axial model over an illustrative
aluminium-like / titanium-like joint at hot and cold extremes, prints the full
state at each extreme, identifies the governing case and shows an area-ratio
sensitivity sweep.

Run with::

    python examples/bimetallic_joint_sanity.py

All material properties are ILLUSTRATIVE MATERIAL INPUT - NOT DESIGN ALLOWABLE.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from thermal_joint import (  # noqa: E402
    AxialMember,
    JointYieldAssessment,
    ThermalJointResult,
    YieldBasis,
    assess_temperature_extremes,
    assess_yield,
    solve_bimetallic_joint,
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
AREA_RATIOS = (0.25, 0.5, 1.0, 2.0, 4.0)
BASE_AREA_MM2 = 100.0


def rule(title: str = "") -> None:
    """Print a section rule with an optional title."""
    if title:
        print()
        print(title)
    print("-" * 78)


def format_margin(margin: float) -> str:
    """Render a margin of safety, including the infinite zero-stress case."""
    return "inf (zero stress, non-governing)" if math.isinf(margin) else f"{margin:+.3f}"


def print_study_basis(member_1: AxialMember, member_2: AxialMember, environment, basis) -> None:
    rule("STUDY BASIS")
    print(f"  {ILLUSTRATIVE_NOTE}")
    print()
    print(f"  member 1: {member_1.label}")
    print(f"  member 2: {member_2.label}")
    print()
    print(f"  {'Property':<34}{'Member 1':>20}{'Member 2':>20}")
    for caption, value_1, value_2 in (
        (
            "Young's modulus E [GPa]",
            member_1.material.elastic_modulus / GPA,
            member_2.material.elastic_modulus / GPA,
        ),
        (
            "CTE alpha [1e-6 / K]",
            member_1.material.thermal_expansion_coefficient / MICROSTRAIN,
            member_2.material.thermal_expansion_coefficient / MICROSTRAIN,
        ),
        (
            "yield strength [MPa]",
            member_1.material.yield_strength / MPA,
            member_2.material.yield_strength / MPA,
        ),
        (
            "cross-sectional area [mm^2]",
            member_1.area / SQUARE_MILLIMETRE,
            member_2.area / SQUARE_MILLIMETRE,
        ),
        (
            "axial stiffness E*A [MN]",
            member_1.axial_stiffness / 1.0e6,
            member_2.axial_stiffness / 1.0e6,
        ),
    ):
        print(f"  {caption:<34}{value_1:>20.4g}{value_2:>20.4g}")

    print()
    print(f"  reference temperature T_ref     {environment.reference_temperature:>10.1f} degC")
    print(f"  cold temperature      T_cold    {environment.cold_temperature:>10.1f} degC")
    print(f"  hot temperature       T_hot     {environment.hot_temperature:>10.1f} degC")
    print(f"  cold excursion        dT_cold   {environment.delta_temperature_cold:>10.1f} K")
    print(f"  hot excursion         dT_hot    {environment.delta_temperature_hot:>10.1f} K")
    print(f"  yield design factor             {basis.design_factor:>10.2f}")


def print_case(
    title: str, result: ThermalJointResult, assessment: JointYieldAssessment
) -> None:
    rule(title)
    print(f"  temperature change dT               {result.delta_temperature:>16.4f} K")
    print(
        "  free thermal strain, member 1       "
        f"{result.member_1_free_thermal_strain / MICROSTRAIN:>16.4f} microstrain"
    )
    print(
        "  free thermal strain, member 2       "
        f"{result.member_2_free_thermal_strain / MICROSTRAIN:>16.4f} microstrain"
    )
    print(
        "  common joint strain                 "
        f"{result.common_strain / MICROSTRAIN:>16.4f} microstrain"
    )
    print(
        f"  stress, member 1                    {result.member_1_stress / MPA:>16.4f} MPa"
        f"  ({'tension' if result.member_1_stress > 0 else 'compression'})"
    )
    print(
        f"  stress, member 2                    {result.member_2_stress / MPA:>16.4f} MPa"
        f"  ({'tension' if result.member_2_stress > 0 else 'compression'})"
    )
    print(f"  internal force, member 1            {result.member_1_force:>16.4f} N")
    print(f"  internal force, member 2            {result.member_2_force:>16.4f} N")
    print(f"  force equilibrium residual N1+N2    {result.force_equilibrium_residual:>16.3e} N")
    print()
    for index, margin in (
        (1, assessment.member_1_margin),
        (2, assessment.member_2_margin),
    ):
        print(
            f"  member {index} allowable                  "
            f"{margin.allowable_stress / MPA:>16.4f} MPa"
        )
        print(
            f"  member {index} preliminary yield margin    "
            f"{format_margin(margin.margin_of_safety):>16}"
            f"   {'PASS' if margin.passes else 'FAIL'}"
        )
    print(
        f"  governing member at this extreme    {assessment.governing_member:>16d}"
        f"   ({assessment.governing_margin.label})"
    )


def print_governing_result(assessment) -> None:
    rule("GOVERNING RESULT")
    governing = assessment.governing_assessment.governing_margin
    print(f"  governing extreme                   {assessment.governing_extreme:>16}")
    print(
        f"  governing member                    {assessment.governing_member:>16d}"
        f"   ({governing.label})"
    )
    print(f"  governing stress                    {governing.stress / MPA:>16.4f} MPa")
    print(f"  governing allowable                 {governing.allowable_stress / MPA:>16.4f} MPa")
    print(
        "  minimum preliminary yield margin    "
        f"{format_margin(assessment.minimum_yield_margin):>16}"
    )
    print(f"  verdict                             {'PASS' if assessment.passes else 'FAIL':>16}")


def print_area_ratio_sensitivity(member_2: AxialMember, delta_temperature: float) -> None:
    rule(f"SENSITIVITY - area ratio A1/A2 at the governing extreme (dT = {delta_temperature:+.1f} K)")
    print(
        f"  {'A1/A2':>7}{'A1 [mm^2]':>12}{'eps_common':>16}{'sigma_1 [MPa]':>16}"
        f"{'sigma_2 [MPa]':>16}{'|N1+N2| [N]':>14}"
    )
    for ratio in AREA_RATIOS:
        member_1 = aluminium_like_member(ratio * BASE_AREA_MM2)
        result = solve_bimetallic_joint(member_1, member_2, delta_temperature)
        print(
            f"  {ratio:>7.2f}{member_1.area / SQUARE_MILLIMETRE:>12.1f}"
            f"{result.common_strain / MICROSTRAIN:>13.2f} ue"
            f"{result.member_1_stress / MPA:>16.4f}{result.member_2_stress / MPA:>16.4f}"
            f"{abs(result.force_equilibrium_residual):>14.2e}"
        )
    print()
    print("  The common strain moves toward the free thermal strain of the axially")
    print("  stiffer (higher E*A) member; the stresses redistribute while the internal")
    print("  forces stay equal and opposite. Stress is not area-independent.")


def main() -> None:
    member_1 = aluminium_like_member(BASE_AREA_MM2)
    member_2 = titanium_like_member(BASE_AREA_MM2)
    environment = radiator_joint_environment()
    basis = YieldBasis(DESIGN_FACTOR)

    print()
    print("=" * 78)
    print("  BIMETALLIC JOINT THERMAL STRESS - MILESTONE 1 SANITY STUDY")
    print("=" * 78)

    print_study_basis(member_1, member_2, environment, basis)

    assessment = assess_temperature_extremes(member_1, member_2, environment, basis)
    print_case("HOT CASE", assessment.hot_result, assessment.hot_assessment)
    print_case("COLD CASE", assessment.cold_result, assessment.cold_assessment)
    print_governing_result(assessment)
    print_area_ratio_sensitivity(member_2, assessment.governing_result.delta_temperature)

    rule()
    print(
        "  This is a perfectly bonded, uniform-temperature axial compatibility model;\n"
        "  interface shear, finite joint length, thermal gradients, plasticity and\n"
        "  fatigue are not included."
    )
    print(f"  {ILLUSTRATIVE_NOTE}.")
    print()


if __name__ == "__main__":
    main()

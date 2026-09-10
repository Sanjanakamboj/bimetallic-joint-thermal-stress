"""Generate the STM-11 portfolio figures.

Every plotted engineering quantity comes from the ``thermal_joint`` package
APIs. Only the canonical study inputs (geometries, sweep ranges, the design
points already selected in the study) are written here as explicit constants.

Run with::

    python figures/make_figures.py

Requires the optional plotting extra::

    pip install -e ".[figures]"

PNG metadata is pinned so repeated runs on the same matplotlib/FreeType build
produce byte-identical files.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from thermal_joint import (  # noqa: E402
    AdhesiveShearBasis,
    AxialRestraint,
    YieldBasis,
    assess_shear_lag_extremes,
    assess_thermal_cycle_fatigue,
    joint_axial_stiffness,
    maximum_allowable_restraint_stiffness,
    solve_restrained_joint,
    solve_shear_lag,
    zero_stress_restraint_stiffness,
)
from thermal_joint.illustrative import (  # noqa: E402
    ADHESIVE_LIKE,
    aluminium_like_member,
    illustrative_cycle_requirement,
    illustrative_fatigue_curves,
    radiator_joint_environment,
    radiator_joint_overlap,
    titanium_like_member,
)

MPA, GPA, MM = 1.0e6, 1.0e9, 1.0e-3
HERE = Path(__file__).resolve().parent
PNG_METADATA = {"Software": None}
DPI = 160

MEMBER_1 = aluminium_like_member(100.0)
MEMBER_2 = titanium_like_member(100.0)
ENVIRONMENT = radiator_joint_environment()
ADHESIVE = ADHESIVE_LIKE
YIELD_BASIS = YieldBasis(1.25)
ADHESIVE_BASIS = AdhesiveShearBasis(1.25)
CURVES = illustrative_fatigue_curves()
REQUIREMENT = illustrative_cycle_requirement()

# The three canonical design points established by the study.
DESIGNS = (
    ("M3 baseline", 20.0, 0.2),
    ("M4 static", 50.0, 1.0),
    ("M5 fatigue", 30.0, 2.0),
)
ALLOWABLE_SHEAR = ADHESIVE_BASIS.allowable_shear_stress(ADHESIVE)

DISCLAIMER = ("Illustrative properties - not design allowables. "
              "Preliminary screening model, not a certified joint.")


def finish(fig, name: str) -> Path:
    """Add the standing disclaimer, save deterministically and report."""
    fig.text(0.5, 0.012, DISCLAIMER, ha="center", va="bottom", fontsize=7.5,
             color="0.35")
    path = HERE / name
    fig.savefig(path, dpi=DPI, metadata=PNG_METADATA)
    plt.close(fig)
    print(f"  wrote {path.name}")
    return path


def figure_1_restraint() -> None:
    """Member stress at the cold extreme against external restraint ratio."""
    stiffness = joint_axial_stiffness(MEMBER_1, MEMBER_2)
    crossing = zero_stress_restraint_stiffness(MEMBER_1, MEMBER_2, 100.0, 2) / stiffness
    boundary = maximum_allowable_restraint_stiffness(
        MEMBER_1, MEMBER_2, ENVIRONMENT, YIELD_BASIS
    ).maximum_stiffness_ratio
    allowable_1 = YIELD_BASIS.allowable_stress(MEMBER_1.material)

    ratios = [i * 15.0 / 600 for i in range(601)]
    cold = ENVIRONMENT.delta_temperature_cold
    results = [
        solve_restrained_joint(
            MEMBER_1, MEMBER_2, cold,
            AxialRestraint.from_stiffness_ratio(r, MEMBER_1, MEMBER_2),
        )
        for r in ratios
    ]

    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ax.plot(ratios, [r.member_1_stress / MPA for r in results],
            color="#c2352b", lw=2.0, label="Member 1, aluminium-like (high CTE)")
    ax.plot(ratios, [r.member_2_stress / MPA for r in results],
            color="#1f5fa8", lw=2.0, label="Member 2, titanium-like (low CTE)")

    ax.axhline(0.0, color="0.4", lw=0.9)
    ax.axhline(allowable_1 / MPA, color="#c2352b", ls=":", lw=1.4,
               label=f"Member 1 yield allowable, {allowable_1 / MPA:.0f} MPa")
    ax.axvline(crossing, color="#1f5fa8", ls="--", lw=1.2)
    ax.axvline(boundary, color="0.3", ls="-.", lw=1.2)

    ax.annotate(f"member 2 stress crosses zero\n$\\eta_r$ = {crossing:.3f}",
                xy=(crossing, 0.0), xytext=(1.9, -46),
                fontsize=8.5, color="#1f5fa8",
                arrowprops=dict(arrowstyle="->", color="#1f5fa8", lw=1.0))
    ax.annotate(f"member 1 reaches its yield\nallowable, $\\eta_r$ = {boundary:.2f}",
                xy=(boundary, allowable_1 / MPA), xytext=(7.4, 128),
                fontsize=8.5, color="0.2",
                arrowprops=dict(arrowstyle="->", color="0.3", lw=1.0))
    ax.plot([0.0], [results[0].member_1_stress / MPA], "o", color="#c2352b", ms=6)
    ax.plot([0.0], [results[0].member_2_stress / MPA], "o", color="#1f5fa8", ms=6)
    ax.annotate("free joint ($\\eta_r$ = 0):\nopposite signs", xy=(0.0, -86.8),
                xytext=(1.0, -108), fontsize=8.5, color="0.2",
                arrowprops=dict(arrowstyle="->", color="0.4", lw=1.0))

    ax.set_xlabel("Restraint stiffness ratio  $\\eta_r = K_r / (E_1A_1 + E_2A_2)$")
    ax.set_ylabel("Member axial stress at the cold extreme [MPa]")
    ax.set_title("External restraint changes the sign and magnitude of member stress\n"
                 "$\\Delta T_{cold}$ = -140 K, free joint through to near-rigid restraint",
                 fontsize=10.5)
    ax.set_xlim(0, 15)
    ax.set_ylim(-130, 260)
    ax.grid(alpha=0.25)
    ax.legend(loc="lower right", fontsize=8.5, framealpha=0.95)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    finish(fig, "fig1_restraint_stress.png")


def figure_2_shear_distribution() -> None:
    """Signed adhesive shear along the baseline overlap, hot and cold."""
    geometry = radiator_joint_overlap()
    length = geometry.overlap_length
    stations = [length * i / 400 for i in range(401)]

    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for label, delta_t, colour in (
        (f"Hot, $\\Delta T$ = +{ENVIRONMENT.delta_temperature_hot:.0f} K",
         ENVIRONMENT.delta_temperature_hot, "#c2352b"),
        (f"Cold, $\\Delta T$ = {ENVIRONMENT.delta_temperature_cold:.0f} K",
         ENVIRONMENT.delta_temperature_cold, "#1f5fa8"),
    ):
        demand = solve_shear_lag(MEMBER_1, MEMBER_2, geometry, ADHESIVE, delta_t)
        ax.plot([x / MM for x in stations],
                [demand.shear_stress(x) / MPA for x in stations],
                color=colour, lw=2.0, label=label)

    cold = solve_shear_lag(MEMBER_1, MEMBER_2, geometry, ADHESIVE,
                           ENVIRONMENT.delta_temperature_cold)
    hot = solve_shear_lag(MEMBER_1, MEMBER_2, geometry, ADHESIVE,
                          ENVIRONMENT.delta_temperature_hot)

    ax.axhline(0.0, color="0.4", lw=0.9)
    for sign in (1.0, -1.0):
        ax.axhline(sign * ALLOWABLE_SHEAR / MPA, color="0.25", ls=":", lw=1.3)
    ax.text(20.5, ALLOWABLE_SHEAR / MPA + 3.0,
            f"$\\pm$ allowable shear, {ALLOWABLE_SHEAR / MPA:.0f} MPa",
            fontsize=8.5, color="0.25")

    ax.axvline(length / MM, color="0.55", ls="--", lw=1.1)
    ax.annotate(f"free edge: cold peak {cold.peak_shear_stress / MPA:.1f} MPa\n"
                f"({cold.peak_shear_stress / ALLOWABLE_SHEAR:.2f}x the allowable)",
                xy=(length / MM, -cold.peak_shear_stress / MPA),
                xytext=(16.0, -64), fontsize=8.5, color="#1f5fa8",
                arrowprops=dict(arrowstyle="->", color="#1f5fa8", lw=1.0))
    ax.annotate(f"hot peak {hot.peak_shear_stress / MPA:.1f} MPa",
                xy=(length / MM, hot.peak_shear_stress / MPA),
                xytext=(26.0, 62), fontsize=8.5, color="#c2352b",
                arrowprops=dict(arrowstyle="->", color="#c2352b", lw=1.0))
    ax.annotate("loaded transfer plane:\njoint interior is\nnearly shear-free",
                xy=(1.2, 1.5), xytext=(3.0, 30), fontsize=8.5, color="0.25",
                arrowprops=dict(arrowstyle="->", color="0.45", lw=1.0))

    ax.set_xlabel("Position along the bonded overlap  $x$ [mm]   "
                  "($x$ = 0 loaded plane, $x = L_b$ free edge)")
    ax.set_ylabel("Signed adhesive shear stress  $\\tau(x)$ [MPa]")
    ax.set_title("Adhesive shear concentrates at the free edge of the overlap\n"
                 "Baseline bondline: 40 x 20 mm, 0.2 mm bondline, "
                 "$G_a$ = 1.0 GPa, transfer length 6.54 mm", fontsize=10.5)
    ax.set_xlim(0, length / MM * 1.03)
    ax.set_ylim(-78, 78)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper left", fontsize=8.5, framealpha=0.95)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    finish(fig, "fig2_shear_distribution.png")


def figure_3_static_trade() -> None:
    """Cold peak shear against bond width for several bondline thicknesses."""
    widths_mm = [10.0 + i * 0.5 for i in range(221)]
    thicknesses_mm = (0.2, 0.5, 1.0, 2.0)
    colours = ("#7a1f6b", "#c2352b", "#d38b16", "#1f5fa8")

    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for thickness_mm, colour in zip(thicknesses_mm, colours):
        peaks = []
        for width_mm in widths_mm:
            geometry = radiator_joint_overlap(40.0, width_mm, thickness_mm)
            peaks.append(
                solve_shear_lag(MEMBER_1, MEMBER_2, geometry, ADHESIVE,
                                ENVIRONMENT.delta_temperature_cold).peak_shear_stress / MPA
            )
        ax.plot(widths_mm, peaks, color=colour, lw=1.9,
                label=f"$t_a$ = {thickness_mm:g} mm")

    ax.axhline(ALLOWABLE_SHEAR / MPA, color="0.2", ls=":", lw=1.5)
    ax.text(78, ALLOWABLE_SHEAR / MPA + 1.6,
            f"static allowable {ALLOWABLE_SHEAR / MPA:.0f} MPa", fontsize=8.5, color="0.2")

    markers = {"M3 baseline": ("s", "#7a1f6b", (28.0, 78.0)),
               "M4 static": ("o", "#d38b16", (58.0, 5.0)),
               "M5 fatigue": ("D", "#1f5fa8", (14.0, 5.0))}
    for label, width_mm, thickness_mm in DESIGNS:
        geometry = radiator_joint_overlap(40.0, width_mm, thickness_mm)
        peak = solve_shear_lag(MEMBER_1, MEMBER_2, geometry, ADHESIVE,
                               ENVIRONMENT.delta_temperature_cold).peak_shear_stress / MPA
        marker, colour, text_at = markers[label]
        ax.plot([width_mm], [peak], marker, color=colour, ms=9,
                markeredgecolor="white", markeredgewidth=1.2, zorder=5)
        ax.annotate(f"{label}\n{width_mm:.0f} mm x {thickness_mm:g} mm",
                    xy=(width_mm, peak), xytext=text_at,
                    fontsize=8.5, color=colour,
                    arrowprops=dict(arrowstyle="->", color=colour, lw=1.0))

    ax.set_xlabel("Bond width  $b$ [mm]")
    ax.set_ylabel("Cold peak adhesive shear  $\\tau_{peak}$ [MPa]")
    ax.set_title("Width and bondline compliance must be traded together\n"
                 "$\\tau_\\infty \\propto b^{-1/2} t_a^{-1/2}$, so no single lever is cheap "
                 "($G_a$ = 1.0 GPa, 40 mm overlap)", fontsize=10.5)
    ax.set_xlim(10, 120)
    ax.set_ylim(0, 100)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", fontsize=8.5, framealpha=0.95)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    finish(fig, "fig3_static_design_trade.png")


def figure_4_static_versus_fatigue() -> None:
    """The central conclusion: passing static shear did not guarantee fatigue life."""
    labels, utilisation, life_ratio = [], [], []
    for label, width_mm, thickness_mm in DESIGNS:
        geometry = radiator_joint_overlap(40.0, width_mm, thickness_mm)
        shear = assess_shear_lag_extremes(
            MEMBER_1, MEMBER_2, ENVIRONMENT, geometry, ADHESIVE, ADHESIVE_BASIS
        )
        fatigue = assess_thermal_cycle_fatigue(
            MEMBER_1, MEMBER_2, ENVIRONMENT, geometry, ADHESIVE, CURVES, REQUIREMENT,
            YIELD_BASIS, ADHESIVE_BASIS,
        )
        labels.append(f"{label}\n{width_mm:.0f} mm x {thickness_mm:g} mm")
        utilisation.append(shear.governing_peak_shear_stress / ALLOWABLE_SHEAR)
        life_ratio.append(fatigue.minimum_life_ratio)

    positions = range(len(labels))
    fig, (left, right) = plt.subplots(1, 2, figsize=(9.6, 4.6))

    bars = left.bar(positions, utilisation,
                    color=["#c2352b" if u > 1.0 else "#3f8f52" for u in utilisation],
                    width=0.58)
    left.axhline(1.0, color="0.2", ls="--", lw=1.4)
    left.text(2.42, 1.06, "limit", fontsize=8.5, color="0.2", ha="right")
    for x, (bar, value) in enumerate(zip(bars, utilisation)):
        left.text(x, value + 0.09, f"{value:.2f}", ha="center", fontsize=9,
                  fontweight="bold")
    left.set_xticks(list(positions))
    left.set_xticklabels(labels, fontsize=8.5)
    left.set_ylabel("Static shear utilisation  $\\tau_{peak}\\,/\\,\\tau_{allow}$")
    left.set_title("Static adhesive shear screen", fontsize=10.5)
    left.set_ylim(0, 4.0)
    left.grid(alpha=0.25, axis="y")

    bars = right.bar(positions, life_ratio,
                     color=["#c2352b" if r < 1.0 else "#3f8f52" for r in life_ratio],
                     width=0.58)
    right.axhline(1.0, color="0.2", ls="--", lw=1.4)
    right.text(-0.42, 1.45, "required life", fontsize=8.5, color="0.2", ha="left")
    right.set_yscale("log")
    for x, (bar, value) in enumerate(zip(bars, life_ratio)):
        right.text(x, value * 1.5, f"{value:.3g}", ha="center", fontsize=9,
                   fontweight="bold")
    right.set_xticks(list(positions))
    right.set_xticklabels(labels, fontsize=8.5)
    right.set_ylabel("Fatigue life ratio  $N_f\\,/\\,N_{required}$  (log)")
    right.set_title(f"Thermal-cycle fatigue screen "
                    f"({REQUIREMENT.required_cycles:,.0f} cycles)", fontsize=10.5)
    right.set_ylim(5e-5, 30)
    right.grid(alpha=0.25, axis="y")

    fig.suptitle("Passing the static adhesive screen did not guarantee fatigue life",
                 fontsize=11.5, y=0.985)
    fig.tight_layout(rect=(0, 0.038, 1, 0.945))
    finish(fig, "fig4_static_versus_fatigue.png")


def main() -> None:
    print("Generating STM-11 portfolio figures")
    figure_1_restraint()
    figure_2_shear_distribution()
    figure_3_static_trade()
    figure_4_static_versus_fatigue()
    print("Done.")


if __name__ == "__main__":
    main()

# Bimetallic Joint Thermal Stress

Thermal-stress analysis of a spacecraft radiator joint made from two dissimilar
metallic members subjected to hot and cold temperature excursions.

**Status: Milestone 1 of an ongoing portfolio project — not complete.**
Milestone 1 delivers only the verified mechanics foundation: the closed-form
axial thermal stress in a perfectly bonded two-member joint under a uniform
temperature change, preliminary elastic yield margins, and hot/cold extreme
evaluation. Later milestones will extend the physics; see
[Limitations](#limitations) for everything deliberately excluded so far.

---

## Objective

> If two dissimilar members are bonded so that they must undergo the same axial
> strain, what common strain, internal force, and thermal stress develop when
> their coefficients of thermal expansion differ?

The primary deliverable of the overall project is a thermal stress calculation
plus a margin at the temperature extremes. Milestone 1 builds and independently
verifies the mechanics that answer the question above.

---

## Assumptions

The Milestone 1 model is a one-dimensional axial compatibility model:

- two members, **perfectly bonded**, so both undergo the same total axial strain
- **uniform temperature** throughout both members
- a **common stress-free reference state** `T_ref` for both members
- **no external axial force** on the joint (self-equilibrating internal forces)
- **linear elasticity**, isotropic, temperature-independent properties
- no bending, no interface slip, no adhesive compliance, no thermal gradient

Because both members share the same reference length, length cancels out of the
compatibility equation and is deliberately not modelled.

---

## Units and sign conventions

SI units are used internally:

| Quantity | Symbol | Unit |
|---|---|---|
| Young's modulus | `E` | Pa |
| Thermal expansion coefficient | `alpha` | 1/K |
| Cross-sectional area | `A` | m² |
| Strain | `eps` | – |
| Stress | `sigma` | Pa |
| Internal force | `N` | N |
| Temperature change | `dT` | K (or °C increment) |

Sign conventions:

- **tensile stress is positive**, compressive stress is negative
- `dT = T - T_ref`; **positive `dT` means heating**
- **positive `alpha` means expansion on heating**, so `alpha * dT > 0` on heating
- `yield_strength` is stored as a **positive magnitude** and compared against `abs(sigma)`

Temperature *differences* may be given in K or in °C: the increments are
numerically identical, so no conversion is performed. Only differences enter the
mechanics, so an entire `ThermalEnvironment` may be expressed on either scale as
long as one scale is used consistently.

---

## Governing equations

**Free (unrestrained) thermal strain** of member `i`:

```
eps_th_i = alpha_i * dT
```

**Compatibility** (perfect bond):

```
eps_1 = eps_2 = eps_common
```

**Constitutive relation** — stress arises only from the difference between the
enforced total strain and the free thermal strain:

```
sigma_i = E_i * (eps_common - alpha_i * dT)
```

**Axial force equilibrium** for a joint with no external axial load:

```
sigma_1 * A_1 + sigma_2 * A_2 = 0
```

Solving the three together gives the **stiffness-weighted common strain**,
implemented directly as an exact closed form (never solved numerically):

```
              E_1 A_1 alpha_1 + E_2 A_2 alpha_2
eps_common =  ---------------------------------  *  dT
                    E_1 A_1 + E_2 A_2
```

**Internal forces** follow as `N_i = sigma_i * A_i`, and the model reports the
equilibrium residual `N_1 + N_2` explicitly rather than discarding it; it is
zero to machine precision.

### Independent closed form (used for verification)

The same problem can be solved for the stresses directly in terms of the CTE
mismatch. These expressions are used in the **tests** as an independent check;
they are never used by the production code path:

```
sigma_1 = [ E_1 E_2 A_2 / (E_1 A_1 + E_2 A_2) ] * (alpha_2 - alpha_1) * dT
sigma_2 = [ E_1 E_2 A_1 / (E_1 A_1 + E_2 A_2) ] * (alpha_1 - alpha_2) * dT
```

### Physical reading

For a high-CTE / low-CTE pair (e.g. aluminium-like vs titanium-like), on
**heating** the high-CTE member wants to expand more, the bond restrains it, and
it goes into **compression** while the low-CTE member goes into **tension**.
**Cooling reverses both signs.** The tests assert these signs from the computed
result rather than hard-coding them.

Note that equilibrium requires equal and opposite *forces*, not equal and
opposite *stresses*. In general:

```
sigma_1 / sigma_2 = -A_2 / A_1
```

which reduces to `-1` only when the areas are equal.

---

## API and result structures

| Object | Purpose |
|---|---|
| `ThermoelasticMaterial` | name, `elastic_modulus`, `thermal_expansion_coefficient`, `yield_strength`, optional `density` |
| `AxialMember` | a `material`, an `area`, an optional `label`; exposes `axial_stiffness` (`E*A`) |
| `common_joint_strain(...)` | the stiffness-weighted common strain |
| `solve_bimetallic_joint(...)` | full solution → `ThermalJointResult` |
| `ThermalJointResult` | `delta_temperature`, `common_strain`, both free thermal strains, both stresses, both internal forces, `force_equilibrium_residual` |
| `YieldBasis` | explicit design basis holding `design_factor` (≥ 1) |
| `assess_yield(...)` | → `JointYieldAssessment` (two `MemberYieldMargin` records, `governing_member`, `minimum_margin`) |
| `ThermalEnvironment` | `reference_temperature`, `cold_temperature`, `hot_temperature`; derives `delta_temperature_cold` / `delta_temperature_hot` |
| `assess_temperature_extremes(...)` | → `TemperatureExtremeAssessment` (both results, both assessments, `governing_extreme`, `governing_member`, `minimum_yield_margin`) |

### Minimal usage

```python
from thermal_joint import (
    AxialMember, ThermalEnvironment, ThermoelasticMaterial,
    YieldBasis, assess_temperature_extremes,
)

aluminium_like = ThermoelasticMaterial("Al-like", 70e9, 23e-6, 270e6)
titanium_like = ThermoelasticMaterial("Ti-like", 110e9, 8.5e-6, 830e6)

joint = (AxialMember(aluminium_like, 100e-6), AxialMember(titanium_like, 100e-6))
environment = ThermalEnvironment(reference_temperature=20.0,
                                 cold_temperature=-120.0,
                                 hot_temperature=120.0)

assessment = assess_temperature_extremes(*joint, environment, YieldBasis(1.25))
print(assessment.governing_extreme, assessment.minimum_yield_margin)
```

---

## Yield-margin convention

Margins are **preliminary elastic yield margins**, not certification margins.

```
sigma_allow_i = yield_strength_i / design_factor
MS_yield_i    = sigma_allow_i / abs(sigma_i) - 1
```

- the **design factor lives in `YieldBasis`, never inside the material**, so
  material data and design policy stay separable; `design_factor` must be ≥ 1,
  and the default of 1.0 leaves the allowable equal to the raw yield strength
- a margin of **exactly zero passes** (`MS >= 0`)
- a member carrying **zero stress** returns `math.inf` — an explicitly
  non-governing margin, never a division by zero
- the **governing member is the one with the smaller margin**, computed rather
  than assumed; the lower-yield material does *not* automatically govern,
  because area and modulus redistribute the stress. An exact tie resolves to
  member 1.

No knockdowns, statistical basis, joint efficiency, plasticity or fatigue
effects are included.

---

## Thermal-extreme logic

`ThermalEnvironment` validates that all three temperatures are finite and that
`cold <= reference <= hot` (equality allowed, which makes that excursion a
zero-`dT` case). `assess_temperature_extremes` then:

1. solves the joint at `dT_cold = T_cold - T_ref` and at `dT_hot = T_hot - T_ref`
2. computes a yield assessment at each extreme
3. selects the **governing extreme from the actual margins** — neither hot nor
   cold is assumed to govern; an exact tie resolves to `"cold"`
4. reports the governing member and the minimum margin over both extremes

With temperature-independent properties and **symmetric** excursions
(`|dT_hot| = |dT_cold|`), the stress magnitudes are identical and only the signs
reverse, so the two margins are equal. With **asymmetric** excursions the larger
`|dT|` governs. Both behaviours are covered by tests.

---

## Representative sanity result

`examples/bimetallic_joint_sanity.py` runs an illustrative equal-area joint.

> **ILLUSTRATIVE MATERIAL INPUT — NOT DESIGN ALLOWABLE.** The values below are
> clean, round, order-of-magnitude figures chosen to make the mechanics legible.
> They are not traceable to any specific alloy, temper, product form or
> statistical basis, so no alloy designation is claimed. The temperatures are an
> illustrative portfolio excursion, not spacecraft qualification temperatures.

| | Member 1 (aluminium-like) | Member 2 (titanium-like) |
|---|---|---|
| `E` | 70 GPa | 110 GPa |
| `alpha` | 23 ×10⁻⁶ /K | 8.5 ×10⁻⁶ /K |
| yield | 270 MPa | 830 MPa |
| area | 100 mm² | 100 mm² |

`T_ref = +20 °C`, `T_cold = -120 °C`, `T_hot = +120 °C` → `dT_cold = -140 K`,
`dT_hot = +100 K`; yield design factor 1.25.

| Quantity | Hot (`dT = +100 K`) | Cold (`dT = -140 K`) |
|---|---|---|
| common strain | +1413.9 µε | −1979.4 µε |
| `sigma_1` | **−62.03 MPa** (compression) | **+86.84 MPa** (tension) |
| `sigma_2` | **+62.03 MPa** (tension) | **−86.84 MPa** (compression) |
| `N_1` / `N_2` | ∓6202.8 N | ±8683.9 N |
| `N_1 + N_2` | ~1e−12 N | ~5e−12 N |
| MS member 1 | +2.482 | **+1.487** |
| MS member 2 | +9.705 | +6.646 |

**Governing case: cold, member 1 (aluminium-like), minimum preliminary elastic
yield margin +1.487 → PASS.** Cold governs because the cold excursion is the
larger one (140 K vs 100 K); member 1 governs because it has by far the lower
allowable, not because of any assumption in the code.

---

## Verification strategy

The implementation is not merely exercised, it is **verified** — the tests
compare it against results derived independently of the production code path:

- **independent closed-form stress equations** — the tests compute
  `sigma_i` directly from the CTE-mismatch form, never by calling the production
  common-strain helper, and compare across a grid of area pairs and `dT` values
- **exact force equilibrium** — `N_1 + N_2` is asserted to vanish to machine
  precision relative to the force magnitude, at every condition tested
- **limiting cases** — identical CTE gives exactly zero stress with the common
  strain equal to the shared free thermal strain; `dT = 0` gives an exactly zero
  state
- **heating/cooling symmetry** — reversing `dT` reverses every sign and
  preserves every magnitude
- **area-ratio identities** — `sigma_1/sigma_2 = -A_2/A_1` in general, reducing
  to `-1` only for equal areas
- **stiffness-weighted strain behaviour** — the common strain is bracketed by
  the two free thermal strains and moves monotonically toward the free strain of
  the higher-`E*A` member as either area or modulus is varied
- **scaling laws** — stress, strain and force are linear in `dT` and in the CTE
  mismatch `|alpha_1 - alpha_2|`
- **exact yield-margin boundaries** — a stress exactly at the allowable gives
  `MS == 0.0` exactly and passes; the zero-stress case returns `inf`
- **hand calculations** — a round-number case (`E` 100/200 GPa, `A` 200/100 mm²,
  `alpha` 20/10 ×10⁻⁶ /K, `dT` +50 K → `eps = 7.5e-4`, `sigma = -25/+50 MPa`,
  `N = ∓5000 N`) is checked against values worked out by hand
- **governing-selection tests** — cases are constructed in which the
  higher-yield member governs, and in which either extreme governs, proving
  neither is hard-coded
- **input validation sweep** — NaN, ±infinity, non-positive modulus/area/yield,
  design factor < 1 and malformed environment ordering are all rejected

Milestone 1 currently has **276 passing tests**.

---

## Limitations

Not modelled in Milestone 1, and **not** to be inferred from these results:

interface shear-lag · adhesive stresses · finite joint length effects ·
bolts and fasteners · contact and slip · plasticity · creep · fatigue ·
thermal gradients through thickness · temperature-dependent properties ·
nonlinear material behaviour · plate bending and warpage · detailed radiator
geometry · optimization · candidate material trade studies · portfolio figures

In short: *this is a perfectly bonded, uniform-temperature axial compatibility
model; interface shear, finite joint length, thermal gradients, plasticity and
fatigue are not included.*

The margins produced are preliminary elastic yield margins for a portfolio
study. They are not certification margins and the material properties are not
design allowables.

---

## Repository structure

```
.
├── pyproject.toml
├── README.md
├── src/thermal_joint/
│   ├── __init__.py          public API
│   ├── _validation.py       shared numeric input validation
│   ├── materials.py         ThermoelasticMaterial, AxialMember
│   ├── joint.py             common strain, stresses, forces, equilibrium residual
│   ├── margins.py           YieldBasis, margins, governing member
│   ├── environment.py       ThermalEnvironment
│   ├── extremes.py          hot/cold assessment and governing extreme
│   └── illustrative.py      illustrative inputs (NOT design allowables)
├── tests/
│   ├── conftest.py                      makes src/ importable without install
│   ├── reference_solution.py            independent closed-form check
│   ├── cases.py                         shared hand-calculation inputs
│   ├── test_materials_and_members.py    tests A–G
│   ├── test_free_thermal_strain.py      tests H–L
│   ├── test_joint_mechanics.py          tests M–U
│   ├── test_limits_and_scaling.py       tests V–AC
│   ├── test_yield_margins.py            tests AD–AK
│   ├── test_thermal_extremes.py         tests AL–AR
│   ├── test_invalid_inputs.py           invalid-input sweep
│   └── test_physical_interpretation.py  expected physics and regression lock
└── examples/
    └── bimetallic_joint_sanity.py
```

The model is pure Python standard library — the closed-form Milestone 1
mechanics does not require NumPy.

---

## Install, test, run

Requires Python 3.10 or newer. From the repository root:

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[test]"
```

Run the test suite:

```bash
pytest
```

Run the sanity study:

```bash
python examples/bimetallic_joint_sanity.py
```

The tests and the example both insert `src/` on `sys.path`, so they also run
directly from a clean checkout without installing the package.

---

## Licensing status

**No license has been chosen for this repository yet.** There is no `LICENSE`
file and no license metadata in `pyproject.toml`. Absent a license, default
copyright applies and no reuse rights are granted. A license will be added only
on an explicit decision by the repository owner.

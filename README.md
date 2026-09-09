# Bimetallic Joint Thermal Stress

Thermal-stress analysis of a spacecraft radiator joint made from two dissimilar
metallic members subjected to hot and cold temperature excursions.

**Status: Milestone 2 of an ongoing portfolio project — not complete.**

- **Milestone 1** — the verified mechanics foundation: closed-form axial thermal
  stress in a perfectly bonded two-member joint under a uniform temperature
  change, preliminary elastic yield margins, hot/cold extreme evaluation.
- **Milestone 2** — finite external axial restraint from the surrounding
  structure, the rigid-restraint limit, restraint sensitivity, and a
  preliminary inverse-design search for the maximum tolerable restraint
  stiffness. The Milestone 1 free-joint solution is preserved *exactly* as the
  zero-restraint limit.

Later milestones will extend the physics; see [Limitations](#limitations) for
everything deliberately excluded so far.

---

## Objective

> **Milestone 1.** If two dissimilar members are bonded so that they must
> undergo the same axial strain, what common strain, internal force, and thermal
> stress develop when their coefficients of thermal expansion differ?

> **Milestone 2.** How much do those stresses increase when the bonded pair is
> not perfectly free to expand, but is attached to surrounding spacecraft
> structure with finite axial stiffness?

The primary deliverable of the overall project is a thermal stress calculation
plus a margin at the temperature extremes. Each milestone builds and
independently verifies the mechanics that answer its question.

---

## Assumptions

The model is one-dimensional axial compatibility:

- two members, **perfectly bonded**, so both undergo the same total axial strain
- **uniform temperature** throughout both members
- a **common stress-free reference state** `T_ref` for both members
- **linear elasticity**, isotropic, temperature-independent properties
- no bending, no interface slip, no adhesive compliance, no thermal gradient

Because both members share the same reference length, length cancels out of the
compatibility equation and is deliberately not modelled.

One Milestone 1 assumption is **relaxed in Milestone 2**:

- Milestone 1 assumed **no external axial force**, so the members were
  self-equilibrating (`N_1 + N_2 = 0`).
- Milestone 2 adds a **single linear external axial restraint** of finite
  stiffness. The restraint is linear and elastic — no contact, slip, backlash
  or nonlinear springs.

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

## Representative Milestone 1 result

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

## Milestone 2 — external axial restraint

### Restraint model

The surrounding structure is represented as a **single linear axial spring**
acting on the common joint strain, with stiffness `K_r` in **newtons [N]** —
force per unit strain, i.e. an equivalent `E*A` for the surrounding load path
(`EA/L` if a reference length is introduced).

`K_r` is deliberately **not** a translational stiffness in N/m, because this
model defines no reference length. The restraint is force-free when
`eps_common == reference_strain`; the canonical Milestone 2 case is
`reference_strain = 0`, meaning the surrounding structure tends to hold the
joint at its reference-temperature length. Nonzero reference strains are
supported and tested.

### Sign convention

Member forces stay tension-positive, exactly as in Milestone 1. The restraint is
a parallel load path carrying its own tension-positive force
`T = K_r (eps_common - eps_ref)`, and the force it exerts **on the joint** is the
reaction to that:

```
N_r = K_r * (reference_strain - eps_common)
```

`restraint_force` in the result is `N_r`, positive when the restraint acts on the
joint in the tensile (positive axial) sense. Joint equilibrium is then:

```
N_1 + N_2 - N_r = 0
```

The equilibrium residual is assembled from these actual signed forces — never
from magnitudes — and is reported, not discarded.

### Closed form

Substituting `sigma_i = E_i (eps - alpha_i dT)` into
`E_1 A_1 (eps - a_1 dT) + E_2 A_2 (eps - a_2 dT) + K_r (eps - eps_ref) = 0`:

```
                 E_1 A_1 alpha_1 dT + E_2 A_2 alpha_2 dT + K_r eps_ref
    eps_common = -----------------------------------------------------
                            E_1 A_1 + E_2 A_2 + K_r
```

The denominator is a **sum** of non-negative stiffnesses, so it is strictly
positive and the solution is unconditionally stable. A denominator of the form
`E_1 A_1 + E_2 A_2 - K_r` would signal an inconsistent force-direction
convention — that trap was audited before the model was coded.

### Limits

| Limit | Behaviour | Enforced by |
|---|---|---|
| `K_r = 0` | Exactly the Milestone 1 free-joint solution | direct regression against the Milestone 1 API |
| `K_r -> inf` | `eps_common -> eps_ref`, `sigma_i -> E_i (eps_ref - alpha_i dT)` | large-stiffness tests against the analytic rigid form |

For `eps_ref = 0` the rigid limit is the classic fully restrained result
`sigma_i = -E_i alpha_i dT`.

### The qualitative transition

This is the central Milestone 2 physics result. The **free** joint is always
self-equilibrating: one member in tension, the other in compression. A
**strongly restrained** joint drives *both* members to the **same sign** —
compression on heating, tension on cooling.

Consequently the low-CTE member must **cross through zero stress** at some
finite restraint level. Setting `eps_common = alpha_i dT` and solving:

```
K_r = dT * (alpha_i * K_joint - S) / (eps_ref - alpha_i * dT)
```

with `K_joint = E_1 A_1 + E_2 A_2` and `S = E_1 A_1 alpha_1 + E_2 A_2 alpha_2`.
For `eps_ref = 0` this collapses to the **dT-independent**
`K_r = S / alpha_i - K_joint`, exposed as
`zero_stress_restraint_stiffness(...)`, which returns `None` when no
admissible `K_r >= 0` exists (as for the high-CTE member, compressive both free
and fully restrained).

A direct consequence, verified rather than assumed: **matching CTEs no longer
means zero stress.** Milestone 1's `alpha_1 = alpha_2 -> sigma = 0` result held
*only* because the joint was externally free. Under restraint both members
develop the same-sign restraint stress.

### Restraint normalisation

```
eta_r = K_r / (E_1 A_1 + E_2 A_2)
```

`eta_r = 0` is free, `eta_r ~ 1` means the surrounding structure is about as
axially stiff as the joint itself, and `eta_r >> 1` is strongly restrained.
Available as `restraint_stiffness_ratio(...)`, as
`AxialRestraint.stiffness_ratio(...)`, and as the constructor
`AxialRestraint.from_stiffness_ratio(...)`.

### Inverse design — maximum tolerable restraint

`maximum_allowable_restraint_stiffness(...)` answers: *what is the largest
restraint stiffness this joint tolerates before its minimum preliminary yield
margin reaches zero?*

Feasibility is resolved **before** any search runs, using the analytic rigid
asymptote, and every outcome is an explicit status:

| Status | Meaning |
|---|---|
| `free_joint_already_fails` | The unrestrained joint already yields |
| `no_finite_limit_within_model` | Even a rigid restraint stays below yield — **no number is invented** |
| `finite_limit` | A boundary exists; returned in both `eta_r` and `K_r` |
| `no_boundary_within_search_bounds` | A boundary exists asymptotically but above the requested bound — reported, never silently expanded |

When a limit exists it is found by deterministic bounded bisection on `eta_r`
(better conditioned than `K_r`), with explicit upper bound, tolerance and
iteration cap.

---

## Representative Milestone 2 result

Same illustrative joint, same environment, design factor 1.25.
`K_joint = E_1 A_1 + E_2 A_2 = 1.80e7 N`; canonical restraint `eta_r = 1.0`
(`K_r = 1.80e7 N`), `eps_ref = 0`.

> **ILLUSTRATIVE** — the restraint is an *illustrative surrounding-structure
> axial restraint*, not a measured spacecraft structural stiffness.

| | Free (`eta_r = 0`) | Canonical (`eta_r = 1`) | Rigid (`eta_r -> inf`) |
|---|---|---|---|
| hot `sigma_1` | −62.03 MPa | −111.51 MPa | −161.00 MPa |
| hot `sigma_2` | **+62.03 MPa** | **−15.74 MPa** | −93.50 MPa |
| cold `sigma_1` | +86.84 MPa | +156.12 MPa | +225.40 MPa |
| cold `sigma_2` | **−86.84 MPa** | **+22.03 MPa** | +130.90 MPa |
| min margin | +1.487 PASS | +0.384 PASS | **−0.042 FAIL** |

The bolded `sigma_2` entries show the sign flip: the titanium-like member
crosses zero stress at **`eta_r = 0.663399`** (`K_r = 1.194e7 N`), independent of
`dT`. The aluminium-like member never changes sign.

Cold governs at every restraint level in this sweep — computed, not assumed.
Member 1 is loaded monotonically harder with restraint, but **member 2 is not
monotonic**: its stress magnitude falls to zero and then grows again with the
opposite sign.

**Inverse design result:** `finite_limit`, maximum tolerable
**`eta_r = 13.7405`** (`K_r = 2.473e8 N`), at which the minimum margin is zero.
This matches an independent analytical derivation of the boundary
(13.740543735) to within the search tolerance. Above that stiffness the joint
no longer meets the preliminary elastic yield criterion at its governing
extreme.

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

Milestone 2 additions (all backward compatible — nothing above changed):

| Object | Purpose |
|---|---|
| `AxialRestraint` | `stiffness` `K_r` [N] ≥ 0, `reference_strain`, optional `label`; `stiffness_ratio(...)` and `from_stiffness_ratio(...)` |
| `joint_axial_stiffness(...)` | `E_1 A_1 + E_2 A_2` [N] |
| `restraint_stiffness_ratio(...)` | `eta_r = K_r / (E_1 A_1 + E_2 A_2)` |
| `solve_restrained_joint(...)` | → `RestrainedThermalJointResult` |
| `RestrainedThermalJointResult` | Milestone 1 fields plus `restraint_stiffness`, `restraint_reference_strain`, `restraint_force`, `free_joint_common_strain`; residual is `N_1 + N_2 - N_r` |
| `rigid_restraint_stresses(...)` | analytic `K_r -> inf` stresses `E_i (eps_ref - alpha_i dT)` |
| `zero_stress_restraint_stiffness(...)` | restraint at which one member's stress crosses zero, or `None` |
| `assess_restrained_temperature_extremes(...)` | → `RestrainedTemperatureExtremeAssessment` |
| `maximum_allowable_restraint_stiffness(...)` | → `MaximumRestraintResult` with a `RestraintLimitStatus` |
| `rigid_restraint_minimum_margin(...)` | minimum margin in the fully restrained asymptote |

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

### Milestone 2 usage

```python
from thermal_joint import (
    AxialRestraint, assess_restrained_temperature_extremes,
    maximum_allowable_restraint_stiffness,
)

# Surrounding structure as stiff as the joint itself, holding it at T_ref length.
restraint = AxialRestraint.from_stiffness_ratio(1.0, *joint, reference_strain=0.0)

restrained = assess_restrained_temperature_extremes(
    *joint, environment, restraint, YieldBasis(1.25)
)
print(restrained.governing_extreme, restrained.minimum_yield_margin)

limit = maximum_allowable_restraint_stiffness(*joint, environment, YieldBasis(1.25))
print(limit.status.value, limit.maximum_stiffness_ratio)
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

Milestone 2 **reuses this machinery unchanged**. `assess_yield` accepts any
result exposing `member_1_stress` and `member_2_stress` (the `MemberStressState`
protocol), so the free and restrained assessments share one definition of a
margin, one zero-stress rule and one governing-member rule.

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

Milestone 2 adds, in the same spirit:

- **zero-restraint regression** — `K_r = 0` is compared directly against the
  Milestone 1 API (`solve_bimetallic_joint`, `assess_temperature_extremes`), not
  against a re-derived copy of the Milestone 1 equations
- **independent back substitution** — the computed strain is substituted into
  `E_1 A_1 (eps - a_1 dT) + E_2 A_2 (eps - a_2 dT) + K_r (eps - eps_ref) = 0`
  using raw scalars, bypassing the result's own force fields
- **full equilibrium under restraint** — `N_1 + N_2 - N_r` vanishes to machine
  precision at every restraint level and reference strain tested
- **rigid-restraint limit** — large-`K_r` numerics are checked against the
  analytic `E_i (eps_ref - alpha_i dT)` form, for zero and nonzero `eps_ref`
- **the qualitative transition** — the free joint is asserted to be
  self-equilibrating and the strongly restrained joint same-sign
- **analytic zero-stress crossing** — the closed-form crossing stiffness is
  checked against the numerically observed sign change, and its
  `dT`-independence for `eps_ref = 0` is verified
- **restraint load sharing** — member forces are shown to stop being equal and
  opposite, while `N_1 + N_2` tracks `N_r`
- **non-monotonicity is tested, not assumed** — member 1's stress magnitude is
  asserted monotone in restraint and member 2's is asserted *not* monotone
- **hand calculations** — a round-number restrained case
  (`K_r = 1.0e7 N`, `eta_r = 0.25`, `dT = +50 K` → `eps = 6.0e-4`,
  `sigma = -40/+20 MPa`, `N_r = -6000 N`) and a nonzero-`eps_ref` variant
- **inverse-design round trip** — the returned boundary matches an independent
  analytical derivation, has a ~zero margin, passes just below and fails just
  above; all four statuses are exercised, including a synthetic joint that
  proves the bisection itself works

The suite currently has **534 passing tests** — the 276 Milestone 1 tests,
unchanged and still green, plus 258 Milestone 2 tests.

---

## Limitations

Not modelled through Milestone 2, and **not** to be inferred from these results:

interface shear-lag · adhesive stresses · finite joint length or overlap ·
bolts and fasteners · contact and slip · nonlinear springs · plasticity ·
creep · fatigue · thermal gradients through thickness ·
temperature-dependent properties · nonlinear material behaviour · plate
bending and warpage · detailed radiator geometry · optimization · candidate
material trade studies · portfolio figures

In short: *this is a perfectly bonded, uniform-temperature axial compatibility
model with a single linear external restraint; interface shear, finite joint
length, thermal gradients, plasticity and fatigue are not included.*

The restraint is one lumped linear spring acting on the common joint strain. It
represents a surrounding load path in aggregate; it is not a model of any
particular bracket, fitting or fastener.

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
│   ├── illustrative.py      illustrative inputs (NOT design allowables)
│   ├── restraint.py         M2: AxialRestraint, restrained closed form, limits
│   ├── restrained_extremes.py  M2: hot/cold assessment under restraint
│   └── inverse.py           M2: maximum tolerable restraint stiffness
├── tests/
│   ├── conftest.py                      makes src/ importable without install
│   ├── reference_solution.py            independent closed-form check
│   ├── cases.py                         shared hand-calculation inputs
│   ├── restrained_cases.py              M2 hand-calculation inputs
│   ├── test_materials_and_members.py    M1 tests A–G
│   ├── test_free_thermal_strain.py      M1 tests H–L
│   ├── test_joint_mechanics.py          M1 tests M–U
│   ├── test_limits_and_scaling.py       M1 tests V–AC
│   ├── test_yield_margins.py            M1 tests AD–AK
│   ├── test_thermal_extremes.py         M1 tests AL–AR
│   ├── test_invalid_inputs.py           M1 invalid-input sweep
│   ├── test_physical_interpretation.py  M1 expected physics and regression lock
│   ├── test_restraint_object.py         M2 tests A–F
│   ├── test_restrained_closed_form.py   M2 tests G–M
│   ├── test_restrained_limits.py        M2 tests N–V
│   ├── test_restraint_sensitivity.py    M2 load sharing, sensitivity, sign flip
│   ├── test_restrained_extremes.py      M2 margins and governing selection
│   └── test_inverse_restraint.py        M2 inverse-design search
└── examples/
    ├── bimetallic_joint_sanity.py       M1 free-joint sanity study
    └── restraint_sensitivity.py         M2 restraint study
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

Run the Milestone 1 sanity study:

```bash
python examples/bimetallic_joint_sanity.py
```

Run the Milestone 2 restraint study:

```bash
python examples/restraint_sensitivity.py
```

The tests and the example both insert `src/` on `sys.path`, so they also run
directly from a clean checkout without installing the package.

---

## Licensing status

**No license has been chosen for this repository yet.** There is no `LICENSE`
file and no license metadata in `pyproject.toml`. Absent a license, default
copyright applies and no reuse rights are granted. A license will be added only
on an explicit decision by the repository owner.

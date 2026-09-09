# Bimetallic Joint Thermal Stress

Thermal-stress analysis of a spacecraft radiator joint made from two dissimilar
metallic members subjected to hot and cold temperature excursions.

**Status: Milestone 3 of an ongoing portfolio project — not complete.**

- **Milestone 1** — the verified mechanics foundation: closed-form axial thermal
  stress in a perfectly bonded two-member joint under a uniform temperature
  change, preliminary elastic yield margins, hot/cold extreme evaluation.
- **Milestone 2** — finite external axial restraint from the surrounding
  structure, the rigid-restraint limit, restraint sensitivity, and a
  preliminary inverse-design search for the maximum tolerable restraint
  stiffness. The Milestone 1 free-joint solution is preserved *exactly* as the
  zero-restraint limit.
- **Milestone 3** — a one-dimensional adhesive shear-lag screen for a finite
  bonded overlap: how the Milestone 1 mismatch force is actually transferred,
  what adhesive shear results, allowable temperature excursion, and
  minimum-overlap inverse design.
- **Milestone 4** — a bounded preliminary design trade over the verified
  mechanics: analytical asymptote scaling, six one-variable sensitivity sweeps,
  inverse sizing for bond width / adhesive thickness / adhesive modulus, a
  bounded width × thickness design map, and a deterministic selection policy.

> **Milestone 3 is a one-dimensional adhesive shear-lag screening model. It does
> not calculate peel stress, edge singularities, adhesive fracture, nonlinear
> response or finite-element joint stresses.**

> **Milestone 4 does not optimize a flight joint. It uses the verified
> one-dimensional shear-lag equations to identify feasible first-order
> combinations of bond width and adhesive compliance within explicitly bounded
> design spaces.**

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

> **Milestone 3.** If the free-joint CTE-mismatch force must be transferred
> through a finite bonded overlap, what axial-force distribution and adhesive
> shear stress develop along that overlap?

> **Milestone 4.** Given that overlap length alone cannot recover the failing
> canonical bondline, which physically meaningful first-order design levers
> *can*?

The primary deliverable of the overall project is a thermal stress calculation
plus a margin at the temperature extremes. Each milestone builds and
independently verifies the mechanics that answer its question.

**Global versus local.** Milestones 1 and 2 are *global*: they determine the
CTE-mismatch force the joint carries, with and without external restraint.
Milestone 3 is *local*: it takes that force as its demand and asks how it is
transferred through a finite bonded overlap. The shear-lag layer imports the
Milestone 1 member force directly from `solve_bimetallic_joint` and never
re-derives it, so the two levels cannot drift apart.

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

## Milestone 3 — adhesive shear-lag over a finite overlap

> **A screening model.** One-dimensional, linear-elastic load transfer and
> adhesive shear only. No peel stress, no edge singularities, no adhesive
> fracture, no nonlinear response, no finite-element joint stresses.

### Assumptions

- two axial adherends, constant overlap length and constant bond width
- thin adhesive layer of constant thickness, constant shear modulus, linear elastic
- adherends carry axial force only
- no bending, no eccentricity, no peel stress, no adhesive normal stress
- no free-edge singularity model, no yielding, no slip or debond
- uniform temperature, uniform properties, perfect bond

### Adhesive data and overlap geometry

`AdhesiveMaterial` carries `shear_modulus`, `shear_strength` and a **mandatory
non-empty `source_note`**, so a screening number can never be mistaken for a
qualified allowable. The shipped record is labelled
`ILLUSTRATIVE ADHESIVE-LIKE INPUT — NOT DESIGN ALLOWABLE`; no commercial
adhesive is named, because the properties are not sourced.

| | Illustrative adhesive | Illustrative overlap |
|---|---|---|
| `G_a` | 1.0 GPa | `L_b` = 40 mm |
| shear strength | 25 MPa | `b` = 20 mm |
| design factor | 1.25 → `tau_allow` = 20 MPa | `t_a` = 0.2 mm |

`bond_area = b * L_b` = 800 mm² is the bonded interface area, kept deliberately
distinct from the adherends' 100 mm² cross-sections. Both values are clean
mid-range picks made *before* the resulting stresses were computed — nothing is
tuned to manufacture a pass or a failure.

### Mismatch strain and the shear-lag parameter

```
d_eps = (alpha_1 - alpha_2) * dT
```

Member 1 minus member 2, used consistently. With `S_i = E_i A_i` and
`C = 1/S_1 + 1/S_2`, the governing equation derived below gives

```
beta^2 = (G_a b / t_a) * ( 1/(E_1 A_1) + 1/(E_2 A_2) )
```

Dimensions check: `G_a b / t_a` is [N/m²], `1/(EA)` is [1/N], so `beta²` is
[1/m²] and `beta` is [1/m]. **Both** adherend compliances enter — a joint
spreads load only as well as its more compliant member allows.

```
transfer length = 1 / beta          lambda = beta * L_b
```

`lambda` (not `beta L_b / 2`) is the single dimensionless overlap metric used
everywhere in this package: the number of transfer lengths the overlap spans.

### Sign convention and derivation

`u_i(x)` is axial displacement, `s(x) = u_1 - u_2` is the relative slip, and the
adhesive is linear: `tau(x) = (G_a / t_a) s(x)`. Adherend forces are
tension-positive, `N_i = E_i A_i (du_i/dx - alpha_i dT)`.

For an axial bar with distributed applied load `q` per unit length,
`dN/dx = -q`. When `s > 0` the adhesive drags adherend 1 backward and adherend 2
forward, so `q_1 = -b tau` and `q_2 = +b tau`:

```
dN_1/dx = +b tau        dN_2/dx = -b tau        d(N_1 + N_2)/dx = 0
```

This is the **sign-mirror of one common textbook ordering**; it is the version
consistent with `tau = +(G_a/t_a)(u_1 - u_2)`, and it is what makes the governing
equation stable rather than oscillatory. With the self-equilibrating free-joint
pair `N_2 = -N_1`:

```
ds/dx   = N_1/S_1 - N_2/S_2 + d_eps = C N_1 + d_eps
d2s/dx2 = C dN_1/dx = C b (G_a/t_a) s

=>  d2s/dx2 - beta^2 s = 0
```

### Boundary conditions

The domain is `x in [0, L_b]`:

- **`x = 0` — the loaded transfer plane.** The inboard edge of the overlap,
  where the free-joint mismatch force pair is fully developed:
  `N_1(0) = N_t` and `N_2(0) = -N_t`, with `N_t` imported from Milestone 1.
- **`x = L_b` — the free edge.** `N_1(L_b) = N_2(L_b) = 0`.

Prescribing the Milestone 1 force at the loaded plane is the **conservative
screening idealisation**: it assumes the mismatch force is fully developed and
requires the overlap to shear-transfer all of it.

### Closed-form distributions

`s = A cosh(beta x) + B sinh(beta x)` with `N_1 = (ds/dx - d_eps)/C`.
`N_1(0) = N_t = -d_eps/C` forces `B = 0`; `N_1(L_b) = 0` then fixes `A`:

```
N_1(x) =  N_t [ 1 - sinh(beta x) / sinh(beta L_b) ]
N_2(x) = -N_1(x)
tau(x) = -(N_t beta / b) cosh(beta x) / sinh(beta L_b)
```

Evaluated in closed form — the ODE is never solved numerically. The hyperbolic
ratios are implemented in an overflow-safe form, so a `lambda = 100` overlap
still evaluates cleanly.

### Peak versus average shear

`|tau|` is proportional to `cosh(beta x)`, whose derivative `beta sinh(beta x)`
is strictly positive for `x > 0`. So `|tau|` increases monotonically across the
overlap and the peak is at the **free edge** `x = L_b` — proven from the
derivative, not assumed.

```
tau_peak = |N_t| beta coth(beta L_b) / b
tau_avg  = |N_t| / (b L_b)                    (whole overlap — no factor of two)
tau_peak / tau_avg = lambda coth(lambda)
```

The concentration factor `lambda coth(lambda)` → 1 as `lambda → 0` (a very
short overlap shears almost uniformly) and → `lambda` as `lambda → ∞`.

The long-overlap asymptote is analytic:

```
tau_inf = |N_t| beta / b
```

and peak shear approaches it **from above**. That makes `tau_inf` a hard lower
bound: no overlap length, however long, can push the peak below it.

### Force-transfer verification

Five identities are enforced by tests: `N_1 + N_2 = 0` at every station;
`b ∫ tau dx` equals the change in adherend force; the loaded-plane force equals
the Milestone 1 demand exactly; the force-transfer residual vanishes to machine
precision; and zero mismatch gives zero force and shear everywhere. The integral
check is done with a Simpson rule **implemented in the test**, not by reusing
the model's own closed-form integral.

### Adhesive margin

```
tau_allow = shear_strength / design_factor
MS_adh    = tau_allow / abs(tau_peak) - 1
```

A **preliminary adhesive shear margin**, not a certification margin. Design
factor ≥ 1, boundary `MS = 0` passes, zero demand returns `math.inf`. The
governing extreme is computed from the two margins, with an exact tie resolving
to `"cold"` as in the earlier milestones.

### Global + local screening

`screen_joint(...)` runs the Milestone 1 yield screen and the Milestone 3
adhesive screen together. The two margins are **reported separately and never
combined numerically** — they measure different failure modes against different
allowables. Only the pass/fail booleans are combined, with a plain boolean AND.

---

## Representative Milestone 3 result

Illustrative joint and environment as before, 40 × 20 mm overlap, 0.2 mm
bondline, `G_a` = 1.0 GPa, adhesive design factor 1.25 (`tau_allow` = 20 MPa).

```
beta = 152.894 /m      transfer length = 6.54 mm      lambda = beta L_b = 6.116
```

| | Hot (`dT = +100 K`) | Cold (`dT = -140 K`) |
|---|---|---|
| transferred force `N_t` | 6202.78 N | 8683.89 N |
| `tau_peak` (free edge) | 47.42 MPa | **66.39 MPa** |
| `tau_avg` | 7.75 MPa | 10.85 MPa |
| peak / average | 6.116 | 6.116 |
| adhesive margin | −0.578 FAIL | **−0.699 FAIL** |

**Cold governs**, at a ratio of exactly 140/100 = 1.4 — peak shear is linear in
`dT`, so the larger excursion wins. Computed from the margins, not assumed.

**Overlap sensitivity** (cold, governing):

| `L_b` [mm] | `lambda` | `tau_peak` [MPa] | `tau_avg` [MPa] | peak/avg |
|---|---|---|---|---|
| 5 | 0.76 | 103.13 | 86.84 | 1.19 |
| 10 | 1.53 | 72.93 | 43.42 | 1.68 |
| 20 | 3.06 | 66.68 | 21.71 | 3.07 |
| 40 | 6.12 | 66.39 | 10.85 | 6.12 |
| 80 | 12.23 | 66.386 | 5.43 | 12.23 |
| 160 | 24.46 | 66.386 | 2.71 | 24.46 |

Average shear falls as `1/L_b` — a 32× reduction across the sweep — while the
peak moves by less than 1% beyond 20 mm, converging on the 66.386 MPa
asymptote. **This is the central Milestone 3 conclusion: overlap length buys
average shear, not peak shear.**

**Adhesive thickness** (0.05 → 1.0 mm) lowers `beta` from 305.8 to 68.4 /m,
lengthens the transfer zone from 3.27 to 14.62 mm and cuts peak shear from
132.8 to 29.9 MPa. **Adhesive modulus** (0.1 → 5 GPa) does the reverse, raising
peak shear from 21.9 to 148.4 MPa.

**Area ratio** `A_Al/A_Ti` (0.25 → 4) is recomputed end to end: the Milestone 1
demand rises from 3065 to 16032 N while `beta` falls from 257 to 113 /m, and
peak shear rises on balance from 39.4 to 90.2 MPa.

**Allowable excursion:** `|dT|_allow = 42.18 K` at this geometry — against study
excursions of −140 K and +100 K.

**Minimum overlap:** status `no_finite_length_within_model`. The allowable
(20 MPa) is below the long-overlap asymptote (66.386 MPa), so **no finite
overlap can pass** and no number is invented. Adding bonded length is the wrong
lever here.

**Combined screen:** member yield **+1.487 PASS**, adhesive shear **−0.699
FAIL**, overall feasibility `False`, failing screen `adhesive shear`. Reported
honestly: this illustrative joint is limited by its bondline, not by its metal.

### Engineering interpretation

- **M1 sets the demand.** The CTE-mismatch force is a global compatibility
  result and does not depend on how the joint is bonded.
- **M3 sets the transfer.** How that force spreads over the overlap depends on
  the adhesive shear stiffness and *both* adherend axial stiffnesses, through
  `beta`.
- **A thicker or more compliant adhesive spreads load over a longer transfer
  region**, lowering `beta` and the peak shear.
- **Longer overlap strongly reduces average shear but eventually gives
  diminishing reduction in peak end shear**, which is floored at `tau_inf`.
- **Cold governs the canonical magnitude** because 140 K > 100 K, for
  temperature-independent properties.
- **Lower predicted peak shear with a softer adhesive is not automatically a
  complete design improvement.** Peel stress, creep, durability, and joint
  deformation are all omitted from this screen, and every one of them tends to
  get worse as the bondline gets softer or thicker.

---

## Milestone 4 — bounded bondline design trade

> **Not an optimization.** This milestone evaluates explicitly bounded grids
> with the verified equations and applies one stated selection rule. It does not
> search a design space, and it makes no claim of global optimality.

### Why overlap length alone could not recover the design

Milestone 3's peak shear is `tau_peak = tau_inf coth(beta L_b)`. Once the
overlap spans more than a few transfer lengths, `coth → 1` and the peak sits on
the long-overlap asymptote `tau_inf`, which does **not** contain `L_b` at all.
At the canonical geometry `lambda = 6.12`, so `coth = 1.00001` — the overlap was
already asymptotic, and the remaining 66.386 MPa is untouchable by more bonded
length.

### Long-overlap asymptotic scaling

Substituting the Milestone 1 demand `N_t = d_eps / C` into `tau_inf = |N_t| beta / b`:

```
tau_inf = |d_eps| * sqrt( G_a / (b * t_a * C) )
```

| lever | exponent on `tau_inf` | cost to halve `tau_inf` |
|---|---|---|
| bond width `b` | **−1/2** (*not* `1/b`) | 4× wider |
| adhesive thickness `t_a` | −1/2 | 4× thicker |
| adhesive modulus `G_a` | +1/2 | 4× softer |
| `|dT|` and `|alpha_1 − alpha_2|` | +1 | 2× smaller |
| adherend compliance `C` | −1/2 | (stiffer adherends make it *worse*) |

Every exponent is verified to 12 decimal places in the tests. The width law is
the important trap: widening the bond also **raises** `beta`, so the naive
`tau_peak ∝ 1/b` is wrong — the true asymptotic law is `b^(-1/2)`.

Because every geometric lever sits under a square root, **halving the peak shear
costs a factor of four** in width, thickness or modulus. Only the thermal and
CTE levers are first order.

### Monotonicity and the uniform-shear floor

Writing `tau_peak = |N_t| f(beta) / b` with `f(u) = u coth(u L_b)`, and using
`sinh(2y)/2 > y`, `f` is **strictly increasing in `beta`**. Hence `tau_peak`
strictly falls with thickness and width, and strictly rises with modulus — all
three established analytically before any bisection runs.

As `beta → 0` (a very thick or very soft bondline) `f(beta) → 1/L_b`, so

```
tau_peak  ->  |N_t| / (b L_b)  =  tau_avg          the uniform-shear FLOOR
```

This floor is the key structural result of Milestone 4. **Thickness and modulus
are both bounded below by it**: if `tau_avg` already exceeds the allowable, no
bondline — however thick or compliant — can pass, and
`required_adhesive_thickness` / `maximum_allowable_adhesive_shear_modulus`
return `no_finite_thickness_within_model` / `no_finite_modulus_within_model`
rather than searching. **Bond width has no floor**, because `tau_avg` itself
falls as `1/b`. That makes width the only unconditionally effective single
lever in this model.

Canonically `tau_avg = 10.855 MPa` against a 20 MPa allowable, so thickness and
modulus do both have finite answers here.

### Sensitivity sweeps

All six sweeps are deterministic, build copies rather than mutating the shipped
records, and compute both extremes at every point (cold governs throughout,
computed from the margins).

| sweep | range | result |
|---|---|---|
| bond width | 10 → 100 mm | 93.92 → 29.69 MPa, **all FAIL** |
| adhesive thickness | 0.05 → 1.5 mm | 132.77 → 24.80 MPa, **all FAIL** |
| adhesive modulus | 5.0 → 0.05 GPa | 148.44 → 16.90 MPa, only 0.05 GPa passes |
| area ratio `A_1/A_2` | 0.25 → 4 | 39.44 → 90.22 MPa, all FAIL |
| excursion `|dT|` | 20 → 160 K | 9.48 → 75.87 MPa, passes to 40 K |
| CTE mismatch scale | 0 → 1.5 | 0 → 99.58 MPa, passes to ×0.25 |

Two findings worth stating plainly:

- **Equal areas are not optimal for adhesive shear.** Shrinking the
  aluminium-like member to `A_1/A_2 = 0.25` cuts peak shear from 66.39 to
  39.44 MPa, because it lowers the Milestone 1 mismatch force faster than it
  raises `beta`. It also cuts the yield margin from +1.487 to +0.762 — the two
  screens pull in opposite directions, which is exactly why they are reported
  separately.
- **The thermal sweep crosses PASS→FAIL between 40 K and 60 K**, bracketing the
  closed-form allowable of 42.18 K from Milestone 3.

### Inverse sizing

Each utility sizes one lever against the governing extreme by deterministic
bounded bisection, after settling feasibility analytically. Bounds, tolerance
and iteration cap are explicit and **never expanded silently**.

| utility | canonical result | status |
|---|---|---|
| `required_bond_width` | **220.35 mm** (11.02× the 20 mm baseline) | `finite_required_width` |
| `required_adhesive_thickness` | **2.499 mm** (12.49× the 0.2 mm baseline) | `finite_required_thickness` |
| `maximum_allowable_adhesive_shear_modulus` | **0.0800 GPa** (12.49× softer) | `finite_maximum_modulus` |

All three demand essentially the same factor on the group `G_a / (b t_a)`,
because the asymptote depends on nothing else. Each is individually extreme — a
220 mm bond, a 2.5 mm bondline, or an 80 MPa-class adhesive. **These are
model-based screening values, not allowables.**

Status vocabularies are explicit: `lower_bound_already_passes`,
`finite_required_*`, `no_boundary_within_search_bounds`, plus
`no_finite_thickness_within_model` / `no_finite_modulus_within_model` where the
uniform-shear floor makes the lever useless, and
`upper_bound_already_passes` / `lower_bound_already_fails` for the modulus
search, whose boundary is a **maximum** rather than a minimum.

### Design-factor sensitivity

Required width scales as `1/tau_allow²`, so the design basis drives the sizing
hard:

| design factor | `tau_allow` | cold margin | required width |
|---|---|---|---|
| 1.0 | 25.00 MPa | −0.623 | 141.03 mm |
| 1.25 | 20.00 MPa | −0.699 | 220.35 mm |
| 1.5 | 16.67 MPa | −0.749 | 317.31 mm |
| 2.0 | 12.50 MPa | −0.812 | 564.11 mm |

Doubling the factor quadruples the bond width. The factor lives in
`AdhesiveShearBasis`, never inside the adhesive record.

### Bounded width × thickness design map

Widths 10–100 mm × thicknesses 0.1–1.0 mm, 42 points, canonical adhesive
(`G_a = 1.0 GPa`, 25 MPa strength, DF 1.25). Cold peak shear [MPa], `*` = passes
both screens:

| `b` \ `t_a` | 0.10 | 0.20 | 0.30 | 0.50 | 0.75 | 1.00 |
|---|---|---|---|---|---|---|
| **10 mm** | 132.77 | 93.92 | 76.79 | 59.88 | 49.61 | 43.78 |
| **20 mm** | 93.88 | 66.39 | 54.21 | 42.02 | 34.41 | 29.94 |
| **30 mm** | 76.66 | 54.20 | 44.26 | 34.29 | 28.02 | 24.30 |
| **40 mm** | 66.39 | 46.94 | 38.33 | 29.69 | 24.25 | 21.01 |
| **50 mm** | 59.38 | 41.99 | 34.28 | 26.56 | 21.68 | **18.78\*** |
| **75 mm** | 48.48 | 34.28 | 27.99 | 21.68 | **17.70\*** | **15.33\*** |
| **100 mm** | 41.99 | 29.69 | 24.24 | **18.78\*** | **15.33\*** | **13.28\*** |

**6 of 42 points are feasible** — co-design works where neither single lever
did. A second illustrative map at `G_a = 0.1 GPa` opens **32 of 42**; the
canonical adhesive is not replaced by that variant, it is shown alongside purely
to demonstrate the effect.

### Preliminary design policy and selected point

Policy `minimum_bond_area`, applied in order:

1. must pass member yield **and** adhesive shear at both extremes
2. smallest bond area (= width × overlap; the overlap is fixed, so this is the
   narrowest bond)
3. tie-break: smaller adhesive thickness
4. tie-break: larger adhesive margin
5. tie-break: earlier position in the grid

| | baseline | selected |
|---|---|---|
| bond width | 20.00 mm | **50.00 mm** |
| adhesive thickness | 0.20 mm | **1.00 mm** |
| bond area | 800 mm² | 2000 mm² |
| hot peak shear | 47.42 MPa | 13.42 MPa |
| cold peak shear | 66.39 MPa | 18.78 MPa |
| hot adhesive margin | −0.578 | **+0.491** |
| cold adhesive margin | −0.699 | **+0.065** |
| member yield margin | +1.487 | +1.487 |
| overall feasible | `False` | **`True`** |

> **The selected point is a preliminary feasible point within the bounded study
> grid, not an optimized flight-joint design.**

A thinner adhesive is not automatically better for manufacturing or durability,
and a wider bond is not automatically worse. Bond area and adhesive volume are
reported as **geometric diagnostics only** — no adhesive density is supplied, so
no adhesive mass is computed or implied, and bond width is not a mass metric.

### Engineering interpretation

- **Why overlap failed.** Past a few transfer lengths peak shear is set by the
  long-overlap asymptote, which does not contain the overlap length.
- **Width.** Widening the bond lowers force-transfer intensity, but it also
  raises `beta`; the net asymptotic law is `b^(-1/2)`, not `1/b`.
- **Adhesive thickness.** A thicker bondline is more compliant, lowers `beta`
  and spreads transfer over a longer distance — down to the uniform-shear floor.
- **Adhesive modulus.** A softer adhesive spreads transfer the same way, but may
  introduce deformation, creep and durability problems this model cannot see.
- **Co-design.** A width or compliance change that is insufficient alone becomes
  effective in combination: neither 50 mm nor 1.0 mm passes by itself, but
  together they do.
- **Metal versus adhesive.** The Milestone 1 metals pass comfortably while the
  Milestone 3 baseline bondline fails. These are distinct margins against
  distinct allowables and must stay separately visible; only the booleans are
  ANDed.
- **Model caution.** The model predicts its largest shear at the free edge —
  exactly where the omitted peel stress and edge singularity matter most. A
  passing screening margin here is necessary, not sufficient.

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

Milestone 3 additions (again purely additive):

| Object | Purpose |
|---|---|
| `AdhesiveMaterial` | `shear_modulus`, `shear_strength`, mandatory `source_note`, optional `notes` |
| `BondedOverlapGeometry` | `overlap_length`, `bond_width`, `adhesive_thickness`; `bond_area`, plus `with_*` copy helpers used by the sweeps |
| `AdhesiveShearBasis` | `design_factor` ≥ 1; `allowable_shear_stress(...)`, `margin_of_safety(...)` |
| `thermal_mismatch_strain(...)` | `(alpha_1 - alpha_2) dT` |
| `adherend_compliance_sum(...)` | `C = 1/(E_1 A_1) + 1/(E_2 A_2)` |
| `shear_lag_parameter(...)` / `transfer_length(...)` / `dimensionless_overlap(...)` | `beta`, `1/beta`, `lambda = beta L_b` |
| `solve_shear_lag(...)` | → `ShearLagDemand` |
| `ShearLagDemand` | scalars plus `shear_stress(x)`, `member_1_force(x)`, `member_2_force(x)`, all range-checked to `[0, L_b]` |
| `long_overlap_peak_shear_stress(...)` | analytic `tau_inf = |N_t| beta / b` |
| `assess_shear_lag_extremes(...)` | → `ShearLagExtremeAssessment` |
| `allowable_temperature_change_for_adhesive_shear(...)` | closed-form `|dT|_allow` [K] |
| `required_overlap_length(...)` | → `RequiredOverlapResult` with an `OverlapLimitStatus` |
| `screen_joint(...)` | → `PreliminaryScreeningResult` (yield AND adhesive) |

Milestone 4 additions (purely additive; no Milestone 1–3 module was touched):

| Object | Purpose |
|---|---|
| `long_overlap_peak_shear(...)` | analytic `tau_inf` for a whole configuration |
| `uniform_shear_floor(...)` | `|N_t| / (b L_b)` — the `beta → 0` floor |
| `JointDesignCandidate` | one evaluated configuration; yield and adhesive margins side by side, plus `bond_area` and `adhesive_volume` |
| `evaluate_joint_design(...)` | → `JointDesignCandidate` at both extremes |
| `SensitivityPoint` | one point of a sweep whose variable is not a bondline dimension |
| `bond_width_sensitivity(...)`, `adhesive_thickness_sensitivity(...)`, `adhesive_modulus_sensitivity(...)` | → tuples of `JointDesignCandidate` |
| `area_ratio_sensitivity(...)`, `thermal_excursion_sensitivity(...)`, `cte_mismatch_sensitivity(...)` | → tuples of `SensitivityPoint` |
| `required_bond_width(...)` | → `RequiredBondWidthResult` with a `BondWidthStatus` |
| `required_adhesive_thickness(...)` | → `RequiredAdhesiveThicknessResult` with an `AdhesiveThicknessStatus` |
| `maximum_allowable_adhesive_shear_modulus(...)` | → `MaximumAdhesiveModulusResult` with an `AdhesiveModulusStatus` |
| `width_thickness_design_map(...)` | bounded grid → tuple of candidates, row-major |
| `select_preliminary_joint_design(...)` | one candidate or `None`, under `DesignSelectionPolicy` |

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

### Milestone 3 usage

```python
from thermal_joint import (
    AdhesiveMaterial, AdhesiveShearBasis, BondedOverlapGeometry,
    assess_shear_lag_extremes, required_overlap_length,
)

adhesive = AdhesiveMaterial(
    "Adhesive-like", shear_modulus=1.0e9, shear_strength=25e6,
    source_note="ILLUSTRATIVE ADHESIVE-LIKE INPUT - NOT DESIGN ALLOWABLE",
)
overlap = BondedOverlapGeometry(
    overlap_length=0.040, bond_width=0.020, adhesive_thickness=0.0002,
)

shear = assess_shear_lag_extremes(
    *joint, environment, overlap, adhesive, AdhesiveShearBasis(1.25)
)
print(shear.governing_extreme, shear.governing_peak_shear_stress, shear.minimum_margin)

sizing = required_overlap_length(
    *joint, environment, overlap, adhesive, AdhesiveShearBasis(1.25)
)
print(sizing.status.value)
```

### Milestone 4 usage

```python
from thermal_joint import (
    required_bond_width, select_preliminary_joint_design, width_thickness_design_map,
)

sizing = required_bond_width(
    *joint, environment, overlap, adhesive, AdhesiveShearBasis(1.25)
)
print(sizing.status.value, sizing.required_bond_width)

grid = width_thickness_design_map(
    *joint, environment, overlap, adhesive,
    bond_widths=[0.010, 0.020, 0.050, 0.100],
    adhesive_thicknesses=[0.0001, 0.0005, 0.0010],
    yield_basis=YieldBasis(1.25), adhesive_basis=AdhesiveShearBasis(1.25),
)
selected = select_preliminary_joint_design(grid)   # a candidate, or None
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

Milestone 3 adds, again in the same spirit:

- **independent numerical integration** — a Simpson rule implemented inside the
  test integrates the production `tau(x)` and reproduces the adherend force
  change, end to end and on sub-intervals; the model's own closed-form integral
  is never used to check itself
- **hand calculations with a round `beta`** — inputs chosen so `beta = 100 /m`
  and `lambda = 2` exactly, giving `tau_avg = 12.5 MPa`, `tau_peak = 25 MPa ×
  coth(2)`, `N_1(10 mm) = -3379.9 N`
- **the peak location is proven, not assumed** — the derivative argument is
  backed by scanning 5001 stations across the overlap
- **self-equilibrium everywhere** — `N_1 + N_2 = 0` checked at 51 stations
- **demand reuse** — the transferred force is asserted *bit-identical* to
  `solve_bimetallic_joint(...).member_1_force`, and separately cross-checked
  against the independent identity `N_t = -d_eps / C`
- **short and long overlap limits** — `lambda coth(lambda) → 1` at
  `lambda < 0.06`, the distribution itself flattening to within 0.1%; and peak
  shear converging monotonically on the analytic asymptote from above with
  strictly diminishing returns per doubling
- **scaling laws** — peak shear exactly linear in `dT`, cold/hot ratio exactly
  1.4 for the canonical environment
- **inverse design** — the allowable `dT` round-trips to `MS = 0` exactly, the
  finite required length matches an independent `atanh` solution, and all four
  overlap statuses are exercised, three of them with synthetic fixtures
- **numerical robustness** — a `lambda = 100` overlap evaluates without
  overflow, thanks to overflow-safe hyperbolic ratios
- **regression on the earlier milestones** — the Milestone 1 stresses and
  margins and the Milestone 2 restrained results, crossing stiffness and
  restraint boundary are all re-asserted from Milestone 3 test files
- **sweeps are non-mutating** — the shipped members, geometry and adhesive are
  asserted unchanged after every sensitivity sweep

Milestone 4 adds:

- **power-law exponents measured, not asserted** — every asymptotic exponent is
  recovered numerically from a factor change and checked to 1e-12
- **the width trap is tested explicitly** — `b^(-1/2)` is confirmed *and* the
  naive `1/b` scaling is asserted to be wrong
- **monotonicity before bisection** — each sweep is asserted monotone in the
  direction the derivation predicts, so every inverse search is well posed
- **the uniform-shear floor** — verified as the `beta → 0` limit with a
  1 kPa-modulus adhesive, and shown to trigger
  `no_finite_thickness_within_model` on a 5 mm-wide bond
- **inverse round trips** — each returned boundary gives `MS ≈ 0`, fails just
  inside and passes just outside; the required width is separately checked
  against a closed-form asymptotic estimate
- **all statuses exercised**, including the two bracket-too-narrow cases
- **selection determinism and order-independence** — the map is shuffled and
  reversed and yields the same pick; the thickness tie-break is tested at equal
  bond area against a candidate with a *larger* margin, proving the stated
  policy order is what actually runs
- **`None` on infeasibility** — a grid with no feasible point returns nothing
  rather than a best-effort pick
- **no mass invented** — the candidate is asserted to expose no mass or density
  attribute, and no combined margin attribute
- **prior milestones re-asserted** from Milestone 4 test files: M1 stresses and
  margins, M2 restrained results, crossing and restraint boundary, M3 baseline
  shear-lag and `no_finite_length_within_model`

The suite currently has **824 passing tests** — the 276 Milestone 1, 258
Milestone 2 and 189 Milestone 3 tests, unchanged and still green, plus 101
Milestone 4 tests.

---

## Limitations

### Milestone 3 (adhesive shear-lag)

1-D shear-lag only · no peel stress · no bending or eccentricity · no free-edge
singularity · no adhesive normal stress · no adhesive plasticity · no
cohesive-zone or fracture model · no debond growth · no creep · no fatigue · no
thermal gradient · no temperature-dependent properties · no fillets or spew
geometry · no surface-preparation effects · no manufacturing defects · no
environmental degradation · illustrative adhesive properties · **no
certification claim**

The prescribed loaded-plane force is a conservative idealisation: it assumes the
mismatch force is fully developed and asks the overlap to carry all of it. Peak
shear at a free edge is exactly where a real bondline is least well represented
by a 1-D model — peel stress and the edge singularity, both excluded here, act
in the same place.

### Milestone 4 (design trade)

Bounded design trade only · no adhesive mass model unless a density with
provenance is supplied · no manufacturability constraint · no minimum or maximum
practical bondline thickness · no adhesive cure or process limits · no
compliance-durability model · no creep · no fatigue · no fracture · no peel · no
edge singularity · no thermal-property variation · no environmental degradation
· illustrative adhesive properties · **inverse results are model-based screening
values, not allowables** · no claim of global optimality

The selected point is a preliminary feasible point within the bounded study
grid, not an optimized flight-joint design. A 2.5 mm bondline or a 220 mm bond
may be perfectly reasonable arithmetic and quite unreasonable hardware; this
model has no way to tell.

### Global model

Not modelled through Milestone 2, and **not** to be inferred from these results:

bolts and fasteners · contact and slip · nonlinear springs · plasticity ·
creep · fatigue · thermal gradients through thickness ·
temperature-dependent properties · nonlinear material behaviour · plate
bending and warpage · detailed radiator geometry · optimization · candidate
material trade studies · portfolio figures

In short: *the global model is a perfectly bonded, uniform-temperature axial
compatibility model with a single linear external restraint, and the local model
is a one-dimensional linear-elastic adhesive shear-lag screen; thermal
gradients, plasticity, fatigue, peel and fracture are not included.*

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
│   ├── inverse.py           M2: maximum tolerable restraint stiffness
│   ├── adhesive.py          M3: adhesive, overlap geometry, shear basis
│   ├── shear_lag.py         M3: beta, closed-form distributions, peak/average
│   ├── shear_lag_extremes.py   M3: hot/cold adhesive shear assessment
│   ├── shear_lag_inverse.py    M3: allowable dT and minimum overlap
│   ├── screening.py         M3: combined yield AND adhesive feasibility
│   ├── design_trade.py      M4: scaling, sweeps, design map, selection policy
│   └── design_inverse.py    M4: width / thickness / modulus inverse sizing
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
│   ├── test_inverse_restraint.py        M2 inverse-design search
│   ├── shear_lag_cases.py               M3 hand-calculation inputs
│   ├── test_adhesive_and_overlap.py     M3 tests A–I
│   ├── test_shear_lag_parameters.py     M3 tests J–Q
│   ├── test_shear_lag_distribution.py   M3 tests R–Y
│   ├── test_shear_lag_limits.py         M3 tests Z–AH
│   ├── test_adhesive_margins.py         M3 tests AI–AO
│   ├── test_shear_lag_inverse.py        M3 tests AP–AW
│   ├── test_milestone_integration.py    M3 tests AX–BC (M1/M2 regression)
│   ├── test_shear_lag_sensitivity.py    M3 tests BD–BI
│   ├── design_cases.py                  M4 canonical fixtures
│   ├── test_design_scaling.py           M4 tests A–F
│   ├── test_design_sensitivity.py       M4 tests G–AE
│   ├── test_design_inverse.py           M4 tests AF–AW
│   └── test_design_map_and_selection.py M4 tests AX–BM
└── examples/
    ├── bimetallic_joint_sanity.py       M1 free-joint sanity study
    ├── restraint_sensitivity.py         M2 restraint study
    ├── adhesive_shear_lag.py            M3 adhesive shear-lag study
    └── joint_design_trade.py            M4 bondline design trade
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

Run the Milestone 3 adhesive shear-lag study:

```bash
python examples/adhesive_shear_lag.py
```

Run the Milestone 4 bondline design trade:

```bash
python examples/joint_design_trade.py
```

The tests and the example both insert `src/` on `sys.path`, so they also run
directly from a clean checkout without installing the package.

---

## Licensing status

**No license has been chosen for this repository yet.** There is no `LICENSE`
file and no license metadata in `pyproject.toml`. Absent a license, default
copyright applies and no reuse rights are granted. A license will be added only
on an explicit decision by the repository owner.

<title>Radiative Heat Transfer in nekRS</title>

# Radiative Heat Transfer in nekRS — the `[RADIATION]` `.par` Card

A gray-diffuse surface-to-surface radiation model, driven by
Monte-Carlo-computed view factors between user-specified boundary-ID
groups. Two things happen when a case turns it on:

1. **Once, at setup** (`Radiation::setup()`, called from `UDF_Setup()`):
   Monte Carlo ray sampling computes the view-factor matrix between
   every pair of radiating surface patches. Pure geometry — independent
   of the flow field, cacheable across runs.
2. **Periodically, during time-stepping** (`Radiation::step(time,
   tstep)`, called from `UDF_ExecuteStep()`): the current wall
   temperatures are gathered, the gray-diffuse radiosity system is
   solved (warm-started from the previous solve), and the resulting
   net radiative flux per patch is written into
   `platform->app->bc->o_usrwrk` for a case's `udfNeumann` to read.

Model: `Jᵢ = εᵢσTᵢ⁴ + (1-εᵢ)Σⱼ Fᵢⱼ Jⱼ`, net flux `qᵢ = εᵢ(σTᵢ⁴ - Σⱼ
Fᵢⱼ Jⱼ)` — the standard diffuse-gray radiosity network, solved via
Gauss-Seidel iteration.

---

## `.par` keys

```ini
[RADIATION]
radiatingBoundaryIDs = 1,2,3
obstructionBoundaryIDs = 1,2,3
nSamples = 16384
seed = 1234
writeMatrix = true
outputFile = myCaseName_radiation
cache = true
emissivity = 0.9,0.9,0.8
stefanBoltzmann = 5.670374419e-8
updateFrequency = 20
radiosityTolerance = 1e-6
radiosityMaxIters = 100
```

| key | default | meaning |
|---|---|---|
| `radiatingBoundaryIDs` | *(required)* | Comma list of raw boundary IDs that participate in the radiative enclosure. Every listed ID must actually appear in the mesh's boundary faces. |
| `obstructionBoundaryIDs` | *(none)* | Comma list of raw boundary IDs that can **block** line-of-sight between radiating patches. Only listed here does a surface occlude — a surface in `radiatingBoundaryIDs` does **not** automatically also occlude (see "Concave enclosures" below). Setting this builds a BVH over the listed surfaces' triangulated faces; leaving it empty skips occlusion testing entirely (cheaper, but wrong for concave enclosures). |
| `nSamples` | `4096` | Requested Monte Carlo samples per patch pair. Internally rounded up to `nStrata²` where `nStrata = ceil(sqrt(nSamples))`, for square stratified sampling. |
| `seed` | `0` | RNG seed for the (stateless, counter-based) sampler. |
| `writeMatrix` | `true` | Write the view-factor matrix to disk (`<outputFile>_viewfactors_groups.csv`, a group-aggregated table, plus `<outputFile>_viewfactors_patches.bin`, the full per-patch matrix). |
| `outputFile` | `<casename>_radiation` | Prefix for all output/cache files. |
| `cache` | `true` | Reuse a previous run's view-factor matrix if a fingerprint of the mesh/config (`<outputFile>.hash`) still matches — skips the expensive Monte Carlo recompute entirely. |
| `emissivity` | `1.0` for every group | Comma list, **one value per unique boundary-ID group among the radiating patches, in ascending boundary-ID order** — not per-patch, and not in whatever order you wrote `radiatingBoundaryIDs`. A length mismatch against the number of distinct groups is a hard error at setup. |
| `stefanBoltzmann` | `5.670374419e-8` | SI value, W/(m²·K⁴). Override only if the case uses a non-SI unit convention throughout. |
| `updateFrequency` | `10` | Radiosity solve and flux update run only every `updateFrequency` steps (`tstep % updateFrequency == 0`); the flux otherwise holds its last computed value. Radiative equilibration is slow relative to CFD timescales, so this is the main cost lever — Yuan et al. report overhead dropping from 620% at `Nr=1` to 12% at `Nr=50` for their pebble-bed case. |
| `radiosityTolerance` | `1e-6` | Gauss-Seidel convergence tolerance for the per-update radiosity solve. |
| `radiosityMaxIters` | `100` | Cap on Gauss-Seidel iterations per update. |

**Patch count ceiling**: the module targets a modest "surface
enclosure" subset of the mesh, not the full CFD boundary. Both the
view-factor matrix storage and the sampling kernel scale as `O(P²)`
in the number of radiating patches `P` (patches = boundary *element
faces*, not GLL nodes); `P` is hard-capped at 8000 (~256MB for the
dense matrix alone). Narrow `radiatingBoundaryIDs` if you hit this.

---

## Wiring it into a case

### `.udf`

```cpp
#include "Radiation.hpp"

void UDF_Setup()
{
  Radiation::setup();
  // ... rest of setup
}

void UDF_ExecuteStep(double time, int tstep)
{
  Radiation::step(time, tstep);
}
```

That's the entire diagnostic-only integration — sufficient if every
radiating wall's temperature is Dirichlet-fixed. `Radiation::step()`
will still solve the radiosity system and print a per-group flux
report every `updateFrequency` steps, useful as a sanity check even
when nothing reads it back.

### Coupling the flux back into the temperature solve

To actually let radiation affect the fluid's energy balance, a wall
needs `udfNeumann` instead of `udfDirichlet`, reading the computed flux
from `bc->usrwrk`:

```cpp
#ifdef __okl__
void udfNeumann(bcData *bc)
{
  if (isField("scalar temperature")) {
    dfloat appliedFlux = 0.0;         // W/m^2, external source at this wall
    if (bc->id == 1) {
      appliedFlux = 0.5;              // e.g. a heater
    }
    // bc->usrwrk holds the net radiative flux LEAVING this wall
    // (positive = net radiator, i.e. hotter than what it sees).
    // The flux actually conducted into the fluid is what's left
    // after radiative losses.
    bc->fluxScalar = appliedFlux - bc->usrwrk[bc->idxVol];
  }
}
#endif
```

**Sign convention**, verified against the actual solver assembly (not
assumed): positive `bc->fluxScalar` heats the fluid. `bc->usrwrk[...]`
is `qᵢ` from the radiosity solve, positive when the patch is a net
radiator. So `appliedFlux - bc->usrwrk[...]` is correct for a wall with
an external heat source: whatever isn't lost to radiation goes into
the fluid. For a purely "adiabatic but still radiatively active" wall
(the paper's floating reradiating-surface treatment), use
`appliedFlux = 0` — the wall's temperature then floats to whatever
value the radiative balance implies, with no other coupling needed.

`boundaryTypeMap` in `[SCALAR TEMPERATURE]` must list `udfNeumann` for
these boundary IDs (and `udfDirichlet` for any genuinely
fixed-temperature walls, `zeroNeumann` for non-radiating symmetry
planes, etc. — see `examples.md` for a full worked case,
`examples/tallCavityNu`).

### `platform->app->bc->o_usrwrk` ownership

The module currently takes **exclusive** ownership of
`platform->app->bc->o_usrwrk` — `setup()` hard-errors if it's already
sized for something else. A case that also needs `o_usrwrk` for an
unrelated boundary condition (e.g. a `gabls1`-style velocity wall
model) isn't supported without extending the module to share offsets.

---

## Boundary-ID bookkeeping

Two issues come up in essentially every real case; see `examples.md`
for the full mechanics and worked fixes.

**Getting boundary IDs into the mesh at all.** nekRS's automatic
CBC-label→numeric-ID import is disabled the moment *any* field gets an
explicit `boundaryTypeMap` — which a radiating case almost always
needs (velocity and temperature need different BC types on the same
raw wall IDs). For `genbox` meshes this means a `usrdat2()` in the
case's `.usr` file assigning `boundaryID` by hand from the text CBC
label (genbox never writes a usable numeric ID for ordinary faces).
`gmsh2nek`-converted meshes don't have this problem — `gmsh2nek` writes
the Physical Surface tag directly into the boundary-ID slot.

**Non-physical boundaries still need to be in the radiative
accounting.** The Monte Carlo ray tracer has no notion of "this
boundary doesn't count" — a symmetry plane, inlet, or outlet excluded
from `radiatingBoundaryIDs` just silently loses whatever energy rays
happen to head toward it, breaking energy conservation for every
*other* surface in the enclosure (rows of the view-factor matrix no
longer sum to ~1). The fix, matching Yuan et al.'s own treatment of
their single-pebble case's inlet/outlet: include the boundary in
`radiatingBoundaryIDs` with `emissivity = 0` — a diffuse mirror,
present for the geometry but contributing zero net flux (`qᵢ =
εᵢ(...) = 0` identically).

**Concave enclosures need `obstructionBoundaryIDs`.** A surface that
can see part of itself around an obstruction (e.g. the outer sphere in
`examples/concentricSpheres`, occluded by the inner one) needs to be
listed in `obstructionBoundaryIDs` — radiating surfaces do not
automatically double as occluders of each other. Symptom if this is
missed: view-factor row sums exceeding 1 (energy leaving a surface
appearing to sum to more than what's physically possible), since
self-view rays that should be blocked aren't.

---

## Near-field sampling (patches close relative to their size)

For two radiating patches whose separation is small relative to their
own size — nearly-touching mesh elements, tightly packed geometry —
naive flat Monte Carlo sampling under-resolves the region that
dominates the view-factor integral (the `cosθᵢcosθⱼ/(πr²)` kernel gets
large there, but occupies a small fraction of each patch's parameter
domain). The kernel detects such pairs (centroid separation less than
roughly two patch diagonals) and biases sampling toward the facing
edge via a Duffy-transform-style reparametrization, with an explicit
Jacobian correction to keep the estimator unbiased. This is internal
to the sampling kernel — no `.par` option controls it, and it should
not require attention in normal use. (Implementation note for anyone
touching `Radiation.okl`: the existing stratified-sampling
decorrelation hash, `wang_hash(k^const) % nStrata`, was validated only
for the original uniform-weight kernel; a naive port of the near-field
bias reused it and introduced a real, confirmed bias from that hash's
low-bit correlation with the new weighting. Close pairs now bypass
stratification entirely for both points instead — see git history on
`Radiation.okl` if extending this further.)

---

## Caching

With `cache = true` (the default), `setup()` computes a fingerprint of
the mesh and `[RADIATION]` configuration and compares it against
`<outputFile>.hash`. On a match, it reloads the view-factor matrix
from `<outputFile>_viewfactors_patches.bin` instead of recomputing —
useful across repeated runs of the same case (e.g. sweeping
emissivity, where the geometry doesn't change). Delete the `.hash`
and `.bin` files (or set `cache = false`) to force a fresh Monte Carlo
recompute after a genuine geometry or sampling-parameter change.

---

## Minimal example (diagnostic only)

```ini
[RADIATION]
radiatingBoundaryIDs = 1,2,3
nSamples = 8192
seed = 1234
emissivity = 0.9,0.9,0.8
updateFrequency = 20
```

```cpp
#include "Radiation.hpp"

void UDF_Setup() { Radiation::setup(); }
void UDF_ExecuteStep(double time, int tstep) { Radiation::step(time, tstep); }
```

With all three boundary IDs Dirichlet-fixed elsewhere in `.par`, this
prints a per-group net-radiative-flux report every 20 steps and writes
`<casename>_radiation_viewfactors_groups.csv` — enough to validate a
new geometry's view factors (row sums ≈1, reciprocity `AᵢFᵢⱼ = AⱼFⱼᵢ`)
before building a coupled case on top of it. See `examples.md` for
what a fully coupled case (`examples/tallCavityNu`) looks like.
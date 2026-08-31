# singlePebble — Case Files, Pipeline, and Assumptions

A coupled radiation-convection case: forced-convection air flow through a
square duct past a single heated spherical pebble, reproducing Section 4.4
("Flow around a single pebble") of Yuan, Dai, Shaver, Gottems, and Merzari,
*"Implementation of the Surface-to-Surface Thermal Radiation Solver in
Spectral Element CFD code NekRS"*.

Unlike `radiationPlates`/`concentricSpheres` (pure radiation diagnostics, no
flow of interest) or `tallCavityNu` (buoyancy-driven natural convection),
this is the module's first forced-convection case: a real inlet/outlet duct
flow with a curved obstacle, combining `concentricSpheres`' gmsh "cubed
sphere" meshing technique (genbox can't do curved geometry) with
`tallCavityNu`'s radiative-convective boundary coupling (`udfNeumann` +
`bc->usrwrk`).

---

## Geometry (Sec 4.4, Table IV)

- Square duct, 8x8 cm² cross-section, 16 cm long (`WX=WY=0.04 m`, `WZ=0.08
  m` half-widths in `make_mesh.py`, i.e. duct spans `x,y in [-0.04,0.04]`,
  `z in [-0.08,0.08]`).
- Spherical pebble, diameter 4 cm (`RP=0.02 m` radius), centered at the
  origin.
- Flow direction: `+z` (inlet at `z=-0.08`, outlet at `z=+0.08`).

### `make_mesh.py` — mesh generation

genbox can't produce curved geometry, so the mesh is built with `gmsh`
(built-in kernel) and converted with `gmsh2nek`, generalizing
`examples/concentricSpheres/make_mesh.py`'s "cubed sphere" technique from
one shell to **two nested shells**:

- **Shell A** (pebble → collar): the pebble surface (sphere, radius `RP`) to
  a small cubed "collar" (flat cube, half-width `RC=0.028 m`). Identical
  technique to `concentricSpheres`: 6 curved panels, each bounded by 2
  great-circle caps (sphere side, exact curvature via `Circle` arcs through
  the origin) and a flat cap (collar side), joined by straight radial
  connectors.
- **Shell B** (collar → duct): the collar cube to the duct's own outer
  walls — a box, not a cube, since the duct cross-section (`WX=WY=0.04`) is
  narrower than its half-length (`WZ=0.08`). Both boundaries of this shell
  are flat, so it's a simpler nested box-in-box mesh using the same 6-panel
  construction, entirely with straight `Line`s. The `+z`/`-z` panels of
  Shell B extend all the way from the collar to the duct's actual inlet/
  outlet planes — no separate duct-extension blocks are needed.

Shell A's outer (collar) caps and Shell B's inner (collar) caps are the
**same geometric surfaces**, built once and referenced by both shells' `Volume`s
— exactly conformal by construction, the same technique
`concentricSpheres.geo` uses to share its 4 radial side surfaces between
adjacent panels.

Run: `python3 make_mesh.py` → `gmsh singlePebble.geo -3 -format msh22 -o
singlePebble.msh -nopopup` → `gmsh2nek` (interactive: `3` / `singlePebble` /
`0` / `0` / `singlePebble` — dimension, fluid mesh name, no solid mesh, no
periodic pairs, output name).

**Resolution** (`N`, `NR_A`, `NR_B` in `make_mesh.py`): `N=5` (4
elements/tangential edge, applied uniformly across both shells — all
panels/shells share edges pairwise, so a single global `N` is required for
conformity, the same constraint `concentricSpheres` has with its one `N`).
`NR_A=3` (2 radial elements, pebble→collar). `NR_B=3` (2 radial elements,
collar→duct) — **uniform across all 6 Shell-B panels**, even though the
`+-x`/`+-y` extension is only 1.2 cm (collar at 2.8 cm to duct wall at 4 cm)
while the `+-z` extension is 5.2 cm (collar to inlet/outlet at 8 cm): each
panel's radial edges are shared with its two neighbors at every corner (the
same O-grid closure constraint `concentricSpheres` has), so independently
grading the short vs. long directions isn't possible without restructuring
the block topology. This is a deliberate, coarse **demo mesh** (384 hex
elements total, matching neither the paper's mesh1/2/3 element counts nor
its convergence study) — good enough to exercise and validate the case
setup, not to reproduce the paper's mesh-converged `Nu_t` numbers. Increase
`N`/`NR_A`/`NR_B` for a production run (watch the `[RADIATION]` module's
`P<=8000` radiating-patch ceiling — see `src/core/plugins/Radiation.cpp` —
well clear of it at this resolution: 192 total participating patches).

Boundary IDs (`Physical Surface`, written directly into the `.re2`'s
`boundaryID` by `gmsh2nek` — no `usrdat2()` fix needed, same as
`concentricSpheres.usr`):

| ID | Name   | Faces | Role |
|----|--------|-------|------|
| 1  | pebble | 96    | heated sphere surface, no-slip, radiating |
| 2  | wall   | 64    | duct side walls (`+-x`,`+-y`), no-slip, adiabatic+radiating |
| 3  | inlet  | 16    | `z=-0.08`, velocity+temperature Dirichlet, non-radiating (mirror) |
| 4  | outlet | 16    | `z=+0.08`, outflow, non-radiating (mirror) |

---

## Case files

### `singlePebble.par`

`[FLUID VELOCITY]`: pebble/wall no-slip (`zeroDirichlet`), inlet
`udfDirichlet` (uniform axial velocity, see below), outlet `zeroNeumann`
(outflow). `rho=0.353`, `viscosity=6.5e-5` (Table IV).

`[SCALAR TEMPERATURE]`: pebble/wall `udfNeumann` (applied flux minus net
radiative loss — real radiative-convective coupling), inlet `udfDirichlet`
(`T=1000K`), outlet `zeroNeumann`. `diffusionCoeff=0.09` (`k`, Table IV —
the paper lists this as `W/m2`, read here as a units typo for `W/(m*K)`,
since `k=0.09` combined with the stated `Pr=0.831` and `mu=6.5e-5` backs out
a physically sensible `Cp~1150 J/(kg*K)` for air at ~1000K; see "Derived
properties" below). `transportCoeff=406.17` (`rho*Cp`).

`[RADIATION]`: all 4 boundary IDs listed in `radiatingBoundaryIDs` (the
Monte Carlo tracer needs a geometrically closed enclosure even where
`emissivity=0`, per the `tallCavityNu.par` front/back-symmetry-plane
lesson). `obstructionBoundaryIDs=1` — the pebble blocks direct wall-to-wall
visibility through its own volume (same reasoning as
`concentricSpheres.par`'s inner-sphere occlusion, generalized from a
concave self-view to a convex in-between obstacle). `emissivity =
0.85,0.85,0.0,0.0` (pebble/wall from Table IV; inlet/outlet=0, "a perfectly
reflecting (mirror) boundary condition", matching the paper's own
treatment). `updateFrequency=10`, a middle-ground default from the paper's
own `Nr` scan (Table VI: 620%/86%/12% overhead at `Nr=1/10/50`).

### `singlePebble.udf`

`udfDirichlet`: inlet (`bID=3`) uniform axial velocity and fixed
temperature. `udfNeumann`: `fluxScalar = appliedFlux - bc->usrwrk[...]`
(pebble `appliedFlux=374 W/m^2`, wall `appliedFlux=0` — a "floating"
reradiating surface, same mechanism as `tallCavityNu`'s top/bottom walls).
`UDF_Setup()`: `Radiation::setup()` + IC (`T=1000K`, `U=0`, matching
`tallCavityNu`'s "uniform IC at the dominant temperature" convention).
`UDF_ExecuteStep()`: `Radiation::step()`.

### `singlePebble.usr`

Plain `zero.usr` template, no `usrdat2()` boundary-ID fix needed (see
mesh section above).

---

## Derived properties

Table IV gives `Re=20`, `Pr=0.831`, but not the inlet velocity or `Cp`
directly:

```
U_in = Re * mu / (rho * D) = 20 * 6.5e-5 / (0.353 * 0.04) = 0.09207 m/s
Cp   = Pr * k / mu          = 0.831 * 0.09 / 6.5e-5       = 1150.6 J/(kg*K)
```

`Cp~1150 J/(kg*K)` is a physically sensible value for air near 1000K, which
is what motivated reading Table IV's `k` units (`W/m2`) as a typo for
`W/(m*K)` rather than taking the value at face value.

---

## Assumptions / choices not fully specified by the paper

- **Duct side-wall thermal BC**: the paper states inlet/outlet
  emissivity=0 explicitly but doesn't state the side-wall thermal BC. Read
  Fig. 8's "adiabatic boundary condition" caption as applying to the duct
  walls (not just the horizontal centerline sampling location) and modeled
  them as zero-applied-flux, radiatively-active surfaces — a floating
  reradiating wall, matching the treatment `tallCavityNu` already validates
  for an analogous "adiabatic but radiating" surface.
- **Inlet velocity profile**: "a uniform velocity is prescribed at the
  inlet" (paper's own words) — implemented as flat/uniform, not a
  developed pipe/duct profile.
- **Outlet BC**: not stated; used a plain `zeroNeumann` outflow (nekRS's
  standard "let it leave" BC), consistent with the paper's inlet/outlet
  emissivity=0 treatment being about radiation only, not the flow BC.
- **Mesh topology/resolution**: the paper's own mesh1/2/3 (11k/52k/208k
  elements, unstructured, presumably tet/poly with boundary-layer
  prisms from STAR-CCM+-style meshing) aren't reproducible with this
  module's structured-hex-only pipeline; built a topologically-equivalent
  coarse structured mesh instead (384 elements) as a demonstration of the
  workflow, not a mesh-convergence match to Table V's `Nu_t` values.
- **Run length**: `endTime=30 s` in the `.par`, estimated from the thermal
  diffusion time scale `L^2/alpha` (`alpha=k/(rho*Cp)~2.2e-4 m^2/s`,
  `L=0.04m` duct half-width, giving `~7.2s`) times a ~4x safety margin to
  approach quasi-steady state — not verified against a converged run in
  this repo (see "Verification" below for what was actually checked).

---

## Verification performed

- `make_mesh.py` → `gmsh` → `gmsh2nek`: clean run, no warnings, exact
  expected element/face counts (384 volume elements; 96/64/16/16
  pebble/wall/inlet/outlet faces; `gmsh2nek`'s own "velocity boundary
  faces: 192" = 96+64+16+16, matching).
- `nekrs --setup singlePebble --build-only 1` (JIT-builds the case without
  executing the time loop) was **attempted but not completed**: the
  local install's nek5000-core build (`nekio.c`, unrelated to this case)
  fails to compile against the system's mpich headers
  (`MPICH_ATTR_TYPE_TAG` / `#define MPI 1` macro collision). Confirmed this
  is a **pre-existing local build-environment issue, not a defect in this
  case's files**, by reproducing the identical failure from a fresh
  `.cache` on `examples/concentricSpheres` (an already-validated case with
  its own working, previously-built cache checked into the repo) in an
  isolated directory. Fixing the local mpich/clang toolchain mismatch is
  outside this case's scope; re-run `--build-only` after that's resolved
  (or on a machine/toolchain where the existing examples build cleanly
  from scratch) to complete this check before a production run.

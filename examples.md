<title>Radiation Examples</title>

# Radiation Examples — What They Are and How They're Built

Four example cases exercise the `[RADIATION]` module
(`src/core/plugins/Radiation.{hpp,cpp}`), in increasing order of
physical coupling:

| example | geometry | walls | coupling |
|---|---|---|---|
| `radiationPlates` | flat parallel plates | all Dirichlet | none (diagnostic only) |
| `concentricSpheres` | curved cubed-sphere shell | all Dirichlet | none (diagnostic only) |
| `tall_cavity` | box enclosure | all Dirichlet | none (diagnostic only, but real NS+energy solve) |
| `tallCavityNu` | box enclosure | 1 flux + 1 Dirichlet + 2 "adiabatic" | full radiative-convective feedback |

"Diagnostic only" means every wall's temperature is Dirichlet-fixed —
`Radiation::step()` still solves the radiosity system and reports the
resulting flux every `updateFrequency` steps, but nothing feeds back
into the temperature solve. This is the cheap, safe way to validate a
new geometry's view factors before trusting it in a coupled run.

---

## Two common problems every example has to solve

Regardless of which example, building a case that combines nekRS's
mesh/BC machinery with `[RADIATION]` runs into the same two structural
issues:

### 1. Getting boundary IDs into the mesh at all

nekRS's automatic CBC-label→numeric-boundary-ID import (`gen_bcmap`)
is disabled the instant **any** field gets an explicit
`boundaryTypeMap` in `.par` — which every one of these cases needs,
since velocity and temperature need *different* BC types on the same
raw wall IDs. Two fixes, depending on how the mesh was built:

- **`genbox` meshes** (`radiationPlates`, `tall_cavity`,
  `tallCavityNu`): `genbox` never writes a usable positive value into
  the `.re2`'s numeric BC-parameter slot for ordinary (non-periodic)
  faces — only the 3-character text CBC label survives. Fix: a
  `usrdat2()` in the case's `.usr` file that manually sets
  `boundaryID(ifc,iel)` by matching `cbc(ifc,iel,1)` against the
  numeric-string labels used in the `.box` file (`'1  '`, `'2  '`,
  etc.). All three `genbox`-based examples carry this same
  boilerplate `usrdat2()`.
- **`gmsh`-generated meshes** (`concentricSpheres`): `gmsh2nek` writes
  the `Physical Surface` tag directly into the `.re2`'s boundary-ID
  slot, so the numeric ID survives the round trip and **no**
  `usrdat2()` fix is needed — the `.usr` file is the plain empty-stub
  template.

### 2. Radiative exchange needs *every* enclosing surface accounted for

The Monte Carlo ray tracer has no notion of symmetry or "this boundary
doesn't matter" — a ray that leaves a patch and heads toward an
unlisted boundary just vanishes from the accounting instead of
returning energy to the enclosure. Two consequences, both first hit
(and fixed) in `tallCavityNu`:

- Any boundary that closes off the domain but isn't a "real" radiating
  surface (a symmetry plane standing in for an infinite/periodic
  extent, an inlet/outlet) must still be listed in
  `radiatingBoundaryIDs`, just with **`emissivity = 0`** — a diffuse
  mirror that neither absorbs nor emits net energy, but is
  geometrically present so the view-factor row sums close correctly.
  This is the same treatment Yuan et al. use for their single-pebble
  case's inlet/outlet.
- Concave surfaces need to self-occlude. `obstructionBoundaryIDs` must
  list any surface that can see itself around an obstruction (e.g. the
  outer sphere in `concentricSpheres`, occluded by the inner one) —
  radiating surfaces do **not** automatically double as occluders of
  each other.

---

## `radiationPlates`

Two `40×40` unit plates separated by a `1`-unit gap (aspect ratio 40),
plus four side walls, all built via `genbox`. The large aspect ratio
approximates the infinite-parallel-plate limit assumed by the
closed-form comparison formulas. `boundaryTypeMap` is `udfDirichlet`
everywhere; `.udf` just calls `Radiation::setup()` and
`Radiation::step()`. No flow of interest — `[FLUID VELOCITY]` uses
plain no-slip walls, since nothing forces the fluid to move.

## `concentricSpheres`

Built with `make_mesh.py`: a Python script that writes a `.geo` file
for `gmsh`'s built-in (not OpenCASCADE) kernel, constructing a "cubed
sphere" — 6 curved hex blocks, one per cube face, each bounded by two
spherical-cap surfaces (exact great-circle arcs via `Circle(id) =
{pStart, origin, pEnd}`, so curvature is exact, not approximated) and
4 shared radial side faces. Meshed with `gmsh -3 -format msh22`, then
converted with `gmsh2nek`. `radiatingBoundaryIDs = 1,2` (inner, outer)
and `obstructionBoundaryIDs = 1,2` (the outer sphere is concave and
needs the inner sphere to occlude its own self-view).

Two gmsh-specific gotchas worth knowing if extending this: (1) a
leftover `Mesh.SubdivisionAlgorithm` setting from an unrelated prior
GUI session, persisted in `~/.gmsh-options`, silently multiplied every
element 8× — fixed by setting `Mesh.SubdivisionAlgorithm = 0;`
explicitly in the generated `.geo` rather than relying on ambient
state; (2) `gmsh2nek` requires a `$PhysicalNames` section in the
`.msh`, which gmsh only writes for *named* physical groups
(`Physical Surface("hot", 1) = {...}`, not the bare-numeric form).

## `tall_cavity`

A tall (aspect ratio 4) rectangular-prism box, `genbox`-built, one hot
wall, one cold wall directly opposite, four side walls at the
hot/cold mean — all Dirichlet, all six walls radiating (the box
interior is convex, so no `obstructionBoundaryIDs` needed). Unlike
`radiationPlates`, this one runs a real Boussinesq natural-convection
Navier–Stokes + energy solve (buoyancy source in `userf`, matching the
`rbc.udf` pattern but with dimensional `g·β` since real
Stefan-Boltzmann flux needs real Kelvin temperatures) — but since
every wall's temperature is still Dirichlet-fixed, `Radiation::step()`
remains a pure diagnostic layered on top of a genuine flow.

## `tallCavityNu`

Reproduces Yuan et al. Sec 4.3: aspect ratio 20 cavity, one wall at a
constant *applied* heat flux, one wall at fixed temperature, top and
bottom "adiabatic" in the paper's sense (no external flux) but still
radiatively active — floating, reradiating surfaces. This requires
the real `udfNeumann`/`o_usrwrk` coupling (see
`radiative_heat_transfer_in_nekRS.md` for the mechanism and sign
convention): the flux actually entering the fluid at the hot wall is
`appliedFlux - radiativeLoss`, and at the two "adiabatic" walls it's
`0 - radiativeLoss` — the same `udfNeumann` function handles both by
branching on `bc->id`. Front/back are symmetry planes standing in for
the paper's quasi-2D domain (a thin-slab extrusion), non-radiating.

`dt` uses `targetCFL=2.0 + max=0.5 + initial=0.01` — an explicit small
initial step, since the case starts from rest (zero velocity), where a
pure CFL target can't sensibly bound anything.
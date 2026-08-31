# Concentric-spheres Nek-VF verification case

Reproduces paper Section 4.2 ("Concentric spheres") using the `nek-vf-case`
skill workflow. Like `parallel_plates_vf`, this is a **radiation-only
verification case** — no forced or buoyant flow — checking Nek-VF's
radiosity solve against the analytic two-concentric-sphere formula (paper
Eq. 12). It is also the case that surfaced and drove the fix of a serious,
previously-undiscovered bug in `RadiativeViewFactor.jl`'s 3-D obstruction
path (see below) — obstruction is not optional here, unlike the other two
new cases: the outer sphere's concave self-view is physically meaningless
without it.

## What the paper specifies vs. what this reproduction derives

The paper (Table II) gives both sphere temperatures (inner/hot = 1000 K,
outer/cold = 500 K) and one global emissivity (0.85), but **does not state
either sphere's radius** — only their ratio is recoverable, and only
indirectly: Eq. 12 combined with the tabulated analytical flux values
implies `A1/A2 = 0.25` (verified: plugging `A1/A2=0.25`, `ε=0.85`, and the
given temperatures into Eq. 12 gives `q₁'' = 4.355e4 W/m²`, matching the
paper's `4.36e4` to rounding). Since `A1/A2 = (r1/r2)²` for spheres, this
mesh uses `r1=0.5 m`, `r2=1.0 m` (`concentric_spheres.geo`) — the absolute
scale is unconstrained by the paper (the flux formula and view factors
depend only on the ratio), so `r2=1.0` was chosen for convenience.

## File-by-file layout

| File | Produced by | Role |
|---|---|---|
| `concentric_spheres.geo` | hand-written | Gmsh "cubed-sphere" script: two concentric shells, `Physical Surface` tags `hot`(1, inner)/`cold`(2, outer) |
| `concentric_spheres.msh` | `gmsh -3 concentric_spheres.geo` | Gmsh mesh (1470 hex elements: 6 curved blocks × 7×7×5 divisions) |
| `concentric_spheres.re2` | `gmsh2nek` | Nek5000 binary mesh + boundary tags |
| `vf_concentric_spheres` | `compute_view_factors.jl` (RadiativeViewFactor.jl), **with obstruction enabled** | Precomputed view factors, NekRS-format (588 walls) |
| `concentric_spheres.par` | hand-written | NekRS case config |
| `concentric_spheres.usr` | hand-written | Legacy Nek5000 backend: `usrdat2` boundary-tag mapping, `userchk` radiosity driver |
| `concentric_spheres.udf` | hand-written | Host-side GPU↔CPU `frad` bridging, `Nr=1` |
| `concentric_spheres.oudf` | hand-written | OKL: `codedFixedValueScalar` (both walls Dirichlet), `fillProp` |
| `view_factors/` | copied verbatim from `tall_cavity_vf_aurora/view_factors/` | Generic radiosity/irradiation solver module |

## How the mesh was generated: a cubed-sphere O-grid

Nek5000/NekRS need hexahedral elements, and Gmsh can't directly
transfinite-mesh a spherical shell — so `concentric_spheres.geo` builds it
as a **cubed sphere**: 6 curved hex blocks, one per face of an imaginary
inscribed cube, each bounded by an inner-sphere patch (4 great-circle
arcs), an outer-sphere patch (4 more arcs, same directions, larger radius),
and 4 straight radial "fin" edges connecting corresponding corners. Each
block is `Transfinite Volume` + `Recombine Volume` meshed; the 8 cube-corner
directions, 12 great-circle arcs (`Circle(id) = {p1, center, p2}`, valid
here because inner/outer corner pairs are always equidistant from the
sphere center), and 24 bounding surfaces were generated programmatically
(a small Python script wrote the `.geo` text) rather than hand-typed, since
the index bookkeeping for 6 shared-edge blocks is mechanical but
error-prone by hand.

```bash
/Applications/Gmsh.app/Contents/MacOS/gmsh -3 concentric_spheres.geo -o concentric_spheres.msh -format msh2
```

Uses the same `Mesh.SubdivisionAlgorithm = 0;` override documented in
`parallel_plates_vf/procedure.md` (stale global Gmsh GUI preference on this
machine).

```bash
export PATH="/Users/coxea3/opentools/Nek5000/bin:$PATH"
printf "3\nconcentric_spheres\n0\n0\nconcentric_spheres\n" | gmsh2nek
```

→ `concentric_spheres.re2` (588 boundary faces, 1 group `MSH` — both `hot`
and `cold` get the same generic non-periodic boundary-condition code from
`gmsh2nek`; see below for why that matters).

## How view factors were computed — and the bug this case found

```bash
julia -t 8 --project=/Users/coxea3/.julia/dev/RadiativeViewFactor.jl \
  ~/.claude/skills/nek-vf-case/scripts/compute_view_factors.jl \
  concentric_spheres.re2 vf_concentric_spheres \
  --backend=cpu --monte-carlo=false --nquad=8 --obstruction=1 --closure-tol=0.02
```

`--obstruction=1` is not optional for this geometry: the outer (cold,
**concave**) sphere partially sees itself around the inner sphere, and the
exact closed-form self-view is `F(outer→outer) = 1 - (r1/r2)² = 0.75` (by
conservation, since `F(outer→inner) = (r1/r2)² = 0.25` by reciprocity from
`F(inner→outer) = 1`, exact for any convex body fully enclosed by another
surface). Computing this correctly requires the inner sphere to be treated
as an *obstructor* of the outer sphere's own self-visibility.

Running this driver first surfaced a row-sum closure error of ~25%
(`F(outer→outer)` computing to ≈1.0 instead of 0.75) — and, diagnostically,
the error got **worse**, not better, when the obstructing inner-sphere mesh
was refined, which is backwards for a discretization/quadrature accuracy
problem and pointed at a logic bug instead. Root-caused via a chord-vs-ball
angular sweep (a chord between two points on the unit sphere enters a
r1=0.5 ball exactly when their angular separation exceeds `2·acos(0.5) =
120°`) down to two real, previously-undiscovered bugs in
`RadiativeViewFactor.jl` itself, both now fixed and covered by a new
`test/obstruction_test.jl`:

1. **`obstruction_groups` couldn't select a single named surface at all**
   on this kind of mesh: `load_re2` groups elements by Nek boundary-code
   (`'W'`/`'P'`/generic), and `gmsh2nek` writes the *same* generic code for
   every non-periodic surface — so `hot` and `cold` were indistinguishable
   via `group` even though `usrdat2`'s numeric `bc_flag` tells them apart.
   Fixed by exposing that discarded tag as `SurfaceElement.phys_tag` and
   adding `split_groups_by_tag`, which the driver script now calls
   automatically (see its `--obstruction=` flag, now numeric-tag-based).
2. **The actual 3-D ray/triangle obstruction test was broken** —
   `BVH.jl`'s `intersect_ray_bvh` (and the GPU equivalent) read each
   triangle's vertex coordinates with the wrong array axes transposed,
   testing three bogus points instead of the real vertices. This affected
   essentially all obstruction queries on any 3-D mesh, on both CPU and
   GPU, independent of this case. See
   `RadiativeViewFactor.jl/v062_changelog.md` (§6-§8) for the full
   root-cause writeup and verification numbers.

After both fixes: row-sum closure error 0.017% (down from ~25%),
reciprocity exact to 3.9e-16, `F(outer→outer)=0.7502`,
`F(outer→inner)=0.2498` (analytic: 0.75, 0.25).

## How the case files are wired together

Same radiation-only-verification pattern as `parallel_plates_vf`: no
flux-type boundary (`codedFixedGradientScalar` unneeded — both boundaries
are Dirichlet), no forced flow, `Nr=1`. The only case-specific difference
is `userchk()` reporting on 2 boundaries instead of 3, and (critically) the
view-factor file having been computed *with* obstruction enabled — get
this wrong for any future obstructed geometry (concave surfaces, pebbles,
internal blockers) and the resulting case will silently under-count
self-view/over-count direct view, exactly as this case did before the fix.

## Build and run (not yet executed)

```bash
$NEKRS_HOME/bin/nekrs --setup concentric_spheres --build-only <ntasks>
$NEKRS_HOME/bin/nekrs --setup concentric_spheres
```
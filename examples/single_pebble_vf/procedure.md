# Single-pebble Nek-VF case

Reproduces paper Section 4.4 ("Flow around a single pebble") using the
`nek-vf-case` skill workflow: a spherical pebble at the center of a square
duct, forced convection plus surface-to-surface radiation. Unlike
`parallel_plates_vf`/`concentric_spheres_vf`, this is a real flow case
(Table IV: Re=20, Pr=0.831), not a radiation-only verification — matching
the paper's own comparison against STAR-CCM+ S2S.

## What the paper specifies vs. what this reproduction derives

The paper's Table IV gives κ, ρ, μ, the pebble's imposed surface heat flux,
inlet temperature, wall emissivity, Re, and Pr directly; duct cross-section
(8×8 cm), length (16 cm), and pebble diameter (4 cm) are given in the text.
Two quantities used by NekRS's `.par`/`.usr` files are **not** given
directly and were derived from Table IV:

- **Inlet velocity** `U0 = Re·μ/(ρ·L) = 20·6.5e-5/(0.353·0.08) ≈ 0.046034 m/s`
  (`L` = duct side length, used as the Reynolds-number length scale — the
  paper doesn't state which length scale it uses, but the duct side is the
  natural/conventional choice for internal flow and is the assumption used
  in `single_pebble.oudf`'s `codedFixedValueVelocity`).
- **Specific heat** `Cp = Pr·κ/μ = 0.831·0.09/6.5e-5 ≈ 1150.6 J/(kg·K)`
  (Table IV gives Pr but not Cp directly), giving
  `rhoCp = ρ·Cp ≈ 406.17 J/(m³·K)` for `single_pebble.par`.

`numSteps=2000`, `dt=1e-3`, and `Nr=10` (`single_pebble.udf`) are
reasonable starting points, not values validated against the paper's own
run — the paper doesn't state its timestep or total run length for this
specific mesh size, only overhead numbers at three much larger mesh
resolutions (paper Table V/VI, 11k-208k elements; this reproduction is
deliberately much coarser, ~2400 elements, for tractability).

## File-by-file layout

| File | Produced by | Role |
|---|---|---|
| `gen_single_pebble_geo.py` | hand-written | Python generator that writes `single_pebble.geo` (see below for why) |
| `single_pebble.geo` | `python3 gen_single_pebble_geo.py` | Gmsh script: two-layer O-grid, `Physical Surface` tags `pebble`(1)/`inlet`(2)/`outlet`(3)/`duct_wall`(4) |
| `single_pebble.msh` | `gmsh -3 single_pebble.geo` | Gmsh mesh (2376 hex elements) |
| `single_pebble.re2` | `gmsh2nek` | Nek5000 binary mesh + boundary tags |
| `vf_single_pebble` | `compute_view_factors.jl`, **with `--obstruction=1`** | Precomputed view factors, NekRS-format (432 walls) |
| `single_pebble.par` | hand-written | NekRS case config |
| `single_pebble.usr` | hand-written | Legacy Nek5000 backend: `usrdat2`, `userchk` radiosity driver |
| `single_pebble.udf` | hand-written | Host-side GPU↔CPU `frad` bridging, `Nr=10` |
| `single_pebble.oudf` | hand-written | OKL: inlet velocity/temperature, pebble/duct-wall flux BCs |
| `view_factors/` | copied verbatim from `tall_cavity_vf_aurora/view_factors/` | Generic radiosity/irradiation solver module |

## How the mesh was generated: a two-layer O-grid

Nek5000/NekRS need hex elements, and a duct with a spherical hole in the
middle isn't directly transfinite-meshable — `gen_single_pebble_geo.py`
builds it as **two nested "cubed sphere" shells**, both using the same
6-block shared-edge topology as `concentric_spheres.geo`:

1. **Layer 1 (curved)**: pebble surface (sphere, R=2cm) ↔ a small
   reference cube (half-side 2cm) circumscribing it, face-centers
   touching the sphere. Corner edges are great-circle arcs on the sphere
   side, straight lines on the cube side, connected by straight radial
   spokes.
2. **Layer 2 (flat)**: the same reference cube ↔ the duct's actual outer
   walls (a non-cubic box, half-extents 4cm/4cm/8cm) — a plain nested-box
   extension reusing Layer 1's outer corners, entirely straight edges.

This gives 12 hex volumes (6 per layer), 2376 elements total. The
geometry was generated **programmatically** (`gen_single_pebble_geo.py`,
kept in this directory) rather than hand-typed in `.geo` syntax directly:
the point/curve/surface/volume index bookkeeping for two layers of 6
shared-edge blocks each (8 points + 12 edges + 6 faces + 12 fins per
layer) is mechanical but highly error-prone by hand.

```bash
/Applications/Gmsh.app/Contents/MacOS/gmsh -3 single_pebble.geo -o single_pebble.msh -format msh2
```

Same `Mesh.SubdivisionAlgorithm = 0;` override as the other two new cases
(stale global Gmsh GUI preference on this machine — see
`parallel_plates_vf/procedure.md`).

**Gmsh reported 96 non-right-hand (inverted-winding) elements** on meshing
(out of 2376) — a known, standard situation for hand-built multi-block
meshes where `Transfinite Volume`'s automatic corner-ordering detection
doesn't consistently pick the same parametrization handedness across
topologically-identical but differently-oriented blocks. `gmsh2nek`
prompts to auto-fix these (a safe, purely combinatorial node-reordering
that doesn't change element geometry); answering "y" resolves it. A
smaller reference cube (0.75×R instead of R, tried to see whether
avoiding the face-center "pinch" — zero shell thickness where the
reference cube face exactly touches the sphere at RC=R — would eliminate
the inversions) made it *worse* (504 inverted elements, not fewer), so the
actual cause is a different, not-fully-root-caused orientation quirk
rather than element degeneracy — kept RC=R since it validates cleanly
(next section) and the fix is standard practice, not a workaround for a
broken mesh.

```bash
export PATH="/Users/coxea3/opentools/Nek5000/bin:$PATH"
printf "3\nsingle_pebble\n0\ny\n0\nsingle_pebble\n" | gmsh2nek
```

(dimension, filename, no solid mesh, **fix left-hand elements: y**, 0
periodic pairs, output basename.)

## How view factors were computed

```bash
julia -t 8 --project=/Users/coxea3/.julia/dev/RadiativeViewFactor.jl \
  ~/.claude/skills/nek-vf-case/scripts/compute_view_factors.jl \
  single_pebble.re2 vf_single_pebble \
  --backend=cpu --monte-carlo=false --nquad=8 --obstruction=1 --closure-tol=0.01
```

`--obstruction=1` (the pebble's Physical Surface tag) is essential: the
pebble sits directly in the line of sight between opposite duct-wall
patches, and without obstruction those patches would incorrectly "see"
each other straight through the pebble. This only works correctly because
of the two `RadiativeViewFactor.jl` fixes documented in
`concentric_spheres_vf/procedure.md` and the package's own
`v062_changelog.md` (§6-§8) — found while building that case, not this
one, but load-bearing here too.

**Plain quadrature, not Monte Carlo, was used for the final file**: an
initial attempt with `--monte-carlo=true` (the skill's usual recommendation
for GPU/large-mesh performance) gave a stubborn ~12.7% row-sum closure
error that did *not* improve with higher `--nquad`/`--n-samples`/`--factor`
— re-running the identical mesh with plain deterministic quadrature
(`--monte-carlo=false --nquad=8`) instead closed to 0.025% immediately,
confirming the mesh geometry itself is sound and pointing at a **separate,
not-yet-root-caused interaction between Monte Carlo sampling (or its
near-pair Duffy patch) and obstruction** as the actual cause — worth
investigating in `RadiativeViewFactor.jl` directly if a future, much
larger obstructed case needs GPU/Monte Carlo performance. For this case's
432-element mesh, CPU deterministic quadrature is fast enough regardless.

## How the case files are wired together

This is the paper's most complete case, combining real forced convection
with radiation:

- **Velocity** (`[VELOCITY] boundaryTypeMap = wall, inlet, outlet, wall`):
  pebble and duct walls no-slip; inlet gets a uniform `codedFixedValueVelocity`
  (`0.046034 m/s` in `+z`, matching the geometry's flow direction — inlet
  is the `z=-Lz` duct face, outlet `z=+Lz`); outlet is a native NekRS
  outflow condition, no coded value needed. **Note**: `single_pebble.par`
  uses `"inlet"`/`"outlet"`/`"wall"`, not `"codedFixedValue"` as
  `tall_cavity_vf_aurora`'s own `.par` does — see the note in
  `parallel_plates_vf/procedure.md`: `"codedFixedValue"` is not a
  recognized `boundaryTypeMap` keyword in this checked-out NekRS source
  (`src/core/bdry/bdryBase.cpp`'s keyword table has no such entry), so all
  three new cases use the verified-correct keywords instead.
- **Temperature** (`boundaryTypeMap = flux, inlet, outlet, flux`): pebble
  gets `codedFixedGradientScalar`'s `374.0 - gflux` (Table IV's imposed
  heat flux plus radiative exchange); duct walls get `0.0 - gflux`
  (adiabatic except for radiation); inlet is Dirichlet at 1000 K; outlet
  is a plain zero-gradient outflow.
- **Zero-emissivity inlet/outlet, without touching `view_factors.f`**: the
  paper models the inlet/outlet as perfectly-reflecting mirrors (ε=0).
  `view_factors.f`'s `indicator_real_wall_ef` — which surface actually
  emits/absorbs versus just reflects — is set purely from each boundary's
  Nek `cbc` code (`1.0` only for `cbc='W  '`, `0.0` for anything else
  non-internal), automatically, with no separate per-boundary-emissivity
  input. Since `usrdat2` already gives inlet/outlet `cbc='v  '`/`'O  '`
  (not `'W  '`) for the flow solve, they automatically become non-emitting
  "mirror" surfaces in the radiosity solve too — exactly the paper's
  intended treatment, with no module changes needed. This resolves the
  concern the `nek-vf-case` `SKILL.md`'s "Known limitations" section
  raises about `vf_eps` being a single global value: that limitation is
  about **wall-to-wall** emissivity variation (e.g. wanting the pebble and
  duct walls at different non-zero values), which this case doesn't need —
  the zero-emissivity inlet/outlet is a completely different, already-
  supported mechanism.

## Build and run (not yet executed)

```bash
$NEKRS_HOME/bin/nekrs --setup single_pebble --build-only <ntasks>
$NEKRS_HOME/bin/nekrs --setup single_pebble
```

Watch for `VF:` log lines from `vf_print_radiation_heat_flux(1)` (pebble
boundary) — the paper's own comparison metric is the pebble's total
Nusselt number (Table V), based on average pebble surface temperature and
volume-averaged fluid temperature, which would need to be computed as a
post-processing step on top of this case's raw `frad`/temperature output.
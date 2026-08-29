# Parallel-plates Nek-VF verification case

Reproduces paper Section 4.1 ("Parallel plates") using the `nek-vf-case`
skill workflow: Gmsh geometry → `gmsh2nek` → `RadiativeViewFactor.jl` view
factors → NekRS case files. This is a **radiation-only verification case** —
there is no forced or buoyant flow, unlike `tall_cavity_vf_aurora`. Its only
purpose is to check that Nek-VF's radiosity solve reproduces the analytical
two-surface radiative-exchange formula (paper Eq. 11) for a closed, diffuse-
gray enclosure with prescribed Dirichlet wall temperatures.

## What the paper specifies vs. what this reproduction assumes

The paper (Table I) gives the hot-plate temperature (1000 K), cold-plate
temperature (500 K), side-wall temperature (750 K), and a single global
emissivity (0.85) for all boundaries — but **does not state the plate
dimensions or plate spacing**. Eq. 11 (the two-infinite-parallel-plate
formula) is only exact in the limit where the plates are large relative to
their separation; the paper's own Nek-VF result (3.92e4 vs. an analytical
3.93e4 W/m2) implies whatever finite geometry they used had a
direct hot→cold view factor close to 1.

Since the exact aspect ratio isn't recoverable from the paper text, this
reproduction uses an assumed 4:1 plate-side-to-spacing ratio (`W=4.0 m`,
`H=1.0 m` in `parallel_plates.geo`) — large enough to approximate the
infinite-plate limit while remaining a genuine closed 6-surface enclosure.
**The resulting Nek-VF flux values will not exactly match Table I** (a
smaller direct view factor than the paper's implied geometry means somewhat
more energy is redistributed to/through the side walls); this is a
documented deviation from the source, not a bug. If the true geometry is
known, only `W` and `H` in `parallel_plates.geo` need to change — everything
downstream (mesh, view factors, case files) regenerates unchanged.

## File-by-file layout

| File | Produced by | Role |
|---|---|---|
| `parallel_plates.geo` | hand-written | Gmsh script: unit-ish box, `Physical Surface` tags `hot`(1)/`cold`(2)/`side`(3) |
| `parallel_plates.msh` | `gmsh -3 parallel_plates.geo` | Gmsh mesh (4400 hex elements, 2nd-order/curved-capable per Nek5000 convention) |
| `parallel_plates.re2` | `gmsh2nek` | Nek5000 binary mesh + boundary tags NekRS reads at startup |
| `vf_parallel_plates` | `compute_view_factors.jl` (RadiativeViewFactor.jl) | Precomputed view factors, NekRS-format (1680 walls) |
| `parallel_plates.par` | hand-written | NekRS case config |
| `parallel_plates.usr` | hand-written | Legacy Nek5000 backend: `usrdat2` boundary-tag mapping, `userchk` radiosity driver |
| `parallel_plates.udf` | hand-written | Host-side GPU↔CPU `frad` bridging, `Nr=1` coupling cadence |
| `parallel_plates.oudf` | hand-written | OKL: `codedFixedValueScalar` (all 3 walls Dirichlet), `fillProp` |
| `view_factors/` | copied verbatim from `tall_cavity_vf_aurora/view_factors/` | Generic radiosity/irradiation solver module (case-independent) |

## How the mesh was generated

```bash
/Applications/Gmsh.app/Contents/MacOS/gmsh -3 parallel_plates.geo -o parallel_plates.msh -format msh2
```

**Gotcha hit and fixed**: a stale global Gmsh GUI preference
(`~/.gmsh-options` had `Mesh.SubdivisionAlgorithm = 2`, "all hexahedra",
left over from unrelated reactor-meshing work on this machine) silently
8×-subdivided every transfinite hex even on command-line invocation,
turning the intended 4400-element mesh into 35200 elements with no error or
warning. Fixed by adding `Mesh.SubdivisionAlgorithm = 0;` at the top of the
`.geo` script itself (overrides the loaded default; does not touch the
user's global preferences file). **Any new `nek-vf-case` geometry on this
machine should include the same override line** — it's now in all three new
cases' `.geo` files.

```bash
export PATH="/Users/coxea3/opentools/Nek5000/bin:$PATH"
printf "3\nparallel_plates\n0\n0\nparallel_plates\n" | gmsh2nek
```

(dimension=3, fluid mesh filename, no solid mesh, 0 periodic pairs, output
basename) — produced `parallel_plates.re2` (1680 boundary faces: 1 group
`MSH`, since none of these faces are periodic/internal).

## How view factors were computed

```bash
julia -t 8 --project=/Users/coxea3/.julia/dev/RadiativeViewFactor.jl \
  ~/.claude/skills/nek-vf-case/scripts/compute_view_factors.jl \
  parallel_plates.re2 vf_parallel_plates \
  --backend=cpu --monte-carlo=true --n-samples=4000 --nquad=6 --closure-tol=0.1
```

Result: row-sum closure error 0.49% max (well within tolerance), reciprocity
exact to machine precision (3.8e-16). CPU/Monte-Carlo was used since 1680
elements is small enough that GPU dispatch overhead isn't worth it; this is
still the same Monte-Carlo-with-near-pair-Duffy-patch kernel the skill
recommends for larger meshes, just run on CPU.

## How the case files are wired together

This case is a stripped-down instance of the same Nek-VF wiring documented
in `tall_cavity_vf_aurora/procedure.md`, with two simplifications specific
to a radiation-only verification case:

1. **No flux-type boundary.** All three walls are `codedFixedValue`
   (Dirichlet) in `[TEMPERATURE] boundaryTypeMap` — there is no
   `codedFixedGradientScalar` in `parallel_plates.oudf` at all, because no
   boundary needs `frad` fed back into the energy equation as a BC. The
   radiative flux (`frad`) is read out purely as a diagnostic via
   `vf_print_radiation_heat_flux`, called three times in `userchk()` (once
   per boundaryID) to match the paper's Table I reporting all three walls,
   rather than the single-call pattern used for `tall_cavity`'s one flux
   wall.
2. **No forced or buoyant flow.** `[VELOCITY] boundaryTypeMap` is `wall`
   (native no-slip) on all three boundaries, with zero initial velocity and
   no `userVelocitySource` — the velocity field stays identically zero for
   all time. `Nr=1` in `parallel_plates.udf` (re-solve radiosity every
   step) since there's no flow-radiation feedback loop to amortize the way
   `tall_cavity`'s `Nr` does; the Dirichlet-prescribed wall temperatures
   converge to their BC values within the first solve, at which point
   `frad` is already the final answer.

The `nekrs_registerPtr('frad', frad)` / `nek::ptr<double>("frad")` bridge,
`vf_read_view_factors`/`vf_calculate_radiation_heat_flux` call sequence, and
`o_usrwrk` upload path are otherwise identical to `tall_cavity_vf_aurora` —
see that case's `procedure.md` for the full call-chain trace through
NekRS's source (`nekInterfaceAdapter.cpp`, `nrs.cpp`).

## Build and run (not yet executed)

```bash
$NEKRS_HOME/bin/nekrs --setup parallel_plates --build-only <ntasks>
$NEKRS_HOME/bin/nekrs --setup parallel_plates
```

Watch for `VF:` log lines from the three `vf_print_radiation_heat_flux`
calls at `istep=0`'s `userchk` — since wall temperatures are Dirichlet and
the mesh/view-factor pipeline above is already closure-validated, the flux
values should stabilize by the first checkpoint (`numSteps=5`).
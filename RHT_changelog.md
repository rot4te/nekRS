<title>RHT Branch Changelog</title>

# `rht` vs `master` — Changelog

Everything below is new relative to `master` (51 files changed, 5608
insertions). It adds a Monte Carlo radiative-heat-transfer module
(`[RADIATION]` in `.par`) plus four example cases exercising it, and one
unrelated one-line fix picked up along the way.

---

## Core module

### `src/core/plugins/Radiation.hpp` / `Radiation.cpp` (new)
The module itself. Two entry points a case calls: `Radiation::setup()`
(once, in `UDF_Setup()`) computes the view-factor matrix between
boundary-ID groups via Monte Carlo ray sampling; `Radiation::step(time,
tstep)` (every step, in `UDF_ExecuteStep()`) solves the gray-diffuse
radiosity equations from current wall temperatures and writes the
resulting net radiative flux into `platform->app->bc->o_usrwrk`, where
a case's `udfNeumann` can read it. See
`radiative_heat_transfer_in_nekRS.md` for full usage.

### `src/core/plugins/RadiationBVH.hpp` / `RadiationBVH.cpp` (new)
A host-built bounding-volume hierarchy over the mesh's triangulated
boundary faces, used by the ray tracer to test for occlusion between
patch pairs when `obstructionBoundaryIDs` is set (needed for concave
enclosures — see below).

### `src/core/plugins/kernels/Radiation.okl` (new)
The device (OCCA) kernel. For every radiating-patch pair, stratified
Monte Carlo samples two points (one per patch, via the patch's own
isoparametric/barycentric Lagrange geometry), accumulates the
`cosθᵢcosθⱼ/(πr²)` view-factor kernel, and optionally tests the segment
against the BVH for occlusion. Includes a Duffy-transform-style
near-field reparametrization for patch pairs flagged geometrically
"close" (see `radiative_heat_transfer_in_nekRS.md` for why, and for the
subtle decorrelation-hash bug found and fixed while implementing it).

## Build / integration plumbing

### `cmake/core.cmake`
Adds the two new `.cpp` files to `CORE_SOURCES` so they compile into
`libnekrs`.

### `src/core/udf/udfMake.hpp`
Registers `Radiation::buildKernel` in the UDF auto-load-plugins table,
so any case whose `.udf` includes `Radiation.hpp` gets its device
kernel built automatically (same mechanism `lpm`, `tavg`, `RANSktau`
etc. already use).

### `src/platform/par/par.cpp` / `src/platform/par/parseRadiation.hpp` (new)
Adds `[radiation]` as a recognized `.par` section and
`parseRadiationSection()`, which validates and forwards its keys
(`radiatingBoundaryIDs`, `obstructionBoundaryIDs`, `nSamples`, `seed`,
`writeMatrix`, `outputFile`, `cache`, `emissivity`, `stefanBoltzmann`,
`updateFrequency`, `radiosityTolerance`, `radiosityMaxIters`) into
`setupAide` options the module reads at runtime.

## Examples (new directories)

### `examples/radiationPlates/`
Two large, closely-spaced parallel plates (hot/cold/sides), all
boundaries Dirichlet-fixed — a pure radiation diagnostic, no flow of
interest. Validates the view-factor solve against the closed-form
Hottel/Feingold parallel-plate formula and, separately, against Yuan et
al.'s Sec 4.1 benchmark. See `radiationPlates/radPlates.md`.

### `examples/concentricSpheres/`
Two concentric spheres (inner hot, outer cold, both Dirichlet), meshed
as a "cubed sphere" via a `gmsh`-driven Python generator
(`make_mesh.py`) since `genbox` can't produce curved geometry.
Validates against Yuan et al. Sec 4.2 and the exact analytic
concentric-sphere view factor. See `concentricSpheres/concentricSpheres.md`.

### `examples/tall_cavity/`
A tall rectangular-prism enclosure, one hot wall, one cold wall
directly opposite, four side walls at the hot/cold mean — all
Dirichlet, all six walls radiating. A real natural-convection
Navier–Stokes + energy solve (not just a radiation diagnostic) coupled
to `Radiation::step()`'s diagnostic flux reporting. See
`tall_cavity/tall_cavity.md`.

### `examples/tallCavityNu/`
The genuinely coupled case: reproduces Yuan et al. Sec 4.3 (tall
cavity, aspect ratio 20, one wall at constant applied heat flux, one
wall at fixed temperature, two "adiabatic" walls that are still
radiatively active). This is the first case to exercise the real
`udfNeumann`/`o_usrwrk` radiative-convective feedback loop end to end
(previous cases only used `Radiation::step()`'s diagnostic output on
Dirichlet walls). See `examples.md` for how it's built and
`radiative_heat_transfer_in_nekRS.md` for the coupling mechanism it
exercises.

## Documentation (new)

- `rvf-dev.md` — file-by-file summary of the original (pre-coupling)
  view-factor-only module.
- `nekRS_rvf.md` — full session transcript from the module's initial
  design/implementation.
- `RHT_changelog.md` (this file), `examples.md`,
  `radiative_heat_transfer_in_nekRS.md` — written after the coupling
  work and all four examples were complete.

## Unrelated

### `3rd_party/adios/examples/simulations/GrayScott.jl/src/analysis/pdfcalc.jl`
One-line syntax fix: `else if` → `elseif` (Julia requires the latter).
Unrelated to the radiation work; picked up incidentally.

### `.gitignore`
Adds `*.code-workspace`.
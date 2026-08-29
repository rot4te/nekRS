# Tall-cavity view-factor radiation case: setup procedure

This document maps the "4.3 Natural convection in tall cavity" case in the
Nek-VF paper (Yuan, Dai, Shaver, Gottems, Merzari — *Implementation of the
Surface-to-Surface Thermal Radiation Solver in Spectral Element CFD code
NekRS*) onto the concrete files in this directory (`tall_cavity_vf_aurora/`),
and lists, in order, every step that was taken to build this case from a
bare template up to a runnable NekRS case on Aurora. Findings below were
checked directly against the file contents and against the NekRS/Nek5000
source (`src/core/nekInterface/`, `src/app/nrs/nrs.cpp`) rather than assumed.

The case-file-level documentation already written into this directory
(`README.md` and `view_factors/README.md`) is the authoritative reference
for *what each routine does*; this document focuses on the *setup workflow*
— what was created, in what order, and why — and on how that maps back to
the paper's Section 2 (numerical method), 3.2 (radiative heat flux
calculation), and 4.3 (tall-cavity verification).

---

## 1. How the paper's model maps onto these files

| Paper element | Equation/Algorithm | File(s) in this case |
| --- | --- | --- |
| Radiosity equation, `J_i = ε_i σ T_i^4 + (1-ε_i) H_i` | Eq. 1 | `view_factors/view_factors.f`, subroutine `vf_calculate_radiation_heat_flux` |
| Irradiation, `H_i = Σ_j F_i→j J_j` | Eq. 2 | same subroutine |
| View factors `F_i→j` (precomputed, not solved in NekRS) | Eq. 6/7 | `vf_tall_cavity` (binary/text view-factor file), read by `vf_read_view_factors` |
| Net radiative flux `q_i = J_i - H_i` | Eq. 4 | `frad_ef`/`frad` arrays, computed in `vf_calculate_radiation_heat_flux`, consumed in `tall_cavity.oudf`'s `codedFixedGradientScalar` |
| Iterative radiosity/irradiation solve (Algorithm 2) | Alg. 2 | `vf_calculate_radiation_heat_flux`'s 5-iteration loop (`niter = 5`), called from `userchk()` in `tall_cavity.usr` |
| Radiosity update frequency `N_r` | Section 3.2 | `tstep % 100 == 0` cadence in `UDF_ExecuteStep` (`tall_cavity.udf`) — this case effectively uses `N_r = 100` |
| Hemi-cube GPU view-factor calculator ("hemicubevf") | Section 3.1 | **not present in this directory or in the NekRS source tree** — it is a separate, standalone OpenGL/GLFW/GLAD tool used offline to produce `vf_tall_cavity`; only its *output* is checked into this case |
| Tall-cavity geometry, H/D = 20, one flux wall + one fixed-T wall | Section 4.3 | `tall_cavity.geo` (`H=1.0`, `L=H/20`, `Physical Surface "hot"`/`"cold"`) |
| Boussinesq buoyancy | Section 4.3 | `userf()` in `tall_cavity.udf` |

One important nuance: the paper's Algorithm 2/Section 3.2 describes the
*generic* Nek-VF radiosity solve; this specific case is built on the
**CPU-side, legacy-Nek5000 Fortran implementation** of that solve
(`view_factors.f`), not the GPU hemi-cube calculator described in Section
3.1. The hemi-cube tool is used **once, offline**, purely to generate the
`vf_tall_cavity` view-factor file consumed at runtime — it plays no further
role once that file exists.

---

## 2. Step-by-step: how this case was built

### Step 1 — Geometry and mesh generation (Gmsh)

`tall_cavity.geo` defines the H×(H/2)×L extruded-rectangle enclosure from
Section 4.3 (`H=1.0`, `L=H/20`, aspect ratio H/D=20 once combined with the
extrusion depth), tagged with the five physical surfaces needed by the
radiation/BC setup (`hot`=1, `cold`=2, `topandbot`=3, `p1`=4, `p2`=5) and one
physical volume (`fluid`). This was meshed with Gmsh:

```
gmsh -3 tall_cavity.geo -o tall_cavity.msh
```

producing `tall_cavity.msh` (present in this directory).

### Step 2 — Convert to Nek5000 mesh format

The Gmsh mesh was converted to Nek5000's native format with the `gmsh2nek`
tool (part of the Nek5000 toolchain, e.g.
`Nek5000/tools/gmsh2nek`), which reads `tall_cavity.msh` and the physical
surface tags and writes:

- `tall_cavity.re2` — the binary mesh/connectivity/boundary-tag file NekRS
  actually reads at startup.
- `tall_cavity.co2` — companion connectivity/curved-side data.

The five physical surfaces (`bc(5,ifc,iel,1)` in `usrdat2`) survive this
conversion as the numeric boundary tags 1–5 that `tall_cavity.usr`'s
`usrdat2()` later maps to Nek's `cbc`/`boundaryID` arrays.

### Step 3 — Case driver files (`.par`, `.udf`, `.oudf`, `.usr`)

These four files were written by hand for this case (they are not
auto-generated); each plays a distinct role in the NekRS/Nek5000 hybrid
architecture:

- **`tall_cavity.par`** — top-level case configuration (INI format), read by
  NekRS at startup. Sets:
  - `polynomialOrder = 5` (this case uses order 5, not the order-7 used for
    the single-pebble case in the paper — the two validation cases use
    different discretizations).
  - `startFrom = tall_cavity0.f00200` — this run restarts from a previous
    checkpoint field file rather than starting cold (see Step 6).
  - `[VELOCITY]`/`[TEMPERATURE] boundaryTypeMap` — indexes into
    `boundaryID` set in `usrdat2`: `wall,wall,wall` for velocity
    (boundaryID 1/2/3 all no-slip) and `flux,codedFixedValue,flux` for
    temperature (boundaryID 1=`hot`→flux, 2=`cold`→fixed value,
    3=`topandbot`→flux). This is the direct realization of the paper's
    "one wall at constant heat flux, one wall at constant temperature"
    tall-cavity setup (Section 4.3, Table III: `f_h`, `T_c`).
  - `rho`/`viscosity`/`rhoCp`/`conductivity` given as negative numbers,
    which is NekRS's shorthand for "store the reciprocal" (i.e.
    `viscosity = -6483.0` means `μ = 1/6483.0`), matching Table III's
    dimensionless fluid properties.
  - `regularization = explicit + nModes=1 + scalingCoeff=0.1` — explicit
    high-pass filter stabilization (relevant for this buoyancy-driven,
    low-Reynolds-number case).

- **`tall_cavity.usr`** — the legacy Nek5000 Fortran "backend" source,
  compiled into a shared object and `dlopen`'d by NekRS
  (`src/core/nekInterface/nekInterfaceAdapter.cpp`'s `userchk()` resolves
  `userchk_ptr` via `dlsym`). This is where the paper's radiosity/
  irradiation algorithm (Algorithm 2) is actually invoked:
  - `include 'view_factors/view_factors.f'` pulls in the whole radiation
    module (documented in `view_factors/README.md`) so its subroutines are
    callable from this file.
  - `usrdat2()` sets `cbc`/`boundaryID` from the mesh's physical-surface
    tags (bridges Step 1/2's geometry tags into Nek5000's BC model).
  - `userchk()` is where, at `istep.eq.0`, `nekrs_registerPtr('frad', frad)`
    exposes the CPU radiative-flux array to the C++/GPU side by name (this
    is a real NekRS/Nek5000 bridging call — confirmed in
    `src/core/nekInterface/nekInterface.f:1109`, which forwards to
    `nekf_registerPtr`), and `vf_read_view_factors('vf_tall_cavity')`
    loads the precomputed view factors (Step 5) once. On every call it
    then runs `vf_calculate_radiation_heat_flux(fac0, option=1)` — the
    actual Eq. 1–4 / Algorithm 2 solve — followed by
    `vf_print_radiation_heat_flux(1)` for monitoring the hot wall.
  - `useric()` sets the initial condition (zero velocity, `T=100`, matching
    `T_c` in Table III).

- **`tall_cavity.udf`** — the C++ host-side glue compiled by NekRS's own
  JIT build step (see Step 4) that runs on the GPU-driver side. It:
  - Registers `userf` (Boussinesq buoyancy source, `factor=0.01` matching
    `β=0.01 K⁻¹` in Table III) and `uservp` (constant properties) as NekRS
    callbacks in `UDF_Setup()`.
  - Allocates `nrs->o_usrwrk` (one field wide) as the device-side channel
    for the radiative flux, and initializes it from the CPU `frad` array
    via `nek::ptr<double>("frad")` — the C++-side counterpart of the
    Fortran `nekrs_registerPtr` call above.
  - In `UDF_ExecuteStep()`, decides the coupling cadence: every 100 steps
    (and on any checkpoint step), it pushes the current GPU temperature
    field back to the Nek5000 CPU arrays (`nrs->copyToNek`, confirmed in
    `src/app/nrs/nrs.cpp:1312`), re-invokes the compiled `.usr`'s
    `userchk()` (`nek::userchk()`, confirmed in
    `src/core/nekInterface/nekInterfaceAdapter.cpp:601`) to re-solve
    radiosity/irradiation with the latest temperatures, then re-uploads
    the refreshed `frad` into `o_usrwrk`. This realizes the paper's
    radiosity-update-frequency parameter `N_r` (Section 3.2) as
    `N_r = 100` for this case.

- **`tall_cavity.oudf`** — OKL (GPU) kernels, `#include`d into the
  `__okl__` block of `tall_cavity.udf`. Implements:
  - No-slip velocity BC on all three wall boundary IDs.
  - Fixed-temperature BC (`s=100.0`) on the `cold` wall.
  - `codedFixedGradientScalar`, which reads `bc->usrwrk[...]` — the
    device buffer populated from `frad` — as the local radiative flux
    `gflux`, and applies net flux `0.5 - gflux` on the `hot` wall (imposed
    heat-generation flux `f_h=0.5` from Table III, minus radiative loss)
    and `-gflux` on `topandbot` (radiation-only, otherwise adiabatic).
    This is the exact point where Eq. 4's `q_i` enters the energy
    equation as a Neumann BC.
  - `fillProp`, filling the constant viscosity/conductivity/ρ/ρCp values
    from Table III as device-side property fields.

### Step 4 — First (offline) run to export wall geometry for the view-factor solver

With `vf_export_all_walls(vf_wall_file_name)` **uncommented** in
`usrdat2()` (it is commented out in the checked-in `tall_cavity.usr`, since
it only needs to run once per mesh and the mesh doesn't change between
production solves), the case was built and briefly run. `vf_export_all_walls`
loops every boundary face that is not internal (`cbc.ne.'E  '` and
`cbc.ne.'   '`) — i.e. every one of the five tagged physical surfaces,
including the periodic pair, since closing the radiation enclosure requires
all boundary faces, not just the three literal `'W  '` walls — and writes
each face's 4 corner vertices, area, and outward normal to a text file. This
produced:

```
quads_file_tall_cavity
```

present in this directory (2400 faces on the header line, one quad record
per face) — a leftover artifact of that one-time preprocessing run, as
documented in `README.md`.

### Step 5 — Offline view-factor computation (hemi-cube tool)

`quads_file_tall_cavity` was handed to the external, GPU-accelerated
hemi-cube view-factor calculator described in Section 3.1 of the paper
("hemicubevf", built on OpenGL/GLFW/GLAD; not part of this repository or
directory). For each of the 2400 quads it rasterizes all other quads onto a
virtual hemi-cube (Algorithm 1) and accumulates the view factor
`F_{i→j}` for every visible pair, satisfying the non-negativity, row-sum
(Eq. 9), and reciprocity (Eq. 10) constraints from Section 2. The output was
written to this directory as:

```
vf_tall_cavity
```

in the exact record format `vf_read_view_factors` expects (header = total
wall count, then per-wall `iw ieg ifc nvwalls` followed by `nvwalls` lines of
`jw jeg jfc F_ij`) — confirmed against `view_factors/README.md`'s documented
file format and the sample records read from the file (e.g. `1 1 1 302`
followed by 302 `jw jeg jfc F_ij` lines).

### Step 6 — Re-comment the export call; restore the case for production runs

Once `vf_tall_cavity` existed, `vf_export_all_walls(...)` was commented back
out in `usrdat2()` (its current state in the checked-in `.usr` file) so
normal runs don't redo the geometry export, and the `.usr`'s `userchk()`
was left pointing at `vf_read_view_factors('vf_tall_cavity')` to load the
precomputed factors at `istep.eq.0` on every subsequent run.

### Step 7 — Build

NekRS's JIT build compiles the `.udf`/`.oudf` (OKL/C++) and the `.usr`
(Fortran, linked into the legacy Nek5000 backend library) for this case.
`cmake.log` in this directory records that build step having been run on
Aurora (Intel GPU cluster at Argonne — `IntelLLVM 2025.2.0`, Cray MPICH,
build directory
`/flare/AdvanceFusionFission/gottems/nekrs_runs/tall_cavity_vf/.cache/udf`),
i.e. this specific case checkout was built and run on Aurora rather than
Frontier, even though a Frontier (`nrsqsub_frontier`, OLCF/Slurm+SLURM
`#SBATCH`, AMD/HIP backend) submission script is also present in the
directory — the latter appears to be a copied-over template/reference
script rather than the one actually used for this Aurora run (Aurora's own
submission mechanism, e.g. `nrsqsub_aurora`/PBS `qsub`, is not checked into
this directory).

### Step 8 — Run / restart

`tall_cavity.par`'s `startFrom = tall_cavity0.f00200` and the presence of
`tall_cavity.nek5000` (`filetemplate: tall_cavity%01d.f%05d`,
`firsttimestep: 0`, `numtimesteps: 166`) show this case directory reflects a
run that was **restarted** from a prior checkpoint field file rather than a
cold start — consistent with the long integration time needed for this
buoyancy + slowly-varying radiative-equilibrium problem
(`numSteps = 4000000`, `dt = 1e-2` in the `.par` file). At runtime the
sequence is:

```
istep = 0:
  userchk() [invoked once via NekRS's normal startup path]
    -> nekrs_registerPtr('frad', frad)
    -> vf_read_view_factors('vf_tall_cavity')     [Step 5's file]
  UDF_Setup()
    -> nrs->o_usrwrk allocated, seeded from frad (still all zero here)

every GPU step:
  tall_cavity.oudf: codedFixedGradientScalar reads bc->usrwrk as gflux

every 100 steps (or on checkpoint):
  UDF_ExecuteStep()
    -> copyToNek(time, tstep)          push GPU T -> CPU
    -> nek::userchk()
         -> vf_calculate_radiation_heat_flux(fac0, option=1)   [Eq.1-4 / Algorithm 2]
         -> vf_print_radiation_heat_flux(1)                    [monitor hot wall]
    -> [every-100-steps branch] re-upload frad -> o_usrwrk
```

### Step 9 — Post-processing / comparison against Section 4.3's Figure 5

Not present in this directory (no post-processing scripts were checked in
here) — the paper's Figure 5 compares NekRS's convective/radiative Nusselt
numbers, decomposed as `Nu_total = Nu_convection + Nu_radiation`, against a
1D model and MOOSE across a range of wall emissivities `ε`. Reproducing that
sweep requires re-running this case at multiple `vf_eps` values (currently
hardcoded to `vf_eps = 1.0` in `tall_cavity.usr`'s `userchk()`) and
extracting/decomposing the Nusselt number from the resulting fields — that
analysis step lives outside this case directory.

---

## 3. Summary of files present vs. their role

| File | Produced in step | Purpose |
| --- | --- | --- |
| `tall_cavity.geo` | 1 | Gmsh geometry/mesh script |
| `tall_cavity.msh` | 1 | Gmsh mesh output |
| `tall_cavity.re2`, `tall_cavity.co2` | 2 | Nek5000-format mesh (via `gmsh2nek`) |
| `tall_cavity.par` | 3 | NekRS case configuration |
| `tall_cavity.usr` | 3 | Legacy Nek5000 Fortran backend (radiation driver) |
| `tall_cavity.udf` | 3 | NekRS host-side C++ (GPU↔CPU radiation coupling) |
| `tall_cavity.oudf` | 3 | OKL GPU kernels (BCs, properties) |
| `view_factors/view_factors.f`, `view_factors/VIEW_FACTORS` | 3 (included) | Radiosity/irradiation solver module |
| `quads_file_tall_cavity` | 4 | Exported wall geometry (offline, one-time) |
| `vf_tall_cavity` | 5 | Precomputed view factors (offline, via hemi-cube tool) |
| `cmake.log` | 7 | Build log from the Aurora build of this case |
| `tall_cavity0.f00200` (referenced, not present) | 8 | Restart checkpoint field file |
| `tall_cavity.nek5000` | 8 | Field-file naming/series metadata for restart/post-processing |
| `nrsqsub_frontier` | — | Reference/template Frontier submission script (not used for this Aurora build, per `cmake.log`) |
| `README.md`, `view_factors/README.md` | — | Case- and module-level documentation (pre-existing in this directory) |

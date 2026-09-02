# tall_cavity — NekRS case: buoyancy-driven convection + view-factor radiation

This case couples NekRS's GPU flow/scalar solver with the CPU-side, legacy-Nek5000
view-factor radiation model implemented in
[`view_factors/view_factors.f`](view_factors/README.md). Three files define the
case-specific physics:

| File | Runs on | Role |
| --- | --- | --- |
| [`tall_cavity.usr`](tall_cavity.usr) | CPU (legacy Nek5000 "backend") | boundary/IC setup, drives the view-factor pipeline each time it's invoked |
| [`tall_cavity.udf`](tall_cavity.udf) | Host (NekRS driver) | registers GPU callbacks, bridges CPU radiation results into GPU device memory |
| [`tall_cavity.oudf`](tall_cavity.oudf) | GPU (OKL kernels) | per-step boundary conditions and material properties, consumes the radiative flux |

The geometry ([`tall_cavity.geo`](tall_cavity.geo)) is a thin extruded 2D
rectangle (`H` x `H/2` x `L`, `L = H/20`) with 5 tagged surfaces:

| Physical Surface | tag | Role in this case |
| --- | --- | --- |
| `hot` | 1 | heated wall |
| `cold` | 2 | fixed-temperature wall |
| `topandbot` | 3 | adiabatic-ish (radiating only) wall |
| `p1` / `p2` | 4 / 5 | periodic pair (the thin extrusion direction) |

---

## `tall_cavity.usr`

```fortran
include 'view_factors/view_factors.f'
```
pulls in every routine documented in
[`view_factors/README.md`](view_factors/README.md) so they're available from
this case's `userchk`/`usrdat*` routines.

### `usrdat2` — boundary condition setup
Maps the mesh's numeric boundary tag (`bc_flag = bc(5,ifc,iel,1)`, i.e. the
`Physical Surface` id from the `.geo` file) to Nek's `boundaryID`/`cbc`
arrays:

- `bc_flag=1` (`hot`), `2` (`cold`), `3` (`topandbot`) → `cbc = 'W  '`
  (wall), `boundaryID = 1/2/3` respectively.
- `bc_flag=4,5` (`p1`/`p2`) → `cbc = 'P  '` (periodic), `boundaryID = 0`.

`boundaryID` here is what `[PRESSURE]/[VELOCITY]/[TEMPERATURE]`'s
`boundaryTypeMap` in [`tall_cavity.par`](tall_cavity.par) index into
(`boundaryTypeMap = flux, codedFixedValue, flux` for TEMPERATURE means
boundaryID 1→flux, 2→codedFixedValue, 3→flux — see the `.oudf` section
below).

The commented-out block at the end,
```fortran
!      call vf_export_all_walls(vf_wall_file_name)
```
is the one-time preprocessing call from
[`view_factors/README.md`](view_factors/README.md#vf_export_all_wallswall_quads_file_name):
it's run once (with this line uncommented) to dump wall geometry to
`quads_file_tall_cavity` — present in this directory as a leftover artifact
of that preprocessing run — for the external view-factor solver. It's
disabled for normal runs since the mesh/geometry doesn't change between
solves.

Note `cbc.eq.'W  '` is only set for surfaces 1–3, matching the
"real wall vs. open-boundary" distinction described in the view-factor
module: for this case *every* non-periodic surface happens to be a `'W  '`
wall, so `indicator_real_wall_ef` is 1.0 everywhere the radiation model
looks — periodic faces are internal (`cbc='P  '`, effectively treated like
`'E  '`/interior) and excluded from the radiation enclosure, so no
open/leaking boundary exists in this particular case (`open_face_option=1`
is used but has no non-wall faces to exercise).

### `useric` — initial conditions
Zero velocity, uniform initial temperature `temp = 100.0`.

### `userchk` — driven every time it's called (see `UDF_ExecuteStep` below)
This is where the view-factor pipeline actually runs at simulation time:

1. `istep.eq.0`: zero `frad` and call
   `nekrs_registerPtr('frad', frad)` so the GPU side (`tall_cavity.udf`) can
   later read this CPU array by name via `nek::ptr("frad")`.
2. `istep.eq.0`: `vf_read_view_factors('vf_tall_cavity')` — loads the
   precomputed view factors from the `vf_tall_cavity` file in this directory
   (produced offline from the `quads_file_tall_cavity` export). The
   partitioned-file alternative
   (`vf_export_partitioned_vf_files`/`vf_read_partitioned_vf_files`, see the
   view-factors doc) is present but commented out.
3. **Every** call (not just `istep.eq.0`): sets up the non-dimensionalization
   factor
   ```fortran
   fac0 = vf_eps*vf_sigma*T0**4.0/f0     ! T0=U0=rho=Cp=1.0 here
   ```
   and calls
   `vf_calculate_radiation_heat_flux(fac0, option=1)` — runs the radiosity/
   irradiation iteration (`open_face_option=1`, ideal diffuse reflection)
   described in the view-factors doc, updating `frad` at every GLL point.
4. `vf_print_radiation_heat_flux(1)` — prints the average radiative flux and
   temperature on sideset `1` (`hot`) each time `userchk` runs, for
   monitoring.

Because `userchk` is CPU-only (legacy Nek5000 arrays: `t`, `cbc`, `bm1`,
etc.), it is **not called every GPU timestep** — see `UDF_ExecuteStep` in
`tall_cavity.udf`, which decides when to invoke it.

---

## `tall_cavity.udf`

Host-side glue between NekRS's GPU solver and the CPU radiation model.

### `uservp(time)` — material properties
Calls the OKL kernel `fillProp` (defined in `tall_cavity.oudf`) once per
call to fill constant viscosity/conductivity/density/rhoCp into
`nrs->o_prop`/`cds->o_prop` (values are hardcoded in the kernel, not
temperature-dependent despite `cds->o_S` being passed in).

### `userf(time)` — Boussinesq-like buoyancy source
```cpp
platform->linAlg->axpby(mesh->Nlocal, factor, o_T, 0.0, o_FUy);
```
adds `factor * T` (`factor = 0.01`) into the **y-momentum** source
(`o_FUy`). Since the cavity is "tall" (extends along `y`, per the `.geo`
file: `H` in `y`, `H/2` in `x`), this is gravity acting in `-y`/`+y` via a
temperature-proportional forcing — the natural-convection driver for this
cavity flow, paired with the radiative heating/cooling applied at the walls.

### `UDF_Setup()` — one-time setup
- Registers `userf`/`uservp` as the velocity-source/property callbacks.
- Allocates `nrs->o_usrwrk` with room for **one** field
  (`1*nrs->fieldOffset`) — this is the device-side buffer the `.oudf`
  boundary conditions will read from.
- Reads the CPU `frad` array via `nek::ptr<double>("frad")` (the pointer
  registered in `userchk`'s `nekrs_registerPtr('frad', frad)` call) and
  copies it into `o_usrwrk` slot 0. At this point `frad` is still all zeros
  (only zeroed, not yet solved — the first real radiosity solve happens in
  `userchk`, invoked later from `UDF_ExecuteStep`).

### `UDF_ExecuteStep(time, tstep)` — coupling cadence
This function decides **how often** the CPU radiation model
(`view_factors.f`'s `vf_calculate_radiation_heat_flux`, via `userchk`) is
re-solved and re-uploaded to the GPU:

- **On every checkpoint step** (`nrs->checkpointStep`): `copyToNek` (push
  current GPU temperature field to the CPU/Nek arrays) then `nek::userchk()`
  — this re-solves radiosity/irradiation with the latest temperature field,
  but note the resulting `frad` is **not** copied back to `o_usrwrk` in this
  branch (only in the branch below), so a checkpoint alone doesn't refresh
  the GPU-visible flux.
- **Every 100 steps** (`tstep % 100 == 0`, `tstep > 0`): same
  `copyToNek` + `nek::userchk()`, **and then** re-reads `frad` from the CPU
  side and re-uploads it into `o_usrwrk` slot 0, so the GPU boundary
  condition (`.oudf`) sees the updated radiative flux going forward.

In other words: the expensive, iterative view-factor radiosity solve is
**not** run every GPU timestep — it's run (and only pushed to the GPU) every
100 steps, treating radiative equilibrium as slowly varying compared to the
convective/diffusive timestep. Between refreshes, the `.oudf` boundary
condition keeps using the last-uploaded `frad` values.

---

## `tall_cavity.oudf`

GPU-side (OKL) boundary conditions and property kernel, `#include`d into
`tall_cavity.udf`'s `__okl__` block.

### `codedFixedValueVelocity` — no-slip
`u = v = w = 0` on every wall (`boundaryTypeMap = wall, wall, wall` in
`tall_cavity.par` applies this to boundaryID 1/2/3, i.e. `hot`, `cold`,
`topandbot`).

### `codedFixedValueScalar` — fixed cold-wall temperature
`s = 100.0`. Applied where `[TEMPERATURE] boundaryTypeMap`'s second entry
(`codedFixedValue`) matches, i.e. `boundaryID = 2` (`cold`). Note this is
the *same* temperature value as the initial condition in `useric`
(`temp = 100.0`) — the cold wall is really "reference temperature",
with the `hot` wall driven by an imposed flux instead (below).

### `codedFixedGradientScalar` — radiative + imposed flux boundary
Applied where `boundaryTypeMap`'s `flux` entries match — `boundaryID = 1`
(`hot`) and `boundaryID = 3` (`topandbot`):
```cpp
gflux = bc->usrwrk[bc->idM + 0 * bc->fieldOffset];  // radiation heat flux from view factor model
bc->flux = 0.0 - gflux;
if (bc->id == 1) bc->flux = 0.5 - gflux;
```
`bc->usrwrk` is exactly the `o_usrwrk` buffer populated in
`tall_cavity.udf` from the CPU `frad` array (`view_factors.f`'s
`vf_calculate_radiation_heat_flux` output) — this is the point where the
view-factor radiation model's result actually enters the flow/scalar
solve. `gflux` is that local face's radiative heat flux; subtracting it
from the imposed flux means **radiative loss reduces the net flux into the
domain** at every `flux`-type wall.
- `topandbot` (`bc->id != 1`): net flux is purely `-gflux` — this wall is
  adiabatic except for radiation exchange.
- `hot` (`bc->id == 1`): net flux is `0.5 - gflux` — an imposed heat
  generation flux of `0.5` on top of the radiative exchange, making this
  the case's actual heat source.

### `fillProp` kernel — constant material properties
```cpp
UPROP[... + 0*uOffset] = 1e-5;    // viscosity
SPROP[... + 0*sOffset] = 0.01;    // conductivity
UPROP[... + 1*uOffset] = 1.0;     // rho
SPROP[... + 1*sOffset] = 1000.0;  // rhoCp
```
Despite computing `rcpTemp = 1.0 / TEMP[id]`, this value is never used —
properties are uniform constants, not actually temperature-dependent
(matches `[VELOCITY]`/`[TEMPERATURE]` `rho`/`viscosity`/`rhoCp`/
`conductivity` in `tall_cavity.par`, which are set to the same values via
negative-viscosity/conductivity shorthand for `1/value`).

---

## End-to-end pipeline

```text
[offline, once per mesh]
  tall_cavity.usr: usrdat2() [with vf_export_all_walls uncommented]
      -> quads_file_tall_cavity  (wall geometry: vertices/area/normal)
      -> [external view-factor solver, not in this repo]
      -> vf_tall_cavity          (view factor file)

[at run start, istep.eq.0]
  tall_cavity.usr: userchk()
      -> nekrs_registerPtr('frad', frad)        makes CPU `frad` visible to UDF
      -> vf_read_view_factors('vf_tall_cavity')  loads view factors (once)

[every GPU timestep]
  tall_cavity.oudf: codedFixedGradientScalar
      -> reads bc->usrwrk (= last-uploaded `frad`) as boundary flux

[every 100 steps, or on checkpoint]
  tall_cavity.udf: UDF_ExecuteStep()
      -> copyToNek(time, tstep)                  push GPU temperature -> CPU
      -> nek::userchk()
           -> vf_calculate_radiation_heat_flux()  re-solve Jr/Gr, recompute frad
           -> vf_print_radiation_heat_flux(1)     report hot-wall flux/T
      -> [every-100-steps branch only] re-upload frad -> o_usrwrk
```

See [`view_factors/README.md`](view_factors/README.md) for the full
subroutine-level documentation of the radiosity solve, view-factor file
formats, and the `VIEW_FACTORS` common-block data model referenced above
(`frad`, `egf_to_iwall`, `view_factor_value`, `vf_eps`, `vf_sigma`, etc.).

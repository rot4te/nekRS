# view_factors.f

Fortran module implementing surface-to-surface radiation heat transfer for
[Nek5000](https://nek5000.mcs.anl.gov/)/NekRS cases via the **view factor
method**. It is `include`d directly into a case's `.usr` file (see
`tall_cavity.usr:1`) and depends on the companion common-block header
`view_factors/VIEW_FACTORS`.

The module is split into two phases that normally run as separate
executables/steps:

1. **Preprocessing** (`vf_export_all_walls`) — run once per mesh to dump every
   boundary face (quad vertices, area, outward normal) to a text file. This
   file is fed to an external view-factor solver (not part of this file) which
   computes face-to-face view factors and produces a "view factor file".
2. **Runtime** (`vf_read_view_factors`, `vf_calculate_radiation_heat_flux`,
   `vf_print_radiation_heat_flux`) — called from `userchk()` every step (or
   once at `istep.eq.0` for the read) to load view factors and compute the
   radiative heat flux applied to wall boundaries.

A partitioned-file variant (`vf_export_partitioned_vf_files` /
`vf_read_partitioned_vf_files`) exists to avoid an MPI all-gather of the
(potentially huge) view factor file — each rank writes/reads only the faces it
owns.

## Data model (`view_factors/VIEW_FACTORS`)

Common block `/VIEWFACTORS/` holds all persistent state:

| Variable | Shape | Meaning |
| --- | --- | --- |
| `farea(nfe,lelt)` | per element/local-face | face area |
| `fnormal(3,nfe,lelt)` | per element/local-face | outward unit normal (area-averaged) |
| `fvertice1..4(3,nfe,lelt)` | per element/local-face | 4 corner vertices of the face quad |
| `fnrmx/y/z(lx1,ly1,lz1,lelt)` | GLL-point level | normal components at every grid point of a face, used as scratch to compute the area-averaged `fnormal` via `surface_int` |
| `egf_to_iwall(6,lelg)` | global element × local face | maps `(ifc, global element)` → global wall index `iwall` |
| `iwall_to_eg(max_walls)`, `iwall_to_f(max_walls)` | indexed by `iwall` | inverse mapping back to `(global element, local face)` |
| `n_visible_walls(6,lelt)` | per element/local-face | number of other faces with a nonzero view factor to this face |
| `jindex_to_jwall(max_visible_walls,6,lelt)` | ditto | global wall index of the `jindex`-th visible face |
| `view_factor_value(max_visible_walls,6,lelt)` | ditto | the view factor `F_{i->j}` values |
| `temp_ef(max_walls)` | per wall | element-face-averaged temperature |
| `indicator_real_wall_ef(max_walls)` | per wall | 1.0 if the boundary condition is `'W  '` (an actual solid wall), 0.0 for other open boundaries (inlet/outlet) that are only included to close the radiation enclosure |
| `Jr(max_walls)` | per wall | radiosity |
| `Gr(max_walls)` | per wall | irradiation |
| `F_sky(max_walls)` | per wall | fraction of a wall's view factor "leaking" to non-wall open boundaries, used as the "sky"/black-hole term in `open_face_option=2` |
| `frad_ef(max_walls)` | per wall | net radiative heat flux at the element-face level |
| `frad(lx1,ly1,lz1,lelt)` | GLL-point level | radiative heat flux mapped back onto the volumetric grid (only nonzero on `'W  '` faces) |
| `vf_crf_firstCalled` | scalar | flag: on first call, zero `Jr`/`Gr` before iterating |
| `vf_eps`, `vf_sigma` | scalar | surface emissivity and Stefan–Boltzmann constant, set by the case's `.usr` file |

Sizing parameters: `max_walls` (4,000,000) is the max total number of
boundary faces across the whole mesh; `max_visible_walls` (10,000) is the max
number of faces visible to any single face; `nfe = 2*ldim` is the number of
local faces per element (6 in 3D).

**Important convention**: throughout this file, "wall" faces are *all*
non-internal, non-empty boundary faces (`cbc.ne.'E  '` and `cbc.ne.'   '`),
not just faces literally tagged `'W  '`. Inlets/outlets are included so the
radiation enclosure is topologically closed ("cloture"), which is required
for view factors from any given face to sum to 1. The `'W  '` tag is used
later to decide which faces are *physical* solid walls (get radiative flux
applied, contribute `indicator_real_wall_ef=1`) versus open boundaries that
only participate in closing the enclosure.

---

## Subroutines

### `vf_export_all_walls(wall_quads_file_name)`

Preprocessing step, run once (e.g. from `usrdat2`, currently commented out in
`tall_cavity.usr:114`) to export every boundary face's geometry for the
external view-factor solver.

1. Loops all local elements/faces; any face where
   `cbc.ne.'E  '` and `cbc.ne.'   '` is counted as a wall face.
2. For each such face:
   - Computes the per-grid-point outward normal via `getSnormal` and stores
     its negation in `fnrmx/y/z` (scratch arrays at GLL-point resolution).
   - Extracts the 4 corner vertices of the face into `fvertice1..4` by
     detecting which logical index (`i`, `j`, or `k`) is degenerate
     (`i0.eq.i1`, etc.) for that face, then reading `xm1/ym1/zm1` at the four
     corner combinations of the other two indices.
   - Computes the area-averaged normal (`fnormal`) and face area (`farea`) via
     `surface_int` on the negated normal components.
3. Counts the global number of wall faces (`iglsum`) and (on rank 0) opens
   `wall_quads_file_name` and writes the total count as a header.
4. Gathers geometry across all MPI ranks in blocks of `lblock=500` global
   elements at a time (via the `scrns_vf` scratch common block and `gop`/`igop`
   reductions), and rank 0 writes, for each wall face: a running index, the
   global element number, the local face index, the four vertices, and the
   normal + area.

This is the only geometry-producing routine; its output feeds an external
view-factor computation tool that is not part of this repository.

### `vf_read_view_factors(view_factor_file)`

Runtime step, called once at `istep.eq.0` (see `tall_cavity.usr:29`) to load
a precomputed view-factor file into memory on every rank (each rank opens and
reads the *entire* file, redundantly).

Format read: a header line with total wall count `nwalls`, followed for each
wall `iw` of a line `iw, ieg, ifc, nvwalls` (global wall index, global element,
local face, number of visible faces), followed by `nvwalls` lines of
`jw, jeg, jfc, view_factor` (only `jw` and `view_factor` are kept).

- Populates the *global* mapping `egf_to_iwall`, `iwall_to_eg`, `iwall_to_f`
  for every wall face regardless of rank ownership.
- Populates the *local* per-rank arrays `n_visible_walls`,
  `view_factor_value`, `jindex_to_jwall` only for faces belonging to elements
  owned by this rank (`gllnid(ieg).eq.nid`); for faces owned by other ranks it
  still reads (and discards) the same number of lines to keep the file
  position in sync.
- Emits a warning (rank 0 only) if `nwalls`/`nvwalls` exceed the
  `max_walls`/`max_visible_walls` array bounds, or if the wall index in the
  file (`iw`) doesn't match the expected sequential index (`iwall`) — a sanity
  check that the file is well-formed.
- Sets `vf_crf_firstCalled = 1` so the next call to
  `vf_calculate_radiation_heat_flux` initializes `Jr`/`Gr` to zero.

### `vf_export_partitioned_vf_files(vf_folder)`

Alternative export path: instead of one shared view-factor file
(`vf_read_view_factors`'s input) that every rank reads in full, this splits
the already-loaded view factors (in `egf_to_iwall`, `n_visible_walls`,
`view_factor_value`, etc.) into one file per MPI rank, `vf_folder/vfp<nid>`,
each rank writing only the wall faces belonging to its own local elements
(`iel = 1,lelt`).

- Rank 0 deletes/recreates `vf_folder` via `execute_command_line('rm -rf ...'
  )` / `mkdir`, followed by `nekgsync()` to barrier before every rank opens
  its own file.
- Same per-face record format as `vf_export_all_walls`'s file but keyed by
  `iwall/ieg/ifc/nvw` header and `jwall/jg/jf/fij` per-visible-face lines
  (fixed-format `71`/`72` instead of list-directed I/O).

Currently unused/commented out in `tall_cavity.usr` — provided as a
performance option for large cases where re-reading the full view-factor file
on every rank (as `vf_read_view_factors` does) is too expensive.

### `vf_read_partitioned_vf_files(vf_folder)`

Counterpart to `vf_export_partitioned_vf_files`: each rank reads only its own
`vf_folder/vfp<nid>` file, populating the same local arrays, then does an
`igop` (integer global sum reduction) on `egf_to_iwall`, `iwall_to_eg`,
`iwall_to_f` to reconstruct the global mappings across ranks (each rank
contributes zeros everywhere except the entries for its own elements, so the
sum recovers the full global map). Also sets `vf_crf_firstCalled = 1`.

### `vf_calculate_radiation_heat_flux(fac0, open_face_option)`

The core per-step (or per-call) radiosity solve. `fac0` is a
non-dimensionalization factor (`eps*sigma*T0**4/f0`, computed by the caller,
e.g. `tall_cavity.usr:48`) and `open_face_option` selects how non-wall
boundaries (inlet/outlet, included only for enclosure closure) are treated:

- **`open_face_option = 1`** — *ideal diffuse reflection*: open faces behave
  like walls with `indicator_real_wall_ef` gating the emitted term, i.e. they
  reflect/re-radiate energy back into the enclosure rather than absorbing it.
- **`open_face_option = 2`** — *black-hole / DOM-matching*: radiation leaving
  through an open face never returns (`Jr_sky = 0.0` by default, i.e. treated
  as leaving to a zero-temperature sky), matching a discrete-ordinates-method
  reference. Only faces literally tagged `'W  '` participate in the
  radiosity/irradiation iteration in this branch; the leakage to open
  boundaries is precomputed once as `F_sky(iwall)` (1 minus the sum of view
  factors to real walls).

Steps:

1. **Map temperature to face level**: for every boundary face (wall or open,
   for enclosure closure), average the volumetric temperature field `t` over
   the face via `surface_int` into `temp_ef(iwall)`, and set
   `indicator_real_wall_ef(iwall)` to 1.0 only if `cbc.eq.'W  '`. Reduced
   across ranks with `gop`; min/max printed for diagnostics.
2. **View factor sum sanity check**: sums `view_factor_value` over all visible
   faces per face (`Fsum`) — should be 1 for a properly closed, normalized
   enclosure. (Diagnostic printing here is currently commented out.)
3. **Radiosity/irradiation iteration**: `eps = vf_eps`. If this is the first
   call since view factors were (re)loaded (`vf_crf_firstCalled.eq.1`),
   zero `Jr`/`Gr`. Then iterate `niter = 5` times:
   - `open_face_option.eq.1`: for every boundary face,
     `Jr(iwall) = fac1*fac0*T^4 + (1 - eps*fac1)*Gr(iwall)` where
     `fac1 = indicator_real_wall_ef(iwall)` (so open faces have
     `Jr = Gr`, i.e. perfect reflection, while walls emit + reflect).
     Then `Gr(iwall) = sum_j F_{iwall->jwall} * Jr(jwall)` over all visible
     faces. Both reduced with `gop` each iteration; min/max of `Jr`/`Gr`
     printed for diagnostics.
   - `open_face_option.eq.2`: same radiosity update but restricted to
     `cbc.eq.'W  '` faces only, without the `fac1` gate on emission
     (`Jr = fac0*T^4 + (1-eps)*Gr`), and `Gr` additionally receives the
     `F_sky(iwall)*Jr_sky` leakage term to account for the black open
     boundary.
4. **Net radiative flux per face**: `frad_ef(iwall) = Jr(iwall) - Gr(iwall)`
   (positive = leaving the wall into the fluid) for every boundary face,
   reduced with `gop`.
5. **Map back to grid level**: for faces tagged `'W  '` only, broadcast the
   scalar `frad_ef(iwall)` to every GLL point on that face (via `facind` to
   get the face's index ranges) into `frad(i,j,k,iel)`. Then apply the
   standard Nek mass-matrix weight/DSSUM/inverse-mass sequence
   (`col2(...,bm1)` → `dssum` → `col2(...,binvm1)`) to average `frad`
   consistently across the corners/edges shared by multiple elements. Min/max
   printed for diagnostics.

The result, `frad`, is the volumetric (GLL-point) field of radiative heat
flux that the case's `.udf`/`.oudf` couples into the energy equation (see
`nekrs_registerPtr('frad', frad)` in `tall_cavity.usr:24`).

### `vf_print_radiation_heat_flux(ss)`

Diagnostic/reporting routine: integrates `frad` and temperature `t` over all
faces belonging to boundary sideset `ss` (matched via `boundaryID(iside,e)`,
distinct from the `cbc` character tag used elsewhere), area-averages them,
and prints the sideset's average radiative heat flux and average temperature
on rank 0. Useful for checking global energy balance / monitoring a
particular wall's heat flux during a run (called each step in
`tall_cavity.usr:53` for sideset 1).

---

## Case-level documentation

For how this module is wired into an actual running case (boundary setup,
GPU/CPU coupling cadence, and where `frad` is consumed as a boundary
condition), see [`../README.md`](../README.md), which documents
`tall_cavity.usr`, `tall_cavity.udf`, and `tall_cavity.oudf`.

## Typical call sequence (see `tall_cavity.usr`)

```text
usrdat2():
    set cbc/boundaryID from mesh tags
    [one-time, normally run separately] vf_export_all_walls(...)
        -> external view-factor solver -> vf_tall_cavity file

userchk(), istep.eq.0:
    vf_read_view_factors('vf_tall_cavity')
    [optional] vf_export_partitioned_vf_files(...) / vf_read_partitioned_vf_files(...)

userchk(), every step:
    vf_calculate_radiation_heat_flux(fac0, option)
    vf_print_radiation_heat_flux(sideset_id)
```

## Notes / gotchas

- `vf_read_view_factors` has every rank read the whole file redundantly —
  fine for small/medium enclosures but a scalability bottleneck for very
  large wall counts, hence the partitioned-file alternative.
- `nwalls`/`nvwalls` bound checks against `max_walls`/`max_visible_walls`
  only *print a warning*, they do not abort — if the file exceeds either
  parameter, arrays will silently overflow. Increase the `parameter`
  values in `VIEW_FACTORS` if the warning appears.
- `open_face_option` must be consistent with how the view factors were
  generated/normalized upstream; the two branches use different subsets of
  faces (`.ne.'E  '/'   '` vs. exactly `'W  '`) and different sky/reflection
  physics.

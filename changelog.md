# Changelog

## 2026-09-01

- Added `examples/concentric_spheres_vf/run_teton.sh`: Slurm batch script to
  run the `concentric_spheres_vf` view-factor case on the INL Teton
  cluster, based on the 2026-06-17 INL Slurm training deck and the
  existing `scripts/nrsmpi` / `tall_cavity_vf_aurora/nrsqsub_frontier`
  launch conventions in this repo.
- Rewrote `examples/concentric_spheres_vf/run_teton.sh` to follow the
  `nrsqsub_frontier`-style generator pattern (source `nrsqsub_utils`, write
  an sbatch script to `s.bin`, submit it), matching the user-supplied
  `nrsqsub_teton-2.sh` reference script verbatim (192 cores/node, CPU
  `serial` backend, `PrgEnv-gnu`/`craype-x86-turin`/`cray-mpich` modules,
  `srun --mpi=cray_shasta`) and defaulting `NEKRS_HOME`/`PATH` from the
  user-supplied `teton_nekrs_env.txt` (`$HOME/.local/nekrs_teton`).
- Added `examples/parallel_plates_vf/run_teton.sh` and
  `examples/single_pebble_vf/run_teton.sh`, copied from the current state
  of `examples/concentric_spheres_vf/run_teton.sh` (same generator
  mechanics, unchanged) with per-case header comments: `parallel_plates`
  mirrors `concentric_spheres_vf`'s resource profile (4400 elements,
  polynomialOrder=5, numSteps=5, trivial radiosity-only verification).
  `single_pebble` is flagged as heavier and unbenchmarked (2376 elements,
  polynomialOrder=7, numSteps=2000, forced flow + conjugate radiative-flux
  BCs) since its procedure.md states it has not been run yet -- its
  `--time` default is called out in-script as an unverified guess.
  Did NOT add a script for `examples/radiationPlates/` -- see chat.
- Fixed `examples/concentric_spheres_vf/concentric_spheres.par`: this
  checked-out nekRS source (v26.0.1, sha 96b3cf9e) does not accept the
  `[VELOCITY]`/`[PRESSURE]`/`[TEMPERATURE]` section names the case was
  written with -- confirmed via `src/platform/par/par.cpp`'s
  `validSections` list and this repo's own `examples/channel/channel.par`
  / `examples/conj_ht/conj_ht.par`. Renamed to `[FLUID VELOCITY]` /
  `[FLUID PRESSURE]`, added `scalars = temperature` to `[GENERAL]`, and
  renamed `[TEMPERATURE]` to `[SCALAR TEMPERATURE]` (temperature is now a
  named scalar, not a dedicated section).
- Applied the identical fix to `examples/parallel_plates_vf/parallel_plates.par`
  and `examples/single_pebble_vf/single_pebble.par` (same
  `[VELOCITY]`/`[PRESSURE]`/`[TEMPERATURE]` -> `[FLUID VELOCITY]`/
  `[FLUID PRESSURE]`/`scalars = temperature` + `[SCALAR TEMPERATURE]`
  rename). Not yet run/verified on Teton -- `concentric_spheres_vf` is the
  one actively being tested; if its `.usr`'s legacy Nek5000 `t(1,1,1,1,1)`
  access still works correctly with the renamed scalar, the same should
  hold for `single_pebble.usr`'s identical pattern, but this hasn't been
  confirmed.
- Removed `maxIterations = 1000` from all three `.par` files' `[FLUID
  PRESSURE]`/`[SCALAR TEMPERATURE]` sections: this checked-out nekRS
  source's `src/platform/par/par.cpp` has no `maxIterations` key in any
  section's valid-key list (`commonKeys`/`ellipticKeys`/`scalarKeys`) --
  it's not renamed, just gone as a per-section knob in this version.
  Confirmed by the run failing with `unknown key: fluid
  pressure::maxiterations` / `unknown key: scalar
  temperature::maxiterations` on `concentric_spheres_vf`.
- Fixed the `[FLUID PRESSURE]` `solver = PFGMRES + nVector = 30` line in
  all three `.par` files: `PFGMRES` isn't a recognized solver token in
  this source (`src/platform/par/parseLinearSolve.cpp`'s validValues list
  has `pgmres`/`flexible`/`gmres` as separate tokens, no combined
  `pfgmres`), and this source's own default pressure solver string
  (`GMRES+FLEXIBLE+NVECTOR=15`, set in the same file) confirms the correct
  replacement is `GMRES+flexible`. Also restored the max-iterations
  behavior removed earlier: `nvector=`/`maxiter=` are extracted via
  `std::regex("nvector=([0-9]+)")` / `std::regex("maxiter=([0-9]+)")`
  (`linearSolverFactory.cpp`, `elliptic.cpp`) with no tolerance for spaces
  around `=`, so `maxIterations` is now set by embedding `+maxIter=1000`
  directly in the solver string, written with no spaces throughout:
  `solver = GMRES+flexible+nVector=30+maxIter=1000`. Checked
  `smootherType = FourthChebyshev+jac` and `initialGuess =
  projection+nVector=30` against the same validation code
  (`checkValidity()` in `par.cpp` is a prefix match, `entry.find(v)==0`,
  so `fourthchebyshev` matching valid-token `fourthcheby` is fine) --
  those two lines needed no changes.
- Fixed a udf.cpp compile error in all three `.udf` files (`nrs.hpp`'s
  `nrs_t` class has no `mesh`/`o_usrwrk` members in this checked-out
  source): `nrs->mesh` -> `nrs->meshV`, `nrs->o_usrwrk` ->
  `platform->app->bc->o_usrwrk` (moved to a `bc` object off
  `platform->app`, per `src/core/bdry/bdryBase.hpp`; same
  `.resize()`/`.copyFrom()` API, confirmed against `src/platform/deviceMemory.hpp`).
  Confirmed the replacement names against this repo's real, currently-built
  examples (`examples/gabls1/gabls.udf`, `examples/turbPipe/turbPipe.udf`),
  not just the compiler's suggestion -- both already use
  `platform->app->bc->o_usrwrk` and `nrs->meshV`. Also checked all three
  cases' `.oudf` files: `single_pebble.oudf`'s `codedFixedGradientScalar`
  already reads `bc->usrwrk[...]` in the current convention, so no changes
  were needed on the OKL/device side, only the host-side `.udf`.
  (Side note: `examples/concentricSpheres/` and `examples/radiationPlates/`
  in this repo are stale build-cache-only directories from an unrelated
  nekRS checkout at `/Users/coxea3/NekCpp/nekRS/` -- not usable as
  in-repo references, despite `examples/concentricSpheres/`'s cached
  `udf.cpp` initially looking like a promising working template.)
- Rewrote all three `.oudf` files: nekRS aborted at runtime with "Cannot
  find required okl function udfDirichlet!" -- traced to
  `src/core/udf/udf.cpp`, which regex-searches the compiled `.oudf` for a
  literal `void udfDirichlet` function; the old per-field function names
  (`codedFixedValueScalar`, `codedFixedValueVelocity`,
  `codedFixedGradientScalar`) don't match, and turned out to be used by
  none of this repo's other examples either (corrects my note two entries
  up -- `single_pebble.oudf`'s `codedFixedGradientScalar` did NOT already
  match the current convention). Converted all three to the current
  convention (one `udfDirichlet`/`udfNeumann` function per case, branching
  with `isField("scalar temperature")` / `isField("fluid velocity")`),
  confirmed against `examples/gabls1/gabls.udf` and
  `examples/turbPipe/turbPipe.udf`. While rewriting, checked the actual
  `bcData` struct (`src/app/nrs/bdry/bcData.h`) instead of assuming field
  names carried over from the old convention -- `bc->s` -> `bc->sScalar`,
  `bc->flux` -> `bc->fluxScalar`, `bc->u`/`bc->v`/`bc->w` ->
  `bc->uxFluid`/`bc->uyFluid`/`bc->uzFluid`, and `bc->idM` (not a real
  struct member) -> `bc->idxVol`, the last one only affecting
  `single_pebble.oudf`'s `usrwrk` indexing.
- Fixed all three `run_teton.sh` scripts: `concentric_spheres_vf` segfaulted
  at runtime (`srun: error: ... Segmentation fault`, ~100 of 192 ranks,
  right after `VF: Reading view factor file: nwalls= 588` in the legacy
  Nek5000 `view_factors.f` module). Root cause: 1470 elements / 192 MPI
  ranks = 7.66 avg, but the build log's `building nekInterface for lx1=6,
  lelt=6, lelg=1470` shows nekRS auto-sized the legacy per-rank element
  buffer (`lelt`) to 6 -- smaller than even the average load, let alone an
  imbalanced partition's worst rank -- so any rank landing >6 local
  elements overflowed that fixed-size Fortran array. `cpu_per_node=192`
  (inherited verbatim from `nrsqsub_teton-2.sh`, sized for a full Turin
  node on production jobs) was simply too many ranks for these small
  meshes. Added a `TASKS_PER_NODE` env var (default 32, comfortably safe
  for all three cases' element counts) so `cpu_per_node` is no longer
  hardcoded at 192; override it explicitly if more parallelism is wanted
  later.

### New Session

- Re-ran `concentric_spheres_vf` at 64 ranks (2 nodes, `TASKS_PER_NODE=32`,
  `lelt=25` -- comfortable headroom over the 23 elements/rank actual max)
  and hit the *same* segfault right after `VF:Reading view factor file:
  nwalls= 588` as the original 192-rank/`lelt=6` run. This rules out
  undersized `lelt` as the actual root cause of that segfault (or at least
  as the only cause) -- the earlier fix was still correct/necessary (a
  real overflow risk existed at 192 ranks), but insufficient. Verified by
  parsing `vf_concentric_spheres` in Python against the exact Fortran read
  logic in `view_factors.f:298-395`: the file itself is well-formed (588
  wall blocks, self-consistent, ends exactly at EOF, max
  visible-walls-per-face 355 vs `max_visible_walls=10000`, max global
  element id 1470 matches `nelg`) -- not a data/format bug. Current
  hypothesis, not yet confirmed: `view_factors.f:326-333` has every MPI
  rank independently `open()` and fully sequentially re-read the *same*
  shared file with no `iostat=`/`err=` checking anywhere in the read loop
  -- a known-fragile pattern under concurrent access on parallel/network
  filesystems, and consistent with the observed failure hitting a subset
  of ranks (not all, not rank 0).
- To isolate the above, tried `TASKS_PER_NODE=1` (single task, ruling out
  any inter-rank race entirely). Got a different, unrelated failure: the
  `udf` build step OOM-killed (`MPICH ERROR ... Abort(1)`, `oom_kill event
  ... task 0: Out Of Memory`) before ever reaching the view-factor read.
  Root cause: `run_teton.sh`'s sbatch header has `--exclusive` but no
  `--mem` request; on this site's Slurm config, `--exclusive` alone does
  not guarantee the job's memory cgroup covers the full node -- without an
  explicit `--mem`, the limit is computed from `ntasks x
  (default-mem-per-cpu)`, so a 1-task job got a ~1/192nd-of-node memory
  cap, well below what compiling `udf.cpp` needs (this wasn't a problem at
  32/64 tasks, where the same per-task default multiplied out to enough
  headroom). Fixed by adding `#SBATCH --mem=0` (Slurm idiom for "all
  memory on the node") to all three `run_teton.sh` scripts. Also lowered
  `concentric_spheres_vf/run_teton.sh`'s `TASKS_PER_NODE` default from 32
  to 8 (avoids both the original segfault's higher rank counts and a
  1-task build, while still leaving ~184 elements/rank, comfortably under
  `lelt`) and added `export GFORTRAN_ERROR_BACKTRACE=1` to its generated
  sbatch script, to get an actual gfortran stack trace/source line the
  next time the `view_factors.f` segfault is hit, instead of only "core
  dumped" with no further detail. Not yet applied to
  `parallel_plates_vf`/`single_pebble_vf`'s `TASKS_PER_NODE` defaults or
  `GFORTRAN_ERROR_BACKTRACE` export -- only `--mem=0` was added to those
  two so far, since `concentric_spheres_vf` is still the one being
  actively debugged.
- Confirmed the race-condition theory: `TASKS_PER_NODE=1 ./run_teton.sh
  concentric_spheres 1 00:15` (with `--mem=0` in place) ran to completion
  cleanly -- `finished with exit code 0`, all 5 steps converged, checkpoint
  file written. This proves the `view_factors.f` segfault only occurs with
  more than one rank, i.e. it's concurrency-dependent, not a fixed
  data/sizing bug -- consistent with every rank independently reopening
  and re-reading the same shared `vf_concentric_spheres` file with no
  `iostat=`/`err=` checking anywhere in the read loop.
  Fixed the race directly in `view_factors.f`'s `vf_read_view_factors`
  (identical file across `examples/concentric_spheres_vf/`,
  `examples/parallel_plates_vf/`, `examples/single_pebble_vf/`, and
  `tall_cavity_vf_aurora/` -- confirmed byte-identical via `diff` before
  patching all four the same way): only rank 0 now opens and reads the
  file (two passes -- first to count total visible-wall records so exact-
  size buffers can be allocated, second to actually store `iw`/`ieg`/
  `ifc`/`nvwalls` headers and the flat `jw`/`view_factor` visible-wall
  data), then `call bcast(...)` (Nek5000's existing MPI_Bcast wrapper,
  `3rd_party/nek5000/core/comm_mpi.f`, already used this way throughout
  nek5000 core, e.g. `ic.f`'s `call bcast(time,wdsize)`) distributes the
  parsed buffers to every rank, using `isize`/`wdsize` (existing Nek5000
  common-block byte-size globals, precision-agnostic) for the byte
  counts. Every rank then applies the exact same `gllnid(ieg).eq.nid`/
  `gllel(ieg)` ownership filter as before, purely from the broadcast
  buffers in memory -- no rank other than 0 ever touches the filesystem.
  Used Fortran `ALLOCATABLE` arrays (not previously used anywhere in
  `3rd_party/nek5000/core`, but standard Fortran 90+ and no `-std=`
  restriction found in `cmake/nek5000.cmake`/`CMakeLists.txt` limiting to
  strict F77) instead of a new hardcoded max constant, since the existing
  `max_walls=4000000`/`max_visible_walls=10000` pattern only bounds
  per-wall counts, not the total flattened visible-wall count across all
  walls (129,108 for this case alone) -- no existing constant to reuse
  safely. `vf_export_partitioned_vf_files`/`vf_read_partitioned_vf_files`
  (the other file-I/O routines in this file) were left untouched: each
  rank there reads its own distinct pre-partitioned file (`vfp<nid>`), so
  there's no shared-file race to begin with, and neither is called from
  any of these cases' `.usr` files anyway.
  Also reverted `concentric_spheres_vf/run_teton.sh`'s `TASKS_PER_NODE`
  default back to 32 (from the temporary 8) now that the actual bug is
  fixed at the source rather than avoided by rank count; updated its
  header comment accordingly. Not yet re-tested on Teton -- this fix
  needs a real multi-rank run to confirm before touching
  `parallel_plates_vf`/`single_pebble_vf`'s `run_teton.sh` defaults or
  declaring the other two cases' `TASKS_PER_NODE`/`GFORTRAN_ERROR_BACKTRACE`
  settings finalized.

## 2026-09-02

- The `view_factors.f` rank-0-reads/bcast rewrite (2026-09-02, above) did
  NOT fix the actual segfault -- a clean rebuild (`rm -rf .cache`) still
  crashed identically. Tried `NEKRS_FFLAGS="-fcheck=all -fbacktrace"`
  (a new env var added to `run_teton.sh`, passed through to the per-case
  nekInterface Fortran build via
  `src/core/nekInterface/nekInterfaceAdapter.cpp`'s `FFLAGS="${NEKRS_FFLAGS}"
  make ...`) to force a bounds-check error instead of a bare segfault --
  this immediately false-positived on `3rd_party/nek5000/core/math.f`'s
  `BLANK(A,N)` (`CHARACTER*1 A(1)` called with `N>1`, standard F77
  storage-association usage that gfortran's `-fcheck=all` can't
  distinguish from a real bug) during nekRS's own MPI bootstrap, before
  any case code runs -- abandoned as unusable here, `NEKRS_FFLAGS`
  reverted to empty/unused in `run_teton.sh`.
  Root-caused properly instead via manual `write(6,*)`+`flush(6)`
  checkpoints (prefixed `VFDBG`) added directly into
  `vf_read_view_factors` at every stage (open, both read passes, every
  `bcast`, allocate, the ownership-filter loop) across all four case
  copies of `view_factors.f`. Result on a 32-rank run: the 12 ranks that
  ultimately segfaulted (exact match against the `srun: error` task
  list) never printed even the unconditional first line of the
  subroutine -- i.e. they crashed *before* `vf_read_view_factors` was
  ever called, and the remaining 20 ranks were left hanging in the
  first `call bcast(nwalls,isize)` (a collective the 12 dead ranks could
  never join), which is what actually produced the observed "stuck at
  bcast" symptom and the eventual Slurm-reported segfaults across all
  32 ranks. This proves the crash was never inside `view_factors.f` in
  the first place -- likely was never the concurrent-file-read race
  either, across every prior test in this saga.
  Actual root cause, found in `concentric_spheres.usr`'s `userchk()`:
  `if (istep.eq.0) call rzero(frad(1,1,1,1),ntot)` is the *first*
  executable statement in the subroutine, and `ntot` is never assigned
  anywhere in `userchk()` -- under Fortran implicit typing it defaults
  to integer, but with no assignment its value is whatever garbage
  happened to be on that rank's stack at entry. `rzero` then zeroes out
  that many elements of `frad` starting from its base address --
  writing out of bounds by an arbitrary, per-rank-random amount whenever
  the garbage value exceeded `frad`'s true size, corrupting adjacent
  memory and segfaulting. This explains every previously-observed
  symptom at once: a random subset of ranks affected (each rank's stack
  garbage differs), rank 0 consistently surviving (coincidental but
  repeatable across many runs -- not something to rely on, just what we
  observed), the crash always appearing to happen "right after" the
  `view_factors.f` output (that output is rank-0-only and printed
  successfully because rank 0 kept surviving; other ranks had no
  per-rank diagnostic output before this VFDBG pass, so nothing showed
  their actual, earlier failure point), and the 1-task run succeeding
  cleanly (no other rank's stack garbage could ever be involved with
  only one rank).
  Confirmed the same bug (`ntot` used, never assigned, as the first line
  of `userchk()`) is present identically in `parallel_plates.usr`,
  `single_pebble.usr`, and `tall_cavity_vf_aurora/tall_cavity.usr`.
  Fixed all four by declaring `ntot` explicitly and computing
  `ntot = lx1*ly1*lz1*lelt` before its first use -- matching the exact
  convention already used a few hundred lines away in the same
  `view_factors.f` (`vf_calculate_radiation_heat_flux` does
  `ntot = lx1*ly1*lz1*lelt` before its own `rzero(frad(1,1,1,1),ntot)`
  call), strongly suggesting this is what the original case author
  intended in `userchk()` too and simply omitted.
  Removed the temporary `VFDBG` write/flush instrumentation from all
  four `view_factors.f` copies, restoring them to the clean
  rank-0-reads/bcast rewrite from 2026-09-02 (kept, since it's a real
  hardening against a genuine, separate concurrent-shared-file-read
  risk pattern -- just not the cause of this particular bug). Reverted
  `run_teton.sh`'s header comment to describe the actual root cause and
  removed the now-unused `NEKRS_FFLAGS` plumbing; left
  `GFORTRAN_ERROR_BACKTRACE=1` in place (harmless, and it may still fire
  usefully for some future, different crash even though it produced
  nothing for this one -- Cray MPICH's own signal handling most likely
  intercepted the SIGSEGV first).
  Not yet re-tested on Teton.
- Fixed `tall_cavity_vf_aurora/tall_cavity.par`/`.udf`/`.oudf` to this
  checked-out nekRS source's actual syntax, the same class of fixes already
  applied to `concentric_spheres_vf`/`parallel_plates_vf`/`single_pebble_vf`
  on 2026-09-01 but not yet done for this fourth case: `.par`'s
  `[PRESSURE]`/`[VELOCITY]`/`[TEMPERATURE]` -> `[FLUID PRESSURE]`/
  `[FLUID VELOCITY]`/`[SCALAR TEMPERATURE]` + `scalars = temperature` in
  `[GENERAL]`, `PFGMRES` -> `GMRES+flexible` with `maxIterations` folded
  into `+maxIter=1000` on the solver line, `boundaryTypeMap`'s
  `codedFixedValue` -> `udfDirichlet`. `.oudf`'s `codedFixedValueVelocity`/
  `codedFixedValueScalar`/`codedFixedGradientScalar` merged into one
  `udfDirichlet`/`udfNeumann` pair gated by `isField(...)`, with
  `bc->u/v/w`->`bc->uxFluid/uyFluid/uzFluid`, `bc->s`->`bc->sScalar`,
  `bc->flux`->`bc->fluxScalar`, `bc->idM`(nonexistent)->`bc->idxVol`.
  `.udf` was more involved than the other three cases since this is the
  only one of the four using `uservp`/`userf` (buoyancy-driven momentum
  source + non-par-file properties) rather than pure par-file properties:
  `nrs->cds`->`nrs->scalar`, `nrs->mesh`->`nrs->meshV`,
  `nrs->o_prop`->`nrs->fluid->o_prop`, `cds->o_prop`->`nrs->scalar->o_prop`,
  `cds->fieldOffset[0]`->`nrs->scalar->fieldOffset()`,
  `nrs->userVelocitySource`->`nrs->userSource` (confirmed via
  `src/app/app.hpp`'s `userProperties_t`/`userSource_t` members -- no
  separate velocity-specific source member exists), and
  `nrs->o_NLT`->`nrs->fluid->o_EXT` for the buoyancy source-injection
  target. The `uservp` rewrite followed `examples/lowMach/lowMach.udf`'s
  exact pattern (`nrs->scalar->fieldOffset()`, `nrs->fluid->o_prop`,
  `nrs->scalar->o_prop`) and the `userf` rewrite followed
  `examples/rbc/rbc.udf`'s (`nrs->fluid->o_EXT`, `nrs->scalar->o_S`,
  `platform->linAlg->axpby`) -- both proven-current, already-building
  examples in this repo, cross-checked since no single existing example
  combines buoyancy *and* radiation to compare against directly.
  **Flagged as unverified**: unlike the other three cases (each confirmed
  running to completion on Teton), `tall_cavity_vf_aurora` has not been
  run at all since these changes -- the `uservp`/`userf` rewrite in
  particular is a syntax-level port only, not something I could test.
  `tall_cavity.usr`'s `ntot` fix (same bug as the other three) was already
  applied earlier the same day, before this entry.
- Mirrored all of the above (the `view_factors.f` rank-0-reads/bcast
  rewrite, the `ntot` fix, and the full `.par`/`.udf`/`.oudf` v26.0.1
  syntax fixes) into the user's global `nek-vf-case` skill
  (`~/.claude/skills/nek-vf-case/`, outside this repo -- generates new
  Nek-VF cases from scratch) so future skill-generated cases start from
  the corrected pattern instead of reproducing this entire debugging
  saga: `assets/view_factors/view_factors.f` replaced with the fixed
  rewrite (its `README.md` updated to describe rank-0-reads/bcast instead
  of "every rank reads redundantly"), `templates/case.usr.template` now
  computes `ntot` explicitly, `templates/case.par.template`/
  `case.udf.template`/`case.oudf.template` updated to the current section
  names/solver syntax/`meshV`/`o_usrwrk` location/unified
  `udfDirichlet`/`udfNeumann` convention, and a new
  `templates/run_teton.sh.template` added (parameterized generalization of
  this repo's `run_teton.sh` scripts, `PROJ_ID` now a hard-required env var
  rather than defaulting to the never-valid `"art"` placeholder). `SKILL.md`
  updated throughout to match (new Step 6 for Teton submission, stale
  `codedFixedValue*`/`nrs->mesh`/`nrs->cds`/`[VELOCITY]`/`[TEMPERATURE]`
  references removed, the `tall_cavity_vf_aurora` "known limitation" note
  rewritten now that the reference case itself is current again rather
  than demonstrating the old convention).
- Compared `parallel_plates_vf`'s converged (step-3-onward) run output
  against the paper's Table I: hot wall 3.7055e4 W/m2 (paper Nek-VF
  3.92e4, analytical 3.93e4 -- ~5.5% low), cold wall -2.9741e4 W/m2
  (paper -3.88e4 -- ~23% low), side wall -0.7246e4 W/m2 (paper -0.99e4 --
  ~27% low). Root-caused to `parallel_plates.geo`'s `W=4.0`/`H=1.0`
  (plate-side/spacing = 4) not being a good enough approximation of Eq.
  11's idealized infinite-parallel-plates limit: the run's own printed
  `totalAA` came out to exactly 16 for all three sidesets (hot-plate area,
  cold-plate area, *and* total side-wall area all equal at this ratio),
  meaning the side walls have as much area as either plate and intercept/
  re-radiate a share of the exchange the analytical formula assumes is
  negligible -- consistent with the hot wall (least side-wall-adjacent)
  deviating least and the cold/side walls (most side-wall-coupled)
  deviating most.
  Regenerated the case's geometry/mesh/view-factor file via the
  `nek-vf-case` skill to test a larger aspect ratio: raised `W` to 20.0 in
  `parallel_plates.geo` (`H` unchanged), documented the reasoning in the
  `.geo` file's own header comment. Remeshed
  (`gmsh -3 parallel_plates.geo -o parallel_plates.msh -format msh2`,
  using the "rvf" conda env's gmsh since `/Applications/Gmsh.app` wasn't
  invoked this time) and reran `gmsh2nek` (`printf
  "3\nparallel_plates\n0\n0\nparallel_plates\n" | gmsh2nek`) --
  confirmed element count unchanged at 4400 hex / 1680 boundary faces,
  since only `W` changed (element *size* grew, not element *count*,
  because `N_W`/`N_H` were left alone). Recomputed view factors with the
  same invocation `procedure.md` documents for the original geometry
  (`compute_view_factors.jl --backend=cpu --monte-carlo=true
  --n-samples=4000 --nquad=6 --closure-tol=0.1`): row-sum closure error
  1.88% max (vs. 0.49% on the original `W=4.0` mesh -- higher but still
  well within the 10% tolerance used, expected given the more elongated
  geometry), reciprocity exact to machine precision (3.05e-16). Overwrote
  `vf_parallel_plates` with the new file; updated
  `examples/parallel_plates_vf/procedure.md`'s "what the paper specifies
  vs. assumes" / mesh-generation / view-factor sections to document the
  change and both sets of results.
  **Not yet re-run on Teton** -- the new geometry/view-factor file has not
  been tested at all; `parallel_plates.par`/`.udf`/`.oudf`/`.usr`/
  `run_teton.sh` were intentionally left untouched (unrelated to the
  geometry change), so the next run should build cleanly, but the actual
  flux values against the paper's Table I are still unverified against
  this new mesh.

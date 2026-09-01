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

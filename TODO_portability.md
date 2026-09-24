# Portability TODO

Ideas to make it easier for other institutes to run this pipeline, raised while working on `ont-bmod-restructuring`. Not urgent, revisit later.

- **Add CI** (e.g. `.github/workflows/test.yml`) running `-profile test,docker` (and ideally `test_annotation`) on every push/PR, across at least two Nextflow versions. Motivated directly by a real bug this session: `publishDir` broke on Nextflow 26.04.6 but not 24.04.4, and only got caught because both versions were tested by hand.
- **Document that SLURM isn't required** — Nextflow supports PBS/Torque, SGE, LSF, HTCondor, cloud batch executors, etc. Adding support for another scheduler only means a new profile block in `nextflow.config` (like the existing `local`/`slurm` ones), no changes to the actual pipeline processes. README currently only documents SLURM.
- **Update README for the newer `--flow` modes** — it still only documents the original single main-flow behavior. No mention of `--flow annotation` (or `--flow dmr` once built), so someone reading the README wouldn't discover they exist.
- **Resolve the `rki-mf1` vs `valegale` remote question** — README tells users to `nextflow pull rki-mf1/ont-methylation`, but this repo's actual remote is `valegale/ONT_methylation`. Is `rki-mf1` the upstream this work should flow back into, or is the README just stale?
- **Container registry robustness (lower priority)** — containers are pulled from a personal Docker Hub namespace (`rkimf1/*`). Docker Hub rate-limits anonymous pulls per IP, which can bite shared HPC clusters pulling through one egress IP. Consider pinning by digest and/or mirroring to a more durable registry (e.g. `ghcr.io`) if this becomes a real problem.

Confirmed already fine, no action needed:
- No hardcoded paths anywhere in the pipeline config (checked via `grep`).
- Test profiles (`test`, `test_meta`, `test_annotation`) already let a new site verify its own Docker/Singularity setup works before running real data.

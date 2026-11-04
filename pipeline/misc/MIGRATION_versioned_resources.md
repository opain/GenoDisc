# Runbook — migrate to the per-reference-epoch resource store

Goal: split the single shared `pipeline_resources` into `resources/{common,v-dense,v-legacy}`
so pinned pipeline versions stop sharing reference-derived resources. Do this in a
supervised window. Everything is **non-destructive** (hardlinks) and **reversible**
(a `.ready` marker arms each epoch; `rm` it to roll back).

All paths assume `GD=/scratch/prj/neurohackpain/GenoDisc`. Run on CREATE as a project
member with write access to the GenoDisc tree.

```bash
GD=/scratch/prj/neurohackpain/GenoDisc
S=$GD/repo/current/pipeline/misc/build_epoch_resdir.sh
```

## Phase 0 — Preconditions (verify, change nothing)
- `resources/` does not exist yet → the web app still uses `pipeline_resources` (unchanged).
- Raw sources staged: `pipeline_resources/data/{ROSMAP_pQTL,ROSMAP-PWAS,Banner-PWAS}`.
- Dev pipeline (`repo/current`) carries the source-wiring + build-not-relocate commits.
- No active jobs need to be drained — hardlinking leaves `pipeline_resources` intact.

## Phase 1 — Choose the copy mode
Hardlinks are instant and use no extra space. Test support:
```bash
: > "$GD/pipeline_resources/.hltest"
ln "$GD/pipeline_resources/.hltest" "$GD/pipeline_resources/.hltest.ln" \
  && echo "HARDLINKS OK -> use --mode link" || echo "NO HARDLINKS -> use --mode copy (slow)"
rm -f "$GD/pipeline_resources/.hltest" "$GD/pipeline_resources/.hltest.ln"
```
Use `--mode link` if OK, else `--mode copy`. (Avoid `--mode move`: it empties
`pipeline_resources` and breaks the fallback + any running jobs.)

## Phase 2 — Build the DENSE store skeleton
```bash
bash "$S" --epoch dense --mode dry-run      # review: common=symlinks, 1kg_hgdp/magma/... relocated,
                                            #         rosmap_smr -> [rebuild] empty, sources -> SHARED
bash "$S" --epoch dense --mode link         # builds resources/common + resources/v-dense
```
This creates `common/` (hardlinks of the reference-independent resources + raw sources),
`v-dense/` (symlinks into common + hardlinked reference-derived dirs), and an **empty**
`v-dense/data/rosmap_smr` to be rebuilt next. `pipeline_resources` is untouched. Nothing
is armed yet.

## Phase 3 — Rebuild rosmap_smr fresh (Step 4), then verify
The script prints the exact command. Dry-run FIRST and confirm it schedules **only**
`format_rosmap_smr_data` (if it wants to rebuild `1kg_hgdp` or anything large, STOP —
mtimes are off):
```bash
source $GD/bin/activate_genodisc.sh
SMK="snakemake prepare_resources --snakefile $GD/repo/current/pipeline/Snakefile \
  --configfile $GD/repo/current/pipeline/prepare_resources.yaml \
  --config resdir=$GD/resources/v-dense --use-conda -c4"
$SMK -n          # REVIEW: expect only rosmap_smr-related jobs
```
Then run it on a compute node:
```bash
sbatch --partition=neurohack_cpu,cpu --mem=8G -c4 --time=02:00:00 \
  --output=$GD/resources/v-dense/logs/build_dense_%j.out \
  --wrap "source $GD/bin/activate_genodisc.sh && $SMK"
```
Verify the rebuilt besd (reads the EUR allele freqs, so ancestry-correct):
```bash
P=$GD/resources/v-dense/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd
ls -la "$P".{besd,esi,epi}                                   # all present, mtime = today
$GD/resources/v-dense/software/smr/smr_linux_x86_64 \
  --beqtl-summary "$P" --query 5.0e-8 --out /tmp/rosmap_check # opens cleanly, returns rows
find $GD/resources/v-dense -type l ! -exec test -e {} \; -print  # (empty = no broken symlinks)
```

## Phase 4 — Go live for DENSE
1. Deploy the Django resdir-by-epoch code to prod (merge `ollie_updates` → `main`, then the
   normal deploy). Until deployed, prod ignores the markers and keeps using
   `pipeline_resources` — so this step is safe to do before or after arming.
2. Arm (atomic go-live — the web app switches dense jobs to v-dense):
   ```bash
   touch $GD/resources/common/.ready
   touch $GD/resources/v-dense/.ready
   ```
3. Smoke-test: submit a small dense job (or `./manage.py shell` →
   `JobService()._build_config_yaml(...)`) and confirm `resdir` is `.../resources/v-dense`.

**Rollback** at any time: `rm $GD/resources/v-dense/.ready` → dense jobs fall straight back
to `pipeline_resources`.

## Phase 5 — LEGACY epoch (best-effort, later)
```bash
bash "$S" --epoch legacy --mode link \
  --pipeline $GD/repo/versions/v1.3.5/pipeline    # a pre-v1.3.6 (legacy-reference) pipeline
```
Caveat: the v1.3.6 migration overwrote some legacy reference-derived resources in place
(e.g. `magma` annotation, `rosmap_smr` esi were rebuilt against the dense reference). Legacy
`data/1kg` is intact, but rebuilding the legacy reference-derived resources needs the OLD
reference and may need old sources that no longer exist. Add the clobbered resources to the
script's `REBUILD_REF` and run the Step-4 build with a legacy pipeline + `resdir=v-legacy`;
accept/flag anything that cannot be faithfully reproduced. Arm only once verified:
`touch $GD/resources/v-legacy/.ready`. Until then, legacy jobs keep using `pipeline_resources`.

## Phase 6 — Finish up
- Push the pipeline `dev` commits and cut a release (`bin/cut_release.sh vX.Y.Z`) so new job
  pins carry the source-wiring, the reference-dep declarations, and the `resource_epoch` tag.
- Push the Django `ollie_updates` commits (deployed in Phase 4).
- Reclaim space ONLY after v-legacy is armed (or you accept no legacy reruns): once nothing
  reads `pipeline_resources`, `rm -rf $GD/pipeline_resources` (hardlinks mean the bytes
  survive in `common/`/`v-dense/`). Keep a short grace period first.

## Quick reference — what each state means
| State | dense job resdir | legacy job resdir |
|---|---|---|
| no markers (now) | pipeline_resources | pipeline_resources |
| common+v-dense armed | **resources/v-dense** | pipeline_resources |
| + v-legacy armed | resources/v-dense | **resources/v-legacy** |

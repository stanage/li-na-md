# High-Throughput MD for the 1000-Electrolyte Oxygen-Only Dataset — Working Plan

**Status:** dataset built, 278/278 solvent force fields verified, 22 ion species
grompp-validated 40/40, all 1000 rows indexed — **campaign submitted**
**Created:** 2026-09-21
**Owner:** eshiemogie@uchicago.edu
**Working dir:** `MD/ht1000/` (all new files live here)

This is a living document. Update it as decisions change — it is the single source of
truth for how the 1000-electrolyte campaign is built and run.

`MD/ht2000/PLAN.md` is the companion document for the sister campaign. It carries the
full derivations that both campaigns share — where the box-and-count math came from, why
the ion force fields are assembled from four sources, the LigParGen failure modes, the
Midway3 cluster gotchas. **This plan does not repeat them**; it records what ht1000 does
*differently* and why, and points at ht2000's plan for the rest.

---

## 0. Quick start

Everything runs with the `moleng` conda python; nothing needs to be installed.

**One-time setup on a fresh checkout.** `runs/` is a symlink to scratch and is not
tracked, because the trajectories do not fit in the pi-chibueze /project quota:

```bash
mkdir -p /scratch/midway3/$USER/ht1000_runs
ln -s /scratch/midway3/$USER/ht1000_runs MD/ht1000/runs
```

```bash
PY=/scratch/midway3/eshiemogie/moleng/bin/python
cd /project/chibueze/stanley/projects/li-na-md/MD/ht1000/scripts

# --- one-time build, in this order (all idempotent, all already done) ---
$PY select_solvent_pool.py     # oxygen-only pool  -> datasets/solvents_oxygen_pool.csv
$PY build_dataset.py           # 1000 rows         -> datasets/electrolyte_dataset_1000.csv
$PY build_ion_ff.py            # 22 ion species    -> ff_ions/
$PY validate_ions.py --cation Li && $PY validate_ions.py --cation Na
$PY import_solvent_ff.py       # copy verified FFs from ht2000's library
$PY fetch_solvent_ff.py --smiles-file ../datasets/solvents_to_fetch.txt   # LOGIN NODE
$PY verify_solvent_ff.py --fix # ALWAYS audit before a campaign
$PY build_run_manifest.py      # index all 1000 rows -> index/run_manifest.csv

# --- THE CAMPAIGN -- one submission, start or continue ---
sbatch ../scripts/submit_campaign.sbatch

# --- status / results / housekeeping at any time ---
$PY campaign_status.py --lanes     # progress, failures, ETA
$PY campaign_status.py --failures  # which runs died and where
$PY collect_results.py             # -> index/results.csv
./maintain.sh                      # refresh indexes + status + push to Box (login node)
./maintain.sh --prune              # ... and free the uploaded trajectories
```

The campaign is **resumable by design**: lanes skip runs that already finished and
`run_md.sh` skips stages that already finished, so the *same* `sbatch` command both
starts the campaign and continues it after a wall-clock kill. Expect to run it ~15 times.
Every other command above also skips work already done.

---

## 1. Goal, and why this campaign exists alongside ht2000

Run the same solvation-MD workflow as ht2000 over a dataset whose solvents can actually
answer the question the analysis asks.

### The problem with ht2000's solvent library

`analyze_solvation.py` measures the **cation–solvent-oxygen RDF**: first peak position,
first minimum, coordination number `cn_solv_O`. Then `cluster_analysis.py` and the
SSIP/CIP/AGG split are all interpreted relative to how well the solvent competes with the
anion for the cation's first shell.

ht2000's 500-solvent library was generated with filters that never mentioned oxygen:

| | ht2000 library | consequence |
|---|---|---|
| contain no oxygen at all | **222 / 500** | `cn_solv_O` is identically 0; the M–O RDF is empty and `rdf_first_peak_A` is fitted to noise |
| pure hydrocarbons (C, H only) | **106 / 500** | not electrolyte solvents in any practical sense — no dipole to solvate an alkali cation |

Those rows are not *wrong*, they are *uninformative*: they consume a full 6.1 ns run each
and return a solvation descriptor that is structurally undefined. Roughly 44% of ht2000's
library is in that state.

### What ht1000 changes

One rule, applied at pool-selection time:

> **Every solvent must contain at least one oxygen atom.**

That removes all 222 O-free molecules and — because a pure hydrocarbon has no oxygen —
all 106 pure hydrocarbons, by construction rather than by a second filter. Everything
else about the campaign is held fixed against ht2000 so the two are directly comparable
and can be pooled into one dataset later.

---

## 2. What the dataset looks like

`datasets/electrolyte_dataset_1000.csv` — 1000 rows, columns
`salt, solvent, concentration, solv_per_ion_pair`.

- **278 unique solvents** in the pool, **274** actually drawn into the dataset (the
  sampler simply never drew the other 4). Every one contains oxygen.
- **40 unique salts** = 20 anions × {Li⁺, Na⁺} — identical to ht2000, same SMILES strings
  character for character, so both campaigns key the same ion force fields.
- **Concentration 0.3–1.5 M** (28-point grid) — identical to ht2000.
- Single solvent per row.
- **MW 88–199** (median 184). Elements present: **C, H, O** always, plus S (49), N (48),
  F (27), Cl (13), P (5). No boron, no silicon.
- **Median 5.82 solvent molecules per ion pair**; 194 rows below 4, **none below 2**.
- 25 rows per salt, drawn without replacement, so no `(salt, solvent)` pair repeats. Each
  solvent appears in 1–10 rows (mean 3.6) against *different* salts — comparable to
  ht2000's 4.1, so neither campaign leans harder on any one molecule.

### 2.1 Where the pool comes from

The pool is the oxygen-bearing **subset of ht2000's existing 500-solvent library**, not a
fresh Electrolyte-GPT run. Generating a new library needs a GPU node; taking the subset
was chosen at the user's direction on 2026-09-21 because the existing library already
contains 278 qualifying molecules with force fields almost entirely in hand.

`select_solvent_pool.py` re-applies every ht2000 filter from scratch rather than trusting
the source CSV, so the pool is a complete statement of what it may contain:

| filter | rule | inherited from ht2000 |
|---|---|---|
| **oxygen** | **≥ 1 O atom** | **no — this is ht1000's addition** |
| allowed elements | `H, C, N, O, F, P, S, Cl, Br, I` (no boron, no silicon) | yes |
| charged / multi-fragment | rejected | yes |
| MW ceiling | ≤ 200 g/mol | yes |
| aldehyde | `[CX3H1](=O)[#6]` — reduces on Li metal | yes |
| protic | `[OX2H,NX3;H1,H2]` — O–H/N–H attacks Li⁰ | yes |
| radical | any unpaired electrons — BOSS cannot type them | yes |
| phosphate/phosphite ester | `[#15;!$([#15]~[#6])]` — no BOSS template | yes |

Re-running the filters over ht2000's 500 rejects exactly 222 molecules, all for
`no_oxygen`. That is a useful check in itself: it confirms the source library already
satisfies every inherited rule, so the oxygen requirement is the only thing separating
the two pools.

### 2.2 Functional-family census

Documented rather than enforced. The pool was fixed before the dataset was drawn, so
quotas could only have been applied at row-sampling time; uniform sampling was chosen
(2026-09-21) and the composition recorded here instead.

| family | molecules (of 278) |
|---|---|
| ether | 126 |
| ketone | 123 |
| ester | 83 |
| nitrile | 30 |
| sulfone | 29 |
| fluorinated | 27 |
| carbonate | 23 |
| chlorinated | 13 |
| amide | 8 |
| phosphoryl | 3 |
| sulfoxide | 2 |

Oxygens per molecule: 1 → 123, 2 → 103, 3 → 34, 4 → 15, 5 → 3.

**The scarce families are the ones worth watching in the results.** Amides (8), sulfoxides
(2) and phosphoryls (3) are thin enough that any conclusion drawn about them rests on a
handful of molecules. Carbonates at 23 are the family closest to real commercial
electrolytes and are adequately but not generously represented.

### 2.3 Known data issues (carried over from ht2000)

1. **Broken nitrate SMILES.** `[O-][N+](=O)=O` has N valence 5 and fails RDKit. Kept
   verbatim in the dataset so both campaigns' `salt` columns match; repaired by the
   `SALT_FIXUPS` map in `scripts/common.py` before anything is simulated, and the
   corrected SMILES is what lands in each run's `summary.json`.
2. **Requested molarity ≠ achieved molarity.** The builder guesses a density
   (ρ_salt = 1.70 g/cm³ flat, ρ_solvent = 1.00 g/cm³) to convert molarity into molecule
   counts, so NPT fixes the density but not the composition. `collect_results.py` records
   `molarity_actual` and `molarity_err_pct` per run and flags anything off by >10%.
   **Report the achieved molarity, not the requested one.** See ht2000 PLAN.md §2.2.
3. **Na⁺ ECC scaling.** q = 0.8 is well established for Li⁺; the same value for Na⁺ is a
   reasonable default but is an assumption worth stating in any writeup.

---

## 3. Isolation from ht2000 — the requirement, and how it is met

The two campaigns run concurrently and must never touch each other's files, in scratch or
in Box. Every writable path resolves through `common.HT`, which is
`Path(__file__).parent.parent` — so it is structurally impossible for an ht1000 script to
write into ht2000.

| | ht2000 | ht1000 |
|---|---|---|
| working dir | `MD/ht2000/` | `MD/ht1000/` |
| dataset | `datasets/` **at the repo root** | `MD/ht1000/datasets/` — inside the campaign dir |
| solvent FF | `MD/ht2000/ff_solvents/` (490) | `MD/ht1000/ff_solvents/` (278) |
| ion FF | `MD/ht2000/ff_ions/` | `MD/ht1000/ff_ions/` — rebuilt from source, not copied |
| run dirs | `/scratch/midway3/eshiemogie/ht2000_runs` | `/scratch/midway3/eshiemogie/ht1000_runs` |
| Box | `.../li-na-md/ht2000/runs` | `.../li-na-md/ht1000/runs` |
| SLURM job | `ht2000`, 3-task array, 1 GPU each | `ht1000`, single job, 1 node, 2 GPUs |
| slurm logs | `logs/slurm/ht2000_%A_%a.out` | `logs/slurm/ht1000_%j.out` |

**`DATASET` deliberately lives under `HT`, not at the repo root.** ht2000 reads
`REPO/datasets/electrolyte_dataset_2000.csv`; putting ht1000's copy inside its own
campaign directory means a dataset path that resolves through `HT` cannot accidentally
point at the other campaign's.

### What the two *do* share, and why that is safe

- **Read-only force-field source data**: `MD/Force Fields/`, `MD/anions/`,
  `MD/ce_solvation_md/`. These are inputs. `build_ion_ff.py` reads them and writes only
  into `ht1000/ff_ions/`, so ht1000's 22 ion species are its own files built from the
  same sources — not copies of ht2000's.
- **`MD/ht2000/ff_solvents/`, read once, at build time only.** `import_solvent_ff.py`
  copied 272 verified solvent force fields out of it (§4). After that one pass nothing in
  ht1000 reads anything under `ht2000/` at any point: deleting ht2000 entirely would not
  affect a single run.
- **The `moleng` conda env and the `gromacs/2025.3` module.** Shared read-only software.

---

## 4. Solvent force fields — 278/278, verified

`ht1000/ff_solvents/` holds `S<12hex>.{itp,gro,pdb,smi,lmp}` for all 278 pool solvents,
keyed by a hash of the canonical SMILES exactly as ht2000 keys them.

| provenance | count | how |
|---|---|---|
| copied from ht2000's verified library | **272** | `import_solvent_ff.py` |
| fetched fresh from LigParGen | **6** | `fetch_solvent_ff.py --smiles-file`, 0.7 min, 6/6 ok |

The 6 fresh ones were never fetched by ht2000 at all — its sampler never drew them into
the 2000-row dataset, so no force field was ever requested. **None of them is a LigParGen
failure**; `index/ff_failures.csv` does not exist for this campaign because nothing has
failed.

### Why copy rather than re-fetch

The 272 are the *same molecules* ht2000 already pulled, so a re-fetch would return the
same files after ~45 minutes of serial LigParGen traffic — and would re-run the single
most dangerous failure mode in this pipeline: **under load the Yale server has been
observed returning another molecule's parameters** (ht2000 PLAN.md §3.1; a `C19H38O`
topology came back for a `C16H32O` request, and it would have simulated the wrong
chemistry silently). Re-asking for 272 already-verified molecules buys nothing and
re-exposes all 272 to that risk.

Copying is not a shortcut around verification. `import_solvent_ff.py` runs
`check_solvent_ff` — full **element composition**, not just atom count, so an isomer swap
cannot pass — against each `.itp` **at the source, before copying**, so a bad file cannot
be laundered into ht1000. `verify_solvent_ff.py` then re-audits the whole library
afterwards: **278 ok, 0 bad.**

---

## 5. Ion force fields — 22 species, 40/40 validated

Built from scratch into `ht1000/ff_ions/` by `build_ion_ff.py`, from the same four
sources ht2000 uses (OPLS-2009IL for 16 anions, fftool/CL&P for TFA and Beti, the old
`ce_solvation_md` runs for FSI, LigParGen for OTs) with the same two corrections applied
— synthesised `[pairs]` from the bond graph, and normalised `[atomtypes]` layouts. Li⁺
and Na⁺ are the Aqvist parameters as tabulated in CL&P `il.ff`. ECC charge scaling
q = 0.8 is applied **to the ions only**.

Every species' charge sums to ±1.000 before scaling. `validate_ions.py` pushes all 40
cation/anion combinations through `gmx grompp`: **20/20 with Li⁺, 20/20 with Na⁺.**

Full derivation of all of this is in ht2000 PLAN.md §4 and is not repeated here.

---

## 6. Simulation protocol — unchanged from ht2000

| stage | integrator | dt | length | notes |
|---|---|---|---|---|
| `em` | steep | — | 5000 steps | emtol 1000 |
| `nvt_heat` | md | 1 fs | 100 ps | annealed 50 K → 298.15 K, v-rescale |
| `npt` | md | 1 fs | 2 ns | C-rescale, 1 bar, τ_p = 2.0 |
| `prod` | md | 1 fs | **4 ns** | NVT, `nstxout-compressed = 5000` |

Common: `cutoff-scheme = Verlet`, `coulombtype = PME`, `rcoulomb = rvdw = 1.2 nm`,
`constraints = h-bonds` (LINCS), `tc-grps = System`, `tau_t = 0.5`, T = 298.15 K.
Box and composition math, `n_salt = 64`: ht2000 PLAN.md §5.1.

Per formulation the campaign produces `prod.xtc` + `prod.tpr`, `solvation.json`,
`clusters.json`, `summary.json` and `density.xvg` — the same artifacts as ht2000 and as
the original `ce_solvation_md` campaign.

### System sizes (all 1000 rows, `index/run_manifest.csv`)

| | atoms |
|---|---|
| min | 5,053 |
| median | **11,640** |
| p90 | 26,218 |
| max | 42,672 |

`L_build` spans 4.78–8.17 nm; total ≈ **14.3M atom-runs**. 0 rows clamped to the
`MIN_SOLVENT = 10` floor. Median is slightly below ht2000's 12,541 — the oxygen-bearing
solvents are marginally heavier on average, which at fixed molarity means slightly fewer
molecules per box.

### Cluster gotchas — all inherited, all already handled

Recorded here only so they are not rediscovered: no `-ntmpi` (Midway3's GROMACS is real
MPI); `gmx_mpi`, never `gmx`; never `module purge` before `module load gromacs/2025.3`;
`run_md.sh` re-execs under `bash -l` because `module` is a login-shell function. Details
in ht2000 PLAN.md §6.

---

## 7. Directory layout

```
MD/ht1000/
├── PLAN.md                   ← this file
├── datasets/
│   ├── solvents_oxygen_pool.csv        278 oxygen-bearing solvents
│   ├── electrolyte_dataset_1000.csv    the 1000 formulations
│   └── solvents_to_fetch.txt           the 6 that needed LigParGen
├── scripts/
│   ├── select_solvent_pool.py ← NEW: the oxygen filter + composition census
│   ├── build_dataset.py       ← NEW: 1000 rows over the pool
│   ├── import_solvent_ff.py   ← NEW: one-time verified copy of shared solvent FFs
│   ├── common.py              ← paths, SALT_FIXUPS, box math, FF integrity checks
│   ├── ligpargen.py           ← LigParGen client (charge-aware)
│   ├── fetch_solvent_ff.py    ← batch solvent fetch (login node, serial)
│   ├── verify_solvent_ff.py   ← audit/repair ff_solvents/ against its SMILES
│   ├── build_ion_ff.py        ← assemble ff_ions/ from 4 sources (one-time)
│   ├── validate_ions.py       ← push every cation/anion pair through grompp
│   ├── setup_run.py           ← build ONE run dir (packmol + top + mdp)
│   ├── run_md.sh              ← the 4-stage GROMACS protocol; stage-resumable
│   ├── analyze_solvation.py   ← solvation.json from prod.xtc
│   ├── cluster_analysis.py    ← ion clusters from prod.xtc
│   ├── build_run_manifest.py  ← index all 1000 rows up front
│   ├── collect_results.py     ← all json → index/results.csv
│   ├── archive_to_box.py      ← push prod.xtc/tpr to Box, optionally prune
│   ├── submit_campaign.sbatch ← THE campaign entry point: 1 node, 2 GPUs, 2 lanes
│   ├── campaign_worker.sh     ← one lane; strides rows, logs, stops before the wall
│   ├── campaign_status.py     ← progress / failures / ETA, read from disk
│   ├── maintain.sh            ← everything the campaign does not do for itself
│   ├── launch_batch.py        ← pilot-era build+submit path (unused here)
│   ├── run_local_batch.sh     ← run a keyfile on an interactive GPU (unused here)
│   └── select_pilot.py        ← pilot selection (unused here)
├── ff_solvents/  278 solvents × 5 files
├── ff_ions/      Li.itp Na.itp <ANION>.itp/.pdb + ions_manifest.json  (22 species)
├── runs/         → /scratch/midway3/eshiemogie/ht1000_runs   (symlink, gitignored)
├── index/        run_manifest.csv, results.csv, solvent_manifest.csv
└── logs/         slurm/  campaign/  + fetch and build logs
```

`launch_batch.py`, `run_local_batch.sh` and `select_pilot.py` are carried over for
parity with ht2000 but are not used — ht1000 goes straight to the campaign driver, since
the pipeline was already validated by ht2000's pilot and there is nothing new to prove.

Run keys: `el<12-hex>` of `canonical(salt)|canonical(solvent)|M` — deterministic, so a row
always maps to the same run dir and the campaign is restartable. All 1000 keys are unique.

---

## 8. Running it: one node, two GPUs, two lanes

### The binding constraint is the node budget, not the GPU count

```
gpu QOS:  MaxSubmitPU = 12    MaxTRESPU = cpu=192, gres/gpu=16, node=4    MaxWall = 1-12:00:00
```

**`node=4` is what shapes this campaign.** ht2000 is a 3-task array of `gpu:1` jobs, and
those three tasks land on three *different* nodes — so it is already holding 3 of the 4
permitted nodes. ht1000 gets the fourth.

| | considered | chosen |
|---|---|---|
| request | 2-task array, `gpu:1` each | **1 job, `--nodes=1 --gres=gpu:2`** |
| nodes consumed | up to 2 (tasks can scatter) | **1** |
| concurrent runs | 2 | 2 |
| node budget left | 0 (over the cap alongside ht2000) | 0 used beyond the fourth |

Same throughput, half the node budget. A 2-task array of single-GPU jobs would have been
the obvious mirror of ht2000, but it can scatter across two nodes and would then put the
user at 5 nodes — over the QOS cap, with the second task pending indefinitely.

### GPU assignment — the one place this differs from ht2000, and it is load-bearing

ht2000's hard-won rule is **never set `CUDA_VISIBLE_DEVICES`** (PLAN.md §6, gotcha 4):
the `gpu` nodes have no cgroup device isolation, `nvidia-smi -L` lists every GPU on the
node including other users', and a 4-GPU concurrency test on a 1-GPU allocation once
trespassed on three GPUs belonging to another user.

Two lanes inside *one* job cannot follow that rule literally — with `CUDA_VISIBLE_DEVICES`
untouched, both lanes would pick the first device and halve each other's throughput. So
`submit_campaign.sbatch` splits the allocation under a stricter rule:

- **Only ever hand a lane a device SLURM already named.** The script splits
  `$CUDA_VISIBLE_DEVICES` on commas and gives lane *L* element *L*. It subsets the
  allocated list and never extends it.
- **Refuse to run if the list is missing**, rather than guessing a device.
- **Refuse to run if SLURM allocated fewer devices than lanes.**

Both checks are fatal, not warnings. Picking a device we were not given is the exact
failure this guards against, and a job that dies at startup is strictly better than one
that quietly runs on someone else's GPU.

### How the work is divided

Each lane claims rows by `row_index % 2 == lane` — **500 rows each**, no shared queue, no
lock file, no master process. `campaign_worker.sh` then:

- **skips any run that already has `clusters.json`** (the last artifact written), so a
  resubmission resumes rather than repeats;
- writes `runs/<key>/pipeline.log` per run with a formulation header and a terminating
  `PIPELINE_ok` / `PIPELINE_FAILED`;
- appends progress to `logs/campaign/lane_NN.csv` and console output to
  `logs/campaign/lane_NN.log`;
- **stops gracefully 45 min before the wall clock** (`RESERVE_S=2700`) rather than being
  killed mid-stage.

Between that and `run_md.sh`'s stage-level resume (it skips any stage whose `.gro` already
exists, and `-cpi`/`-maxh` make a wall-clock kill recoverable mid-stage), the same command
is both "start" and "continue".

### Throughput estimate

6.1 ns of MD per run (0.1 heat + 2 NPT + 4 prod) at a median 11.6k atoms, on Quadro
RTX 6000s. Budget roughly **30–60 min per median run**, so 1000 runs / 2 concurrent
≈ **250–500 h wall**, call it **~15 resubmissions** at the 36 h cap. The heavy tail (the
42k-atom rows) is not split out — a run that does not finish simply resumes next time.

Packmol is ~290 s per system and runs inline on the node's CPUs while that lane's GPU
idles between runs. 1000 × 290 s ≈ 80 CPU-hours over 2 lanes, small next to the MD.

### Disk

Trajectories dominate. At 4 ns / 5000-step output and a median 11.6k atoms, budget
roughly **40–60 GB** of `prod.xtc` — which does **not** fit in the pi-chibueze /project
headroom alongside ht2000's. Hence `runs/` on scratch (2 TB) plus `archive_to_box.py`
pushing `prod.xtc` + `prod.tpr` to Box and, with `--prune`, freeing the local copy after
`rclone check` verifies the remote by hash.

**Box archiving must run from a login node** — compute nodes have no outbound internet,
the same constraint that keeps LigParGen fetching off them. `./maintain.sh` wraps the
whole housekeeping pass (refresh manifest → refresh results → status → archive) behind a
`flock` so two copies cannot tear `index/results.csv`.

---

## 9. Execution stages

- [x] **S0 — Scaffolding.** `MD/ht1000/` tree, scripts copied from ht2000 and re-pointed,
      `runs/` → scratch, Box destination created.
- [x] **S1 — Solvent pool.** `select_solvent_pool.py` → 278 oxygen-bearing solvents from
      ht2000's 500; 222 rejected, every one of them for `no_oxygen`.
- [x] **S2 — Dataset.** `build_dataset.py` → 1000 rows, 274 unique solvents, 40 salts,
      0.3–1.5 M, median 5.82 solvent per ion pair, none below 2.
- [x] **S3 — Ion FF library.** 2 cations + 20 anions in `ff_ions/`, all charge-checked,
      **40/40 grompp-validated**.
- [x] **S4 — Solvent FF library.** 272 copied + verified, 6 fetched from LigParGen (6/6
      ok, 0 failures). `verify_solvent_ff.py`: **278/278 composition-verified.**
- [x] **S5 — Index.** `build_run_manifest.py` → 1000/1000 rows keyed and runnable, 1000
      unique keys, 0 clamped, 0 missing force fields.
- [x] **S6 — Campaign driver.** `submit_campaign.sbatch` rewritten for one node / two
      GPUs / two lanes with the device-subsetting rule above.
- [x] **S7 — Launch.** `sbatch submit_campaign.sbatch`.
- [ ] **S8 — Run to completion.** Resubmit until `campaign_status.py` reports complete.
- [ ] **S9 — Collect.** `collect_results.py` → single results table.

**No pilot was run, deliberately.** ht2000 already validated the protocol, the force
fields and the analysis end to end on this exact software stack; ht1000 changes the
solvent *selection*, not the physics or the code path. The first submission is itself the
trial — it is resumable and costs nothing to cancel.

---

## 10. Open questions / decisions to revisit

1. **The pool is a subset, so ht1000 and ht2000 share solvents.** All 278 of ht1000's
   solvents also appear in ht2000's library, and 272 of them in its dataset. The
   *formulations* are disjoint (different concentration draws and salt pairings), but if
   the two campaigns are ever pooled, deduplicate on `(salt, solvent, concentration)`
   rather than assuming independence. Generating a genuinely new library was the original
   plan and was dropped for lack of a GPU node at build time — worth revisiting if a
   larger or more diverse pool is ever wanted.
2. **Thin functional families.** Amide (8), sulfoxide (2), phosphoryl (3). Any
   family-level conclusion about these rests on a handful of molecules. Stratified row
   sampling was considered and rejected in favour of uniform sampling (2026-09-21); the
   alternative is to regenerate a library with explicit quotas.
3. **Oxygen presence is not oxygen *availability*.** The filter requires an O atom, not a
   *coordinating* one — a sterically buried ether oxygen counts. If `cn_solv_O` comes back
   near zero for a subset of rows, that is the likely cause, and it is a real chemical
   result rather than a data defect.
4. **Nitrile nitrogen is not counted.** 30 pool solvents contain a nitrile, which
   coordinates alkali cations through N, but `analyze_solvation.py` only measures M–O.
   Their true first shell is under-counted. Extending the analysis to a general
   donor-atom set would fix it for both campaigns.
5. **`n_salt = 64` kept**, as in ht2000 and the original campaign, giving an 8× spread in
   system size (5k–43k atoms).
6. **4 ns production** is enough for solvation structure but not for converged transport
   (σ, η). A subset would need much longer runs if conductivity becomes a target.
7. **`prod` is NVT at the NPT-equilibrated box.** Density comes from `npt.edr`, not from
   production.

---

## 11. Change log

- **2026-09-21** — Campaign created and launched.
  - Replicated ht2000's workflow into `MD/ht1000/` with every writable path re-pointed;
    verified no ht1000 script writes outside its own tree, and that the only shared paths
    are read-only force-field sources.
  - **Oxygen requirement added** as the one chemical difference (§1): ht2000's library was
    44% O-free and 21% pure hydrocarbon, which makes the M–O RDF — the campaign's primary
    descriptor — undefined for those rows.
  - Pool taken as a subset of ht2000's library rather than generated afresh, at the user's
    direction (no GPU node available for Electrolyte-GPT at build time). Re-applying every
    inherited filter rejected exactly 222 molecules, all for `no_oxygen`, confirming the
    source library already satisfied the rest.
  - Dataset sized at 1000 rows to match the two GPUs available.
  - Solvent force fields **copied and re-verified** rather than re-fetched (§4), at the
    user's direction: the 272 shared molecules are already composition-verified, and
    re-asking LigParGen would re-expose all of them to the wrong-molecule failure mode for
    no gain. 6 genuinely new solvents fetched, 6/6 ok.
  - Campaign driver changed from ht2000's 3-task array to **one job on one node with two
    GPUs** (§8), because the QOS `node=4` cap — not the GPU count — is the binding
    constraint while ht2000 holds three nodes.
  - ht2000's "never set `CUDA_VISIBLE_DEVICES`" rule replaced with a stricter one that
    permits the two-lane split: subset the SLURM-provided device list, never extend it,
    and fail hard if it is missing or too short.

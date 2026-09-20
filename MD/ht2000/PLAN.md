# High-Throughput MD for the 2000-Electrolyte Dataset — Working Plan

**Status:** pilot validated (2/2 runs); force fields complete; ready for the full campaign
**Created:** 2026-09-19
**Owner:** eshiemogie@uchicago.edu
**Working dir:** `MD/ht2000/` (all new files live here)

This is a living document. Update it as decisions change — it is the single source of
truth for how the 2000-electrolyte high-throughput MD campaign is built and run.

---

## 0. Quick start

Everything runs with the `moleng` conda python; nothing needs to be installed.

```bash
PY=/scratch/midway3/eshiemogie/moleng/bin/python
cd /project/chibueze/stanley/projects/li-na-md/MD/ht2000/scripts

# one-time: ion force fields (22 species) + check them through grompp
$PY build_ion_ff.py
$PY validate_ions.py --cation Li && $PY validate_ions.py --cation Na

# solvent force fields -- LOGIN NODE ONLY (needs internet), serial by necessity
$PY fetch_solvent_ff.py                    # all 490, ~65 min (DONE: 487/490)
$PY fetch_solvent_ff.py --limit 10         # or just a slice

# pilot (already done; select_pilot.py picks coverage-maximising formulations)
$PY select_pilot.py -n 2
$PY launch_batch.py --from-csv ../index/pilot_rows.csv --build-only -j 4
$PY launch_batch.py --from-csv ../index/pilot_rows.csv --submit-only
# ...or, if you already hold an interactive GPU, skip the queue entirely:
bash run_local_batch.sh ../index/pilot_order.txt 1 16

# full campaign
$PY verify_solvent_ff.py --fix                        # ALWAYS audit first
$PY launch_batch.py --build-array --tag full          # packing -> caslake CPU nodes
$PY launch_batch.py --submit-only --array --throttle 12 --tag full   # MD -> gpu

# status / results at any time
$PY collect_results.py
squeue -u $USER
```

Re-running any of these is safe — every stage skips work that is already done.

---

## 1. Goal

Replicate the high-throughput solvation-MD workflow that was previously run on ~800 CE
electrolytes (archived in `MD/ce_solvation_md/`) for the new **2000-entry Li/Na
electrolyte dataset** at `datasets/electrolyte_dataset_2000.csv`.

Per formulation, produce the same artifacts the old campaign produced:

| artifact | meaning |
|---|---|
| `prod.xtc` + `prod.tpr` | 4 ns production trajectory |
| `solvation.json` | CN to solvent O, CN to anion, SSIP/CIP/AGG fractions, RDF first peak |
| `summary.json` | formulation spec + counts + box + atom count |
| `density.xvg` | NPT-average density |

**Order of work:** force fields → run builder → pilot → full 2000.
Status: force fields done and verified, pilot run and validated, full campaign not yet launched.

---

## 2. What the dataset looks like

`datasets/electrolyte_dataset_2000.csv` — 2000 rows, columns `salt, solvent, concentration`.

- **490 unique solvents**, all GPT-generated (`generated_500_organic_solvents.csv`).
  All parse in RDKit. 8–25 heavy atoms, MW 113–359. Elements: C,H,N,O,S,F,Cl,P,B.
- **40 unique salts** = 20 anions × {Li⁺, Na⁺}.
- **Concentration** 0.3–3.0 M (28-point grid).
- Single solvent per row (the old campaign had up to 5 — simpler here).

### 2.1 Known data issues (carry forward, do not silently ignore)

1. **Broken nitrate SMILES.** `[O-][N+](=O)=O` has N valence 5 and fails RDKit. Affects
   **100 rows** (LiNO3 + NaNO3). Correct form is `[O-][N+](=O)[O-]`.
   → Handled by a `SALT_FIXUPS` map in `scripts/common.py`. The source CSV is *not*
   edited (it lives outside `MD/`); the corrected SMILES is what gets simulated and is
   recorded in each run's `summary.json`.
2. **Chemically implausible solvents.** The generative model emitted species that are not
   viable battery solvents — e.g. acid chlorides (`CC(=O)Cl`), peroxides (`...OO1...`),
   aldehydes. The old workflow had `filter_li_metal.py` (non-PFAS / aprotic / no
   peroxide) for exactly this. **Decision: do not filter.** The task is to simulate the
   2000 as given; filtering is a separate scientific call for you to make. If you decide
   to drop them, the old `filter_li_metal.py` can be run over
   `index/solvent_manifest.csv` to produce an exclusion list — no filtering hook is
   built into the pipeline today.
3. **Na⁺ is new.** The old campaign was Li-only. Na⁺ parameters are taken from the same
   CL&P table as the old Li⁺ (see §4.2), so Li and Na are treated consistently.

4. **⚠ Most of the dataset is not in a dilute-electrolyte regime.** This is the most
   important thing on this page and it is a property of the dataset, not of the
   pipeline. The generated solvents are heavy (mean MW **242**, vs ~90 for DME or EC),
   but the concentrations span 0.3–3.0 M, a range calibrated for light carbonates and
   ethers. Combining the two leaves very little solvent per ion pair:

   | solvent molecules per salt formula unit | rows |
   |---|---|
   | < 1 (more salt than solvent) | ~110 |
   | < 2 (solvate ionic-liquid regime) | **849** |
   | < 4 (cannot fill a 4-coordinate first shell with solvent alone) | **1463 / 2000** |
   | median | **2.3** |

   By comparison the old campaign's typical 1 M LiFSI/DME run had **7.5**. Concretely,
   the pilot's 3.0 M NaFSI row packs 54 solvent molecules against 64 ion pairs.

   Consequence: for most rows `cn_solv_O` will be small and `agg` will be ≈1 — the
   solvation descriptors saturate and lose the ability to discriminate between
   formulations, which is presumably what they are wanted for. The runs are still
   valid MD; they are just mostly sampling concentrated/solvate electrolytes.

   **Not acted on** — the brief is to simulate the 2000 as given, and the
   high-concentration regime is itself relevant to Li-metal anodes (LHCEs). But the
   choice belongs to you; see §10.2 for the options.

---

## 3. The answer to "where do the .itp / .gro files come from?"

**The old campaign fetched them from the Yale LigParGen web service**, one HTTP round-trip
per unique SMILES. The scraper is preserved at

```
MD/ce_solvation_md/simulation_workflow-main/simulation_workflow-main/scripts/ligpargen_fetch.py
```

It POSTs a SMILES to `https://zarbi.chem.yale.edu/cgi-bin/results_lpg.py`, then downloads
the generated `.itp`, `.gro`, `.pdb`, `.lmp` (OPLS-AA topology with **CM1A-LBCC** charges).
You can see the fingerprint in every old run — `sol1.itp` literally begins:

```
; GENERATED BY LigParGen Server
; Jorgensen Lab @ Yale University
```

**Verified working from a Midway3 login node on 2026-09-19.** Test fetches succeeded for
DME plus 5 randomly drawn dataset solvents (5/5), ~15 s each, all four file types returned.
This is the mechanism we reuse.

Two consequences:

- **Zero reuse from the old library.** The 490 new solvents were checked against the 59
  solvent SMILES recoverable from the old run dirs: **overlap = 0**. Every solvent must be
  fetched fresh.
- **Fetching must happen on a login node.** Compute nodes have no outbound internet.
  So FF generation is a separate, pre-MD stage — it is not part of the SLURM job.

### 3.1 Solvent FF stage (`scripts/fetch_solvent_ff.py`)

- Key each solvent by `S<12-hex>` of the canonical SMILES → `ff_solvents/<key>.{itp,gro,pdb,smi}`.
- Idempotent: skips keys already complete, so it can be re-run after failures.
- **Strictly serial (`--workers 1`).** This is not politeness, it is correctness: the
  LigParGen backend runs BOSS in a shared server-side working directory, so concurrent
  submissions clobber each other. Measured on the pilot set:
  **7/10 spurious failures at 2 workers, 0/10 at 1.** The failing SMILES all succeeded
  immediately on a serial retry.

> ### ⚠ The failure mode that matters: LigParGen can return the *wrong molecule*
>
> Concurrency does not only cause spurious errors. Of the 3 requests that *appeared to
> succeed* during the 2-worker batch, **one came back with a different molecule's
> parameters** — a topology of `C19H38O` was returned for `CCCCC(C)C(CCCC(C)=O)CC(C)C`
> (`C16H32O`). The file is valid GROMACS input and grompp, mdrun and the analysis would
> all have accepted it without a murmur; it would simply have simulated the wrong
> chemistry. This was caught only because `setup_run.py` cross-checks the packed atom
> count against RDKit.
>
> Mitigation, now mandatory at three points:
> 1. `fetch_solvent_ff.py` verifies **element composition** (not just atom count — that
>    would miss an isomer swap) against the requested SMILES before keeping a download,
>    discards and retries on mismatch.
> 2. `verify_solvent_ff.py` re-audits the whole library; `--fix` deletes and re-fetches
>    bad entries. **Run this once after the full 490-solvent fetch, before S7.**
> 3. `setup_run.py` refuses to build a run dir from an .itp that does not match its
>    SMILES.
>
> Current library status: 10/10 verified (1 corrupted entry found and repaired).
- Fallback chain already in the scraper: `cm1abcc/opt=0` → `cm1a/opt=1` → `cm1a/opt=2`.
- Failures are written to `index/ff_failures.csv` with the returned HTML for diagnosis.
  Expect a handful (boron- and phosphorus-containing species are the likely casualties:
  2 B atoms and 7 P atoms across the whole set).
- **Measured rate: ~9 s/solvent serially** (10/10 of the pilot solvents succeeded), so
  490 solvents ≈ **75 min**, not the 2–3 h originally estimated.
- Charge matters: the original scraper hardcoded the server's `dropcharge` field to 0,
  which is why it only ever handled neutral solvents. `scripts/ligpargen.py` passes the
  real formal charge, which is what makes anions fetchable at all (see §4.2).

---

## 4. Force field composition

Matching the old campaign: **OPLS-AA (LigParGen/CM1A-LBCC) solvents + CL&P ions**, with
**ECC charge scaling q = 0.8 applied to the ions only** (solvent charges untouched — this
is standard electronic-continuum correction practice and is exactly what the old runs did:
`Li.itp` carries charge `+0.800`, and FSI⁻ in `an.itp` sums to `−0.800`).

### 4.1 The constraint that drove the design

**LigParGen cannot parameterise hypervalent inorganic anions.** Tested directly: BF4⁻,
FSI⁻, TFA⁻ and Beti⁻ are all rejected by the Yale server ("Problem found in the file
format") even with the correct net charge, because BOSS has no templates for them.
Organic anions (acetate, tosylate) *do* work, but only with `chargetype=cm1a` — the
default `cm1abcc` (CM1A-LBCC) is parameterised for neutral molecules only.

So no single source covers all 20 anions, and the library is assembled from four.

| source | path | provides |
|---|---|---|
| **OPLS-2009IL (Acevedo)** | `MD/Force Fields/2009IL FF Orlando/2009IL/` | 16 of the 20 anions, with matching `_atomtypes.itp` and PDB geometries |
| CL&P + fftool | `MD/Force Fields/clandp-master/` + vendored `tools/fftool` | TFA, Beti |
| Old runs | `ce_solvation_md/solvation/runs/*/an.itp` | FSI, exactly as previously simulated |
| LigParGen | web | OTs |

### 4.2 Per-anion sourcing (all 20 resolved)

| source | anions |
|---|---|
| **2009IL** | OAc(ACE), BF4, BNZ, Br, Cl, ClO4, DCA, HCOO, MS, NO3, TFSI(NTF2), PF6, PROP, SCN, TCM, OTf(TFO) |
| **fftool/CL&P** | TFA, Beti |
| **old runs** | FSI |
| **LigParGen** | OTs |

**Cations:** Li⁺ σ=2.126 Å, ε=0.07648 kJ/mol; Na⁺ σ=3.330 Å, ε=0.01160 kJ/mol — the
Aqvist parameters as tabulated in CL&P `il.ff`. Li⁺ reproduces the old campaign's
`Li.itp` exactly, so Li and Na are treated on the same footing.

### 4.3 Two corrections applied while assembling (`build_ion_ff.py`)

1. **Missing 1-4 `[pairs]`.** The 2009IL ITPs ship with `moleculetype/atoms/bonds/
   angles/dihedrals` and *no* `[pairs]` section. GROMACS does not synthesise the 1-4
   pair list from bonds — it only generates *parameters* for pairs that are listed — so
   using these files as-is would silently drop the OPLS 1-4 scaled LJ/Coulomb terms for
   every anion big enough to have them (OAc, PROP, BNZ, MS, TFSI, OTf, DCA, TCM).
   `[pairs]` is therefore generated from the bond graph (all atom pairs exactly 3 bonds
   apart). The small/monatomic anions are unaffected — they have no 1-4 pairs.
2. **Mixed `[atomtypes]` layouts.** 2009IL uses a 6-column form, LigParGen a 7-column
   form whose second field is a *bonded type* rather than an atomic number. Ion rows are
   normalised to the canonical `name at.nr mass charge ptype sigma epsilon`; LigParGen
   solvent rows are left verbatim with only the type name prefixed, which is exactly what
   the old campaign's `system.top` contains.

**FSI needed un-scaling.** The archived `an.itp` was stored *already* ECC-scaled (charges
summing to −0.8). It is divided back out so the whole library is uniformly unscaled; the
recovered values (F −0.13, S +1.02, N −0.66, O −0.53) are exactly the published CL&P FSI
charges, which confirms the correction.

### 4.4 Validation — done

`scripts/validate_ions.py` packs one cation + one anion and runs `gmx grompp` for every
combination. **40/40 pass** (20 anions × Li⁺ and Na⁺), and every ion's charge sums to
±1.000 before scaling. This is a one-time cost of 22 species, cached in `ff_ions/` and
reused by all 2000 runs.

---

## 5. Simulation protocol (reproduced from the old campaign)

Read directly off the archived `*.mdp` and `submit.sh`. Unchanged except for the cluster
adaptation in §6.

| stage | integrator | dt | length | notes |
|---|---|---|---|---|
| `em` | steep | — | 5000 steps | emtol 1000 |
| `nvt_heat` | md | 1 fs | 100 ps | annealed 50 K → 298.15 K, v-rescale |
| `npt` | md | 1 fs | 2 ns | C-rescale, 1 bar, τ_p = 2.0 |
| `prod` | md | 1 fs | **4 ns** | NVT, `nstxout-compressed = 5000` |

Common: `cutoff-scheme = Verlet`, `coulombtype = PME`, `rcoulomb = rvdw = 1.2 nm`,
`constraints = h-bonds` (LINCS), `tc-grps = System`, `tau_t = 0.5`, T = 298.15 K.

### 5.1 Box and composition math — reverse-engineered and verified

The original builder (`elytefm`) was **not** backed up, so the sizing rule was recovered
numerically from the archived `summary.json` files:

```
n_salt    = 64                                    (fixed)
V_target  = n_salt / (M · 0.6022)                 [nm³]   — molarity sets the final volume
L_build   = (V_target / 0.65)^(1/3)               [nm]    — pack at 65% of target density
V_salt    = n_salt · MW_salt / (ρ_salt · N_A)     ρ_salt = 1.70 g/cm³
n_solvent = (V_target − V_salt) · ρ_solv · N_A / MW_solv    ρ_solv = 1.00 g/cm³
```

**Validation against the archive:** `L_build` reproduces `box_build_nm` **exactly on 100%**
of single-solvent runs; `n_solvent` matches on **48/53**. The 5 misses are all degenerate
super-concentrated cases (M = 5.5 and 10.0) where the formula goes negative and the old
code clamped to a floor — outside our 0.3–3.0 M range, so immaterial here. A floor of 10
solvent molecules is kept for safety.

Packmol then packs into the (deliberately expanded) `L_build` box and NPT compresses to the
true density. Packmol's "ENDED WITHOUT PERFECT PACKING" exit is non-fatal by design.

### 5.2 Projected system sizes for the 2000

With `n_salt = 64`:

| | atoms |
|---|---|
| min | 2,134 |
| median | 6,400 |
| p90 | 20,394 |
| max | 46,337 |

`L_build` spans 3.79–8.17 nm; total ≈ 19M atom-runs. Median is smaller than the old
campaign (~11k atoms) but the tail is heavier, because low molarity + heavy solvents
inflate the box. **`n_salt` is exposed as a CLI flag** — dropping to 32 halves everything
(median 3,201 atoms, max 23,145) if the tail proves too expensive. Default stays at 64 for
fidelity with the old work.

---

## 6. Cluster adaptation (Midway3 — this is NOT the old cluster)

The archived `submit.sh` targets a different machine entirely (`/net/scratch2/qinanh/...`,
a self-built GROMACS, a `general` partition, `rclone`-to-Box). All of that is replaced:

| | old campaign | here (Midway3) |
|---|---|---|
| GROMACS | self-built 2025 at `/net/scratch2/qinanh/gromacs_gpu` | `module load gromacs/2025.3` — **binary is `gmx_mpi`**, CUDA 12.2, AVX-512 |
| GPU | A100 80 GB | **Quadro RTX 6000, 24 GB (Turing, sm_75)** |
| partition | `general`, `--gres=gpu:1` | `gpu` (11 nodes × 4 GPUs, 48 cores/node) |
| account | — | `pi-chibueze` |
| Python | `/home/qinanh/.conda/envs/byteff2` | `/scratch/midway3/eshiemogie/moleng` (rdkit 
+ MDAnalysis 2.10 + **packmol 21.2.1**) |
| packmol | conda `lammps` env | `/scratch/midway3/eshiemogie/moleng/bin/packmol` |
| archiving | `rclone` → Box, then delete xtc | keep local under `MD/ht2000/runs/` |

mdrun flags carried over: `-nb gpu -bonded gpu -pme gpu -update gpu -ntomp 8`
(EM stays on `-bonded cpu -pme cpu`, as in the original). Note the deliberate
absence of `-ntmpi` — see gotcha 0 below.

### Running without the queue

This session holds an interactive allocation (`midway3-0281`, **Tesla V100 16 GB**,
24 cores), and the `gpu` partition was backed up behind other users, so the pilot was
run there directly via `scripts/run_local_batch.sh <keyfile> <concurrency> <omp>`
rather than through `sbatch`. `run_md.sh` is stage-resumable and identical in both paths,
so work started locally can be finished by the queue and vice versa.

**Do not share a GPU between runs.** Measured on the V100 with ~3k-atom systems:
one run alone gives **818 ns/day**; three concurrent runs give **142 + 349 + 148 =
639 ns/day aggregate** — slower in total, with the 24-core node oversubscribed to load
30 (3 x 8 OpenMP threads plus the already-running job). A single small system already
saturates the GPU, so `run_local_batch.sh` should be used with concurrency **1**. This
does not affect the sbatch path, where each job gets its own GPU.

### Four cluster gotchas, all found the hard way

0. **Drop `-ntmpi 1` from every mdrun call.** Midway3's GROMACS is built against
   *real MPI*, not thread-MPI, and `-ntmpi` is a hard error there:
   `Setting the number of thread-MPI ranks is only supported with thread-MPI`.
   The old campaign's scripts all carry `-ntmpi 1` because that build was
   thread-MPI. Running `gmx_mpi` without `mpirun` already gives the single rank we
   want; `-ntomp` sets the OpenMP width. This would have failed all 2000 runs
   identically at the first mdrun.

1. **`gmx` does not exist on Midway3 — only `gmx_mpi`**, invoked without `mpirun` for
   these single-rank GPU runs.
2. **Never `module purge` before `module load gromacs/2025.3`.** The site module defines
   `SOFTPATH`, which the gromacs modulefile reads; purging makes the load die with
   `eval set [array get env SOFTPATH] / wrong # args`. This cost a probe job to find and
   is exactly the kind of thing that would have failed 2000 times silently.
3. **`module` is a shell function from the login profile.** `run_md.sh` re-execs itself
   under `bash -l` if it is not defined, so it works whether invoked from a batch script,
   a job array, or by hand.

> **Performance caveat:** the old campaign's timings (~680 ns/day, ~15 min/run) were on
> A100s. These are Quadro RTX 6000s — roughly 3–4× slower for GROMACS — so expect
> ~40–60 min per median run rather than ~15. The pilot exists partly to measure this;
> §9's estimate is revised from the pilot's actual `ns/day`.

---

## 7. Directory layout

```
MD/ht2000/
├── PLAN.md                  ← this file
├── scripts/
│   ├── common.py             ← SMILES/keys, SALT_FIXUPS, box math, FF integrity checks
│   ├── ligpargen.py          ← LigParGen client (charge-aware; fixes the old scraper)
│   ├── fetch_solvent_ff.py   ← batch solvent fetch, verified (login node, serial)
│   ├── verify_solvent_ff.py  ← audit/repair ff_solvents/ against its SMILES
│   ├── build_ion_ff.py       ← assemble ff_ions/ from 4 sources (one-time)
│   ├── validate_ions.py      ← push every cation/anion pair through grompp
│   ├── select_pilot.py       ← deterministic, coverage-driven pilot selection
│   ├── setup_run.py          ← build ONE run dir (packmol + top + mdp + submit.sh)
│   ├── run_md.sh             ← the 4-stage GROMACS protocol; stage-resumable
│   ├── launch_batch.py       ← build + submit (individually or as a job array)
│   ├── analyze_solvation.py  ← reproduce solvation.json from prod.xtc
│   └── collect_results.py    ← gather all summary/solvation json → results table
├── ff_solvents/  S<hash>.{itp,gro,pdb,smi}
├── ff_ions/      Li.itp Na.itp <ANION>.itp/.pdb + ions_manifest.json
├── runs/         el<hash>/ ... (one dir per formulation, mirrors old layout)
├── index/        formulation_index.csv, ff_failures.csv, batch manifests
└── logs/
```

Run keys: `el<12-hex>` of `canonical(salt)|canonical(solvent)|M` — deterministic, so a row
always maps to the same run dir and the campaign is restartable.

---

## 8. Execution stages

- [x] **S0 — Scaffolding.** Dirs + `common.py`.
- [x] **S1 — Ion FF library.** 2 cations + 20 anions in `ff_ions/`, all charge-checked,
      **40/40 grompp-validated**.
- [x] **S2 — Solvent FF for the pilot.** 10/10 fetched.
- [x] **S3 — Run builder.** `setup_run.py`; 10 pilot dirs built, all with matching atom
      counts and net charge exactly 0.000.
- [x] **S4 — Pilot batch.** Scaled to 2 runs at the user's request; both completed all
      four stages on the local V100. Results below.
- [x] **S5 — Analysis.** `analyze_solvation.py` reproduces the old schema; Li/Na RDF
      peaks and coordination numbers are physically correct.
- [x] **S6 — Full solvent FF fetch.** 487/490 fetched (65 min, serial); 3 failures are
      boron/phosphate esters LigParGen cannot do.
- [x] **S6b — Audit the library.** 487/487 composition-verified.
- [ ] **S7 — Full campaign.** 2000 runs, one throttled job array.
- [ ] **S8 — Collect.** `collect_results.py` → single results table.

### Build cost (measured)

Packmol dominates run-dir construction: **~290 s per system** (the old campaign's logs
show ~256 s for the same protocol, so this is inherent, not a regression). For 2000 runs
that is ≈**160 CPU-hours of packing**, which is why `launch_batch.py` parallelises the
build (`-j`) and separates `--build-only` from `--submit-only`. Run the full build as a
CPU-partition array job rather than on a login node.

Topology generation is decoupled from packing: rebuilding a run dir reuses an existing
`mixture.pdb` whose atom count still matches, so force-field fixes do not re-pay the
packmol cost.

### Pilot success criteria (S4)

1. All 4 GROMACS stages exit 0.
2. `prod.xtc` exists and has ≥ 100 frames.
3. NPT density is physically sane (roughly 0.7–2.0 g/cm³).
4. System net charge ≈ 0 (|q| < 0.05 e) — grompp warning check.
5. `solvation.json` produced with Li⁺/Na⁺ CN in a believable range.
6. Wall time per run recorded, to extrapolate the 2000-run cost.

### Pilot results

Scaled down from 10 runs to **2** on 2026-09-19 at the user's request; the other eight
run dirs were deleted. The two kept still span both cations, both concentration
extremes and both box-size extremes:

**Both runs completed all four stages plus analysis.**

| | NaFSI 3.0 M | LiFSI 0.3 M |
|---|---|---|
| run key | `el7378a10bbeb1` | `el4f3bbcd28325` |
| atoms | 2,584 | 31,132 |
| prod ns/day (V100) | 646 | 232 |
| wall time, 6.1 ns | ~13 min | ~40 min |
| density (g/cm³) | 1.397 | 0.983 |
| RDF M–O peak (Å) | **2.43** | **2.11** |
| shell cutoff (Å) | 3.25 | 3.27 |
| cn_solv_O | 1.09 | 3.65 |
| cn_anion | 4.32 | 1.31 |
| ssip / cip / agg | 0.000 / 0.035 / **0.965** | 0.228 / 0.626 / 0.146 |

Three things this establishes:

1. **The physics is right.** The Na–O first peak (2.43 Å) is correctly larger than the
   Li–O one (2.11 Å) — the old campaign measured 1.96 Å for Li–O with *ether* oxygens,
   and 2.11 Å is the expected small lengthening for carbonyl O under ECC q = 0.8.
   Densities are sensible. Nothing here was tuned to produce these numbers.
2. **The coordination numbers are robust to the cutoff.** The detected cutoff (~3.26 Å
   for both) sits on a genuine plateau of the running CN: moving it from 2.6 → 3.3 Å
   changes CN by only 3.43 → 3.63 (Li) and 1.00 → 1.08 (Na). So the automatic
   first-minimum choice is not doing anything load-bearing, which is what we want when
   the same code has to run unattended over 2000 systems.
3. **The concentration problem from §2.1 is real and visible.** The dilute LiFSI run
   gives a well-spread SSIP/CIP/AGG distribution (0.23 / 0.63 / 0.15) — informative.
   The 3 M NaFSI run, with 54 solvent molecules to 64 ion pairs, is 96.5% aggregates
   with SSIP = 0.000 — saturated, and indistinguishable from any other concentrated
   row. Expect the latter behaviour for roughly 1500 of the 2000 formulations.

### Solvent force-field library — complete

477 fetched in 65 min (serial), plus 10 from the pilot. **487 of 490 solvents available,
and `verify_solvent_ff.py` reports 487/487 composition-verified.**

The 3 failures are exactly the predicted casualties — LigParGen/BOSS has no templates
for boron esters or hypervalent phosphate esters:

| solvent | rows affected |
|---|---|
| `CC(C)OP(=O)(OC(C)C)OC(C)C` (triisopropyl phosphate) | 5 |
| `CCC1CCC(B(OC)OC)C1` (boronate ester) | 5 |
| `CCCCCCCC1COB(CC)O1` (dioxaborolane) | 4 |

**14 of 2000 rows (0.7%) are blocked; 1986 are runnable.** If those 14 matter, the
options are GAFF/antechamber or OpenFF for those three species; otherwise drop them.

### Pilot selection

10 rows chosen to *stress the pipeline*, not at random: both cations, a spread of anion
types (perfluorinated imide, simple monatomic, carboxylate, nitrile-based, oxoanion), and
both concentration extremes — so failures surface now rather than at run 1500.

---

## 9. Scaling to 2000

### The binding constraint: the `gpu` QOS

```
MaxJobsPU = 12      MaxTRESPU = cpu=192, gres/gpu=16, node=4      MaxWall = 1-12:00:00
```

**At most 12 of our jobs run at once**, no matter how many are queued. That decides the
shape of the campaign:

- Submit **one throttled job array** (`launch_batch.py --array --throttle 12`), not 2000
  independent `sbatch` calls — the latter would park ~1990 pending jobs in a shared
  queue for no gain. The array driver skips any task whose `solvation.json` already
  exists, so resubmitting the same array resumes the campaign.
- Throughput estimate: 6.1 ns of MD per run (0.1 heat + 2 NPT + 4 prod). At a median
  ~6.4k atoms that is roughly 10–25 min per run on one GPU, so
  **2000 runs / 12 concurrent × ~20 min ≈ 2.5–4 days** of wall time, longer with the
  heavy tail and queue contention.
- Nodes are `midway3-[0277-0286,0294]`, 4 GPUs and 48 cores each; 8 cores per run leaves
  the node comfortably shareable.

### Practical notes

- **Split the heavy tail.** The 46k-atom, low-molarity systems should go in their own
  array with a longer `--hours`; the default 8 h suits the median case. `run_md.sh` is
  stage-resumable, so a task killed at the wall clock resumes from its last finished
  stage when the array is resubmitted.
- **Build first, separately.** `--build-only -j 16` on a CPU node; packing 2000 systems
  is ~160 CPU-hours and has no business running on a login node or holding a GPU.
- Disk: trajectories dominate. At 4 ns / 5000-step output and median 6.4k atoms, budget
  roughly 10–50 MB per `prod.xtc` → **~20–100 GB total**. Check the quota on
  `/project/chibueze` before S7; if tight, thin `nstxout-compressed`.
- Unlike the old campaign, nothing is pushed to Box and no trajectories are deleted —
  `prod.xtc` stays next to its `prod.tpr` so the analysis can be re-run.

---

## 10. Open questions / decisions to revisit

1. **Keep `n_salt = 64`?** Faithful to the old work, but produces a 20× spread in system
   size. Alternative: fix the box (~5 nm) and let `n_salt` float with concentration.
   *Current decision: keep 64, revisit after the pilot timing.*
2. **What to do about the concentration/MW mismatch (§2.1 item 4)?** Options, roughly in
   increasing order of intervention:
   - **(a) Nothing.** Run all 2000 as specified and report the solvent:salt ratio
     alongside each result, so the saturated points can be excluded at analysis time.
     *This is the current default.*
   - **(b) Analyse the dilute subset separately.** ~537 rows have ≥4 solvent per salt
     and behave like conventional electrolytes; treat the rest as a second population.
   - **(c) Re-express concentration.** Regenerate the dataset in molality or salt mole
     fraction rather than molarity, so heavy and light solvents are compared at matched
     stoichiometry rather than matched mol/L.
   - **(d) Cap molarity per solvent** so that every row keeps at least ~4 solvents per
     ion pair.

   (c) and (d) mean changing `datasets/electrolyte_dataset_2000.csv`, which is outside
   `MD/` and outside this task's scope — flagging, not doing.
3. **Solvent density is assumed 1.0 g/cm³ for everything** (as the old campaign did).
   For MW-240 branched esters/ketones the true value is ~0.85–1.05, so the initial box
   is off by up to ~15%. NPT fixes the density, but it shifts the *composition* slightly
   because `n_solvent` is derived from it. Low priority; fixable with a
   group-contribution density estimate if it ever matters.
4. **Filter implausible solvents?** Currently flagged but not dropped. Your call.
5. **Na⁺ ECC scaling.** q = 0.8 is well-established for Li⁺ in carbonates/ethers; for Na⁺
   the same 0.8 is a reasonable default but is an assumption worth stating in any writeup.
6. **4 ns production** is enough for solvation structure (what we want) but *not* for
   converged transport (σ, η). If conductivity is a downstream target, a subset will need
   much longer runs.
7. **`prod` is NVT at the NPT-equilibrated box** (as in the old work). Density therefore
   comes from `npt.edr`, not from production.

---

## 11. Change log

- **2026-09-19 (a)** — Initial plan. Explored `ce_solvation_md`; confirmed LigParGen as
  the .itp/.gro source and verified it works from Midway3; reverse-engineered and
  validated the box/count math; inventoried ion FF sources; found the nitrate SMILES bug;
  confirmed zero solvent overlap with the old library; mapped the Midway3 software stack
  (`gmx_mpi`, `moleng` env, `gpu` partition, `pi-chibueze` account).
- **2026-09-19 (b)** — Force fields built and validated. Findings that changed the plan:
  - LigParGen **rejects hypervalent inorganic anions**, so the ion library is sourced
    from four places rather than mostly LigParGen (§4.1–4.2). Discovered the complete
    OPLS-2009IL library already vendored under `MD/Force Fields/2009IL FF Orlando/`,
    which covers 16 of the 20 anions with matching atomtypes and PDBs.
  - The original scraper hardcoded net charge to 0; fixed in `scripts/ligpargen.py`.
  - LigParGen **must be called serially** — concurrency produces spurious failures.
  - 2009IL ITPs omit `[pairs]`; generated from the bond graph (§4.3).
  - Archived FSI was already ECC-scaled; un-scaled for the library.
  - Per-molecule charge neutralisation added: LigParGen's 4-decimal rounding summed to
    −0.092 e over a 924-molecule box and tripped grompp's Ewald-with-net-charge warning.
  - Packmol measured at ~290 s/system → build must be parallelised (§8).
  - **LigParGen can silently return another molecule's parameters** under concurrency
    (1 of 10 pilot solvents). Composition verification is now mandatory at fetch,
    audit and build time — see the callout in §3.1. This is the single most dangerous
    thing found so far, because nothing downstream would have flagged it.
  - Quantified the concentration/molecular-weight mismatch in the dataset (§2.1 item 4):
    1463 of 2000 rows have fewer than 4 solvent molecules per ion pair.
  - Cluster gotchas found with a probe job: `module purge` breaks the gromacs
    modulefile; GPUs are Quadro RTX 6000, not A100 (§6).
  - `gpu` QOS caps this user at 12 running jobs / 16 GPUs → the full campaign must be
    one throttled job array (§9).

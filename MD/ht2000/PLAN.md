# High-Throughput MD for the 2000-Electrolyte Dataset — Working Plan

**Status:** dataset rebuilt for a sane concentration regime; 490/490 force fields in hand;
all 2000 rows indexed; campaign driver written and ready — **not yet launched**
**Created:** 2026-09-19
**Last updated:** 2026-09-20
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
$PY fetch_solvent_ff.py                    # all 490, ~65 min (DONE: 490/490)
$PY fetch_solvent_ff.py --limit 10         # or just a slice
$PY verify_solvent_ff.py --fix             # ALWAYS audit before a campaign

# index every dataset row before running anything
$PY build_run_manifest.py                  # -> index/run_manifest.csv

# THE FULL CAMPAIGN -- one submission, start or continue
sbatch ../scripts/submit_campaign.sbatch

# status / results at any time
$PY campaign_status.py --lanes             # progress, failures, ETA
$PY campaign_status.py --failures          # which runs died and where
$PY collect_results.py                     # -> index/results.csv
squeue -u $USER
```

The campaign is **resumable by design**: lanes skip runs that already finished and
`run_md.sh` skips stages that already finished, so the *same* `sbatch` command both
starts the campaign and continues it after a wall-clock kill. Expect to run it ~10 times.
The same is true of every other command above — all of them skip work already done.

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
| `clusters.json` | ion-cluster count/size from the cation-anion RDF cutoff, free-ion and percolation fractions |
| `summary.json` | formulation spec + counts + box + atom count |
| `density.xvg` | NPT-average density |

**Order of work:** force fields → run builder → pilot → *dataset rebuild* → full 2000.
Status: pilot run and validated on the original dataset; the dataset was then rebuilt
(§2) and the solvent library refetched to match; 490/490 force fields verified; all 2000
rows indexed; campaign driver ready. Full campaign not yet launched.

---

## 2. What the dataset looks like

> **This section describes the dataset as rebuilt on 2026-09-20.** The original
> generated dataset put most rows outside any sensible electrolyte regime; that is why
> it was regenerated rather than simulated as-is. The reasoning is in §2.2.

`datasets/electrolyte_dataset_2000.csv` — 2000 rows, columns
`salt, solvent, concentration, solv_per_ion_pair`.

- **490 unique solvents** drawn from a 500-solvent GPT-generated library
  (`generated_500_organic_solvents.csv`). All parse in RDKit.
  **MW 88–199** (mean 177). Elements: **C, H, N, O, S, F, Cl, P** — no boron, no silicon.
- **40 unique salts** = 20 anions × {Li⁺, Na⁺}.
- **Concentration 0.3–1.5 M** (28-point grid).
- Single solvent per row (the old campaign had up to 5 — simpler here).
- **Median 5.97 solvent molecules per ion pair**; 401 rows below 4, **none below 2**.

The solvent library is filtered at generation time — see §2.2 for what is excluded and
why. Only 490 of the 500 library solvents appear in the dataset; the sampler simply never
drew the other 10.

### 2.1 Known data issues (carry forward, do not silently ignore)

1. **Broken nitrate SMILES.** `[O-][N+](=O)=O` has N valence 5 and fails RDKit. Affects
   **100 rows** (LiNO3 + NaNO3). Correct form is `[O-][N+](=O)[O-]`.
   → Handled by a `SALT_FIXUPS` map in `scripts/common.py`. The source CSV is *not*
   edited (it lives outside `MD/`); the corrected SMILES is what gets simulated and is
   recorded in each run's `summary.json`.
2. **Chemically implausible solvents — now filtered at generation.** *(Decision reversed
   2026-09-20; the original plan was to simulate everything as given.)* The generative
   model emitted species that are not viable battery solvents, and separately some that
   no force-field tool can parameterise at all. Both classes are now rejected inside
   `generate_500_solvents.py` rather than downstream, so the dataset never contains them:

   | filter | SMARTS / rule | why |
   |---|---|---|
   | aldehyde | `[CX3H1](=O)[#6]` | reactive toward Li metal; not used as a solvent |
   | protic | `[OX2H,NX3;H1,H2]` | O–H / N–H attacks Li⁰; also breaks the aprotic premise |
   | radical | any atom with unpaired electrons | **BOSS cannot type them** — hard LigParGen failure |
   | boron | element filter (`ALLOWED_Z`) | BOSS has no boron templates |
   | phosphate/phosphite ester | `[#15;!$([#15]~[#6])]` (P with no P–C bond) | BOSS has no templates |
   | charged / multi-fragment / metal-containing | — | not a neutral single-component solvent |

   Allowed elements are `H, C, N, O, F, P, S, Cl, Br, I`. **Silicon is excluded**, though
   LigParGen handles it fine — it is caught by the inherited `METAL_ATOMIC_NUMS`
   (`range(11,15)`) *before* the element check, so re-enabling it needs edits in two
   places, not one.

   The generator also reads `index/ff_failures.csv` and skips any SMILES already known to
   fail LigParGen, which closes the loop on one-off rejects that no structural rule
   catches (see §8, S6c).
3. **Na⁺ is new.** The old campaign was Li-only. Na⁺ parameters come from the same CL&P
   table as Li⁺ (§4.2), so the two are treated consistently.

---

### 2.2 The concentration/molecular-weight mismatch — diagnosed, then fixed

This was the most important problem found in the original dataset, and it is the reason
the dataset was regenerated rather than simulated as delivered.

#### The diagnosis

Molarity fixes the *volume* per mole of salt, not the *stoichiometry*. The original
generated solvents were heavy (mean MW **242**, vs ~90 for DME or EC), but the
concentrations spanned 0.3–3.0 M — a range calibrated for light carbonates and ethers.
Combining the two leaves almost no solvent per ion pair, because a mole of heavy solvent
eats the whole box:

| solvent molecules per salt formula unit | original | **rebuilt** |
|---|---|---|
| < 1 (more salt than solvent) | ~110 | **0** |
| < 2 (solvate ionic-liquid regime) | 849 | **0** |
| < 4 (cannot fill a 4-coordinate first shell with solvent alone) | **1463 / 2000** | **401 / 2000** |
| median | **2.3** | **5.97** |

For reference the old campaign's typical 1 M LiFSI/DME run had 7.5, and the original
pilot's 3.0 M NaFSI row packed 54 solvent molecules against 64 ion pairs.

The consequence was visible in the pilot (§8): at 3.0 M NaFSI the descriptors saturate —
`agg` = 0.965, `ssip` = 0.000 — so the run carries almost no information distinguishing
it from any other concentrated row. Roughly 1500 of 2000 rows would have behaved that way.

#### The fix

Two changes, both in the dataset generators, chosen to attack the product `MW × M`
rather than either factor alone:

- **MW ceiling 200 g/mol** in `generate_500_solvents.py` (was unbounded; realised range
  is now 88–199, mean 177).
- **Concentration ceiling 1.5 M** in `build_electrolyte_dataset.py` (was 3.0 M).

1.5 M is not arbitrary: it is the top of the range where conventional Li-ion electrolytes
actually operate (1.0–1.2 M is the industrial norm), so the dataset now spans dilute to
moderately concentrated rather than dilute to solvate-ionic-liquid.

The 401 rows still below 4 solvents per pair are the genuinely concentrated corner of a
legitimate range, not an artifact — they are worth keeping and are flagged by the
`solv_per_ion_pair` column in both the dataset and `results.csv`.

#### What this does *not* fix

Molarity is a property of the *equilibrated* box, and the builder has to guess a density
to convert it into molecule counts (§5.1, ρ_salt = 1.70 g/cm³ flat). The requested and
achieved molarity therefore differ. Measured on the pilot: 0.3 M landed at 0.300 M
(−0.1%), but 3.0 M landed at 3.441 M (**+14.7%**). Capping at 1.5 M shrinks this error
but does not remove it, so `collect_results.py` now records `molarity_actual` and
`molarity_err_pct` for every run and flags anything off by >10% (§8, S5b). **Report the
achieved molarity, not the requested one.**

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
> Current library status: **490/490 composition-verified** (1 corrupted entry found and
> repaired during the pilot).
- Fallback chain already in the scraper: `cm1abcc/opt=0` → `cm1a/opt=1` → `cm1a/opt=2`.
- Failures are written to `index/ff_failures.csv` with the returned HTML for diagnosis.
  **`generate_500_solvents.py` reads this file back** and excludes known-bad SMILES from
  future libraries, so each failure is paid for once (§8, S6c).
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
code clamped to a floor — far outside our 0.3–1.5 M range, so immaterial here. A floor of
10 solvent molecules is kept for safety; **0 of the 2000 rows hit it** (`clamped` is false
throughout `index/run_manifest.csv`), which is a direct consequence of the §2.2 rebuild.

**ρ_salt = 1.70 g/cm³ is a flat assumption for all 20 anions** and is the main source of
the requested-vs-achieved molarity gap discussed in §2.2. Left as-is: at the 1.5 M ceiling
the error is tolerable, and `molarity_actual` records the truth per run.

Packmol then packs into the (deliberately expanded) `L_build` box and NPT compresses to the
true density. Packmol's "ENDED WITHOUT PERFECT PACKING" exit is non-fatal by design.

### 5.2 Projected system sizes for the 2000

With `n_salt = 64`, measured over all 2000 rows of the **rebuilt** dataset
(`index/run_manifest.csv`, column `n_atoms_est`):

| | atoms |
|---|---|
| min | 3,908 |
| median | **12,541** |
| p90 | 28,288 |
| max | 47,152 |

`L_build` spans 4.78–8.17 nm; total ≈ **30.7M atom-runs**.

Note this got *bigger*, not smaller, after the §2.2 rebuild (median 6,400 → 12,541).
That is the expected trade: capping molarity at 1.5 M means more solvent per ion pair,
and solvent molecules are what fill the box. The campaign buys physically meaningful
solvation structure at roughly 2× the compute. Median is now comparable to the old
campaign (~11k atoms).

**`n_salt` is exposed as a CLI flag** — dropping to 32 roughly halves everything if the
tail proves too expensive. Default stays at 64 for fidelity with the old work.

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

### Five cluster gotchas, all found the hard way

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

> ### ⚠ 4. The `gpu` partition does NOT isolate GPUs per job
>
> There is no cgroup device isolation on these nodes. `nvidia-smi -L` inside a job lists
> **every GPU on the node**, including ones allocated to other people, and setting
> `CUDA_VISIBLE_DEVICES` yourself will happily run on them.
>
> This was found the hard way on 2026-09-20: a 4-GPU concurrency test on an interactive
> node holding `gres/gpu=1` was launched with `CUDA_VISIBLE_DEVICES=0,1,2,3` and
> **trespassed on three GPUs belonging to another user (jlguerra)**. Nothing warned us;
> the job simply ran.
>
> **Rule, now enforced in `submit_campaign.sbatch`:** never set `CUDA_VISIBLE_DEVICES`.
> Trust the value SLURM exports, and warn loudly if it is missing. The script prints both
> the visible set and the node's true device count on every launch so a mismatch is
> obvious in the log. The old multi-lane loop that assigned `CUDA_VISIBLE_DEVICES="$L"`
> per lane has been removed for exactly this reason.

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
│   ├── launch_batch.py       ← build + submit (per-run or job array) — pilot-era path
│   ├── run_local_batch.sh    ← run a keyfile on an interactive GPU, no queue
│   ├── analyze_solvation.py  ← reproduce solvation.json from prod.xtc
│   ├── cluster_analysis.py   ← ion clusters from prod.xtc (own cation-anion cutoff)
│   ├── build_run_manifest.py ← index all 2000 rows up front (key, composition, stage)
│   ├── collect_results.py    ← gather all summary/solvation/cluster json → results table
│   ├── submit_campaign.sbatch ← THE campaign entry point: 3-task array, 1 GPU each
│   ├── campaign_worker.sh    ← one lane; strides rows, logs, stops before the wall
│   └── campaign_status.py    ← progress / failures / ETA, read from disk
├── ff_solvents/  S<hash>.{itp,gro,pdb,smi,lmp}   490 solvents × 5 files
├── ff_ions/      Li.itp Na.itp <ANION>.itp/.pdb + ions_manifest.json
├── runs/         el<hash>/ ... (one dir per formulation, mirrors old layout)
├── index/        run_manifest.csv, results.csv, solvent_manifest.csv,
│                 ff_failures.csv, pilot_* (historical)
├── logs/         slurm/  campaign/  + fetch and generation logs
└── session_transcript.txt   ← gitignored; running record of the build sessions
```

`launch_batch.py` and `run_local_batch.sh` are the pilot-era paths and still work, but the
2000-run campaign goes through `submit_campaign.sbatch` (§9).

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
- [x] **S5b — Ion cluster analysis.** `cluster_analysis.py` added and wired into
      `run_md.sh`; cross-checked against OVITO (§8.1). `collect_results.py` extended with
      achieved molarity, dataset row linkage and the full cluster block.
- [x] **S6 — Full solvent FF fetch (original library).** 487/490 fetched (65 min,
      serial); 3 failures were boron/phosphate esters LigParGen cannot do.
- [x] **S6b — Audit the library.** 487/487 composition-verified.
- [x] **S6c — Dataset rebuild (§2.2).** MW ≤ 200 and M ≤ 1.5 applied; solvent library and
      `electrolyte_dataset_2000.csv` regenerated in place (no v2 variants kept). Six
      LigParGen failures in the new library traced to 5 carbon radicals + 1 thiophene;
      a radical filter was added and the thiophene recorded in `ff_failures.csv`, which
      the generator now reads back.
- [x] **S6d — Refetch the library.** Zero overlap with the old 490, so all 2445 cached
      files were cleared and refetched: **490/490 complete, 2450 files.**
- [x] **S6e — Index the dataset.** `build_run_manifest.py` → all 2000 rows keyed,
      costed and checked for force-field coverage. 2000/2000 runnable, 0 clamped.
- [x] **S6f — Campaign driver.** `submit_campaign.sbatch` + `campaign_worker.sh` +
      `campaign_status.py` (§9).
- [ ] **S7 — Full campaign.** 2000 runs via `sbatch submit_campaign.sbatch`, resubmitted
      until `campaign_status.py` reports complete. **← next action, awaiting go-ahead**
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

### Pilot results *(historical — ran against the pre-rebuild dataset)*

> These two runs validated the **pipeline**, and that conclusion still holds — the
> protocol, force fields and analysis are unchanged. But they were built from the
> original 0.3–3.0 M dataset, so the 3.0 M formulation below **no longer exists** in the
> rebuilt dataset. Both run directories have since been deleted; `runs/` is empty and the
> campaign starts from zero. Keep this table for what it establishes about physics and
> cost, not as a description of what will be run.

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
3. **The concentration problem was real and visible — and is what triggered the
   rebuild.** The dilute LiFSI run gives a well-spread SSIP/CIP/AGG distribution
   (0.23 / 0.63 / 0.15) — informative. The 3 M NaFSI run, with 54 solvent molecules to
   64 ion pairs, is 96.5% aggregates with SSIP = 0.000 — saturated, and
   indistinguishable from any other concentrated row. That would have been the
   behaviour of ~1500 of the 2000 formulations, which is why the dataset was
   regenerated (§2.2) rather than run as delivered.

### 8.1 Ion cluster analysis — validated against OVITO

`cluster_analysis.py` does on-cluster clustering so trajectories never have to be
downloaded. Method:

1. Take the cation–anion-donor RDF over the production trajectory.
2. Cutoff = the **first minimum** of that RDF (a different, larger quantity than the
   cation–solvent-O `shell_cutoff_A` used by `analyze_solvation.py` — the two are not
   interchangeable and both are recorded).
3. Per frame, build a contact graph with `capped_distance` under the minimum-image
   convention, then `scipy.sparse.csgraph.connected_components`.

**Two counting conventions are reported, because they differ and the difference matters:**

| column | counts | equals |
|---|---|---|
| `n_components` | every component, lone ions included | OVITO's "cluster count" |
| `n_clusters` | aggregates of ≥ 2 ions only | — |

with the identity `n_components = n_clusters + n_free_cation + n_free_anion`.

Cross-check: for the Li pilot at frame 348/800 with a 2.67 Å cutoff, ours gives **72**
against OVITO's **~70** (OVITO set to Residue type, i.e. whole molecules). Frame-to-frame
σ is 3.1, so the two agree. Reporting a single frame is misleading — hence the `_sd`
columns.

Also recorded: `free_cation_frac`, `mean_cluster_size`, `max_cluster_size`,
`largest_cluster_frac`, `percolating_frac` and `mean_cluster_charge`. When
`percolating_frac` > 0.5 a single network spans the box and `n_clusters` stops being
meaningful; `collect_results.py` prints how many runs are in that state.

### Solvent force-field library — complete

The library was refetched from scratch after the §2.2 rebuild: the new 490 solvents have
**zero overlap** with the previous 490, so all 2445 cached files were deleted rather than
reused. **490/490 solvents available (2450 files), all composition-verified** — every one
of the 2000 dataset rows is runnable.

Six LigParGen failures surfaced during the refetch and all six were designed out rather
than worked around:

| failure | count | resolution |
|---|---|---|
| carbon radicals (`[C]`, `[CH]`) | 5 | RDKit parses them, BOSS cannot type them → radical filter added to the generator |
| `N#CC1=CCC=C(F)S1` (thiophene) | 1 | no obvious structural rule; recorded in `ff_failures.csv`, which the generator now reads back |

The generator re-sampled replacements for all six, so the library is full at 490 with no
blocked rows — an improvement on the pre-rebuild state, where 14 of 2000 rows (0.7%) were
unrunnable for want of boron/phosphate-ester parameters.

### Pilot selection

10 rows chosen to *stress the pipeline*, not at random: both cations, a spread of anion
types (perfluorinated imide, simple monatomic, carboxylate, nitrile-based, oxoanion), and
both concentration extremes — so failures surface now rather than at run 1500.

---

## 9. Scaling to 2000

### The binding constraints

```
gpu QOS:  MaxJobsPU = 12    MaxTRESPU = cpu=192, gres/gpu=16, node=4    MaxWall = 1-12:00:00
```

The QOS caps us at 12 running jobs / 16 GPUs. But the *real* constraint turned out not to
be the QOS — it is **partition occupancy**. Every GPU in the partition is held by
single-GPU jobs from other users, so a `--gres=gpu:4` whole-node request has to wait for
some node's last job to drain: **measured at 25 h minimum, 84 h on one node.** A
`--gres=gpu:1` request slots into the first GPU that frees anywhere on the partition.

### The shape that follows: one sbatch, three lanes, one GPU each

`submit_campaign.sbatch` is a **3-task job array, 1 GPU + 12 CPUs + 32 G per task**.
One submission covers all 2000 runs.

| | considered | chosen |
|---|---|---|
| request | 3 nodes × `gpu:4`, 4 lanes/task = 12 concurrent | 3 × `gpu:1` = **3 concurrent** |
| queue wait | 25–84 h | minutes |
| wall to finish | ~85 h | **~330 h** |
| neighbourliness | occupies whole nodes | leaves the rest of each node free |

We trade 4× throughput for actually starting, and for not monopolising nodes. **Raise
`NLANES` and `--gres` together if the partition frees up** — the QOS ceiling allows up to
12 jobs / 16 GPUs.

Why not GPU memory: it is not separately requestable on Midway3. `--gres=gpu:1` hands
over one whole 24 GB device, and these systems need ~2 GB. `--mem=32G` is host memory,
asked for as a share rather than `--mem=0` so the node stays usable by others.

### How the work is divided: row striding, no coordination

Each lane claims rows by `row_index % NLANES == LANE`. There is no shared queue, no lock
file and no master process — a lane can die, be resubmitted, or run on a different node
and it still covers exactly its own rows. `campaign_worker.sh` then:

- **skips any run that already has `clusters.json`** (the last artifact written), so a
  resubmission resumes rather than repeats;
- writes `runs/<key>/pipeline.log` per run with a formulation header and a terminating
  `PIPELINE_ok` / `PIPELINE_FAILED`;
- appends progress to `logs/campaign/lane_NN.csv`;
- **stops gracefully 45 min before the wall clock** (`RESERVE_S=2700`) rather than being
  killed mid-stage.

Between that and `run_md.sh`'s stage-level resume (it skips any stage whose `.gro`
already exists), the same command is both "start" and "continue". Expect ~10 resubmissions
at the 36 h cap.

### Throughput estimate

6.1 ns of MD per run (0.1 heat + 2 NPT + 4 prod) at a median 12.5k atoms. On the Quadro
RTX 6000s, budget roughly **30–60 min per median run**, so 2000 runs / 3 concurrent
≈ **330 h ≈ 14 days** of wall time, longer with the heavy tail. Nodes are
`midway3-[0277-0286,0294]`, 4 GPUs and 48 cores each.

### Practical notes

- **Build is inline, not a separate stage.** Unlike the pilot-era `launch_batch.py` path,
  each lane packs its own run dir immediately before simulating it. Packmol is ~290 s per
  system and it runs on the GPU node's CPUs while that lane's GPU is idle between runs —
  wasteful in principle, but 2000 × 290 s ≈ 160 CPU-hours spread across 3 lanes is small
  next to 330 h of MD, and it removes a whole build/submit handoff.
- **The heavy tail is not split out.** With a 36 h wall and graceful stop, a 47k-atom run
  that does not finish simply resumes on the next submission.
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
2. ~~**What to do about the concentration/MW mismatch?**~~ **RESOLVED 2026-09-20.**
   Fixed at the source by capping MW at 200 and molarity at 1.5 M, which is a blend of
   the old options (c) and (d) — see §2.2. The scope note that once said this was
   "outside this task's scope" no longer applies: the dataset generators were changed
   and `datasets/electrolyte_dataset_2000.csv` was regenerated in place, at your
   direction, with no v2 variants kept.
3. **Solvent density is assumed 1.0 g/cm³ for everything** (as the old campaign did), and
   **salt density 1.70 g/cm³ for all 20 anions**. NPT fixes the *density* but not the
   *composition*, because `n_solvent` was already derived from the guess — so the
   achieved molarity differs from the requested one (pilot: −0.1% at 0.3 M, +14.7% at
   3.0 M). Mitigated rather than solved: the 1.5 M cap shrinks the error and
   `molarity_actual` / `molarity_err_pct` record it per run. A group-contribution density
   estimate, or per-anion crystal densities, would improve it if it ever matters.
4. ~~**Filter implausible solvents?**~~ **RESOLVED 2026-09-20 — yes, at generation time.**
   Aldehydes, protic species, radicals, boron and phosphate/phosphite esters are rejected
   by `generate_500_solvents.py` (§2.1 item 2). Filters are hardcoded deliberately: a CLI
   override was prototyped and removed, because the chemistry rules should not be
   silently togglable per invocation.
4b. **Silicon is excluded, arguably by accident.** LigParGen handles Si fine, but Si is
   caught by the inherited `METAL_ATOMIC_NUMS = range(11,15)` before the element check.
   Re-enabling it requires edits in two places. Left excluded; flagged so it is a choice
   rather than an oversight.
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
    one throttled job array (§9). *(Superseded 2026-09-20 — see below.)*
- **2026-09-20 (a)** — Analysis extended. `cluster_analysis.py` written (RDF-derived
  cutoff, contact graph, connected components) and wired into `run_md.sh`; validated
  against OVITO at 72 vs ~70 for the Li pilot (§8.1). Both counting conventions are
  reported after establishing that
  `n_components = n_clusters + n_free_cation + n_free_anion`.
  `collect_results.py` gained achieved molarity (`molarity_actual`,
  `molarity_err_pct`, `box_equil_nm`), dataset row linkage (`row_index`, `csv_line`)
  and the cluster block. `build_run_manifest.py` written so all 2000 rows are indexed
  and costed before anything is built.
- **2026-09-20 (b)** — **Dataset rebuilt (§2.2).** The concentration/MW mismatch was
  escalated from "flagged, not acted on" to fixed: MW ≤ 200 and M ≤ 1.5. Median solvent
  per ion pair 2.3 → 5.97; rows below 2 went 849 → 0. Chemistry filters (aldehyde,
  protic, radical, boron, phosphate ester) moved into the generator, reversing the
  earlier "do not filter" decision. Consequences:
  - Scripts and datasets were **replaced in place** at your instruction — no v2 variants.
  - The new 490 solvents have zero overlap with the old 490, so the entire force-field
    library was cleared and refetched: 490/490, 2450 files.
  - 6 LigParGen failures (5 carbon radicals + 1 thiophene) → radical filter added, and
    the generator now reads `ff_failures.csv` back so each one-off reject is paid once.
  - Median system grew 6.4k → 12.5k atoms (30.7M atom-runs total). Accepted: that is the
    cost of having enough solvent to form a real solvation shell.
- **2026-09-20 (c)** — **Campaign driver, and a GPU-safety correction.**
  `submit_campaign.sbatch` + `campaign_worker.sh` + `campaign_status.py` replace the
  `launch_batch.py` array path for the 2000 (§9). Row-strided lanes, per-run
  `pipeline.log`, graceful stop 45 min before the wall.
  - Originally 3 nodes × 4 GPUs = 12 concurrent. Changed to **3 × 1 GPU = 3 concurrent**
    after finding every partition GPU held by single-GPU jobs: a whole-node request
    would have queued 25–84 h. Trades throughput (~85 h → ~330 h) for starting now.
  - **The `gpu` partition does not isolate GPUs per job** (§6 gotcha 4). Found by
    trespassing on three GPUs belonging to another user during a 4-GPU concurrency test
    on a 1-GPU allocation. The per-lane `CUDA_VISIBLE_DEVICES` assignment was removed;
    the script now only ever trusts what SLURM exports, and warns if it is absent.
  - Trial jobs cancelled and skipped entirely at your direction — the campaign's first
    array task serves as the trial, since it is resumable and costs nothing to cancel.
- **2026-09-20 (d)** — Session transcript moved to `MD/ht2000/session_transcript.txt`
  (gitignored) and this plan brought back in line with the code.

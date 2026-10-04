# Figure 5 — host adiposity, stromal composition, and the compound screen

Code for **panel a** and for the pipeline that nominated the compounds screened in **panel b** of Figure 5 of
*A neurotrophic secretome from adipose-progenitor-like stromal cells promotes prostate cancer lineage plasticity*.

Everything here is the version used in the manuscript. Two inputs were replaced during revision, and the reasons
are stated in the scripts themselves and summarised under **What changed** below.

---

## Layout

```
R/common.R                        paths (environment variables), helpers, palette, reference gene sets
R/fig5a_01_bmi_composition.R      panel a, analysis: stromal composition against BMI in the Beltran cohort
R/fig5a_02_panel.R                panel a, the print-ready panel
R/fig5b_01_signature_reversal.R   panel b, branch B: LINCS L1000 reversal in the Beltran NEPC-like biopsies
R/fig5b_02_compound_pipeline.R    panel b, branch A + convergence: the compound list that was screened
R/fig5b_03_schematic.R            the full pipeline schematic (supplementary)
R/fig5b_04_panel_outline.R        the condensed pipeline strip used as a main-figure panel
data/                             derived tables and gene signatures (see below)
```

Run in order:

```bash
Rscript R/fig5a_01_bmi_composition.R     # writes the per-sample scores and the statistics
Rscript R/fig5a_02_panel.R               # panel a
Rscript R/fig5b_01_signature_reversal.R  # branch B (queries the live L1000FWD API)
Rscript R/fig5b_02_compound_pipeline.R   # branch A, drug retrieval, cross-evidence ranking
Rscript R/fig5b_03_schematic.R
Rscript R/fig5b_04_panel_outline.R
```

## Paths

No path is hard-coded. `R/common.R` reads these environment variables and falls back to the defaults used on the
machine the figure was made on:

| Variable | What it points at |
|---|---|
| `F5_REPO` | this repository (default: the parent of `R/`) |
| `F5_DATA` | working directory for intermediate results |
| `F5_RR1384` | Beltran RR1384 bulk RNA-seq: counts matrix and sample key |
| `F5_CLINICAL` | the WCM clinical table (BMI, therapy, histology), `.xlsx` |
| `F5_OUT` | where finished panels are written (default `figures/`) |

## Panel a — what it tests and why it is a contrast

Every stromal signature scores high wherever there is more stroma, so the PI16+ and POSTN+ scores are collinear
in bulk and each tracks stromal content. A raw "the PI16+ score rises with BMI" result is therefore
indistinguishable from "stromal content rises with BMI". The primary outcome is a **contrast**,

```
delta = z(PI16+ specific) − z(POSTN+ specific)
```

which is invariant to a shared multiplicative content term: a biopsy with twice the stroma raises both scores and
leaves `delta` unchanged. The question becomes whether adiposity changes the **composition** of the stroma rather
than its **amount**. The total stromal score is reported alongside and entered as a covariate. The three gene sets
are made mutually exclusive before scoring and the disjointness is asserted at run time.

The analysis is **exploratory**: with roughly 25 samples per category a null result is uninformative, and the
script says so. Repeat biopsies per patient are reported through a mixed model with `(1|record_id)`.

## Panel b — how the screened compounds were nominated

Two branches converge on one approved-drug list.

**Branch A — PI16+ / POSTN+ stromal-state transcription factors.** Within the fibroblast lineage of the 15
single-cell datasets used in Figures 1 and 2 (11 human, 4 mouse; 10,322 PI16+, 6,356 POSTN+, 25,394 rest-of-lineage
cells), a transcription factor enters the list when it is differentially expressed in one of the two states, in the
same direction, in at least three datasets of a species (**tier 1**), or when one of the supporting analyses run on
the same cells returns it (**tier 2**). That gives **138 state transcription factors** — 50 PI16+, 89 POSTN+,
37 tier 1 — and **100 upstream kinases and receptors** for 74 of them from OmniPath.

**Branch B — ligand co-expression and signature reversal.** The top 250 positively and negatively correlated genes
with *MDK* and with *TNFSF12* in the Beltran NEPC-like biopsies are submitted to LINCS L1000FWD; the compounds whose
signature opposes that programme are collected. The AR-altered biopsies of the same cohort are run as a contrast,
so a compound that reverses the programme in neuroendocrine disease specifically can be told apart from one that
reverses it in any metastatic prostate tissue.

**Convergence.** Approved drugs against either target set come from Open Targets, kept where maximum clinical stage
is approval. Compounds are ranked on three independent axes — acts on a state TF, acts on an upstream kinase or
receptor, reverses the ligand programme — and a compound supported on more than one axis ranks above one supported
on a single axis.

**Result: 343 compounds.** 318 FDA-approved from drug-target mapping (163 acting on a state TF, 155 on an upstream
kinase or receptor); 26 further compounds from signature reversal, 19 of which do not reverse the programme in the
AR-altered contrast. Perhexiline is the only compound supported on two axes.

Of the ten compounds screened in panel b, **five are recovered by this pipeline**, all through branch A:
fulvestrant and tamoxifen (ESR1), betamethasone (NR3C1), sorafenib (KIT, PDGFRB, RAF1, RET) and trametinib
(MAP2K2). Atorvastatin, mirabegron, paroxetine, salbutamol and sildenafil are not recovered. Four of those five
had been nominated by earlier runs on cohorts that are not neuroendocrine disease, each on a single L1000 signature
in a single unrelated cell line; mirabegron was never nominated by any version. `data/ten_tested_compounds.csv`
records this per compound.

## What changed, and why

**The GSE235036 route was removed.** The original branch A step 3 intersected a promoter-derived transcription-factor
list with `cluster.markers.sig.csv %>% filter(cluster == "c5")`. That file is the **whole-tissue** clustering of the
single mouse dataset GSE235036, in whose numbering **c5 is the entire stromal compartment** — not a state. Its
markers separate stroma from epithelium, in one mouse, and cannot distinguish PI16+ from POSTN+ because the two were
never resolved. It is also unrelated to *atlas* c5, the POSTN-high mesenchymal state, with which it is easily
confused. For the record: 1,800 c5 markers, 96 of them transcription factors, only 15 of which survive the state
analysis.

**The MDK-promoter motif step was removed.** A HOCOMOCO scan of the *MDK* promoter with an active-chromatin filter
discarded transcription factors on motif presence at one of the two ligand genes — *TNFSF12* was never scanned the
same way — and does not discriminate regardless: 1,116 transcription factors have a hit in that window.

**The bulk cohort was swapped.** Branch B previously ran on TCGA PRAD, SU2C and an ARCHS4 adipose compendium. TCGA
PRAD is primary, hormone-naive, overwhelmingly adenocarcinoma prostate tissue — the wrong cohort for a question
about neuroendocrine disease. It now runs on the Beltran metastatic cohort only.

**Declared substitutions.** Upstream kinase–substrate edges come from OmniPath because the full PhosphoSitePlus
dataset is licence-restricted and cannot be redistributed; approved drugs come from Open Targets because Selleck
Chemicals publishes no programmatic catalogue. Every compound in the list can be looked up in the Selleck catalogue
by name.

## Data in this repository

| File | What it is |
|---|---|
| `data/state_TFs_fig1_2.csv` | the 138 PI16+/POSTN+ state transcription factors, with tier and supporting evidence |
| `data/target_genes.csv` | the 234 genes queried against Open Targets (state TFs plus upstream regulators) |
| `data/compound_list.csv` | the 343 nominated compounds with their cross-evidence scores |
| `data/ten_tested_compounds.csv` | the account of the ten compounds screened in panel b |
| `data/l1000_opposite_beltran.csv` | the L1000 signatures opposing the ligand programme, by disease group |
| `data/signatures/PI16_PROGENITOR_200.txt` | the 200-gene progenitor signature used to call PI16+ cells |
| `data/signatures/FIBROBLAST_STROMA.txt` | the 200-gene pan-fibroblast signature |

### Inputs that are NOT in this repository

- **`data/signatures/POSTN_HIGH_100.txt`** — the 100-gene POSTN-high signature. Derived from the mouse mesenchymal
  atlas; it lives on the analysis volume and must be supplied before `R/fig5a_01_bmi_composition.R` will run.
- **Patient-level data.** The RR1384 counts matrix, the sample key and the clinical table (BMI, therapy, histology)
  are controlled-access and are deliberately not published here. Point `F5_RR1384` and `F5_CLINICAL` at your own
  copies.
- **Upstream dependencies of branch A.** `R/fig5b_02_compound_pipeline.R` reads the per-state differential-expression
  tables produced by the Figure 1 and 2 analysis under `$F5_DATA`: `state_tf_network/`, `nepc_tf_replication/`,
  `collectri_per_dataset/`, `nepc_tf_consistency/`, `human_tf_strategy/` and `nepc_tf_expression/`. The finished
  product of that chain is shipped as `data/state_TFs_fig1_2.csv`, so the rest of the pipeline can be reproduced
  from this repository alone.
- **OmniPath** is downloaded on first run into `$F5_DATA/omnipath/`.

## Requirements

R ≥ 4.4, with `dplyr`, `tidyr`, `tibble`, `ggplot2`, `readxl`, `writexl`, `data.table`, `httr`, `jsonlite`,
`org.Hs.eg.db`, `AnnotationDbi`; optionally `lmerTest` for the mixed model.

Branch B queries the live LINCS L1000FWD API and the drug step queries the live Open Targets GraphQL API; both are
cached locally after the first run.

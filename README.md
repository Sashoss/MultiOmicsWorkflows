# Multi‑omics Regulatory Landscape Workflows

> **Integrated RNA‑Seq, ATAC‑Seq & ChIP‑Seq pipelines for transcriptional network discovery**

<hr>
<br>

## Overview

This repository hosts **validated, end‑to‑end Jupyter‑based workflows** that jointly analyze transcriptomic (RNA‑Seq), chromatin accessibility (ATAC‑Seq) and transcription‑factor binding (ChIP‑Seq) data.  Together, they enable researchers to:

1. **Map transcription‑factor (TF) binding sites** via ChIP‑Seq peak calling and motif analysis.
2. **Interrogate chromatin openness** at those sites using ATAC‑Seq to retain only peaks residing in accessible regulatory regions.
3. **Quantify differential gene expression** with RNA‑Seq to confirm functional consequences of TF binding.
4. **Integrate the three layers** to build high‑confidence TF➜target gene regulatory networks, prioritize candidate enhancers, and generate mechanistic hypotheses for further validation.

Each workflow has been benchmarked on human hg38, mouse mm10 and custom genomes, and is readily adaptable to additional organisms.

<hr>
<br>

## Repository Layout

```text
repo_root/
├── Notebook/                # Interactive, modular notebooks
│   ├── RNASeq/              # RNA‑Seq preprocessing & DGE
│   ├── ATACSeq/             # ATAC‑Seq QC, peak calling, footprinting
│   └── ChIPSeq/             # ChIP‑Seq QC, peak calling, motif scan
│
├── envs/                    # Conda YAMLs for each assay
└── README.md                # <‑‑ you are here
```

<hr>
<br>


## Scientific Rationale

| Layer        | Key Output                                   | Purpose in Integration                                                                       |
| ------------ | -------------------------------------------- | -------------------------------------------------------------------------------------------- |
| **ChIP‑Seq** | High‑confidence TF peaks & motifs            | Identify candidate regulatory regions bound by a TF of interest.                             |
| **ATAC‑Seq** | Accessible chromatin peaks; footprint scores | Filter ChIP peaks to those occurring in open chromatin and predict co‑bound factors.         |
| **RNA‑Seq**  | Differentially expressed (DE) genes          | Validate regulatory impact by linking accessible ChIP peaks to genes showing coordinated DE. |

Combining these layers sharpens specificity: **a gene whose promoter harbors an accessible ChIP peak for the TF *and* is significantly up‑ or down‑regulated constitutes a high‑confidence direct target.**

<hr>
<br>

## Applications

* Discovery of TF regulatory circuits in development & disease
* Prioritisation of enhancer–gene pairs for CRISPR perturbation screens
* Cross‑species comparison of conserved regulatory logic
* Drug‑target identification by coupling TF binding with expression changes

---

<hr>
<br>
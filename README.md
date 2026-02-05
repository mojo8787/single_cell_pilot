# Single-Cell Biofilm Heterogeneity Pilot

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18499914.svg)](https://doi.org/10.5281/zenodo.18499914)

## Overview

Bacterial populations exhibit substantial transcriptional heterogeneity at the single-cell level, with individual cells adopting distinct phenotypes—including bet-hedging strategies, prophage induction, and phage-infected states—that are invisible to bulk sequencing [1]. M3-Seq (Massively-parallel Multiplexed Microbial sequencing) addresses this by combining combinatorial cell indexing with post hoc rRNA depletion, enabling profiling of hundreds of thousands of bacterial cells and revealing rare subpopulations in *E. coli* and *B. subtilis* [1,2]. This pilot reanalyses publicly available M3-Seq data from *E. coli* MG1655 (GSE231935) [3] to demonstrate population heterogeneity and cluster-specific gene expression patterns relevant to biofilm regulation and mobile genetic element (MGE)–host interactions.

The analysis pipeline uses [Scanpy](https://scanpy.readthedocs.io/) [4] for quality control, normalization, dimensionality reduction (PCA/UMAP), clustering, and differential expression. Findings support the concept that *E. coli* populations harbour distinct transcriptional states—consistent with curli ON/OFF bistability, stress responses, and MGE-mediated rewiring of regulatory networks [1,5]—and validate the approach for future work on biofilm evolution and rare subpopulation characterization.

> **Cite this work:** See [CITATION.cff](CITATION.cff) or the *How to cite* section below.

## Dataset

- **Source:** NCBI GEO [GSE231935](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231935) [3]
- **Original study:** Wang et al., *Nature Microbiology* 2023 [1]
- **Content:** M3-Seq — *E. coli* MG1655, 7,681 cells × 4,396 genes (post-QC: 3,525 cells × 2,619 genes)
- **Relevance:** Rare subpopulations (bet-hedging, prophage induction, phage-infected cells) [1,2]

## Setup

```bash
cd single_cell_pilot
python3 -m venv .venv
source .venv/bin/activate   # On Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Optional: `pip install leidenalg` for Leiden clustering (otherwise KMeans on PCA is used).

## Usage

### 1. Download data

```bash
python scripts/01_download_data.py
```

Fetches GSE231935_RAW.tar (~42 Mb) from GEO and extracts CSV files to `data/`.

### 2. Run analysis

```bash
python scripts/02_analyze.py
```

Pipeline: load → QC → normalize → PCA/UMAP → Leiden clustering → differential expression → figures.

## Outputs

| File | Description |
|------|-------------|
| `output/umap_clusters.png` | UMAP colored by cluster — transcriptional heterogeneity |
| `output/violin_total_counts.png` | Total counts per cluster |
| `output/heatmap_top5.png` | Top 5 DE genes per cluster |
| `output/dotplot__top3.png` | Dotplot of top 3 DE genes per cluster |
| `output/summary_stats.csv` | Cell counts, gene counts, cluster sizes |

## Results

### UMAP clusters (8 transcriptional states)

![UMAP clusters](https://raw.githubusercontent.com/mojo8787/single_cell_pilot/main/output/umap_clusters.png)

### Total counts per cluster

![Violin plot](https://raw.githubusercontent.com/mojo8787/single_cell_pilot/main/output/violin_total_counts.png)

### Top 5 differentially expressed genes per cluster

![Heatmap](https://raw.githubusercontent.com/mojo8787/single_cell_pilot/main/output/heatmap_top5.png)

### Top 3 DE genes per cluster (dotplot)

![Dotplot](https://raw.githubusercontent.com/mojo8787/single_cell_pilot/main/output/dotplot__top3.png)

## Findings & Future Work

- **UMAP clusters:** *E. coli* populations exhibit distinct transcriptional states at single-cell resolution, supporting population heterogeneity (e.g., curli ON/OFF bistability, MGE-induced states) [1,5].
- **Differential expression:** Cluster-specific genes include stress response, prophage-related, and biofilm-associated loci—consistent with MGE-mediated rewiring of regulatory networks [1,5].
- **Future directions:** Findings validate the approach for follow-up studies on biofilm regulation, MGE-host interactions, and rare subpopulation characterization. The pipeline is reusable for similar bacterial single-cell datasets.

## Zenodo

1. Push this repo to GitHub (e.g. `mojo8787/single_cell_pilot`).
2. Go to [Zenodo](https://zenodo.org) → Log in → GitHub → Enable the repo.
3. Create a release (e.g. `v1.0.0`) — Zenodo will mint a DOI.

Metadata is in [CITATION.cff](CITATION.cff) and [.zenodo.json](.zenodo.json).

## Project structure

```
single_cell_pilot/
├── data/                    # Downloaded GSE231935 data
├── scripts/
│   ├── 01_download_data.py  # Fetch from GEO
│   └── 02_analyze.py        # Scanpy pipeline
├── output/                  # Figures and stats
├── .zenodo.json             # Zenodo metadata
├── CITATION.cff             # Citation metadata (Zenodo, GitHub)
├── LICENSE                  # MIT
├── requirements.txt
└── README.md
```

## Author

**Almotasem Bellah Younis, PhD**  
Division of Microbial Ecology (DOME), Centre for Microbiology and Environmental Systems Science (CeMESS), University of Vienna  
[motasem.youniss@gmail.com](mailto:motasem.youniss@gmail.com) · [motasemyounis.com](https://motasemyounis.com) · [ORCID 0000-0003-2070-2811](https://orcid.org/0000-0003-2070-2811)

## How to cite

```bibtex
@software{younis_single_cell_pilot_2026,
  author = {Younis, Almotasem Bellah},
  title = {Single-Cell Biofilm Heterogeneity Pilot},
  year = {2026},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.18499914},
  url = {https://github.com/mojo8787/single_cell_pilot},
  version = {1.0.1}
}
```

Or use the [CITATION.cff](CITATION.cff) file. Zenodo DOI: [10.5281/zenodo.18499914](https://doi.org/10.5281/zenodo.18499914)

## References

1. Wang B, Lin AE, Yuan J, Novak KE, Koch MD, Wingreen NS, Adamson B, Gitai Z. Single-cell massively-parallel multiplexed microbial sequencing (M3-seq) identifies rare bacterial populations and profiles phage infection. *Nature Microbiology* 2023;8(10):1846–1862. [DOI:10.1038/s41564-023-01462-3](https://doi.org/10.1038/s41564-023-01462-3) · [PMID:37653008](https://pubmed.ncbi.nlm.nih.gov/37653008/)

2. Wang B, Lin AE, Yuan J, Koch MD, Adamson BS, Wingreen N, Gitai Z. Massively-parallel Microbial mRNA Sequencing (M3-Seq) reveals heterogeneous behaviors in bacteria at single-cell resolution. *bioRxiv* 2022. [DOI:10.1101/2022.09.21.508688](https://doi.org/10.1101/2022.09.21.508688)

3. Wang B, Lin AE, Yuan J, Koch MD, Adamson BS, Wingreen N, Gitai Z. Massively-parallel Microbial mRNA Sequencing (M3-Seq) reveals heterogenous behaviors in bacteria at single-cell resolution. *Gene Expression Omnibus* 2023;GSE231935. [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231935](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231935)

4. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* 2018;19:15. [DOI:10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)

5. Veening J-W, Smits WK, Kuipers OP. Bistability, epigenetics, and bet-hedging in bacteria. *Annual Review of Microbiology* 2008;62:193–210. [DOI:10.1146/annurev.micro.62.081307.163002](https://doi.org/10.1146/annurev.micro.62.081307.163002)

## License

Code: MIT. Data: see GEO/NCBI terms. See [LICENSE](LICENSE).

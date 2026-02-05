# Single-Cell Biofilm Heterogeneity Pilot

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Bioinformatics reanalysis of bacterial single-cell RNA-seq data (GSE231935, M3-Seq) demonstrating transcriptional heterogeneity in *E. coli* populations. Findings support future work on biofilm/MGE-related projects and population heterogeneity.

> **Cite this work:** See [CITATION.cff](CITATION.cff) or the *How to cite* section below.

## Dataset

- **Source:** NCBI GEO [GSE231935](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231935)
- **Citation:** Wang B et al., Nature Microbiology 2023 (PMID: 37653008)
- **Content:** M3-Seq (Massively-parallel Microbial mRNA Sequencing) — *E. coli* MG1655, 7,681 cells × 4,396 genes
- **Relevance:** Identifies rare subpopulations (bet-hedging, prophage induction, phage-infected cells), aligning with MGE/biofilm heterogeneity

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

## Findings & Future Work

- **UMAP clusters:** *E. coli* populations exhibit distinct transcriptional states at single-cell resolution, supporting population heterogeneity (e.g., curli ON/OFF bistability, MGE-induced states).
- **Differential expression:** Cluster-specific genes include stress response, prophage-related, and biofilm-associated loci — consistent with MGE-mediated rewiring of regulatory networks.
- **Future directions:** Findings validate the approach for follow-up studies on biofilm regulation, MGE-host interactions, and rare subpopulation characterization. The pipeline is reusable for similar bacterial single-cell datasets.

## Zenodo

1. Push this repo to GitHub (e.g. `mojo8787/single_cell_pilot`).
2. Go to [Zenodo](https://zenodo.org) → Log in → GitHub → Enable the repo.
3. Create a release (e.g. `v1.0.0`) — Zenodo will mint a DOI.
4. Replace `XXXXXXX` in the DOI badge above with your Zenodo record ID.

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
  url = {https://github.com/mojo8787/single_cell_pilot},
  version = {1.0.0}
}
```

Or use the [CITATION.cff](CITATION.cff) file. After depositing on Zenodo, replace the DOI badge above and add the Zenodo DOI to your citations.

## License

Code: MIT. Data: see GEO/NCBI terms. See [LICENSE](LICENSE).

# Cardiac Hypertension Transcriptomics

Chamber-specific RNA-seq differential expression analysis of the hypertensive female human heart. Re-analysis of GSE282618 (Milstone et al., 2025).

*Co-authored with Micaiah McDonald, M.S. — Meharry Medical College, MSBD501*

<img width="600" height="450" alt="pca_pooled_degs" src="https://github.com/user-attachments/assets/d15b8ff7-5df4-44ee-ae46-e50dfc9dab43" />

<img width="600" height="450" alt="volcano_pooled" src="https://github.com/user-attachments/assets/4128d3f9-dc22-454c-a3b0-0860e4f56223" />

## Overview
Hypertension affects roughly 1.4 billion adults worldwide and is a leading driver of heart failure — but most transcriptomic studies pool left and right ventricle (LV/RV) samples into a single model, which can obscure chamber-specific biology. This project re-analyzes public postmortem cardiac RNA-seq data with three separate DESeq2 models (pooled, LV-only, RV-only) to test whether hypertensive remodeling is shared across ventricles or chamber-specific.
 
**Research questions:**
1. How does the hypertensive transcriptome differ from control in female ventricular tissue?
2. Are hypertension-associated transcriptional changes shared between LV and RV, or chamber-specific?
3. Does modeling choice (pooled vs. chamber-stratified) meaningfully change which genes are detected?
## Key Findings
- **213 DEGs** identified in the pooled model (padj < 0.05, |LFC| > 1), with downregulation predominating (114 down vs. 99 up)
- Chamber-stratified analysis found **58 LV-specific and 71 RV-specific DEGs** — but only **11 were shared between chambers (~14–19% overlap)**, indicating hypertensive remodeling is largely chamber-specific, not global
- Top LV hit: **YTHDC1** (padj ≈ 10⁻¹¹), an m6A RNA methylation reader previously linked to dilated cardiomyopathy
- Top RV hit: **TNC** (padj ≈ 10⁻³⁴), an extracellular matrix glycoprotein and established marker of cardiac fibrosis — the single most significant gene across the entire study
- **10 genes were significant across all three models** (pooled, LV, RV), representing the highest-confidence, chamber-independent hypertension targets: *HLX-AS1, HECW1, LSM14A, CD5L, SERPINA5, CATSPERT, DDX52, DHFR*, and two unannotated Ensembl IDs
- A PCA restricted to the 213 significant DEGs cleanly separated HTN from Control samples on PC1 (70% variance) — even though global PCA on all genes showed no condition separation, since donor-level variability dominated
**Takeaway:** pooled models can detect a hypertension signal, but chamber-stratified analysis reveals substantially more — and more specific — biology than pooling alone. Given the small cohort (n=5 donors), findings are hypothesis-generating rather than confirmatory.
 
## Dataset
- **Source:** [GEO GSE282618](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282618) (Milstone et al., 2025, *Physiological Reports*)
- **Samples:** 5 postmortem female donors (2 control, 3 hypertensive), each contributing one LV and one RV sample → 10 total samples
- **Genes tested:** 34,737, after filtering genes with <10 counts in at least 2 samples
- **Note:** hypertensive donors were on average 29 years older than controls, and comorbidities (T2DM, pulmonary hypertension) were present only in the HTN group — a real confound the report discusses directly rather than glossing over
## Pipeline
 
**1. Preprocessing** (`rnaseq_pipeline_refined.sh`) — raw FASTQ → aligned, quantified counts
```
SRA download → FastQC/MultiQC (raw) → Trimmomatic (paired-end trim)
→ FastQC/MultiQC (post-trim) → HISAT2 alignment (hg38) → samtools sort/index
→ featureCounts quantification
```
Run a single-sample test with `bash scripts/preprocessing/rnaseq_pipeline_refined.sh --test`, or the full 10-sample run with `bash scripts/preprocessing/rnaseq_pipeline_refined.sh`.
 
**2. Differential expression** (`DESeq2_pipeline_refined.Rmd`) — three DESeq2 models from the same count matrix:
- **Pooled:** `~ ventricle + condition`, all 10 samples, ventricle as covariate
- **LV-only:** `~ condition`, 5 samples
- **RV-only:** `~ condition`, 5 samples
All three apply apeglm shrinkage to log2 fold changes and classify DEGs at padj < 0.05, |LFC| > 1. See the rendered [`DESeq2_pipeline_refined.html`](notebooks/deg/DESeq2_pipeline_refined.html) for full output with all figures inline.
 
## Quality Control & Preprocessing Results
- **Pre-trim QC:** no adapter contamination detected; per-base quality, per-sequence quality, and N-content all passed across all 10 samples. Some overrepresented sequences and elevated duplication levels were flagged (expected for RNA-seq, where highly-expressed transcripts naturally duplicate).
- **Trimming:** Trimmomatic (paired-end; LEADING:3, TRAILING:3, MINLEN:36) retained **93.7–99.9% of read pairs** across samples — see [pre-trim](quality_check/multiqc/pre_trim_multiqc_report.html) and [post-trim](quality_check/multiqc/post_trim_multiqc_report.html) MultiQC reports.
- **Post-trim QC:** consistent improvement in base quality and a substantial reduction in overrepresented sequences; residual duplication warnings are expected and not a data quality concern.
- **Raw count matrix:** 78,899 genes quantified across 10 samples via featureCounts — this is the unfiltered matrix; the DESeq2 analysis stage (above) filters this down to 34,737 genes retained for testing (≥10 counts in ≥2 samples).
## Sample Metadata
| Donor | Chamber | Condition | SRA ID | GEO ID |
|---|---|---|---|---|
| 1 | RV | HTN | SRR31443256 | GSM8647323 |
| 1 | LV | HTN | SRR31443257 | GSM8647322 |
| 3 | RV | HTN | SRR31443254 | GSM8647325 |
| 3 | LV | HTN | SRR31443255 | GSM8647324 |
| 4 | RV | Control | SRR31443252 | GSM8647327 |
| 4 | LV | Control | SRR31443253 | GSM8647326 |
| 5 | RV | HTN | SRR31443250 | GSM8647329 |
| 5 | LV | HTN | SRR31443251 | GSM8647328 |
| 7 | RV | Control | SRR31443248 | GSM8647331 |
| 7 | LV | Control | SRR31443249 | GSM8647330 |
 
Full table: [`Sample_Table.xlsx`](data/Sample_Table.xlsx)
 
 
## Repository Structure
```
├── scripts/
│   └── preprocessing/
│       └── rnaseq_pipeline_refined.sh      # Bash: SRA → aligned counts
├── notebooks/
│   └── deg/
│       ├── DESeq2_pipeline_refined.Rmd     # R: 3 DESeq2 models + all figures
│       └── DESeq2_pipeline_refined.html    # Rendered analysis notebook
├── quality_check/
│   └── multiqc/
│       ├── pre_trim_multiqc_report.html
│       └── post_trim_multiqc_report.html
├── data/
│   └── Sample_Table.xlsx                   # Donor/chamber/condition/GEO map
├── report/
│   └── Final_Team_Report.docx
├── presentation/
│   └── Final_Team_Presentation.pptx
└── plots/                                  # Exported figures (PCA, volcano, MA, heatmaps)
```
 
## Tech Stack
**Preprocessing:** Bash · SRA Toolkit · FastQC · MultiQC · Trimmomatic · HISAT2 · SAMtools · featureCounts (Subread)
**Analysis:** R · DESeq2 · apeglm · AnnotationDbi / org.Hs.eg.db · limma · ggplot2 · pheatmap · ggrepel · dplyr
 
## Limitations
- Small cohort (5 donors, 10 samples total) — chamber-stratified models are underpowered; results are hypothesis-generating, not confirmatory
- Age is confounded with condition (HTN donors ~29 years older on average)
- Comorbidities (T2DM, pulmonary hypertension) present only in the HTN group add unmodeled transcriptional noise
- Individual donor-level covariates (cause of death, race/ethnicity) were unavailable and could not be statistically controlled for
## Authors
**Destiny Pounds** — [portfolio](https://dpounds24.github.io/portfolio/) · [LinkedIn](https://linkedin.com/in/destiny-pounds)
**Micaiah McDonald**, M.S. — co-author
 
Preprocessing pipeline refinement, QC, DESeq2 modeling, and report/presentation writing were shared between both authors; see the presentation's individual-contributions slide for the full breakdown.
 
## Acknowledgments
- Original dataset: Milstone, Z. J., Moreira-Bouchard, J. D., et al. (2025). Multiomics investigation of the female hypertensive human heart. *Physiological Reports*, 13, e70586.
- AI tools (Copilot/ChatGPT) were used to assist with code commenting and error handling in the preprocessing script; all logic was reviewed by the authors.
 

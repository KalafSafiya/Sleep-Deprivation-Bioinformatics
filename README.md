# Analysis of Gene Expression Changes During Sleep Deprivation

## 1. About the Project

**What happens inside the human body when we stay awake for too long?**

Sleep is important for maintaining normal body functions. When a person stays awake for an extended period, it may affect different biological processes at the molecular level.

This project investigates how **prolonged sleep deprivation affects gene expression in human blood**.

Gene expression is the process through which information stored in genes is used by cells to produce functional products. By studying changes in gene expression, we can explore how the body responds to sleep deprivation and what changes occur during recovery.

The project uses publicly available human gene-expression data and bioinformatics techniques to identify changes in individual genes and groups of biologically related genes.

## 2. Research Question

**How does prolonged sleep deprivation affect gene expression in human blood, and how do these changes differ between baseline, sleep deprivation, and recovery periods?**

## 3. Dataset Used

The project uses the publicly available dataset **GSE98564**, obtained from the NCBI Gene Expression Omnibus (GEO).

| Dataset Information | Details |
|---|---|
| GEO Accession | GSE98564 |
| Organism | Human |
| Biological Sample | Blood |
| Technology | Affymetrix Human Gene 1.0 ST Microarray |
| Participants | 17 |
| Control Participants | 6 |
| Sleep-Deprived Participants | 11 |
| Total Samples | 199 |
| Transcript Clusters | 8,454 |

The dataset contains samples collected during three study phases:

- **Baseline:** Measurements before the experimental sleep-deprivation period.
- **Experimental:** Measurements during the experimental period.
- **Recovery:** Measurements during the recovery period.

The dataset includes repeated measurements from participants across the study phases.

The expression data used in the analysis were already processed using the RMA (Robust Multi-array Average) normalization method. Therefore, the analysis used the processed GEO expression data rather than applying normalization again.

**Dataset:** [GSE98564 – NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98564)

## 4. Tools and Technologies

The following tools and packages were used:

| Tool / Package | Purpose |
|---|---|
| R | Main programming language for the analysis |
| RStudio | Environment used to write and execute R scripts |
| GEOquery | Retrieve the gene-expression dataset from GEO |
| limma | Perform differential gene-expression analysis |
| AnnotationDbi | Map transcript identifiers to gene information |
| hugene10sttranscriptcluster.db | Annotate transcript clusters |
| org.Hs.eg.db | Obtain human gene annotations |
| clusterProfiler | Perform Gene Set Enrichment Analysis |
| enrichplot | Visualize pathway enrichment results |
| ggplot2 | Create data visualizations |
| pheatmap | Generate the gene-expression heatmap |
| Git and GitHub | Manage and store project files and analysis code |

## 5. What Was Done?

The project was completed through several stages.

### 5.1 Data Collection and Preparation

- Retrieved the GSE98564 dataset from NCBI GEO.
- Imported the processed gene-expression matrix and sample information into R.
- Prepared sample metadata, including participant IDs, group information, and study phases.
- Checked the consistency between expression data and sample metadata.
- Examined missing values and expression-value distributions.

### 5.2 Exploratory Data Analysis (EDA)

Exploratory Data Analysis was performed to understand the structure and characteristics of the dataset.

This included:

- Checking the number of samples and transcript clusters.
- Examining the distribution of gene-expression values.
- Checking missing values.
- Examining the number of participants in each group and study phase.

### 5.3 Differential Gene Expression Analysis

Differential expression analysis was performed using the `limma` package.

The purpose was to identify transcript clusters whose expression levels differed between study phases.

The analysis included six comparisons:

1. Sleep-deprived group: Experimental vs Baseline.
2. Sleep-deprived group: Recovery vs Experimental.
3. Sleep-deprived group: Recovery vs Baseline.
4. Control group: Experimental vs Baseline.
5. Experimental difference: Sleep-deprived group compared with the control group.
6. Recovery difference: Sleep-deprived group compared with the control group.

Repeated measurements from the same participants were accounted for using `duplicateCorrelation` and a participant-level blocking structure.

Multiple testing was addressed using the Benjamini-Hochberg (BH) adjustment method.

### 5.4 Gene Annotation and TAF4 Analysis

The significant transcript cluster identified in the sleep-deprived group's Recovery vs Baseline comparison was further investigated.

The transcript cluster was mapped to its corresponding gene using annotation databases.

The identified gene was:

**TAF4 – TATA-box binding protein associated factor 4**

Its expression levels were examined across the different study phases and groups.

### 5.5 Principal Component Analysis (PCA)

Principal Component Analysis was performed to explore patterns and similarities in gene-expression profiles.

The analysis used the 500 most variable transcript clusters.

A PCA plot was generated to visualize the samples based on their expression profiles.

### 5.6 Hierarchical Clustering and Heatmap

Hierarchical clustering was used to explore similarities between gene-expression patterns.

A heatmap was generated using the 30 most variable transcript clusters.

The heatmap helps visualize expression patterns across samples and study groups.

### 5.7 Gene Set Enrichment Analysis (GSEA)

Gene Set Enrichment Analysis was performed using the `clusterProfiler` package.

Unlike individual gene analysis, GSEA examines groups of related genes to identify biological pathways that show coordinated expression patterns.

The analysis used:

- Ranked genes based on moderated t-statistics.
- KEGG biological pathway information.
- Human organism identifier (`hsa`).
- Benjamini-Hochberg multiple-testing correction.

The analysis identified significantly enriched pathways using an adjusted p-value threshold of 0.05.

### 5.8 Investigation of Important Biological Pathways

Two enriched pathways were investigated in greater detail.

#### Ribosome Pathway

The Ribosome pathway is associated with the cellular machinery responsible for protein synthesis.

The analysis included:

- Identifying the leading-edge genes contributing to pathway enrichment.
- Mapping gene identifiers to gene symbols and names.
- Extracting gene-expression statistics.
- Examining the direction of the moderated t-statistics.
- Generating a pathway enrichment plot.

#### Oxidative Phosphorylation Pathway

Oxidative phosphorylation is a biological process involved in cellular energy production within mitochondria.

The analysis included:

- Identifying leading-edge genes.
- Annotating the genes using human gene databases.
- Extracting their expression statistics.
- Examining the direction of their moderated t-statistics.

## 6. Key Findings

The following are the main findings obtained from the analysis.

| Analysis | Main Finding |
|---|---|
| Dataset | 199 blood samples from 17 participants |
| Differential Expression | No individually significant transcript clusters in the primary comparison |
| TAF4 | One significant transcript cluster in the sleep-deprived group's Recovery vs Baseline comparison |
| KEGG GSEA | 62 significantly enriched pathways out of 310 |
| Ribosome | NES = 2.743; adjusted p-value = 4.24 × 10⁻¹² |
| Oxidative Phosphorylation | NES = 1.904; adjusted p-value = 0.000677 |

### Important Observations

**1. Individual Gene Analysis**

The primary comparison identified no transcript clusters that remained statistically significant after multiple-testing correction.

**2. TAF4**

One significant transcript cluster corresponding to TAF4 was identified in the sleep-deprived group's Recovery vs Baseline comparison.

Its estimated expression was approximately 14.3% lower during recovery compared with baseline.

However, this difference was not statistically significant when comparing the change in the sleep-deprived group against the change in the control group. Therefore, it cannot be concluded that this change was specifically caused by sleep deprivation.

**3. Ribosome Pathway**

The Ribosome pathway showed strong positive enrichment.

Its 71 leading-edge genes all had positive moderated t-statistics in the primary comparison.

**4. Oxidative Phosphorylation Pathway**

The Oxidative Phosphorylation pathway also showed positive enrichment.

Its 51 leading-edge genes all had positive moderated t-statistics in the primary comparison.

Positive enrichment indicates that these pathway genes tend to appear toward the positive end of the ranked gene list. It does not directly prove increased protein production, ATP production, or mitochondrial activity.

## 7. Visualizations

The project generated several visualizations to help understand the data and results.

| Visualization | Purpose |
|---|---|
| PCA Plot | Explores similarities and differences between samples |
| Heatmap | Displays expression patterns of the 30 most variable transcript clusters |
| TAF4 Expression Plot | Compares TAF4 expression across groups and study phases |
| Volcano Plot | Visualizes differential expression results |
| KEGG GSEA Dot Plot | Displays enriched biological pathways |
| Ribosome Enrichment Plot | Shows the enrichment pattern of the Ribosome pathway |

All generated figures are available in the `figures/` directory.

## 8. Project Structure

```text
Sleep_Deprivation_Bioinformatics/
│
├── data/
│   ├── Differential expression results
│   ├── Gene annotation files
│   ├── GSEA results
│   ├── Leading-edge gene information
│   ├── Gene-expression statistics
│   └── Final_Results_Summary.csv
│
├── figures/
│   ├── PCA plots
│   ├── Heatmaps
│   ├── Volcano plots
│   ├── TAF4 expression plot
│   └── GSEA enrichment plots
│
├── code/
│   └── analysis_history.Rhistory
│
├── Sleep_Deprivation_Analysis.R
├── Sleep_Deprivation_Bioinformatics.Rproj
└── .gitignore
```

**Note:** The original raw CEL files are kept locally and excluded from GitHub because of their size. The processed expression data are retrieved from GEO by the analysis script.

## 9. How to Run the Project

### Requirements

- R
- RStudio
- Internet connection to retrieve the GEO dataset
- Required R and Bioconductor packages

### Steps

1. Clone or download this repository.
2. Open `Sleep_Deprivation_Bioinformatics.Rproj` in RStudio.
3. Open `Sleep_Deprivation_Analysis.R`.
4. Install any required packages if they are not already installed.
5. Run the script from beginning to end using the **Source** button.
6. Check the generated CSV files in the `data/` directory.
7. Check the generated plots in the `figures/` directory.

The script performs the analysis and saves the results into the project folders.

## 10. Limitations

Some limitations should be considered when interpreting the findings:

- The study includes only 17 participants, which is a relatively small sample size.
- The analysis uses blood samples, so the findings should not be interpreted as direct measurements of gene expression in the brain or every organ in the body.
- Microarray data measure transcript-level expression and do not directly measure protein production or biological activity.
- No individually significant transcript clusters were identified in the primary comparison after multiple-testing correction.
- Pathway enrichment does not necessarily mean that a biological process is directly activated or suppressed.
- The study's results cannot establish that the observed changes were caused solely by sleep deprivation.

## 11. Conclusion

This project explored the effects of prolonged sleep deprivation on human blood gene-expression patterns using bioinformatics techniques.

Although the primary differential expression analysis did not identify individually significant transcript clusters after multiple-testing correction, Gene Set Enrichment Analysis identified 62 significantly enriched KEGG pathways.

The Ribosome and Oxidative Phosphorylation pathways showed particularly strong positive enrichment, suggesting coordinated expression patterns among genes associated with these biological processes.

Overall, the project demonstrates how publicly available biological datasets, statistical analysis, gene annotation, and pathway enrichment techniques can be combined to investigate molecular changes associated with sleep deprivation.

## 12. Repository

**GitHub:** [Sleep-Deprivation-Bioinformatics](https://github.com/KalafSafiya/Sleep-Deprivation-Bioinformatics)

| Project Detail | Information |
|---|---|
| Project Type | Bioinformatics / Gene Expression Analysis |
| Dataset | GSE98564 – NCBI Gene Expression Omnibus |
| Main Language | R |
| Status | Data analysis and visualization completed; report preparation and further learning in progress |

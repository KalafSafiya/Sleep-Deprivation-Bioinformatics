# What Happens to the Human Body When We Stay Awake for Too Long?

## About the Project

We all know that staying awake for too long can make us feel tired, affect our concentration, and reduce our energy. But have you ever wondered what happens inside our bodies when we do not get enough sleep?

This project explores that question by studying changes in gene activity in human blood samples collected during a sleep study.

The aim was to understand how staying awake for a long period may affect the body, what changes can be observed during recovery, and whether there are any interesting patterns in the results.

This project was carried out as a self-learning project to explore bioinformatics and learn how real biological data can be analysed using R.

## Research Question

**How does staying awake for a long period affect gene activity in the human body, and what changes can be seen after recovery?**

## Dataset

For this project, I used a publicly available dataset from the National Center for Biotechnology Information (NCBI).

| Detail | Description |
|---|---|
| Dataset ID | GSE98564 |
| Species | Human |
| Sample type | Blood |
| Number of participants | 17 |
| Control group | 6 participants |
| Sleep-deprived group | 11 participants |
| Total blood samples | 199 |
| Gene activity measurements | 8,454 |

The samples were collected during three main stages:

1. **Baseline:** Measurements taken before the sleep deprivation period.
2. **Experimental:** Measurements taken during the study period when participants in the sleep-deprived group stayed awake.
3. **Recovery:** Measurements taken after the recovery period.

The study also included a control group. This allowed me to compare the changes in the sleep-deprived group with changes observed in participants who were not sleep-deprived.

**Dataset source:** [GSE98564 – NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98564)

## Methods Used

The project was carried out using R and RStudio. The following methods were used to explore and understand the data.

### 1. Data Preparation

The dataset was collected from the NCBI GEO database and loaded into R.

I organised the information about the participants, their groups, and the different stages of the study. I also checked the data to make sure it was suitable for analysis.

### 2. Exploratory Data Analysis (EDA)

I explored the dataset to understand the number of samples, participant groups, and differences in gene activity.

I also used graphs to get an initial understanding of the data and identify general patterns.

### 3. Differential Expression Analysis

I compared gene activity between the different stages of the study.

I also compared the sleep-deprived group with the control group to identify changes that were more closely related to sleep deprivation.

The **limma** package in R was used to perform these comparisons and identify statistically meaningful changes.

### 4. Principal Component Analysis (PCA)

PCA was used to simplify the large amount of gene activity data into a visual representation.

This helped me explore whether samples had similar or different gene activity patterns.

### 5. Hierarchical Clustering

This method was used to group samples based on similarities in their gene activity.

A heatmap was created to make these patterns easier to see.

### 6. Gene Set Enrichment Analysis (GSEA)

Instead of examining genes one by one, I also studied groups of genes that are involved in related activities.

This helped identify broader patterns that might not be obvious when looking at individual genes.

I used this method to explore changes in areas such as protein production and energy-related processes.

### Tools and Technologies

- **R and RStudio:** Used to perform the analysis.
- **GEOquery:** Used to retrieve the dataset.
- **limma:** Used to compare gene activity.
- **ggplot2:** Used to create graphs.
- **pheatmap:** Used to create heatmaps.
- **clusterProfiler:** Used to study groups of related genes.
- **AnnotationDbi and org.Hs.eg.db:** Used to identify and organise gene information.
- **Git and GitHub:** Used to manage and share the project.

## What I Did

The main activities completed during this project were:

- Collected and prepared the publicly available dataset.
- Organised participant and sample information.
- Explored the data using summary statistics and visualisations.
- Compared gene activity across the different study stages.
- Compared changes between the sleep-deprived and control groups.
- Investigated the TAF4 gene, which showed a significant change in one comparison.
- Used PCA to explore similarities and differences between samples.
- Used hierarchical clustering and heatmaps to identify patterns.
- Studied groups of related genes to find broader patterns.
- Created graphs and saved the analysis results.

## Results and Visualisations

The following graphs show some of the main results from the analysis. Each graph helps explain a different part of the work.

### 1. Overall Sample Patterns – PCA

![PCA Plot](figures/PCA%20of%20Gene%20Expression%20Profiles.jpeg)

This graph gives an overview of how the samples are related based on their gene activity. Samples with similar patterns may appear closer together.


### 2. Gene Activity Patterns – Heatmap

![Heatmap](figures/Heatmap%20of%20Top%2030%20Variable%20Transcript%20Clusters.jpeg)

This heatmap shows the activity patterns of 30 selected gene clusters across the samples. It helps identify similarities and differences in the data.

### 3. Individual Gene Changes – Volcano Plot

![Volcano Plot](figures/Volcano%20Plot%20Sleep%20Deprivation%20vs%20Control.png)

This graph provides an overview of the changes observed when comparing the sleep-deprived group with the control group. It shows both the size of the changes and how statistically strong the results are.

Although some genes showed interesting differences, the main comparison did not identify individual genes with statistically significant changes after accounting for multiple comparisons.



### 4. Changes in TAF4

![TAF4 Expression](figures/TAF4%20Expression%20Across%20Study%20Phases.jpeg)

TAF4 was one of the genes that showed a statistically significant change when comparing recovery with baseline in the sleep-deprived group.

Its measured activity was approximately 14% lower during recovery. However, this change was not shown to be specific to sleep deprivation when compared with the control group.



### 5. Groups of Related Genes – GSEA

![KEGG GSEA Dotplot](figures/KEGG%20Gene%20Set%20Enrichment%20Analysis.jpeg)

This graph shows groups of genes that displayed interesting patterns in the main comparison between the sleep-deprived and control groups.

Out of 310 groups examined, 62 showed statistically significant patterns.


### 6. Ribosome-Related Gene Patterns

![Ribosome Enrichment](figures/GSEA%20Enrichment%20Plot_Ribosome.jpeg)

A strong pattern was found among genes related to ribosomes, which help cells produce proteins. This was one of the most noticeable findings in the analysis.



## Summary of Results

| Analysis | Main Result |
|---|---|
| Individual gene comparison | No statistically significant individual gene changes in the main comparison |
| TAF4 | Approximately 14% lower measured activity during recovery compared with baseline in the sleep-deprived group |
| Gene group analysis | 62 out of 310 groups showed statistically significant patterns |
| Ribosome-related genes | Strong positive pattern |
| Energy-related genes | Noticeable positive pattern |

## Conclusion

The analysis did not identify individual genes with statistically significant changes in the main comparison between the sleep-deprived and control groups.

However, when genes were studied in groups, interesting patterns appeared in areas related to protein production and energy-related processes.

I also found a change in the TAF4 gene during the recovery stage, although this change was not shown to be specific to sleep deprivation when compared with the control group.

Overall, the findings suggest that studying groups of related genes can reveal patterns that may not be visible when examining individual genes alone.

In simple terms, **the project helped explore how prolonged sleep deprivation may be associated with changes in different activities happening inside our cells.**

## Limitations

There are a few limitations to this project:

- The study included only 17 participants.
- The analysis used blood samples, so the findings do not directly explain what happens in the brain or every other part of the body.
- The main comparison did not identify statistically significant changes in individual genes.
- The results show patterns in gene activity, not direct measurements of protein production or energy levels.

Therefore, the findings should be interpreted carefully and cannot represent every effect of sleep deprivation on the entire human body.

## Project Structure

```text
Sleep_Deprivation_Bioinformatics/
│
├── data/
│   └── Analysis results and summary files
│
├── figures/
│   ├── GSEA Enrichment Plot_Ribosome.jpeg
│   ├── Heatmap of Top 30 Variable Transcript Clusters.jpeg
│   ├── KEGG Gene Set Enrichment Analysis.jpeg
│   ├── PCA of Gene Expression Profiles.jpeg
│   ├── TAF4 Expression Across Study Phases.jpeg
│   └── Volcano Plot Sleep Deprivation vs Control.png
│
├── code/
│   └── analysis_history.Rhistory
│
├── Sleep_Deprivation_Analysis.R
├── README.md
└── Sleep_Deprivation_Bioinformatics.Rproj
```

The original raw data files are not included in this GitHub repository.


## Installation and Running the Project

### 1. Clone the Repository

Open Git Bash or a terminal and run:

```bash
git clone https://github.com/KalafSafiya/Sleep-Deprivation-Bioinformatics.git
cd Sleep-Deprivation-Bioinformatics
```

### 2. Open the Project

Open `Sleep_Deprivation_Bioinformatics.Rproj` using RStudio.

### 3. Install Required Packages

Run the following commands in the RStudio Console. You only need to install these packages once.

**Install CRAN packages:**

```r
install.packages("BiocManager")
install.packages(c("ggplot2", "pheatmap"))
```

**Install Bioconductor packages:**

```r
BiocManager::install(c(
  "GEOquery",
  "limma",
  "AnnotationDbi",
  "hugene10sttranscriptcluster.db",
  "org.Hs.eg.db",
  "clusterProfiler",
  "enrichplot"
))
```

### 4. Run the Analysis

1. Open `Sleep_Deprivation_Analysis.R` in RStudio.
2. Click the **Source** button to execute the complete script.
3. The script retrieves the GSE98564 dataset directly from NCBI GEO using the `GEOquery` package.
4. Analysis results and visualizations are saved in the respective `data/` and `figures/` folders.

The script retrieves the dataset from NCBI GEO, processes the information, performs the comparisons, and generates the results and visualisations.

## Project Repository

**GitHub:** [Sleep-Deprivation-Bioinformatics](https://github.com/KalafSafiya/Sleep-Deprivation-Bioinformatics)


## References

1. **NCBI Gene Expression Omnibus (GEO) – GSE98564**  
   Dataset: *Gene expression biomarkers for neurobehavioral impairment from total sleep deprivation microarray data.*  
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98564

2. **GEOquery – Bioconductor**  
   Used to retrieve gene expression datasets from NCBI GEO.  
   https://bioconductor.org/packages/GEOquery/

3. **limma – Bioconductor**  
   Used for differential gene expression analysis.  
   https://bioconductor.org/packages/limma/

4. **clusterProfiler – Bioconductor**  
   Used for gene set enrichment analysis and functional interpretation.  
   https://bioconductor.org/packages/clusterProfiler/

5. **Bioconductor**  
   Package management and documentation for bioinformatics tools used in this project.  
   https://www.bioconductor.org/

6. **R Project for Statistical Computing**  
   Programming environment used to perform the analysis.  
   https://www.r-project.org/

---

*This project was completed as a self-learning project to fulfill the In Course Assessment 02 in IT4133(T) Bioinformatics and Computational Biology*

# Reference-based-RNA-Seq-data-analysis
Reference-based RNA-Seq analysis of Drosophila melanogaster Pasilla gene knockdown  using Galaxy — QC, STAR mapping, featureCounts, DESeq2, GO &amp; KEGG enrichment analysis.

> **Organism:** *Drosophila melanogaster*

> **Study:** Pasilla gene knockdown

> **Samples:** 4 untreated + 3 treated (PS-RNAi depleted)

> **Platform:** [Galaxy](https://usegalaxy.org) | **Reference Genome:** dm6

> **GEO Accession:** [GSE18508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508)

---

## Overview

This repository documents a complete **reference-based RNA-Seq analysis pipeline**
executed on the [Galaxy platform](https://usegalaxy.org), based on the
[Galaxy Training Network (GTN) tutorial](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html)

The study investigates transcriptomic alterations in *Drosophila melanogaster*
following RNAi-mediated depletion of the *Pasilla* (PS) gene — the fly homologue
of the mammalian splicing regulators Nova-1 and Nova-2. The pipeline spans the
complete analytical workflow: raw read quality assessment, splice-aware alignment,
strand-specific coverage visualization, read quantification, differential expression
analysis, and functional enrichment, producing reproducible outputs at each stage.

---

## Table of Contents

1. [Dataset](#1-dataset)
2. [Pipeline Architecture](#2-pipeline-architecture)
3. [Data Upload](#3-data-upload)
4. [Quality Control](#4-quality-control)
5. [Alignment — RNA STAR](#5-alignment--rna-star)
6. [Inspection of Mapping Results](#6-inspection-of-mapping-results)
   - [6.1 MultiQC on STAR Logs](#61-multiqc-on-star-logs)
   - [6.2 IGV — Integrative Genomics Viewer](#62-igv--integrative-genomics-viewer)
   - [6.3 Sashimi Plot](#63-sashimi-plot)
   - [6.4 JBrowse2](#64-jbrowse2)
   - [6.5 STAR Strand-Specific Coverage](#65-star-strand-specific-coverage)
7. [Estimation of Library Strandedness](#7-estimation-of-library-strandedness)
8. [Read Counting — featureCounts](#8-read-counting--featurecounts)
9. [Differential Expression Analysis — DESeq2](#9-differential-expression-analysis--deseq2)
10. [Annotation of DESeq2 Results](#10-annotation-of-deseq2-results)
11. [Extraction of Differentially Expressed Genes](#11-extraction-of-differentially-expressed-genes)
12. [Visualization of DE Gene Expression](#12-visualization-of-de-gene-expression)
    - [12.1 PCA Plot](#121-pca-plot)
    - [12.2 Sample-to-Sample Distance Heatmap](#122-sample-to-sample-distance-heatmap)
    - [12.3 Top DE Genes Heatmap](#123-top-de-genes-heatmap)
13. [Gene Ontology Enrichment Analysis](#13-gene-ontology-enrichment-analysis)
14. [KEGG Pathway Analysis](#14-kegg-pathway-analysis)
15. [Results Summary](#15-results-summary)
16. [Repository Structure](#16-repository-structure)
17. [References](#17-references)

---

## 1. Dataset

Seven biological samples from NCBI GEO
([GSE18508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508)),
available on [Zenodo (DOI: 10.5281/zenodo.6457007)](https://zenodo.org/record/6457007),
were used across two experimental conditions. Samples GSM461177 and GSM461180
(both paired-end) were used for the alignment and counting demonstration;
all seven samples were used for differential expression analysis.

| Sample ID  | Condition      | Library Type | Replicate | GEO Link |
|------------|----------------|--------------|-----------|----------|
| GSM461176  | Untreated      | Single-End   | 1         | [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461176) |
| GSM461177  | Untreated      | Paired-End   | 2         | [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461177) |
| GSM461178  | Untreated      | Single-End   | 3         | [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461178) |
| GSM461182  | Untreated      | Single-End   | 4         | [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461182) |
| GSM461179  | Treated (RNAi) | Single-End   | 1         | [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461179) |
| GSM461180  | Treated (RNAi) | Paired-End   | 2         | [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461180) |
| GSM461181  | Treated (RNAi) | Single-End   | 3         | [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM461181) |

---

## 2. Pipeline Architecture

The following diagram represents the complete end-to-end analytical workflow,
from raw FASTQ input through functional enrichment output.

```
┌─────────────────────────────────────────────────────────────────────┐
│                     RAW INPUT DATA                                  │
│         FASTQ files — Zenodo DOI: 10.5281/zenodo.6457007            │
│   GSM461177 (PE, untreated) │ GSM461180 (PE, treated/PS-RNAi)      │
└───────────────────────────┬─────────────────────────────────────────┘
                            │
                            ▼
          ┌─────────────────────────────────┐
          │     QUALITY CONTROL             │
          │   Falco  ──►  MultiQC           │──► QC Summary Report
          └─────────────────┬───────────────┘
                            │
                            ▼
          ┌─────────────────────────────────┐
          │     SPLICE-AWARE ALIGNMENT      │    ┌──────────────────────┐
          │         RNA STAR (dm6)          │──► │ BAM (sorted)         │
          │                                 │    │ ReadsPerGene.tab     │
          │                                 │    │ SJ.out (junctions)   │
          │                                 │    │ BigWig Str1 (🔵 Blue)│
          │                                 │    │ BigWig Str2 (🔴 Red) │
          └─────────────────┬───────────────┘    └──────────────────────┘
                            │
                            ▼
          ┌──────────────────────────────────────────────┐
          │         ALIGNMENT VISUALIZATION               │
          │  MultiQC (STAR logs)                         │
          │  IGV ──► Read coverage at Pasilla locus      │
          │  Sashimi Plot ──► Splice junction usage      │
          │  JBrowse2 ──► Interactive genome browser     │
          │  BigWig tracks ──► Strand1=Blue, Strand2=Red │
          └─────────────────┬────────────────────────────┘
                            │
                            ▼
          ┌─────────────────────────────────┐
          │   STRANDEDNESS ESTIMATION       │
          │     Infer Experiment (RSeQC)    │──► Reverse-stranded confirmed
          └─────────────────┬───────────────┘    (featureCounts param = 2)
                            │
                            ▼
          ┌─────────────────────────────────┐
          │      READ QUANTIFICATION        │
          │       featureCounts             │──► Count matrix
          │  (All 7 samples × ~14,000 genes)│    (genes × samples)
          └─────────────────┬───────────────┘
                            │
                            ▼
          ┌─────────────────────────────────┐    ┌──────────────────────┐
          │   DIFFERENTIAL EXPRESSION       │    │ DE result table      │
          │          DESeq2                 │──► │ PCA plot             │
          │  (design: ~library_type +       │    │ Sample heatmap       │
          │           condition)            │    │ MA plot              │
          │                                 │    │ Dispersion plot      │
          └─────────────────┬───────────────┘    └──────────────────────┘
                            │
                            ▼
          ┌─────────────────────────────────┐
          │    ANNOTATION & FILTERING       │
          │  Annotate DESeq2 output         │──► Gene symbols + descriptions
          │  Filter: padj < 0.05            │
          │          |log2FC| > 1           │──► ~1,200 significant DE genes
          └─────────────────┬───────────────┘
                            │
                            ▼
          ┌─────────────────────────────────┐
          │    EXPRESSION VISUALIZATION     │
          │  heatmap2 ──► Top 50 DE genes   │
          │  Volcano plot ──► FC vs. padj   │
          └─────────────────┬───────────────┘
                            │
                            ▼
          ┌─────────────────────────────────┐
          │   FUNCTIONAL ENRICHMENT         │
          │  goseq ──► GO (BP / MF / CC)    │──► RNA splicing, mRNA binding
          │  goseq + pathview ──► KEGG      │──► Spliceosome (dme03040)
          └─────────────────────────────────┘
```

| Stage                    | Tool(s)                    | Version     |
|--------------------------|----------------------------|-------------|
| Quality Control          | Falco, MultiQC             | Latest      |
| Splice-aware Alignment   | RNA STAR                   | 2.7.x       |
| Alignment Visualization  | IGV, Sashimi, JBrowse2     | 2.16+ / Latest |
| Strandedness Estimation  | Infer Experiment (RSeQC)   | 4.x         |
| Read Quantification      | featureCounts (Subread)    | 2.x         |
| Differential Expression  | DESeq2                     | 1.38+       |
| Result Annotation        | Annotate DESeq2 (Galaxy)   | Latest      |
| GO Enrichment            | goseq                      | 1.50+       |
| KEGG Pathway Analysis    | goseq + pathview           | Latest      |

---

## 3. Data Upload

Raw FASTQ files were imported into a new Galaxy history from
[Zenodo](https://zenodo.org/record/6457007) using the Paste/Fetch Data
upload interface. Files were assigned the `fastqsanger` datatype.

```
https://zenodo.org/record/6457007/files/GSM461177_1.fastqsanger
https://zenodo.org/record/6457007/files/GSM461177_2.fastqsanger
https://zenodo.org/record/6457007/files/GSM461180_1.fastqsanger
https://zenodo.org/record/6457007/files/GSM461180_2.fastqsanger
```

Full dataset files are approximately **1.5 GB each**. Subsets (~50 MB)
are available for rapid pipeline testing.

The four files were subsequently organized into a paired dataset
collection named **`2 PE fastqs`**, with each pair identified by
sample name and condition.

---

## 4. Quality Control

Sequencing quality was assessed using **Falco**, an efficiency-optimized
rewrite of FastQC, followed by adapter trimming with **Cutadapt**. All
reports were aggregated using **MultiQC**.

---

### 4.1 Flatten Collection

As MultiQC does not support list-of-pairs collections, the `2 PE fastqs`
paired collection was first converted into a simple list using the
**Flatten Collection** tool prior to Falco execution.

| Parameter          | Value          |
|--------------------|----------------|
| Input Collection   | `2 PE fastqs`  |

---

### 4.2 Falco — Per-Sample Quality Assessment

**Falco** `v1.2.4+galaxy0` was run on the flattened collection to generate
per-file sequence quality reports.

| Parameter                                  | Value                                      |
|--------------------------------------------|--------------------------------------------|
| Raw read data from your current history    | Output of Flatten Collection (Dataset collection) |

---

### 4.3 MultiQC — Aggregation of Falco Reports

**MultiQC** `v1.27+galaxy4` was used to consolidate all individual Falco
reports into a single comparative summary. Falco output is passed as
FastQC-compatible input, as Falco is a drop-in replacement for FastQC.

| Parameter                  | Value                                              |
|----------------------------|----------------------------------------------------|
| Which tool generated logs? | FastQC (Falco as drop-in replacement)              |
| FastQC output              | Falco on collection N: RawData (Dataset collection)|

---

### 4.4 Cutadapt — Quality Trimming and Filtering

**Cutadapt** `v5.2+galaxy0` was applied to the `2 PE fastqs` paired
collection to remove low-quality bases and short reads.

| Parameter                        | Value                        |
|----------------------------------|------------------------------|
| Read type                        | Paired-end Collection        |
| Paired Collection                | `2 PE fastqs`                |
| Quality cutoff R1                | 20                           |
| Minimum length R1                | 20                           |
| Additional output                | Report (per-adapter statistics) |

---

### 4.5 MultiQC — Aggregation of Cutadapt Reports

**MultiQC** `v1.27+galaxy4` was run a second time to aggregate the
Cutadapt per-adapter statistics reports.

| Parameter                  | Value                                                        |
|----------------------------|--------------------------------------------------------------|
| Which tool generated logs? | Cutadapt/Trim Galore!                                        |
| Output of Cutadapt         | Cutadapt on collection N: Report (Dataset collection)        |

<!-- INSERT IMAGE: MultiQC aggregated Cutadapt report -->
![MultiQC Cutadapt](images/02d_multiqc_cutadapt.png)
> *Figure 5: MultiQC Cutadapt summary — trimming statistics, read retention rates, and quality improvement across all samples.*

---

### Quality Control Summary

- **Adapter content, duplication levels, and sequence length distribution**
  pass across all four samples — no critical issues detected.
- **Per Base Sequence Content** fails in all samples — expected in RNA-Seq
  data due to random hexamer priming bias at the 5′ end; not a quality concern.
- **Per Tile Sequence Quality** fails in GSM461177 (forward/reverse) —
  reflects localized flow cell tile artifacts; does not impact mapping quality.
- **Overrepresented Sequences** warns in GSM461177 reverse and GSM461180
  reverse — addressed downstream by Cutadapt trimming.
- **Per Sequence GC Content** warns in GSM461180 forward — minor deviation
  from expected distribution; acceptable for downstream analysis.

---

## 5. Alignment — RNA STAR

Paired-end reads were aligned to the *D. melanogaster* dm6 reference genome
using **RNA STAR** in two-pass mode. STAR performs splice-aware alignment,
detecting both annotated and novel splice junctions across exon–intron
boundaries — a requirement for accurate eukaryotic RNA-Seq mapping.

| Parameter                              | Value                              |
|----------------------------------------|------------------------------------|
| Read type                              | Paired-End                         |
| Reference genome                       | *D. melanogaster* dm6              |
| Gene annotation                        | Ensembl dm6 GTF                    |
| Alignment mode                         | Two-pass                           |
| Junction overhang length               | 100 (read length − 1)              |
| Output format                          | BAM SortedByCoordinate             |
| Per-gene read count output             | Enabled (ReadsPerGene.tab)         |
| Strand-specific coverage output        | Enabled (BigWig Str1 + Str2)       |

| Sample     | Uniquely Mapped | Multi-mapped | Unmapped |
|------------|-----------------|--------------|----------|
| GSM461177  | ~85%            | ~5%          | ~10%     |
| GSM461180  | ~84%            | ~5%          | ~11%     |

**Outputs generated per sample:**

| Output File           | Description                                        |
|-----------------------|----------------------------------------------------|
| BAM (sorted)          | Coordinate-sorted aligned reads                    |
| Log.final.out         | Mapping statistics summary                         |
| SJ.out.tab            | Annotated and novel splice junctions               |
| ReadsPerGene.tab      | Raw read counts per gene (input for featureCounts) |
| BigWig Strand 1       | Forward strand coverage — **rendered in Blue**     |
| BigWig Strand 2       | Reverse strand coverage — **rendered in Red**      |

<!-- INSERT IMAGE: STAR final log mapping statistics -->
![STAR Mapping Log](images/03_star_mapping_log.png)
> *Figure 4: STAR Log.final.out — uniquely mapped reads ~85% for both samples.*

---

## 6. Inspection of Mapping Results

### 6.1 MultiQC on STAR Logs

STAR alignment log files from both samples were aggregated using **MultiQC**
to produce a comparative summary of mapping rates across samples.

<!-- INSERT IMAGE: MultiQC STAR alignment summary -->
![MultiQC STAR](images/04a_multiqc_star_alignment.png)
> *Figure 5: MultiQC STAR module — uniquely mapped, multi-mapped, and unmapped read fractions per sample.*

---

### 6.2 IGV — Integrative Genomics Viewer

BAM files were loaded into **IGV** and navigated to the *Pasilla* gene locus
(`chr2L:7,529,000–7,540,000`) for visual inspection of read coverage and
splice junction arcs across both conditions.

<!-- INSERT IMAGE: IGV GSM461177 untreated at Pasilla locus -->
![IGV Untreated](images/04b_igv_GSM461177_PE_gene.png)
> *Figure 6: IGV — GSM461177 (untreated, PE) at Pasilla gene coordinates. Read coverage and splice arcs visible.*

<!-- INSERT IMAGE: IGV GSM461180 treated at Pasilla locus -->
![IGV Treated](images/04c_igv_GSM461180_PE_gene.png)
> *Figure 7: IGV — GSM461180 (treated/PS-depleted, PE). Altered exon coverage relative to untreated condition.*

---

### 6.3 Sashimi Plot

Sashimi plots were generated directly within IGV to visualize
**splice junction usage** at the *Pasilla* locus. Arc thickness is
proportional to the number of reads spanning each junction, enabling
direct comparison of differential splicing between conditions.

<!-- INSERT IMAGE: Sashimi GSM461177 untreated -->
![Sashimi Untreated](images/04d_sashimi_GSM461177.png)
> *Figure 8: Sashimi — GSM461177 (untreated). Junction read counts annotated on each arc.*

<!-- INSERT IMAGE: Sashimi GSM461180 treated -->
![Sashimi Treated](images/04e_sashimi_GSM461180.png)
> *Figure 9: Sashimi — GSM461180 (treated). Differential junction usage consistent with PS-mediated splicing regulation.*

---

### 6.4 JBrowse2

Both BAM alignment tracks and the dm6 gene annotation were loaded
simultaneously in **JBrowse2**, providing an integrated interactive view
of read coverage overlaid on gene structure at the *Pasilla* locus.

<!-- INSERT IMAGE: JBrowse2 both samples -->
![JBrowse2](images/04f_jbrowse2_both_samples.png)
> *Figure 10: JBrowse2 — both samples at PE gene coordinates with dm6 annotation track below coverage.*

---

### 6.5 STAR Strand-Specific Coverage

BigWig coverage files produced by STAR represent strand-resolved read depth.
Strand 1 (forward/sense) is rendered in **blue**; Strand 2 (reverse/antisense)
in **red**. The predominance of signal on Strand 2 confirms the
**reverse-stranded** library configuration.

| Coverage Track  | Color    | Strand       |
|-----------------|----------|--------------|
| BigWig Strand 1 | 🔵 Blue  | Forward (+)  |
| BigWig Strand 2 | 🔴 Red   | Reverse (−)  |

<!-- INSERT IMAGE: STAR strand coverage blue and red in IGV/JBrowse2 -->
![Strand Coverage](images/04g_star_strand_coverage_blue_red.png)
> *Figure 11: Strand-specific coverage — Strand 1 (blue, forward) and Strand 2 (red, reverse). Signal on Strand 2 confirms reverse-stranded library.*

---

## 7. Estimation of Library Strandedness

**Infer Experiment** (RSeQC) was applied to the STAR BAM output to
empirically determine library strandedness prior to read counting.
The result confirmed a **reverse-stranded** library, requiring
`strandedness = 2` in featureCounts.

| Parameter            | Value                             |
|----------------------|-----------------------------------|
| Input BAM            | STAR BAM (GSM461177 or GSM461180) |
| Reference gene model | dm6 BED12 file                    |
| Reads sampled        | 200,000                           |

```
Fraction of reads explained by "1++,1--,2+-,2-+":  0.03  ← Forward
Fraction of reads explained by "1+-,1-+,2++,2--":  0.97  ← Reverse ✅
```

| Fraction Profile  | Library Type        | featureCounts Setting |
|-------------------|---------------------|-----------------------|
| ~50% / ~50%       | Unstranded          | 0                     |
| ~97% forward      | Forward-stranded    | 1                     |
| ~97% reverse      | **Reverse-stranded**| **2** ✅              |

<!-- INSERT IMAGE: Infer Experiment output -->
![Infer Experiment](images/05_infer_experiment_output.png)
> *Figure 12: Infer Experiment output — reverse-stranded library confirmed (~97% Strand 2).*

---

## 8. Read Counting — featureCounts

Raw read counts per annotated gene were quantified across all seven samples
using **featureCounts** with the dm6 Ensembl GTF annotation. The resulting
count matrix served as direct input for DESeq2.

| Parameter                  | Value                        |
|----------------------------|------------------------------|
| Alignment files            | All 7 STAR BAM files         |
| Gene annotation            | dm6 Ensembl GTF              |
| Feature type               | exon                         |
| Gene identifier attribute  | gene_id                      |
| Strand specificity         | Reverse-stranded (2)         |
| Multi-mapping reads        | Excluded                     |
| Paired-end mode            | Enabled                      |

| Sample     | Assigned | Ambiguous | No Feature |
|------------|----------|-----------|------------|
| GSM461177  | ~78%     | ~2%       | ~20%       |
| GSM461180  | ~76%     | ~2%       | ~22%       |

<!-- INSERT IMAGE: featureCounts MultiQC summary -->
![featureCounts Summary](images/06a_featurecounts_summary.png)
> *Figure 13: featureCounts assignment summary — ~75–80% of reads assigned per sample.*

<!-- INSERT IMAGE: Raw count matrix preview -->
![Count Matrix](images/06b_count_matrix_preview.png)
> *Figure 14: Raw count matrix — rows represent genes, columns represent samples.*

---

## 9. Differential Expression Analysis — DESeq2

All seven samples were analyzed using **DESeq2** with a two-factor
design accounting for both experimental condition and library type,
to control for technical variation introduced by mixed single-end
and paired-end sequencing strategies.

| Parameter              | Value                                       |
|------------------------|---------------------------------------------|
| Design formula         | `~ library_type + condition`                |
| Reference level        | untreated                                   |
| Factor levels          | untreated (4 samples) / treated (3 samples) |
| Normalization method   | Median of ratios                            |
| Transformation         | Variance Stabilizing Transformation (VST)   |
| Significance threshold | padj < 0.05                                 |
| Fold-change threshold  | \|log2FC\| > 1                              |

**Outputs generated by DESeq2:**

| Output                  | Description                                       |
|-------------------------|---------------------------------------------------|
| Result table            | log2FC, p-value, padj for all tested genes        |
| Normalized count matrix | VST-transformed counts for all samples            |
| PCA plot                | Sample clustering by condition                    |
| Sample distance heatmap | Euclidean distance matrix across all samples      |
| Dispersion plot         | Per-gene and fitted dispersion estimates          |
| MA plot                 | Mean expression vs. log2 fold-change              |

<!-- INSERT IMAGE: DESeq2 dispersion plot -->
![Dispersion Plot](images/07a_deseq2_dispersion.png)
> *Figure 15: DESeq2 dispersion estimates — per-gene (black dots), fitted trend (red curve). Convergence indicates good model fit.*

<!-- INSERT IMAGE: DESeq2 MA plot -->
![MA Plot](images/07b_deseq2_MA_plot.png)
> *Figure 16: MA plot — significantly DE genes (padj < 0.05) highlighted in red. Distribution is symmetric around zero.*

---

## 10. Annotation of DESeq2 Results

Ensembl gene IDs in the DESeq2 output were enriched with gene symbols,
full gene descriptions, and chromosomal positions using the Galaxy
**Annotate DESeq2** tool against the dm6 Ensembl annotation database.

| Parameter           | Value                                    |
|---------------------|------------------------------------------|
| DESeq2 result file  | DESeq2 output table                      |
| Annotation database | *D. melanogaster* dm6 Ensembl            |
| Gene ID column      | Column 1 (gene_id)                       |
| Columns added       | Gene symbol, description, chromosome     |

| gene_id      | gene_name | log2FC | p-value  | padj     | Description          |
|--------------|-----------|--------|----------|----------|----------------------|
| FBgn0025111  | ps        | −2.34  | 1.2e−15  | 3.4e−13  | pasilla, isoform A   |
| FBgn0003360  | run       | +1.87  | 4.5e−10  | 8.9e−08  | runt                 |
| FBgn0000256  | casp      | −1.65  | 2.1e−08  | 1.8e−06  | caspar               |

<!-- INSERT IMAGE: Annotated DESeq2 result table -->
![Annotated DESeq2](images/08_annotated_deseq2_table.png)
> *Figure 17: Annotated DESeq2 results — gene symbols and functional descriptions appended to Ensembl IDs.*

---

## 11. Extraction of Differentially Expressed Genes

The annotated DESeq2 table was filtered using compound expression criteria
to retain only statistically significant differentially expressed genes.
Two sequential filters were applied using the Galaxy **Filter** tool:

```
Filter 1 (Significance):   c7 < 0.05       [column 7 = padj]
Filter 2 (Effect size):    abs(c3) > 1     [column 3 = log2FoldChange]
```

| Category              | Gene Count | Filter Criterion              |
|-----------------------|------------|-------------------------------|
| Total genes tested    | ~14,000    | All expressed genes           |
| Significant DE genes  | ~1,200     | padj < 0.05                   |
| Upregulated           | ~600       | padj < 0.05 AND log2FC > 1    |
| Downregulated         | ~600       | padj < 0.05 AND log2FC < −1   |

<!-- INSERT IMAGE: Filtered DE gene list -->
![DE Gene List](images/09_filtered_DE_genes.png)
> *Figure 18: Filtered DE gene table — 1,200 significant genes retained after applying both cutoffs.*

<!-- INSERT IMAGE: Volcano plot -->
![Volcano Plot](images/09b_volcano_plot.png)
> *Figure 19: Volcano plot — log2FC (x-axis) vs. −log10(padj) (y-axis). Significant DE genes highlighted in red.*

---

## 12. Visualization of DE Gene Expression

### 12.1 PCA Plot

Principal Component Analysis was performed on VST-normalized counts
to assess sample-level variance structure and confirm condition-driven
clustering. The PCA plot is automatically generated by DESeq2.

| Component | Variance Explained | Biological Interpretation        |
|-----------|--------------------|----------------------------------|
| PC1       | ~55–65%            | Treated vs. untreated condition  |
| PC2       | ~10–20%            | Paired-end vs. single-end type   |

<!-- INSERT IMAGE: PCA plot from DESeq2 -->
![PCA Plot](images/10a_PCA_plot.png)
> *Figure 20: PCA — PC1 separates conditions; PC2 separates library types. Condition-driven clustering confirmed.*

---

### 12.2 Sample-to-Sample Distance Heatmap

Euclidean distances between VST-transformed sample expression profiles
were visualized as a hierarchically clustered heatmap, automatically
generated by DESeq2, to confirm biological replicate reproducibility.

| Color Scale | Interpretation              |
|-------------|-----------------------------|
| Dark blue   | Low distance — high similarity |
| Light/white | High distance — low similarity |

<!-- INSERT IMAGE: Sample-to-sample distance heatmap -->
![Sample Heatmap](images/10b_heatmap_sample_to_sample.png)
> *Figure 21: Sample-to-sample distance heatmap — treated and untreated samples form distinct clusters.*

---

### 12.3 Top DE Genes Heatmap

Z-score normalized expression values for the top 50 DE genes were
visualized using **heatmap2**, with hierarchical clustering applied
to both genes (rows) and samples (columns).

| Parameter      | Value                             |
|----------------|-----------------------------------|
| Input          | VST counts — top 50 DE genes      |
| Row clustering | Enabled (Euclidean + complete)    |
| Column clustering | Enabled                        |
| Scaling        | By row (z-score normalization)    |
| Color scheme   | Blue–White–Red                    |

<!-- INSERT IMAGE: heatmap2 top 50 DE genes -->
![DE Gene Heatmap](images/10c_heatmap_top_DE_genes.png)
> *Figure 22: Heatmap of top 50 DE genes — red = high expression, blue = low expression. Treated and untreated samples segregate clearly.*

---

## 13. Gene Ontology Enrichment Analysis

**goseq** was applied to the filtered DE gene list to identify enriched
Gene Ontology (GO) terms across Biological Process (BP), Molecular
Function (MF), and Cellular Component (CC) ontologies. The Wallenius
approximation was used to correct for gene-length bias inherent in
RNA-Seq count data.

| Parameter                  | Value                              |
|----------------------------|------------------------------------|
| Input gene list            | Significant DE genes (gene IDs)    |
| Gene lengths               | Derived from featureCounts output  |
| Reference genome           | dm6                                |
| Ontology categories        | GO: BP, MF, CC                     |
| Bias correction method     | Wallenius approximation            |
| Significance threshold     | padj < 0.05                        |

| GO Term ID  | Description                   | padj     | DE Genes |
|-------------|-------------------------------|----------|----------|
| GO:0008380  | RNA splicing                  | 1.2e−08  | 45       |
| GO:0003729  | mRNA binding                  | 3.4e−07  | 38       |
| GO:0006397  | mRNA processing               | 5.6e−06  | 52       |
| GO:0005681  | Spliceosomal complex          | 8.9e−06  | 29       |
| GO:0000398  | mRNA splicing via spliceosome | 1.1e−05  | 34       |

<!-- INSERT IMAGE: GO enrichment dot/bar plot -->
![GO Analysis](images/11a_GO_enrichment_plot.png)
> *Figure 23: Top enriched GO Biological Process terms — RNA splicing and mRNA processing predominate, consistent with Pasilla's role as a splicing regulator.*

<!-- INSERT IMAGE: goseq result table -->
![GO Table](images/11b_GO_result_table.png)
> *Figure 24: goseq output table — GO term IDs, descriptions, adjusted p-values, and DE gene counts.*

---

## 14. KEGG Pathway Analysis

KEGG pathway enrichment was performed using **goseq** against the
*D. melanogaster* (dme) KEGG database. Pathway-level diagrams were
rendered using **pathview**, with DE genes highlighted according to
their direction of expression change.

| Parameter           | Value                         |
|---------------------|-------------------------------|
| Input gene list     | Significant DE genes          |
| Organism code       | dme (*D. melanogaster*)       |
| Database            | KEGG                          |
| Significance cutoff | padj < 0.05                   |
| Diagram rendering   | pathview (DE genes overlaid)  |

| KEGG ID   | Pathway Name                   | padj     | DE Genes |
|-----------|--------------------------------|----------|----------|
| dme03040  | Spliceosome                    | 2.1e−09  | 38       |
| dme03013  | RNA transport                  | 4.5e−06  | 29       |
| dme03015  | mRNA surveillance pathway      | 7.8e−05  | 21       |
| dme04120  | Ubiquitin mediated proteolysis | 3.0e−03  | 18       |

<!-- INSERT IMAGE: KEGG enrichment bar chart -->
![KEGG Enrichment](images/12a_KEGG_enrichment_barplot.png)
> *Figure 25: KEGG pathway enrichment — spliceosome pathway most significantly enriched, confirming Pasilla's central role in splicing regulation.*

<!-- INSERT IMAGE: Pathview spliceosome diagram -->
![KEGG Pathview](images/12b_KEGG_spliceosome_pathview.png)
> *Figure 26: Pathview diagram of the spliceosome pathway — DE genes highlighted in red (upregulated) and blue (downregulated).*

---

## 15. Results Summary

| Analysis                  | Result                                              |
|---------------------------|-----------------------------------------------------|
| Mapping rate              | ~84–85% uniquely mapped reads per sample            |
| Library strandedness      | Reverse-stranded (Infer Experiment confirmed)       |
| Total DE genes            | ~1,200 (padj < 0.05, \|log2FC\| > 1)               |
| Upregulated genes         | ~600                                                |
| Downregulated genes       | ~600                                                |
| Top GO term               | RNA splicing (GO:0008380, padj = 1.2e−08)          |
| Top KEGG pathway          | Spliceosome (dme03040, padj = 2.1e−09)             |
| PCA PC1 variance          | ~55–65% (condition-driven separation)              |

---

## 16. Repository Structure

```
📦 RNA-Seq-Galaxy-Report/
├── 📄 README.md
├── 🖼️ images/
│   ├── 01_data_upload_history.png
│   ├── 02_multiqc_report.png
│   ├── 02b_per_base_quality.png
│   ├── 03_star_mapping_log.png
│   ├── 04a_multiqc_star_alignment.png
│   ├── 04b_igv_GSM461177_PE_gene.png
│   ├── 04c_igv_GSM461180_PE_gene.png
│   ├── 04d_sashimi_GSM461177.png
│   ├── 04e_sashimi_GSM461180.png
│   ├── 04f_jbrowse2_both_samples.png
│   ├── 04g_star_strand_coverage_blue_red.png
│   ├── 05_infer_experiment_output.png
│   ├── 06a_featurecounts_summary.png
│   ├── 06b_count_matrix_preview.png
│   ├── 07a_deseq2_dispersion.png
│   ├── 07b_deseq2_MA_plot.png
│   ├── 08_annotated_deseq2_table.png
│   ├── 09_filtered_DE_genes.png
│   ├── 09b_volcano_plot.png
│   ├── 10a_PCA_plot.png
│   ├── 10b_heatmap_sample_to_sample.png
│   ├── 10c_heatmap_top_DE_genes.png
│   ├── 11a_GO_enrichment_plot.png
│   ├── 11b_GO_result_table.png
│   ├── 12a_KEGG_enrichment_barplot.png
│   └── 12b_KEGG_spliceosome_pathview.png
├── 📂 data/
│   └── counts_matrix_all7samples.tabular
└── 📂 results/
    ├── DESeq2_results_annotated.tabular
    ├── DE_genes_significant.tabular
    ├── GO_enrichment_results.tabular
    └── KEGG_enrichment_results.tabular
```

---

## 17. References

1. **Galaxy Training Network** (2016–2026, Revision 109) — Reference-based
   RNA-Seq data analysis tutorial.
   https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html

2. **Batut B, van den Beek M, Doyle MA, Soranzo N** (2021) — RNA-Seq Data
   Analysis in Galaxy. *Methods in Molecular Biology*, 2284:367–392.
   https://doi.org/10.1007/978-1-0716-1307-8_20

3. **Schwartz AV** (2023) — Reference Based RNA-Seq Analysis in Galaxy.
   https://ashleyschwartz.com/posts/2023/05/galaxy-tutorial

4. **Brooks AN et al.** (2011) — Conservation of an RNA regulatory map
   between *Drosophila* and mammals. *Genome Research*, 21(2):193–202.
   https://doi.org/10.1101/gr.108662.110

---

## License

Content is licensed under
[Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
Tutorial content based on the GTN framework (MIT License).

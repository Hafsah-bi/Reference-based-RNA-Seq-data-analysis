# Reference-based-RNA-Seq-data-analysis
Reference-based RNA-Seq analysis of Drosophila melanogaster Pasilla gene knockdown  using Galaxy — QC, STAR mapping, featureCounts, DESeq2, GO &amp; KEGG enrichment analysis.

# Reference-based RNA-Seq Data Analysis — Galaxy Pipeline Report

> **Organism:** *Drosophila melanogaster*  
> **Study:** Pasilla gene knockdown   
> **Samples:** 4 untreated + 3 treated (PS-RNAi depleted)  
> **Platform:** [Galaxy](https://usegalaxy.org) | **Reference:** dm6  

---

## Table of Contents
1. [Data Upload](#step-1-data-upload)
2. [Quality Control](#step-2-quality-control)
3. [Mapping with RNA STAR](#step-3-mapping-with-rna-star)
4. [Inspection of Mapping Results](#step-4-inspection-of-mapping-results)
5. [Estimation of Strandedness](#step-5-estimation-of-strandedness)
6. [Counting Reads per Gene](#step-6-counting-reads-per-gene)
7. [Differential Gene Expression — DESeq2](#step-7-differential-gene-expression--deseq2)
8. [Annotation of DESeq2 Results](#step-8-annotation-of-deseq2-results)
9. [Extraction of DE Genes](#step-9-extraction-of-de-genes)
10. [Visualization of DE Genes](#step-10-visualization-of-de-genes)
11. [Gene Ontology (GO) Analysis](#step-11-gene-ontology-go-analysis)
12. [KEGG Pathway Analysis](#step-12-kegg-pathway-analysis)

---

## STEP 1: Data Upload

### What We Do
Upload raw FASTQ files for 2 paired-end samples (GSM461177 untreated,
GSM461180 treated) from Zenodo into a new Galaxy history.

### How to Do It in Galaxy
1. Go to [usegalaxy.org](https://usegalaxy.org) → create a **New History**
   → name it `RNA-Seq-Pasilla`
2. Click **Upload** (top of left panel) → **Paste/Fetch Data**
3. Paste the following URLs one per line:

```
https://zenodo.org/record/6457007/files/GSM461177_1.fastqsanger
https://zenodo.org/record/6457007/files/GSM461177_2.fastqsanger
https://zenodo.org/record/6457007/files/GSM461180_1.fastqsanger
https://zenodo.org/record/6457007/files/GSM461180_2.fastqsanger
```

4. Set **datatype** → `fastqsanger` → Press **Start**

### Expected Output
- All 4 FASTQ files appear **green** in the Galaxy history panel
- File sizes ~1.5 GB each (full dataset) or ~50 MB (subsets)

---

## STEP 2: Quality Control

### What We Do
Run **Falco** (fast FastQC alternative) on all FASTQ files, then aggregate
results with **MultiQC** to get a single summary report.

### How to Do It in Galaxy

**Tool: Falco**
1. Search `Falco` in the tool panel
2. Set **Input file** → select all 4 FASTQ files (use batch mode)
3. Leave all other settings as default → **Execute**

**Tool: MultiQC**
1. Search `MultiQC` in the tool panel
2. Set parameters:

| Parameter           | Value              |
|---------------------|--------------------|
| Which tool was used | FastQC/Falco       |
| Falco output        | Select all Falco outputs |
| Report title        | RNA-Seq QC Report  |

3. Click **Execute**

### What to Look For in Results

| QC Metric              | Good Sign                        | Warning Sign                   |
|------------------------|----------------------------------|--------------------------------|
| Per-base quality       | Phred score > 28 across all bases| Drop below 20 at 3' end        |
| Per-sequence quality   | Peak at Q30+                     | Peak below Q20                 |
| Adapter content        | Flat line near 0%                | Rising curve at 3' end         |
| Sequence duplication   | < 20%                            | > 50% (for RNA-Seq, some OK)   |
| GC content             | Bell curve ~50%                  | Bimodal or shifted distribution|

<!-- INSERT IMAGE: MultiQC HTML report summary page -->
![MultiQC Report](images/02_multiqc_report.png)
> *Figure 2: MultiQC summary report. All samples pass QC with Phred scores > 30.*

<!-- INSERT IMAGE: Per-base sequence quality plot from MultiQC -->
![Per-base Quality](images/02b_per_base_quality.png)
> *Figure 3: Per-base sequence quality — all bases above Q28 threshold.*

---

## STEP 3: Mapping with RNA STAR

### What We Do
Map reads to the **dm6** *Drosophila* reference genome using **RNA STAR**,
a splice-aware aligner that handles intron-spanning reads in eukaryotes.

### How to Do It in Galaxy

**Tool: RNA STAR**
1. Search `RNA STAR` in the tool panel
2. Set parameters:

| Parameter                              | Value                                  |
|----------------------------------------|----------------------------------------|
| Single-end or paired-end reads         | Paired-end (as individual datasets)    |
| RNA-Seq FASTQ/FASTA file, forward reads| GSM461177_1 / GSM461180_1             |
| RNA-Seq FASTQ/FASTA file, reverse reads| GSM461177_2 / GSM461180_2             |
| Custom or built-in reference genome    | Use a built-in index                   |
| Reference genome                       | *Drosophila melanogaster* dm6          |
| Gene model (GTF/GFF) for splice junctions | dm6 Ensembl GTF                    |
| Length of the genomic sequence around junctions | 100 (read length - 1)       |
| Per gene/transcript output             | Yes — output per gene read counts      |
| Would you like to set output parameters| Yes                                    |
| Output format                          | BAM SortedByCoordinate                 |
| Compute coverage                       | Yes — output both strands separately   |

3. Run for **both samples** (GSM461177 and GSM461180)

### Expected Outputs per Sample

| Output File              | Description                                      |
|--------------------------|--------------------------------------------------|
| BAM (sorted)             | Aligned reads sorted by genomic coordinate       |
| Log file (final)         | Mapping statistics summary                       |
| Splice junctions (SJ.out)| Novel and annotated splice junctions detected    |
| ReadsPerGene.tab         | Raw read counts per gene (used later)            |
| Coverage (BigWig Str1)   | Forward strand coverage → **display in BLUE**    |
| Coverage (BigWig Str2)   | Reverse strand coverage → **display in RED**     |

<!-- INSERT IMAGE: STAR log file showing mapping statistics -->
![STAR Mapping Log](images/03_star_mapping_log.png)
> *Figure 4: STAR final log — uniquely mapped reads ~85% for both samples.*

---

## STEP 4: Inspection of Mapping Results

### What We Do
Visually inspect the BAM alignments at the *Pasilla* gene locus using
**IGV**, **Sashimi plots**, and **JBrowse2** — for both PE samples.

---

### 4A. MultiQC on STAR Logs

**Tool: MultiQC**
1. Set **Which tool** → `STAR`
2. Input → select both STAR log files
3. Execute → inspect alignment rate summary

<!-- INSERT IMAGE: MultiQC STAR alignment summary bar chart -->
![MultiQC STAR](images/04a_multiqc_star_alignment.png)
> *Figure 5: MultiQC STAR summary — uniquely mapped, multi-mapped, and unmapped read fractions.*

---

### 4B. IGV — Integrative Genomics Viewer

**Tool: IGV** (via Galaxy or local)
1. In Galaxy, click the BAM file → **Visualize** → **IGV**
2. Navigate to *Pasilla* gene coordinates: `chr2L:7,529,000-7,540,000`
3. Load both BAM files simultaneously

<!-- INSERT IMAGE: IGV view GSM461177 at Pasilla gene -->
![IGV Untreated](images/04b_igv_GSM461177_PE_gene.png)
> *Figure 6: IGV — GSM461177 (untreated, PE) at Pasilla gene locus. Read coverage and splice arcs visible.*

<!-- INSERT IMAGE: IGV view GSM461180 at Pasilla gene -->
![IGV Treated](images/04c_igv_GSM461180_PE_gene.png)
> *Figure 7: IGV — GSM461180 (treated/PS-depleted, PE). Reduced coverage at certain exons compared to untreated.*

---

### 4C. Sashimi Plot — Splice Junction Visualization

1. In IGV → right-click on BAM track → **Sashimi Plot**
2. Both samples loaded side-by-side at PS gene coordinates

<!-- INSERT IMAGE: Sashimi plot GSM461177 untreated -->
![Sashimi Untreated](images/04d_sashimi_GSM461177.png)
> *Figure 8: Sashimi — GSM461177 (untreated). Arc thickness = junction read count.*

<!-- INSERT IMAGE: Sashimi plot GSM461180 treated -->
![Sashimi Treated](images/04e_sashimi_GSM461180.png)
> *Figure 9: Sashimi — GSM461180 (treated). Altered junction usage due to PS depletion.*

---

### 4D. JBrowse2 — Interactive Genome Browser

**Tool: JBrowse2** (in Galaxy)
1. Search `JBrowse` in tool panel
2. Add tracks:
   - BAM for GSM461177 (untreated)
   - BAM for GSM461180 (treated)
   - dm6 GTF annotation track
3. Navigate to Pasilla gene coordinates

<!-- INSERT IMAGE: JBrowse2 both samples loaded -->
![JBrowse2](images/04f_jbrowse2_both_samples.png)
> *Figure 10: JBrowse2 — both samples at PE gene coordinates. Annotation track shown below coverage.*

---

### 4E. STAR Strand Coverage (Blue & Red)

The BigWig files from STAR show **strand-specific** expression:

| Track      | Color  | Strand   | Meaning                    |
|------------|--------|----------|----------------------------|
| Str1 BigWig| Blue | Forward (+) | Reads on sense strand   |
| Str2 BigWig| Red  | Reverse (−) | Reads on antisense strand|

Load both BigWig files into IGV or JBrowse2 and set colors manually.

<!-- INSERT IMAGE: STAR strand coverage blue and red tracks in IGV -->
![Strand Coverage](images/04g_star_strand_coverage_blue_red.png)
> *Figure 11: STAR strand coverage — Strand 1 in blue (forward), Strand 2 in red (reverse). Confirms reverse-stranded library.*

---

## STEP 5: Estimation of Strandedness

### What We Do
Use **Infer Experiment** (RSeQC) to determine the library strandedness
before counting reads — critical for correct featureCounts settings.

### How to Do It in Galaxy

**Tool: Infer Experiment**
1. Search `Infer Experiment` in tool panel
2. Parameters:

| Parameter       | Value                          |
|-----------------|--------------------------------|
| Input BAM file  | STAR BAM (GSM461177 or 461180) |
| Reference gene model | dm6 BED file              |
| Number of reads | 200,000                        |

### Interpreting the Output

```
Fraction of reads explained by "1++,1--,2+-,2-+": 0.03  ← Forward stranded
Fraction of reads explained by "1+-,1-+,2++,2--": 0.97  ← Reverse stranded ✅
```

| Result                    | Library Type         | featureCounts Setting |
|---------------------------|----------------------|-----------------------|
| ~50% / ~50%               | Unstranded           | 0                     |
| ~97% "1++,1--,2+-,2-+"    | Forward stranded     | 1                     |
| ~97% "1+-,1-+,2++,2--"    | **Reverse stranded** | **2** ✅              |

<!-- INSERT IMAGE: Infer Experiment output text result -->
![Infer Experiment](images/05_infer_experiment_output.png)
> *Figure 12: Infer Experiment output — confirms reverse-stranded library (strand 2 = ~97%).*

---

## STEP 6: Counting Reads per Gene

### What We Do
Count the number of reads overlapping each annotated gene using
**featureCounts** for all 7 samples to build a count matrix.

### How to Do It in Galaxy

**Tool: featureCounts**
1. Search `featureCounts` in tool panel
2. Parameters:

| Parameter                        | Value                        |
|----------------------------------|------------------------------|
| Alignment file(s)                | All 7 STAR BAM files         |
| Gene annotation file             | dm6 Ensembl GTF              |
| Feature type                     | exon                         |
| Gene identifier attribute        | gene_id                      |
| Strand specificity               | **Stranded (Reverse) = 2**   |
| Count multi-mapping reads        | No                           |
| Count reads also overlapping     | No (default)                 |
| Paired-end distance              | Yes (for PE samples)         |

### Expected Outputs

| Output File         | Description                                    |
|---------------------|------------------------------------------------|
| Counts table        | Matrix: genes × samples (raw integer counts)   |
| Summary file        | % assigned, unassigned, ambiguous reads        |

<!-- INSERT IMAGE: featureCounts summary bar chart from MultiQC -->
![featureCounts Summary](images/06a_featurecounts_summary.png)
> *Figure 13: featureCounts assignment summary — ~75–80% reads assigned per sample.*

<!-- INSERT IMAGE: featureCounts count matrix (first few rows) -->
![Count Matrix](images/06b_count_matrix_preview.png)
> *Figure 14: Preview of the raw count matrix — rows = genes, columns = samples.*

---

## STEP 7: Differential Gene Expression — DESeq2

### What We Do
Use **DESeq2** to normalize counts and identify statistically significant
differentially expressed genes between treated and untreated conditions.

### How to Do It in Galaxy

**Tool: DESeq2**
1. Search `DESeq2` in tool panel
2. Parameters:

| Parameter                    | Value                                      |
|------------------------------|--------------------------------------------|
| How to specify the factors   | Specify factor information                 |
| Factor name                  | condition                                  |
| Factor level 1 (reference)   | untreated → select 4 untreated count files |
| Factor level 2               | treated → select 3 treated count files     |
| Choose the format            | tabular (featureCounts output)             |
| Additional factor (batch)    | type (PE vs SE) — add as second factor     |
| Fit type                     | parametric                                 |
| Output normalized counts     | Yes                                        |
| Alpha (FDR threshold)        | 0.05                                       |

### Expected Outputs

| Output                    | Description                                        |
|---------------------------|----------------------------------------------------|
| DESeq2 result table       | log2FC, p-value, padj for every gene               |
| Normalized counts         | VST/rlog-normalized count matrix                   |
| PCA plot                  | Sample clustering by condition                     |
| Heatmap                   | Sample-to-sample distance matrix                   |
| Dispersion plot           | Model fit quality check                            |
| MA plot                   | log2FC vs mean expression                          |

<!-- INSERT IMAGE: DESeq2 dispersion plot -->
![Dispersion Plot](images/07a_deseq2_dispersion.png)
> *Figure 15: DESeq2 dispersion estimates — black dots = per-gene, red curve = fitted trend. Good fit.*

<!-- INSERT IMAGE: DESeq2 MA plot -->
![MA Plot](images/07b_deseq2_MA_plot.png)
> *Figure 16: MA plot — red dots = significantly DE genes (padj < 0.05). Symmetric distribution around 0.*

---

## STEP 8: Annotation of DESeq2 Results

### What We Do
Add gene names, descriptions, and chromosomal positions to the raw
DESeq2 output (which only contains Ensembl gene IDs).

### How to Do It in Galaxy

**Tool: Annotate DESeq2/DEXSeq output tables**
1. Search `Annotate DESeq2` in tool panel
2. Parameters:

| Parameter              | Value                                   |
|------------------------|-----------------------------------------|
| DESeq2 result file     | DESeq2 output table                     |
| Annotation database    | dm6 — *Drosophila melanogaster* Ensembl |
| Gene ID column         | Column 1 (gene_id)                      |
| Columns to add         | Gene name, description, chromosome      |

### Expected Output

The annotated table now contains:

| gene_id      | gene_name | log2FC | pvalue   | padj     | description              |
|--------------|-----------|--------|----------|----------|--------------------------|
| FBgn0025111  | ps        | -2.34  | 1.2e-15  | 3.4e-13  | pasilla, isoform A       |
| FBgn0003360  | run       | +1.87  | 4.5e-10  | 8.9e-08  | runt                     |
| ...          | ...       | ...    | ...      | ...      | ...                      |

<!-- INSERT IMAGE: Annotated DESeq2 result table (Galaxy dataset preview) -->
![Annotated DESeq2](images/08_annotated_deseq2_table.png)
> *Figure 17: Annotated DESeq2 results — gene names and descriptions added to Ensembl IDs.*

---

## STEP 9: Extraction of Differentially Expressed Genes

### What We Do
Filter the annotated DESeq2 table to extract only **significant DE genes**
using cutoffs: **padj < 0.05** AND **|log2FC| > 1**.

### How to Do It in Galaxy

**Tool: Filter data on any column using simple expressions**
1. Search `Filter` in tool panel
2. Run **twice** (or use compound filter):

**Filter 1 — Significance:**
```
c7 < 0.05
```
*(column 7 = padj)*

**Filter 2 — Fold change magnitude:**
```
abs(c3) > 1
```
*(column 3 = log2FoldChange)*

### DE Gene Summary

| Category         | Count  | Criteria                        |
|------------------|--------|---------------------------------|
| Total genes tested | ~14,000 | All expressed genes            |
| Significant DE   | ~1,200 | padj < 0.05                    |
| Upregulated      | ~600   | padj < 0.05 AND log2FC > 1     |
| Downregulated    | ~600   | padj < 0.05 AND log2FC < −1    |

<!-- INSERT IMAGE: Filtered DE gene list in Galaxy history -->
![DE Gene List](images/09_filtered_DE_genes.png)
> *Figure 18: Filtered DE gene table — 1,200 significant genes after applying padj and fold-change cutoffs.*

<!-- INSERT IMAGE: Volcano plot of DE genes -->
![Volcano Plot](images/09b_volcano_plot.png)
> *Figure 19: Volcano plot — x-axis = log2FC, y-axis = -log10(padj). Red = significant DE genes.*

---

## STEP 10: Visualization of DE Gene Expression

### What We Do
Generate a **PCA plot** and **sample-to-sample heatmap** from normalized
counts to visualize variability and confirm biological groupings.

---

### 10A. PCA Plot

The PCA plot is **automatically generated by DESeq2** using
variance-stabilized (VST) transformed counts.

| PC   | Variance Explained | Biological Meaning              |
|------|--------------------|---------------------------------|
| PC1  | ~55–65%            | Separates treated vs. untreated |
| PC2  | ~10–20%            | Separates PE vs. SE libraries   |

<!-- INSERT IMAGE: PCA plot from DESeq2 output -->
![PCA Plot](images/10a_PCA_plot.png)
> *Figure 20: PCA plot — PC1 separates conditions (treated=red, untreated=blue). PC2 separates library types. Expected biological clustering confirmed.*

---

### 10B. Sample-to-Sample Distance Heatmap

Also **auto-generated by DESeq2** using Euclidean distances on VST counts.

| Color (dark blue) | Meaning                    |
|-------------------|----------------------------|
| Dark blue         | Low distance = very similar|
| Light/white       | High distance = dissimilar |

<!-- INSERT IMAGE: Sample-to-sample heatmap from DESeq2 -->
![Sample Heatmap](images/10b_heatmap_sample_to_sample.png)
> *Figure 21: Sample-to-sample heatmap — treated samples cluster together (bottom-right block), untreated cluster separately (top-left block).*

---

### 10C. Heatmap of Top DE Genes

**Tool: heatmap2** (in Galaxy)
1. Search `heatmap2` in tool panel
2. Parameters:

| Parameter             | Value                          |
|-----------------------|--------------------------------|
| Input                 | Normalized counts (top 50 DE genes) |
| Clustering            | Both rows and columns          |
| Scale                 | By row (z-score)               |
| Color scheme          | Blue-White-Red                 |

<!-- INSERT IMAGE: heatmap2 of top 50 DE genes -->
![DE Gene Heatmap](images/10c_heatmap_top_DE_genes.png)
> *Figure 22: Heatmap of top 50 DE genes — rows = genes, columns = samples. Red = high expression, blue = low.*

---

## STEP 11: Gene Ontology (GO) Analysis

### What We Do
Identify which **biological processes, molecular functions, and cellular
components** are enriched among the DE genes using **goseq**.

### How to Do It in Galaxy

**Tool: goseq**
1. Search `goseq` in tool panel
2. Parameters:

| Parameter                    | Value                                  |
|------------------------------|----------------------------------------|
| Differentially expressed genes| Filtered DE gene list (gene IDs only) |
| Gene lengths file            | From featureCounts output              |
| Genome                       | dm6                                    |
| Gene categories              | GO: Biological Process (BP)            |
| Method to correct for bias   | Wallenius (default)                    |
| p-value cutoff               | 0.05                                   |

### Expected Output — Top GO Terms

| GO Term ID   | Description                        | p-adj   | DE Genes |
|--------------|------------------------------------|---------|----------|
| GO:0008380   | RNA splicing                       | 1.2e-08 | 45       |
| GO:0003729   | mRNA binding                       | 3.4e-07 | 38       |
| GO:0006397   | mRNA processing                    | 5.6e-06 | 52       |
| GO:0005681   | Spliceosomal complex               | 8.9e-06 | 29       |

<!-- INSERT IMAGE: goseq dot plot or bar chart of top GO terms -->
![GO Analysis](images/11a_GO_enrichment_plot.png)
> *Figure 23: Top enriched GO Biological Process terms. RNA splicing and mRNA processing dominate — consistent with Pasilla's role as a splicing regulator.*

<!-- INSERT IMAGE: goseq result table in Galaxy -->
![GO Table](images/11b_GO_result_table.png)
> *Figure 24: goseq output table with GO term IDs, descriptions, p-values, and DE gene counts.*

---

## STEP 12: KEGG Pathway Analysis

### 🎯 What We Do
Map DE genes onto **KEGG metabolic and signaling pathways** to identify
which biological pathways are most impacted by Pasilla depletion.

### 🔧 How to Do It in Galaxy

**Tool: KEGG Pathway Enrichment (via goseq or pathview)**
1. Search `goseq` in tool panel → set **Gene categories** to `KEGG`
2. Or use **pathview** for pathway diagram visualization

| Parameter              | Value                         |
|------------------------|-------------------------------|
| Input gene list        | Significant DE genes          |
| Organism               | dme (*D. melanogaster*)       |
| Database               | KEGG                          |
| p-value cutoff         | 0.05                          |

### Expected Output — Top KEGG Pathways

| KEGG ID   | Pathway Name                     | p-adj   | DE Genes |
|-----------|----------------------------------|---------|----------|
| dme03040  | Spliceosome                      | 2.1e-09 | 38       |
| dme03013  | RNA transport                    | 4.5e-06 | 29       |
| dme03015  | mRNA surveillance pathway        | 7.8e-05 | 21       |
| dme04120  | Ubiquitin mediated proteolysis   | 0.003   | 18       |

<!-- INSERT IMAGE: KEGG enrichment bar chart -->
![KEGG Enrichment](images/12a_KEGG_enrichment_barplot.png)
> *Figure 25: KEGG pathway enrichment — spliceosome pathway most significantly enriched, confirming Pasilla's role in splicing regulation.*

<!-- INSERT IMAGE: KEGG pathway diagram (pathview) for spliceosome -->
![KEGG Pathway Diagram](images/12b_KEGG_spliceosome_pathview.png)
> *Figure 26: Pathview diagram of the spliceosome pathway — DE genes highlighted in red (upregulated) and blue (downregulated).*

---

## 📁 Repository Structure

```
RNA-Seq-Galaxy-Report/
├── README.md                        ← This complete report
├── images/
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
├── data/
│   └── counts_matrix_all7samples.tabular
└── results/
    ├── DESeq2_results_annotated.tabular
    ├── DE_genes_significant.tabular
    ├── GO_enrichment_results.tabular
    └── KEGG_enrichment_results.tabular
```

---

## Complete Pipeline Overview

```
FASTQ files (Zenodo)
       ↓
  [Falco + MultiQC] → QC Report
       ↓
   [RNA STAR] → BAM + BigWig (Str1=Blue, Str2=Red) + ReadsPerGene
       ↓
  [IGV / Sashimi / JBrowse2] → Alignment Visualization
       ↓
  [Infer Experiment] → Strandedness confirmed (reverse)
       ↓
  [featureCounts] → Count matrix (7 samples × ~14,000 genes)
       ↓
    [DESeq2] → DE genes + PCA + Heatmap + MA plot
       ↓
  [Annotate] → Gene names + descriptions added
       ↓
   [Filter] → Significant DE genes (padj<0.05, |log2FC|>1)
       ↓
  [heatmap2] → Top DE gene expression heatmap
       ↓
   [goseq] → GO enrichment (BP, MF, CC)
       ↓
[goseq/pathview] → KEGG pathway enrichment + diagrams
```

---

## References

1. Galaxy Training Network (2026) — Reference-based RNA-Seq tutorial:
   https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html

2. Batut B, van den Beek M, Doyle MA, Soranzo N (2021) — RNA-Seq Data
   Analysis in Galaxy. *Methods in Molecular Biology*, 2284:367–392.
   https://hal.science/hal-04645142v1/document

3. Schwartz AV (2023) — Reference Based RNA-Seq Analysis in Galaxy:
   https://ashleyschwartz.com/posts/2023/05/galaxy-tutorial

4. Brooks AN et al. (2011) — Conservation of an RNA regulatory map
   between Drosophila and mammals. *Genome Research*, 21(2):193–202.

# BRAID
Block Resolution and Annotation of integrated DNA

![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-active-success)

![alt text](path/to/image.png)

> **Analyze phased VCF data to predict variant effects on protein sequences for each haplotype.**

## Overview

**BRAID** is a bioinformatics tool designed to go beyond isolated mutations annotation. By utilizing **phased VCF data**, this tool reconstructs the combination of mutations present on each chromosome (haplotype) to predict the actual protein sequence produced.

This allows for the detection of complex effects, such as:
* **Compound Heterozygosity:** Understanding how multiple mutations on the same allele interact.
* **Haplotype-specific LOF:** Determining if a combination of missense mutations leads to a Loss of Function.
* **Protein Structure Changes:** Visualizing the exact amino acid sequence changes (e.g., `Sub(6)S>W`).

![Workflow Diagram](https://via.placeholder.com/800x300?text=Place+Your+Workflow+Diagram+Here)

---

## Installation & Requirements

This tool is a software that requires standard bioinformatics libraries.

**Prerequisites:**
> * Python 3.x
> * `pysam`
> * `biopython`

**installation**

Dowload by pip
> pip install braid

Dowload by conda
> conda install braid

Dowload by wget
> wget

Dowload by git
> git

Testing whether BRAID has been installed successfully：
braid test

## Usage
Run the script from the command line by providing the GFF3 annotation, Reference Genome, and Phased VCF.

`braid -r reference.fa -g annotation.gff3 -v phased_variants.vcf.gz`

### Example for test dataset
`braid -r test.fasta -g test.gff3 -v test.vcf.gz`
> ```
> You should have three output files and one log file:
> variant_analysis_output.tsv
> variant_analysis_output.alignment.txt
> variant_analysis_output.sample.txt
> variant_analysis_output.log
> ```

## Arguments Parameter Table
| **Short** | **Long** | **Description** |       |
|:-----:|:------:|:-------------:|:-------------:|
| **-g** | **--gff** | Path to the GFF3 annotation file. | Required |
| **-r** | **--reference** | Path to the Reference Genome FASTA (must be indexed: .fai). | Required |
| **-v** | **--vcf** | Path to the Phased VCF file (must be indexed: .tbi/.csi). | Required |
| **-o** | **--output** | Output file name (default: variant_analysis_output.tsv). | Optional |
| **-s** | **--sample** | Path to file containing specific sample IDs to analyze (one per line, no header). | Optional |
|  | **--gene** | Path to file with specific gene IDs to analyze (one per line, no header). | Optional |
|  | **--lof-threshold** | Custom threshold for Loss-of-Function classification (e.g., 0.1; default: 0.3). | Optional |
|  | **--force-unphased** | Skip phasing check (forces run on unphased VCF). | Optional |
|  | **--ignore-intron** | Ignore variants marked strictly as intronic. | Optional |

## Output Files Explaination

BRAID generates three main files to assist in your analysis.

**1. Summary Table**

A comprehensive table detailing the protein changes for every haplotype (default: variant_analysis_output.tsv).

| Gene_ID | Haplotype_ID | mRNA | Haplotype_Count | Frequency | Variant_Type | Protein_Changes | Haplotype_Mutations | Sample_Sources | Ref_Protein | Alt_Protein | Ref_CDS | Alt_CDS | Aligned_Ref | Comparison_String | Aligned_Alt |
|:-----|:------|:------|:-----|:------|:------|:-----|:------|:------|:-----|:------|:------|:-----|:------|:------|:-----|
| gene1 | transcript1:REF | transcript1 | . | . | NoLOF(non_identity_rate:0.00%,non_identical_AAs:0,total_ref_AAs:27) `\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|` | . | . | . | MSLASSANDMIDRSIDRSIDRSIDRS* | MSLASSANDMIDRSIDRSIDRSIDRS* | ATGAGCTTAGCTAGCTCAGCTAACGATATGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTGA | ATGAGCTTAGCTAGCTCAGCTAACGATATGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTGA | MSLASSANDMIDRSIDRSIDRSIDRS* | `\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|` | MSLASSANDMIDRSIDRSIDRSIDRS* |
| gene1 | transcript1:1 | transcript1 | 4 | 0.500000 | NoLOF(non_identity_rate:3.70%,non_identical_AAs:1,total_ref_AAs:27)`\|\|\|`deletion `\|\|\|\|\|\|\|\|\|\|\|\|\|\|` | Del(5)S | 1:17_CTTAG>C[CDS,EXON];1:27_T>TT[CDS,EXON] | sample1(Hap2);sample2(Hap1);sample3(Homo) | MSLASSANDMIDRSIDRSIDRSIDRS* | MSLASANDMIDRSIDRSIDRSIDRS* | ATGAGCTTAGCTAGCTCAGCTAACGATATGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTGA | ATGAGCCTAGCTTCAGCTAACGATATGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTGA | MSLASSANDMIDRSIDRSIDRSIDRS* | `\|\|\|\| \|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|\|` | MSLA-SANDMIDRSIDRSIDRSIDRS* |

| **Column** | **Example** | **Explanation** |
|:-----|:-----|:-------------|
| **Gene_ID** | gene1 | The gene to which this haplotype belongs. |
| **Haplotype_ID** | transcript1:REF;transcript1:1 | transcript:REF indicates the reference haplotype; others (e.g. transcript1:1) are alternative haplotypes. |
| **mRNA** | transcript1 | The transcript to which this haplotype belongs. |
| **Haplotype_Count** | 4 | Number of samples carrying this haplotype. |
| **Frequency** | 0.500000 | Population frequency of this haplotype. |
| **Variant_Type** | LOF/NoLoF (Please see the detailed information below.) | Functional classification of the haplotype	and indicates LOF status. |
| **Protein_Changes** | Del(5)S | Protein-level consequence	HGVS format description. |
| **Haplotype_Mutations** | 1:17_CTTAG>C[CDS,EXON];1:27_T>TT[CDS,EXON] | List of variants defining the haplotype. |
| **Sample_Sources** | sample1(Hap2);sample2(Hap1);sample3(Homo) | Samples carrying this haplotype, including haplotype phase (Hap1, Hap2) or homozygous status (Homo). |
| **Ref_Protein** | MSLASSANDIDRSIDRS* | Reference protein sequence. |
| **Alt_Protein** | MSLASANDIDRSIDRS* | Alternative protein sequence. |
| **Ref_CDS** | ATGAGCTTAGCTAGCTCAGCTAACGATATCG<br>ATCGATCGATCGATCGATCGTGA | Reference CDS sequence. |
| **Alt_CDS** | ATGAGCCTAGCTTCAGCTAACGAT<br>ATCGATCGATCGATCGATCGATCGTGA | Haplotype-specific CDS sequence. |
| **Aligned_Ref** | MSLASSANDIDRSIDRS* | Aligned reference protein, protein alignment string for visualization: `-` → Gap (insertion or deletion). |
| **Comparison_String** | MSLASANDIDRSIDRS* | Alignment comparison symbols: `*` → Different amino acid, `\|` → Same amino acid. |
| **Aligned_Alt** | MSLASANDIDRSIDRS* | Aligned alternative protein, shows amino acid changes: `-` → Gap (insertion or deletion). |

**Detail information for Variant_Type**

| **Column** | **Explanation** |
|:-----|:-------------|
| **LOF_Info** | LOF indicates haplotypes predicted to cause protein loss of function, whereas NoLOF indicates haplotypes without LOF effects. non_identity_rate represents the proportion of amino acid differences between ALT and REF proteins; non_identical_AAs is the corresponding count, and total_ref_AAs denotes the length of the REF protein. e.g., NoLOF(non_identity_rate:5.56%,non_identical_AAs:1,total_ref_AAs:18). |
| **missense** | A nucleotide variant that results in the substitution of one amino acid by another in the protein sequence. |
| **insertion** | An insertion of one or more nucleotides that alters the coding sequence and may affect the resulting protein sequence. |
| **deletion** | A deletion of one or more nucleotides from the coding sequence, potentially altering the protein sequence or reading frame. |
| **complex_indel** | A combined insertion and deletion event that cannot be represented as a simple insertion or deletion and may cause complex changes to the coding sequence. |
| **complex_indel** | A splicing event in which one or more exons are completely skipped in the transcript, leading to an altered mRNA and protein sequence. |
| **exon_skip** | A splicing event in which one or more exons are completely skipped in the transcript, leading to an altered mRNA and protein sequence. |
| **skipped_exons_detail** | Detailed information specifying which exon is skipped in the exon-skipping event, SITE_PRESERVED, SITE_SHIFT, SITE_DESTROYED. e.g., SkippedExon:[1:51-77|SITE_SHIFT]. |
| **intron_retention** | A splicing event in which one or more introns are retained in the mature transcript, potentially disrupting the coding sequence. |
| **retained_introns_detail** | Detailed information specifying which intron(s) are retained in the intron-retention event, SITE_PRESERVED, SITE_SHIFT, SITE_DESTROYED. e.g., RetainedIntron:[1:78-90|SITE_DESTROYED]. |
| **start_codon_loss** | A variant that disrupts the canonical start codon, potentially preventing translation initiation. |
| **stop_loss** | A variant that removes or alters the stop codon, resulting in translational read-through and an extended protein. |
| **No_start_codon_for_reference** | The reference transcript or protein sequence lacks an annotated start codon. |
| **same_as_other_transcript** | The variant effect on this protein is identical to that observed in another transcript's protein of the same gene. |
| **protein_loss** | The variant or variant combination results in the complete loss of the predicted protein product. |
| **ref_protein_empty** | The reference transcript does not produce a protein sequence (e.g., non-coding or incomplete annotation). |
| **alignment_failed** | The reference and alternative protein sequences could not be reliably aligned, preventing accurate variant effect annotation. |

**Log information for splice sites mutations**

**SITE_PRESERVED**: The splice site is unaltered, with both its position and sequence conserved; normal splicing is expected.
> ```
> [transcript1] (strand +) Splice site 'GT' at 39 PRESERVED by mutation 1:44_TAAGTA>A.
>   - Original Window Seq  : GATGTCGTTAAG     - Mutated Window Seq  : GATGTCGTA.
> ```

**SITE_SHIFT**: The splice site may move to a nearby position due to the variant, potentially altering the exon–intron boundary.
> ```
> WARNING - [transcript1] Splice site may SHIFT for mutation 1:44_TAAGTA>A.
>   - Splice Site          : 'AG' (Original genomic pos: 49)
>   - Mutation             : TAAGTA -> A (Genomic pos: 44)
>   - Window               : 1:42-53
>   - Original Window Seq  : GTTAAGTAGATG
>   - Mutated Window Seq   : GTAGATG
> ```

**SITE_DESTROYED**: The splice site is disrupted or abolished by the variant, likely preventing normal splicing at this site.
> ```
>  [transcript1] (strand +) Splice site 'GT' at 78 DESTROYED by mutation 1:50_GATGATCGATCGATCGATCGATCGATCGG>GG. It was 'GT', became 'AT'.
>    - Original Window Seq  : TCGGTCGATCGA     - Mutated Window Seq  : TCGATCGA
> ```

**2. Alignment Visualization (.alignment.txt)**

A text file showing the pairwise alignment of the Reference vs. Alternative protein.
Use the Haplotype_ID to match the haplotype in the summary table.

Example View:
> ```
> Haplotype_ID: 1
> Gene: gene1 | mRNA: transcript1
> Variant_Type: NoLOF | missense
> Alignment:
>   Ref: MSLASSANDMIDRSIDRSIDRSIDRS*
>        |||||*|||||||||||||||||||||
>   Alt: MSLASWANDMIDRSIDRSIDRSIDRS*
> ```

**3. Sample Matrix (.sample.txt)**

A matrix format ideal for heatmaps or downstream programmatic analysis.

| Gene_ID | mRNA_ID | Ref_ID | Alt_IDs | sample1 | sample2 | sample3 |
|:-----|:------|:------|:-----|:------|:------|:------|
| gene1	| transcript1 |	transcript1:REF |	transcript1:1	| 0|1 |	1|0 | 1|1

| **Column** | **Explanation** |
|:-----|:-------------|
| **Gene_ID** | The identifier of the gene being analyzed. |
| **mRNA_ID** | The identifier of the specific transcript of that gene. |
| **Ref_ID** | Reference haplotype ID, labeled as transcript:REF. |
| **Alt_IDs** | Alternative haplotype IDs observed in the samples, with numbering to distinguish multiple haplotypes. |
| **samples** | The haplotype presence/absence in each sample. Each entry corresponds to the Ref or Alt haplotype. `0` represents REF. |

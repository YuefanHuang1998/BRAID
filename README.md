# BRAID
Block Resolution and Annotation of Integrated DNA

![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-active-success)

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
> * `pysam` (for VCF/BAM handling)
> * `biopython` (for sequence manipulation)
> * `pandas` (for data tables)

## installation 

### method 1
`pip install braid`

### Method 2
`conda install braid`

### Method 3
`Wget`

### Method 4
`git `

### To test installation:
braid test

## Usage
Run the script from the command line by providing the GFF3 annotation, Reference Genome, and Phased VCF.

> braid -g annotation.gff3 -r reference.fa -v phased_variants.vcf.gz 

### Example for test dataset
`braid -r test.fasta -g test.gff3 -v test.vcf.gz`
> You should have three output files: 
> variant_analysis_output.tsv
> variant_analysis_output.alignment.txt
> variant_analysis_output.sample.txt
> variant_analysis_output.log

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
The tool generates three main files to assist in your analysis.

1. Summary Table

A comprehensive table detailing the protein changes for every haplotype (default: variant_analysis_output.tsv).

| Gene_ID | Haplotype_ID | mRNA | Haplotype_Count | Frequency | Variant_Type | Protein_Changes | Haplotype_Mutations | Sample_Sources | Ref_Protein | Alt_Protein | Ref_CDS | Alt_CDS | Aligned_Ref | Comparison_String | Aligned_Alt |
|:-----:|:------:|:------:|:-----:|:------:|:------:|:-----:|:------:|:------:|:-----:|:------:|:------:|:-----:|:------:|:------:|:-----:|
| gene1 | transcript1:REF | transcript1 | . | . | NoLOF(non_identity_rate:0.00%,non_identical_AAs:0,total_ref_AAs:27)`|||||||||||||||||` | . | . | . | MSLASSANDMIDRSIDRSIDRSIDRS* | MSLASSANDMIDRSIDRSIDRSIDRS* | ATGAGCTTAGCTAGCTCAGCTAACGATATGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTGA | ATGAGCTTAGCTAGCTCAGCTAACGATATGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTGA | MSLASSANDMIDRSIDRSIDRSIDRS* | `|||||||||||||||||||||||||||` | MSLASSANDMIDRSIDRSIDRSIDRS* |
| gene1 | transcript1:1 | transcript1 | 4 | 0.500000 | NoLOF(non_identity_rate:3.70%,non_identical_AAs:1,total_ref_AAs:27)|||deletion`||||||||||||||` | Del(5)S | 1:17_CTTAG>C[CDS,EXON];1:27_T>TT[CDS,EXON] | sample1(Hap2);sample2(Hap1);sample3(Homo) | MSLASSANDMIDRSIDRSIDRSIDRS* | MSLASANDMIDRSIDRSIDRSIDRS* | ATGAGCTTAGCTAGCTCAGCTAACGATATGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTGA | ATGAGCCTAGCTTCAGCTAACGATATGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTGA | MSLASSANDMIDRSIDRSIDRSIDRS* | `|||| ||||||||||||||||||||||` | MSLA-SANDMIDRSIDRSIDRSIDRS* |

| **Column** | **Example** | **Explanation** |
|:-----:|:------:|:-------------:|
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
| **Ref_CDS** | ATGAGCTTAGCTAGCTCAGCTAACGATATCGATCGATCGATCGATCGATCGTGA | Reference CDS sequence. |
| **Alt_CDS** | ATGAGCCTAGCTTCAGCTAACGATATCGATCGATCGATCGATCGATCGTGA | Haplotype-specific CDS sequence. |
| **Aligned_Ref** | MSLASSANDIDRSIDRS* | Aligned reference protein, protein alignment string for visualization: `-` → Gap (insertion or deletion). |
| **Comparison_String** | MSLASANDIDRSIDRS* | Alignment comparison symbols: `*` → Different amino acid, `|` → Same amino acid. |
| **Aligned_Alt** | MSLASANDIDRSIDRS* | Aligned alternative protein, shows amino acid changes: `-` → Gap (insertion or deletion). |

Detail information for Variant_Type 

| **Column** | **Explanation** |
|:-----:|:-------------:|
| **LOF_Info** | LOF indicates haplotypes predicted to cause protein loss of function, whereas NoLOF indicates haplotypes without LOF effects. non_identity_rate represents the proportion of amino acid differences between ALT and REF proteins; non_identical_AAs is the corresponding count, and total_ref_AAs denotes the length of the REF protein. e.g., NoLOF(non_identity_rate:5.56%,non_identical_AAs:1,total_ref_AAs:18). |
| **missense** | A nucleotide variant that results in the substitution of one amino acid by another in the protein sequence. |
| **insertion** | An insertion of one or more nucleotides that alters the coding sequence and may affect the resulting protein sequence. |
| **deletion** | A deletion of one or more nucleotides from the coding sequence, potentially altering the protein sequence or reading frame. |
| **complex_indel** | A combined insertion and deletion event that cannot be represented as a simple insertion or deletion and may cause complex changes to the coding sequence. |

        "missense",
        "insertion",
        "deletion",
        "complex_indel",
        "exon_skip",
        "skipped_exons_detail",
        "intron_retention",
        "retained_introns_detail",
        "try_to_exon_skip",
        "try_to_intron_retention",
        "start_codon_loss",
        "stop_loss",
        "No_start_codon_for_reference",
        "same_as_other_transcript",
        "protein_loss",
        "ref_protein_empty",
        "alignment_failed"


2. Alignment Visualization (.alignment.txt)
A human-readable text file showing the pairwise alignment of the Reference vs. Alternative protein.

Example View:

Plaintext
Haplotype_Analysis_ID: 1
Gene: gene1 | mRNA: transcript1
Variant_Type: NoLOF | missense
Alignment:
  Ref: MSLASSANDMIDRSIDRSIDRSIDRS*
       |||||*|||||||||||||||||||||
  Alt: MSLASWANDMIDRSIDRSIDRSIDRS*
3. Sample Matrix (.sample.txt)
A matrix format ideal for heatmaps or downstream programmatic analysis.

Gene_ID	mRNA_ID	Ref_ID	Alt_IDs	sample1	sample2	sample3
gene1	transcript1	transcript1:REF	transcript1:1	0|1	1|0	1|1

# nf-ampgenomecov

**nf-ampgenomecov** is a modular and scalable [Nextflow](https://www.nextflow.io/) pipeline designed to generate **genome and amplicon-level coverage profiles and visualizations** from BAM alignment file inputs.  It is particularly suited for **tiled amplicon sequencing data** — such as viral panels or multiplexed PCR assays — and enables high-resolution coverage diagnostics across target regions.

---

## ✨ Key Features

- Accepts a directory contains **aligned BAM files** as input
- Computes **per-genome** and **per-amplicon region** coverage and depth
- Supports **tiled amplicon schemes**
- Customizable reference metadata and amplicon definitions
- Schema-driven configuration with validation
- Compatible with Docker, Singularity, or Conda environments

---

## 🧬 Typical Use Case

This pipeline is useful in:
- **Pathogen genomics** (e.g., SARS-CoV-2, MPXV, viral panels)
- **Amplicon-based sequencing** workflows
- Public health genomic surveillance
- Microbial genome completeness and coverage evaluation

---

## 🧾 Input Requirements

### 1. **BAM files**  
Pre-aligned BAM files with reads mapped to a reference genome.  
All BAMs should be coordinate-sorted and indexed.

```bash
/path/to/sample1.primertrimmed.rg.sorted.bam
/path/to/sample2.primertrimmed.rg.sorted.bam
...
```

example command

```bash
nextflow run ../main.nf --bam_dir ./bam --bam_file_pattern *.primertrimmed*.bam -profile singularity --outdir results/ --input bam --ref reference/reference.fasta --primer_bed bed/amplicon.bed 
```


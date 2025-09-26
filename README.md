


# nf-ampgenomecov

**nf-ampgenomecov** bioinformatics pipeline, takes aligned andsorted bam files, primer bed file, fasta format reference file as input to generate **genome and amplicon-level depth and coverage profiles and visualizations**.  It is particularly suited for **tiled amplicon sequencing data** — such as viral panels or multiplexed PCR assays — and enables high-resolution coverage diagnostics across target regions.

## Pipeline summary
The pipeline uses [mosdepth](https://github.com/brentp/mosdepth) to calculate per-base and windowed sequencing depth, then summarizes and plots these coverage profiles 

## 🧬 Typical Use Case
This pipeline is useful in:
- **Pathogen genomics** (e.g., SARS-CoV-2, MPXV, viral panels)
- **Amplicon-based sequencing** workflows
- Public health genomic surveillance
- Microbial genome completeness and coverage evaluation

## Quick start

>If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with -profile test before running the workflow on actual data.

### Check workflow options
You can clone or download the nf-ampgenomecov from github to local computer or you can directly run the pipeline from github. To check the pipeline command line options:

```{r df-drop-ok, class.source="bg-success"}
# running directly from github without downloading or cloning
nextflow run xiaoli-dong/nf-ampgenomecov -r revision_number(e.g:e34cdfb) --help
```

### required inputs
1. bam directory" a directory contains the sorted bam files with or without index files
2. bam file naming pattern
3. bed files
4. primer suffix
5. reference fasta file
### Run the pipeline:
Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash

nextflow run ../main.nf \
  --bam_dir ./bam \
  --bam_file_pattern *.primertrimmed*.bam \
  -profile singularity \
  --outdir results \
  --input bam \
  --ref reference/reference.fasta \
  --primer_bed bed/amplicon.bed \
  --primer_left_suffix _LEFT \
  --primer_right_suffix _RIGHT \
  -resume 

```  
>* Notes: Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits
xiaoli-dong/nf-ampgenomecov was originally written by Xiaoli Dong.


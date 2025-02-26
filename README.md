# General Description

This script computes the Distribution of Spatial Variation in SNP Frequencies (DSVSF) and the expected of heterozygosity either from an output of PhyloGeoSim (v2.0.beta_210319) (`--input_mode pgs`) or from a VCF file (`--input_mode vcf`).

# Authors

Chedly Kastally & Patrick Mardulyn

# Copyright

Copyright © 2024 Chedly Kastally & Patrick Mardulyn.

All rights reserved.

# Usage

## Usage: `pgs` input mode

Example:

```.bash

Rscript ./DSVSF_computation.R \
    --input_file <pgs_output_file> \
    --input_mode pgs \
    --threshold <threshold_value> \
    --output_prefix <output_prefix>

```

In `pgs` input mode, four output files are produced:

- \*dsvsf_output_table.tsv: a tab separated table with the DSVSF patterns and their count in the input file
- \*global_he_output.tsv:   a tab separated table with the loci ID and their He (at the entire sample level)
- \*global_Fst_output.txt:  a text file with the value of the overall Fst (Weir & Cockerham 1984, Evolution)
- \*all_pwFst_output.tsv:   a table with pairwise Fst for each population pair (Weir & Cockerham 1984, Evolution)

## Usage: `vcf` input mode

Example:

```.bash

Rscript ./DSVSF_computation.R \
    --input_file <vcf_input> \
    --input_mode vcf \
    --population_file <population_file> \
    --posterior_n_file <posterior_n_file> \
    --vcf_missing_data proj \
    --threshold <threshold_value> \
    --output_prefix <output_prefix>

```

In `vcf` input mode, one output is produced: 

- \*dsvsf_output_table.tsv: a tab separated table with the DSVSF patterns and their count in the input file

In `vcf` input mode, two additional input files are required:

- a population file, which indicates how the samples are grouped together. It needs to be a tab-separated values file (`.tsv`) with two columns: `sample_id` and `population_id`
- posterior sample size file, which indicates the projected sample size for each population. It needs to be a `.tsv` file with two columns: `population_id` and `posterior_n`

### Handling missing data

Further, a method to handle missing data in the VCF has to be provided with `--vcf_missing_data`:

- `randomSample`: at each locus, and for each population, a random set of samples are picked to compute the allele count (the size of the set is defined in the `posterior_n_file`. This is faster to compute, but add randomness in the computation of patterns.
- `proj`: at each locus, and for each population, a projection using the hypergeometric distribution is used to determine the patterns observed. This should recover the true set of patterns, but is cpu-intensive and requires a small simplification during the projection: patterns with probability $< 0.01$ (at the locus level) are excluded from the computation.

# Citation

Kastally, C., Dellicour, S., Hardy, O. J., Gilbert, M., & Mardulyn, P. (2021).
_Estimating Migration of Gonioctena quinquepunctata (Coleoptera: Chrysomelidae) Inside a Mountain Range in a Spatially Explicit Context._
Insect Systematics and Diversity, 5(5). https://doi.org/10.1093/isd/ixab019

# Usage

Usage: ./DSVSF_computation.R -i input_file [-t 0.1] -m pgs [-o output_prefix]

Options:

	-i INPUT_FILE, --input_file=INPUT_FILE
		Input file [Default: example/out_simulation_SNPs_1.txt]

	-m INPUT_MODE, --input_mode=INPUT_MODE
		Input mode: pgs or vcf [Default: pgs]

	-n POSTERIOR_N_FILE, --posterior_n_file=POSTERIOR_N_FILE
		File with posterior sample sizes, required if using a vcf [Default: example/posterior_n_per_pop_file_2.tsv]

	-p POPULATION_FILE, --population_file=POPULATION_FILE
		population file, required if using a vcf [Default: example/pop_file_2.tsv]

	-s VCF_MISSING_DATA, --vcf_missing_data=VCF_MISSING_DATA
		Method to handle missing data: randomSample or proj [Default: randomSample]

	-t THRESHOLD, --threshold=THRESHOLD
		Threshold used to compute the DSVSF [Default: 0.1]

	-o OUTPUT_PREFIX, --output_prefix=OUTPUT_PREFIX
		Output prefix [Default: Output_]

	-h, --help
		Show this help message and exit

# Notes about using VCF as input

Make sure of the following:

## VCF contains only bi-allelic loci.

To make sure, you can use `bcftools` to keep only bi-allelic loci

```{.bash}

bcftools view -m2 -M2 input.vcf > output.vcf

```

Note that currently the script checks that only the following genotypes are present: "0/0"; "0/1"; or "1/1".

This is more restrictive than processing only bi-allelic loci since a bi-allelic locus with genotypes: 1/1; 1/2 and 2/2 would be bi-allelic in the data but excluded by this analysis. These genotypes may be more common if the reference genome used is of a different species than the samples from which sequencing reads were obtained.

Reach out if this is an issue.

## Posterior sample size and ploidy

The minimum sample size at a locus to make a projection of the SFS at a (haploid) size $n_t$ is $n_t$.

When using a VCF, the script expect diploid samples, but the posterior sample sizes must be indicated in a number of (haploid) gene copies.

Loci with observed size (after removing missing data) below $n_t$ for any given population are dropped during the projection (a warning message is issued for those loci).

If the total number of samples is below the posterior sample size, the analysis stops early since no loci in the data could ever have a sample size big enough for a projection.

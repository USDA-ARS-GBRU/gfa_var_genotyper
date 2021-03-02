# gfa_var_genotyper:

Usage:

    gfa_var_genotyper.pl [options]

Options:

  general options:

     -v --var     gfa variants file (required)

     -p --pack    vg pack table file(s)

     --packlist   text file containing list of pack file paths
                    1 file per line

     --ploidy     1 (haploid) or 2 (diploid) currently supported
                    default: 1

  genotyping options:

    Genotyping thresholds can be set manually or dynamically by creating
    coverage count histograms that are provide to GenomeScope for modeling.
    (modeling is not recommended for GBS data) When modeling is enabled, any
    samples that fail to generate a valid model will fallback to using the
    static parameters.

    Static parameter genotyping sets a minimum total coverage (all alleles)
    threshold. The remaining variants are processed as either low or high
    total coverage. Low coverage variants are genotyped using minimum allele
    counts. High coverage variants are genotyped using minimum allele
    percentages.

   modeling:

     -m --model   use genomescope model to define genotyping thresholds
                    default: disabled

     --modeldir   output directory for modeling files
                    default: gfa_var_genotyper_models

     --gs         path to GenomeScope executable

   static parameters:

     --min_tot_cov               minimum total coverage (all alleles)
                                   default: 3

     --max_low_cov_tot_cov       maximum "low coverage" total coverage
                                   default: 9

     --min_low_cov_allele_count  minimum "low coverage" allele count
                                   default: 3

     --min_high_cov_allele_pct   minimum "high coverage" allele percent
                                   default: 10

  help:

     -h --help    display help menu

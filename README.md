# gfa_var_genotyper:

Options:

  general options:

     -v --var       gfa variants file (required)

     -p --pack      vg pack table or segment coverage file(s)
                      ex: -p sample.1.pack.table sample.2.pack.table

     --packlist     text file containing list of pack file paths
                      1 file per line

    1 or more pack (table or segment coverage) files may be specified using
    -p/--pack and/or --packlist. Pack table files are generated using the
    `vg pack -d` command. Pack segment coverage files are generated using
    pack_table_to_seg_cov.pl, which is part of the gfa_var_genotyper
    project. (https://github.com/brianabernathy/gfa_var_genotyper) Pack
    files may be uncompressed or compressed using either gzip or bzip2. (.gz
    or .bz2 file extension)

     --ploidy       1 (haploid) or 2 (diploid) currently supported
                      default: 1

     --rm_inv_head  remove variants with inverted head node
                      default: disabled

    Inversion variants can result in a negative sign suffix in the head node
    (POS) field. Conventionally, this field represents the position in the
    reference genome and negative values may cause issues with tools that
    use vcfs. To remove such variants from vcf ouput, use --rm_inv_head.
    Note, the reciprocal variant of the inversion should be called
    regardless, so the variant information is still retained for most
    practical purposes.

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
                    default: autodetect in $PATH (if available)

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

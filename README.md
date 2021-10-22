# gfa_var_genotyper

### gfa_var_genotyper.pl

Usage:

    gfa_var_genotyper.pl [options] > genotypes

Options:

  general options:

     -v --var       gfa variants file (required)

     -p --pack      vg pack table or segment coverage file(s)
                      ex: -p sample.1.pack.table sample.2.pack.table

     --packlist     text file containing list of pack or segment coverage
                      file paths (1 file per line)

    1 or more pack (table or segment coverage) files may be specified using
    -p/--pack and/or --packlist. Pack table files are generated using the
    `vg pack -d` command. Pack segment coverage files are generated using
    pack_table_to_seg_cov.pl, which is part of the gfa_var_genotyper
    project. (https://github.com/brianabernathy/gfa_var_genotyper) Pack
    files may be uncompressed or compressed using either gzip or bzip2. (.gz
    or .bz2 file extension)

     --edge_cov     use coverage from only first variant node position
                    (adjacent to head node) when calculating allele
                    depth (AD) values
                      default: use average coverage of first variant node

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

    *Currently dynamic, model-based thresholds have not been fully
    implemented. At the moment, selecting -m/--model will result in coverage
    histograms and genomescope plots being generated, but all calls will
    still be based on the static parameters.*

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

     --low_cov                   call GT based solely on allele with highest
                                   coverage (--min_tot_cov still applies)
                                   default: disabled

     --min_tot_cov               minimum total coverage (all alleles)
                                   default: 3

     --max_tot_cov               maximum total coverage (all alleles)
                                   default: no maxiumum

     --max_low_cov_tot_cov       maximum "low coverage" total coverage
                                   default: 9

     --min_low_cov_allele_count  minimum "low coverage" allele count
                                   default: 3

     --min_high_cov_allele_pct   minimum "high coverage" allele percent
                                   default: 10

  help:

     -h --help    display help menu

---

### pack_table_to_seg_cov.pl

Usage:

    pack_table_to_seg_cov.pl [options] > pack.seg.cov

    pack_table_to_seg_cov.pl [options] | gzip > pack.seg.cov.gz

Options:

     -p --pack    vg pack table file

     -s --sort    sort by node (increases memory consumption and output file
                  compressibility)
                    default: disabled

     -h --help    display help menu

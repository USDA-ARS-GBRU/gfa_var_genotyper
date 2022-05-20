# gfa_var_genotyper

### gfa_var_genotyper.pl

Usage:

    gfa_var_genotyper.pl [options] > genotypes

Options:

  general options:

     -v --var       gfa variants file (required)

     -p --pack      vg pack edge table(s)
                      ex: -p sample.1.pack.edge.table sample.2.pack.edge.table.gz

     --packlist     text file containing list of pack edge tables
                      (1 file per line)

    1 or more vg (https://github.com/vgteam/vg) pack edge tables may be
    specified using -p/--pack and/or --packlist. Pack edge tables are
    generated using the `vg pack -D` command. Pack edge tables may be
    uncompressed or compressed using either gzip or bzip2. (.gz or .bz2 file
    extension)

     --ploidy       1 (haploid) or 2 (diploid) currently supported
                      default: 1

     --rm_inv_head  remove variants with inverted head node
                      default: disabled

    Inversion variants can result in a negative sign prefix in the 'POS'
    field head node id. Conventionally, this field represents the position
    in the reference genome and negative values may cause issues with tools
    that use vcfs. To remove such variants from vcf ouput, use
    --rm_inv_head. Note, the reciprocal variant of the inversion should be
    called regardless, so the variant information is still retained for most
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

### gfa_variants.pl

Usage:

    gfa_variants.pl -g graph.gfa [options] > graph.vcf

Description:

    gfa_variants.pl writes node-based GFA graph variants to a VCF-like file

    Variants are based on branch positions in the graph. gfa_variants.pl
    understands when a reverse-oriented variant is the same as a forward
    variant and so produces biologically appropriate variants within long
    inversions that are ortholgous despite their structural difference. This
    behavior constrasts to available tools that we know of.

    Output is VCF-like, but instead of the 'POS' field being expressed in
    linear coordinates they are instead ordered graph node ids. If you have
    used xmfa_tools.pl with sorting enabled
    (https://github.com/brianabernathy/xmfa_tools), the node order will
    reflect the primary sort reference that you specified (then secondary,
    etc.). In other words, they are collinear with the genome sequence of
    that reference. Though these node ids do not represent useful physical
    distances, they can be determined using gfa_nodes_to_linear_coords.pl
    and vcf_node_to_linear_coords.pl to convert graph-based vcfs to linear
    coordinates.

    The gfa_variants.pl output POS field refers to the node id corresponding
    to the 'head node'. 'Allelic nodes', which serve as the REF and ALT
    fields, are the two possible branches extending from the head node. If
    these two branches return to the same node over after a single allelic
    node, node sequences are used to determine the explicit base change.
    Otherwise, variants are encoded as 'long', 'dense' or 'multipath' with
    the relevant node id suffix.

    Simple indels are encoded differently than conventional vcf format.
    Instead of giving the reference base and the insertion (or vice versus
    for deletions), gfa_variants.pl uses the '-' character and the indel
    sequence. We find this more intuitive and in keeping with the variation
    graph concept, but this difference may break tools expecting a
    traditional vcf.

    Inversion variants can result in a negative sign prefix in the 'POS'
    field head node id. Conventionally, this field represents the position
    in the reference genome and negative values may cause issues with tools
    that use vcfs. Note, the reciprocal variant of the inversion should be
    called regardless, so the variant information is still retained for most
    practical purposes.

Options:

  general options:

     -g --gfa        genome gfa file, vg-based (required)

     -c --chr_delim  pattern used to split genotype from chromosome in path names
                       default: '\.'

     -h --help       display help menu

---

### gfa_nodes_to_linear_coords.pl

Usage:

    gfa_nodes_to_linear_coords.pl [options] > gfa_node_coords.out

Description: 

    gfa_nodes_to_linear_coords.pl uses an input GFA file to generate
    linear reference coordinates for each graph node

Options:

  general options:

     -g --gfa       genome gfa file, vg-based (required)

     -h --help      display help menu


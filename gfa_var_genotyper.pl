#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $gfa_var_file;
my $deconstruct_format;
my $pack_file;
my $pack_label;
my @pack_files = ();
my %pack_labels = ();
my $pack_list_file;
my %node_lens = ();
my %vcf_nodes = ();
my $use_low_cov;
my $min_tot_cov = 3;
my $max_tot_cov;
my $max_low_cov_tot_cov = 9;
my $min_low_cov_allele_count = 3;
my $min_high_cov_allele_pct = 10;
my $ploidy = 1;
my $rm_inv_head = 0;
my $use_model = 0;
my $model_dir = 'gfa_var_genotyper_models';
my $gs_path;
my $kmer_size = 31;
my $help;

my %rec_gt_nodes = ();
my %rec_id_allele_counts = ();
my %pack_covs = ();
my @pack_labels = ();
my %pack_label_counts = ();
my %pack_file_counts = ();

parse_args();
parse_gfa_var_file();

foreach my $pack_file (@pack_files) {
	parse_pack_file($pack_file);
}

print_gfa_var_gts();

exit(0);


sub print_gfa_var_gts {
	my $gfa_var_fh;

	if ($gfa_var_file =~ /\.gz$/) {
		open($gfa_var_fh, '-|', "gzip -dc $gfa_var_file") or error("can't open $gfa_var_file: $!");
	}

	elsif ($gfa_var_file =~ /\.bz2$/) {
		open($gfa_var_fh, '-|', "bzip2 -dc $gfa_var_file") or error("can't open $gfa_var_file: $!");
	}

	else {
		open($gfa_var_fh, '<', $gfa_var_file) or error("can't open $gfa_var_file: $!");
	}

	while (my $line = <$gfa_var_fh>) {
		chomp($line);

		if ($line =~ /^##/) {
			print(STDOUT "$line\n");

			next();
		}

		my ($chr, $pos, $id, $ref, $alts, $qual, $filter, $info, $format, @gts) = split(/\t/, $line);

		if (! defined($id)) {
			next();
		}

		if ($chr eq '#CHROM') {
			print(STDOUT join("\t", $chr, $pos, $id, $ref, $alts, $qual, $filter, $info, $format, @gts, @pack_labels), "\n");

			next();
		}

		if ($pos =~ /\-/) {
			if ($rm_inv_head == 1) {
				next();
			}

			$pos =~ s/^\-//;
			$pos =~ s/\-$//;

			my $warning = 'inverted head node found, see --rm_inv_head option for more information';
			my $info_str = 'inverted_head_node';

			if ($info eq '.') {
				$info = "$info_str";
			}

			else {
				$info .= ";$info_str";
			}

			warning("$warning\n\tchr: $chr pos: $pos");
		}

		my $rec_id = "$chr\t$pos\t$id";
		my $allele_count = $rec_id_allele_counts{$rec_id};

		if (! defined($allele_count)) {
			next();
		}

		$format .= ':AD';

		print(STDOUT join("\t", $chr, $pos, $id, $ref, $alts, $qual, $filter, $info, $format, @gts));

		foreach my $pack_label (@pack_labels) {
			my %gt_cov = ();
			my %gt_cov_order = ();
			my $tot_cov = 0;
			my @ads = ();

			foreach my $gt (0..$allele_count - 1) {
				my $cov = 0;

				if (exists($rec_gt_nodes{$rec_id}{$gt})) {
					my $link_nodes = $rec_gt_nodes{$rec_id}{$gt};
            		my ($head_node, $var_node) = split(/\t/, $link_nodes);

					if (exists($pack_covs{$pack_label}{$head_node}{$var_node})) {
						$cov = $pack_covs{$pack_label}{$head_node}{$var_node};
					}

					elsif (exists($pack_covs{$pack_label}{$var_node}{$head_node})) {
						$cov = $pack_covs{$pack_label}{$var_node}{$head_node};
					}

					else {
						#warning("pack edge node pair not found: pack_label: $pack_label\tnode1: $head_node\tnode2: $var_node");
					}
				}

				$gt_cov{$gt} = $cov;
				$gt_cov_order{$cov}{$gt}++;
				$tot_cov += $cov;

				push(@ads, $cov);
			}

			my $ref_gt = 0;
			my $alt_gt = 1;
			my $gt_rec = '.';

			if ($ploidy == 2) {
				$gt_rec .= '/.';
			}

			if ($allele_count > 2) {
				# multiallelic: select top 2 cov alleles to represent ref/alt
				my $count = 0;

				foreach my $cov (reverse sort { $a <=> $b } keys %gt_cov_order) {
					foreach my $gt (sort { $a <=> $b } keys %{$gt_cov_order{$cov}}) {
						if ($count == 0) {
							$ref_gt = $gt;
						}

						elsif ($count == 1) {
							$alt_gt = $gt;
						}

						else {
							last();
						}

						$count++;
					}
				}
			}

			if ($tot_cov < $min_tot_cov) {
				# skip: low cov
			}

			elsif (defined($max_tot_cov) && $tot_cov > $max_tot_cov) {
				# skip: high cov
			}

			else {
				my $ref_cov = 0;
				my $ref_cov_pct = 0;
				my $alt_cov = 0;
				my $alt_cov_pct = 0;

				if (exists($gt_cov{$ref_gt})) {
					$ref_cov = $gt_cov{$ref_gt};
					$ref_cov_pct = $ref_cov / $tot_cov * 100;
				}

				if (exists($gt_cov{$alt_gt})) {
					$alt_cov = $gt_cov{$alt_gt};
					$alt_cov_pct = $alt_cov / $tot_cov * 100;
				}


				if (defined($use_low_cov)) {
					if ($ref_cov > $alt_cov) {
						$gt_rec = "$ref_gt";
					}

					elsif ($alt_cov > $ref_cov) {
						$gt_rec = "$alt_gt";
					}
				}

				elsif ($tot_cov <= $max_low_cov_tot_cov) {
					if ($ref_cov >= $min_low_cov_allele_count && $alt_cov >= $min_low_cov_allele_count) {
						if ($ploidy >= 2) {
							$gt_rec = "$ref_gt/$alt_gt";
						}
					}

					elsif ($ref_cov >= $min_low_cov_allele_count) {
						$gt_rec = "$ref_gt";

						if ($ploidy == 2) {
							$gt_rec = "./$ref_gt";
						}
					}

					else {
						$gt_rec = "$alt_gt";

						if ($ploidy == 2) {
							$gt_rec = "./$alt_gt";
						}
					}
				}

				else {
					if ($ref_cov_pct >= $min_high_cov_allele_pct && $alt_cov_pct >= $min_high_cov_allele_pct) {
						if ($ploidy >= 2) {
							$gt_rec = "$ref_gt/$alt_gt";
						}
					}

					elsif ($ref_cov_pct < $min_high_cov_allele_pct) {
						$gt_rec = "$alt_gt";

						if ($ploidy == 2) {
							$gt_rec = "$alt_gt/$alt_gt";
						}
					}

					else {
						$gt_rec = "$ref_gt";

						if ($ploidy == 2) {
							$gt_rec = "$ref_gt/$ref_gt";
						}
					}
				}
			}

			print(STDOUT "\t$gt_rec:", join(',', @ads));
		}

		print(STDOUT "\n");
	}

	close($gfa_var_fh);

	return(0);
}


sub parse_pack_file {
	my $pack_file = shift();

	my %cov_histo = ();

	my $pack_fh;

	if ($pack_file =~ /\.gz$/) {
		open($pack_fh, '-|', "gzip -dc $pack_file") or error("can't open $pack_file: $!");
	}

	elsif ($pack_file =~ /\.bz2$/) {
		open($pack_fh, '-|', "bzip2 -dc $pack_file") or error("can't open $pack_file: $!");
	}

	else {
		open($pack_fh, '<', $pack_file) or error("can't open $pack_file: $!");
	}

	my $pack_label;

	if (exists($pack_labels{$pack_file})) {
		$pack_label = $pack_labels{$pack_file};
	}

	else {
		$pack_label = $pack_file;
		$pack_label =~ s/^.*\///;
		$pack_label =~ s/\t/_/g;
		$pack_label_counts{$pack_label}++;

		if (scalar($pack_label_counts{$pack_label}) > 1) {
			warning("duplicate label found: $pack_label");
		}
	}

	push(@pack_labels, $pack_label);

	print(STDERR "processing pack file: $pack_file\tlabel: $pack_label\n");

	while (my $line = <$pack_fh>) {
		chomp($line);

		if ($line =~ /^from\.id/) {
			next();
		}

		my ($node1_id, $node1_start, $node2_id, $node2_stop, $cov) = split(/\t/, $line);

		if ((defined($node1_id) && defined($node1_start) && defined($node2_id) && defined($node2_stop) && defined($cov)) == 0) {
			error("invalid format found in pack edge table: $pack_file\n\t$line");
		}

		if ($use_model == 1) {
			if ($cov > 0) {
				$cov_histo{$cov}++;
			}
		}

		if (! exists($vcf_nodes{$node1_id}) || ! exists($vcf_nodes{$node2_id})) {
			next();
		}

		$pack_covs{$pack_label}{$node1_id}{$node2_id} = $cov;
	}

	close($pack_fh);


	if ($use_model == 1) {
		my $hist_dir = "$model_dir/cov_hist";

		if (! -e $hist_dir) {
			mkdir($hist_dir) or error("$!");
		}

		my $cov_histo_file = "$hist_dir/$pack_label.cov.histo";

		open(COV, '>', $cov_histo_file) or error("can't open cov histo file: $!");

		foreach my $cov (sort { $a <=> $b } keys %cov_histo) {
			print(COV "$cov\t$cov_histo{$cov}\n");
		}

		close(COV);


		my $gs_dir = "$model_dir/genomescope";

		if (! -e $gs_dir) {
			mkdir($gs_dir) or error("$!");
		}

		my $gs_base = "$gs_dir/$pack_label";

		my $gs_cmd = "$gs_path -i $cov_histo_file -o $gs_dir -k $kmer_size -p $ploidy -n $pack_label 1> $gs_base.gs.stdout 2> $gs_base.gs.stderr";

		system($gs_cmd);
	}

	return(0);
}


sub parse_gfa_var_file {
	my $gfa_var_fh;

	if ($gfa_var_file =~ /\.gz$/) {
		open($gfa_var_fh, '-|', "gzip -dc $gfa_var_file") or error("can't open $gfa_var_file: $!");
	}

	elsif ($gfa_var_file =~ /\.bz2$/) {
		open($gfa_var_fh, '-|', "bzip2 -dc $gfa_var_file") or error("can't open $gfa_var_file: $!");
	}

	else {
		open($gfa_var_fh, '<', $gfa_var_file) or error("can't open $gfa_var_file: $!");
	}

	while (my $line = <$gfa_var_fh>) {
		chomp($line);

		if ($line =~ /^#/) {
			next();
		}

		my ($chr, $pos, $id, $ref, $alts, $qual, $filter, $info, $format, @gts) = split(/\t/, $line);

		if (! defined($id)) {
			next();
		}

		my $rec_id = "$chr\t$pos\t$id";
		my @link_nodes = ();

		if (! defined($deconstruct_format)) {
			if ($info =~ /AT\=/) {
				$deconstruct_format = 1;
			}

			else {
				$deconstruct_format = 0;
			}
		}

		# PanPipes style graph variants
		if ($deconstruct_format == 0) {
			if ($pos =~ /\-/) {
				if ($rm_inv_head == 1) {
					next();
				}

				$pos =~ s/\-$//;
			}

			my $type;

			if ($ref =~ s/^complex\.//) {
				$type = 'complex';
			}

			elsif ($ref =~ s/^multiPath\.//i) {
				$type = 'multipath';
			}

			else {
				$type = 'snp';

				($ref, $alts) = split('-', $id);
			}

			my @alts = split(',', $alts);
			@link_nodes = ($ref, @alts);

			my $head_node = $pos;

			foreach my $index (0..$#link_nodes) {
				$link_nodes[$index] = "$head_node\t$link_nodes[$index]";
			}
		}

		# vg deconstruct/minigraph-cactus style graph variants
		elsif ($deconstruct_format == 1) {
			foreach my $info_field (split(';', $info)) {
				if ($info_field =~ /^AT\=/) {
					$info_field =~ s/^AT\=//;

					my %node_pairs = ();

					foreach my $allele_node_str (split(',', $info_field)) {
						$allele_node_str =~ s/[<>]/,/g;
						$allele_node_str =~ s/^,//;

						my @allele_nodes = split(',', $allele_node_str);

						foreach my $trailing_index (1..$#allele_nodes) {
							$node_pairs{$allele_nodes[$trailing_index - 1]}{$allele_nodes[$trailing_index]}++;
						}
					}

					foreach my $allele_node_str (split(',', $info_field)) {
						$allele_node_str =~ s/[<>]/,/g;
						$allele_node_str =~ s/^,//;

						my @allele_nodes = split(',', $allele_node_str);

						foreach my $trailing_index (1..$#allele_nodes) {
							my $leading_node = $allele_nodes[$trailing_index - 1];
							my $trailing_node = $allele_nodes[$trailing_index];

							if ($node_pairs{$leading_node}{$trailing_node} > 1) {
								next();
							}

							push(@link_nodes, "$leading_node\t$trailing_node");
							last();
						}
					}
				}
			}
		}

		$rec_id_allele_counts{$rec_id} = $#link_nodes + 1;

		foreach my $gt (0..$#link_nodes) {
			my $link_nodes = $link_nodes[$gt];
			my ($leading_node, $trailing_node) = split(/\t/, $link_nodes);

			$vcf_nodes{$leading_node}++;
			$vcf_nodes{$trailing_node}++;
			$rec_gt_nodes{$rec_id}{$gt} = $link_nodes;
		}
	}

	close($gfa_var_fh);

	return(0);
}


sub warning {
	my $msg = shift();

	if (defined($msg)) {
		print(STDERR "warning: $msg\n");
	}

	return(0);
}


sub error {
	my $msg = shift();

	if (defined($msg)) {
		print(STDERR "error: $msg\n");
	}

	exit(0);
}


sub arg_error {
	my $msg = shift();

	if (defined($msg)) {
		print(STDERR "error: $msg\n\n");
	}

	print(STDERR "use option -h or --help to display help menu\n");

	exit(0);
}


sub parse_args {
	if ($#ARGV == -1) {
		pod2usage(-exitval => 0, -verbose => 1);
	}

	GetOptions ('v|var=s' => \$gfa_var_file,
				'p|pack=s' => \$pack_file,
				'l|label=s' => \$pack_label,
				'packlist=s' => \$pack_list_file,
				'ploidy=i' => \$ploidy,
				'rm_inv_head' => \$rm_inv_head,
				'm|model' => \$use_model,
				'modeldir=s' => \$model_dir,
				'gs=s' => \$gs_path,
				'low_cov' => \$use_low_cov,
				'min_tot_cov=i' => \$min_tot_cov,
				'max_tot_cov=i' => \$max_tot_cov,
				'max_low_cov_tot_cov=i' => \$max_low_cov_tot_cov,
				'min_low_cov_allele_count=i' => \$min_low_cov_allele_count,
				'min_high_cov_allele_pct=i' => \$min_high_cov_allele_pct,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -verbose => 1);
	}

	if (! defined($gfa_var_file)) {
		arg_error('gfa variants file required');
	}

	if (defined($pack_file)) {
		push(@pack_files, $pack_file);
		$pack_file_counts{$pack_file}++;

		if (defined($pack_label)) {
			$pack_label =~ s/\t/_/g;
			$pack_labels{$pack_file} = $pack_label;
			$pack_label_counts{$pack_label}++;
		}
	}

	if (defined($pack_list_file)) {
		open(PACKLIST, '<', $pack_list_file) or error("can't open pack list file: $!");

		foreach my $line (<PACKLIST>) {
			chomp($line);

			my ($file, $label) = split(/\t/, $line);

			if (exists($pack_file_counts{$file})) {
				warning("duplicate pack file found: $file, only the first occurence will be processed");

				next();
			}

			push(@pack_files, $file);
			$pack_file_counts{$file}++;

			if (defined($label)) {
				$label =~ s/\t/_/g;
				$pack_labels{$file} = $label;
				$pack_label_counts{$label}++;
			}
		}

		close(PACKLIST);
	}

	if (! @pack_files) {
		arg_error('at least 1 pack table file required');
	}

	foreach my $label (keys %pack_label_counts) {
		if (scalar($pack_label_counts{$label}) > 1) {
			warning("duplicate pack label found: $label");
		}
	}

	if ($ploidy != 1 && $ploidy != 2) {
		arg_error('ploidy must be either 1 or 2');
	}

	if (defined($use_low_cov) && $ploidy > 1) {
		arg_error('ploidy must be 1 when --lowcov is used');
	}

	if ($use_model == 1) {
		$model_dir =~ s/\/$//;

		if (! -e $model_dir) {
			my $cmd = "mkdir -p $model_dir";

			system($cmd);
		}

		if (defined($gs_path)) {
			if (! -e $gs_path) {
				error("genomescope path: $gs_path does not exist");
			}
		}

		else {
			$gs_path = qx(which genomescope.R);

			chomp($gs_path);

			if (! defined($gs_path) || $gs_path eq '') {
				error("genomescope not found in \$PATH, specify using --gs option");
			}

			print(STDERR "using genomescope found at $gs_path\n");
		}
	}

	return(0);
}


__END__

=head1 NAME

gfa_var_genotyper.pl

=head1 SYNOPSIS

gfa_var_genotyper.pl [options] > genotypes

=head1 DESCRIPTION

gfa_var_genotyper.pl genotypes variants called directly from GFA files

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

=head2 general options

 -v --var       graph variants file (required)
                  can parse either vg deconstruct/minigraph-cactus
                  format or gfa_variants.pl (PanPipes) format 

1 or more vg (https://github.com/vgteam/vg) pack edge tables may be
specified using -p/--pack and/or --packlist. Pack edge tables are
generated using the `vg pack -D` command. Pack edge tables may be
uncompressed or compressed using either gzip or bzip2. (.gz or .bz2
file extension)

For a single pack table and label use -p/--pack and -l/--label. For
multiple pack table files and labels use --packlist.

 -p --pack      vg pack edge table
                  ex: -p sample.1.pack.edge.table
                  ex: -p sample.2.pack.edge.table.gz

 -l --label     vg pack edge table label displayed in vcf header

 --packlist     text file containing list of pack edge table files
                  and labels (labels optional but recommended)
                  (1 file and label per line, tab-delimited)

 --ploidy       1 (haploid) or 2 (diploid) currently supported
                  default: 1

Variants in gfa_variants.pl (PanPipes) format can contain inversion
variants, which can result in a negative sign prefix in the 'POS'
field head node id. Conventionally, this field represents the position
in the reference genome and negative values may cause issues with tools
that use vcfs. To remove such variants from vcf ouput, use --rm_inv_head.
Note, the reciprocal variant of the inversion should be called regardless,
so the variant information is still retained for most practical purposes. 

 --rm_inv_head  remove variants with inverted head node
                only applies to gfa_variants.pl (PanPipes) format
                  default: disabled

=head2 genotyping options

*Currently dynamic, model-based thresholds have not been fully
implemented. At the moment, selecting -m/--model will result in
coverage histograms and genomescope plots being generated, but
all calls will still be based on the static parameters.*

Genotyping thresholds can be set manually or dynamically by creating
coverage count histograms that are provide to GenomeScope for modeling.
(modeling is not recommended for GBS data) When modeling is enabled,
any samples that fail to generate a valid model will fallback to
using the static parameters. 

Static parameter genotyping sets a minimum total coverage (all alleles)
threshold. The remaining variants are processed as either low or high
total coverage. Low coverage variants are genotyped using minimum
allele counts. High coverage variants are genotyped using minimum
allele percentages.

=head3 modeling

 -m --model   use genomescope model to define genotyping thresholds
                default: disabled

 --modeldir   output directory for modeling files
                default: gfa_var_genotyper_models

 --gs         path to GenomeScope executable
                default: autodetect in $PATH (if available)

=head3 static parameters

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

=head2 help

 -h --help    display help menu

=cut

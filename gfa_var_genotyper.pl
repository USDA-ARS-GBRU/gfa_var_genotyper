#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $gfa_var_file;
my @pack_files = ();
my $pack_list_file;
my %node_lens = ();
my %gfa_nodes = ();
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

my %rec_node_gts = ();
my %rec_gt_nodes = ();
my %rec_id_allele_counts = ();
my %pack_covs = ();
my @file_bases = ();

parse_args();
parse_gfa_var_file();

foreach my $pack_file (@pack_files) {
	parse_pack_file($pack_file);
}

print_gfa_var_gts();

exit(0);


sub print_gfa_var_gts {
	open(GFA_VARS, '<', $gfa_var_file) or error("can't open gfa var file: $!");

	while (my $line = <GFA_VARS>) {
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
			print(STDOUT join("\t", $chr, $pos, $id, $ref, $alts, $qual, $filter, $info, $format, @gts, @file_bases), "\n");

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

		my $head_node = $pos;

		foreach my $file_base (@file_bases) {
			my %gt_cov = ();
			my %gt_cov_order = ();
			my $tot_cov = 0;
			my @ads = ();

			foreach my $gt (0..$allele_count - 1) {
				my $cov = 0;

				if (exists($rec_gt_nodes{$rec_id}{$gt})) {
					my $var_node = $rec_gt_nodes{$rec_id}{$gt};
					my $node1 = $head_node;
					my $node2 = $var_node;

					if ($node1 > $node2) {
						my $tmp = $node1;

						$node1 = $node2;
						$node2 = $tmp;
					}

					if (exists($pack_covs{$file_base}{$node1}{$node2})) {
						$cov = $pack_covs{$file_base}{$node1}{$node2};
					}

					else {
						warning("pack edge node not found: file_base: $file_base\tnode1: $node1\tnode2: $node2");
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

	return(0);
}


sub parse_pack_file {
	my $pack_file = shift();

	my %cov_histo = ();
	my $file_base = $pack_file;

	my $pack_fh;

	if ($pack_file =~ /\.gz$/) {
		open($pack_fh, '-|', "gzip -dc $pack_file") or error("can't open pack file: $!");
	}

	elsif ($pack_file =~ /\.bz2$/) {
		open($pack_fh, '-|', "bzip2 -dc $pack_file") or error("can't open pack file: $!");
	}

	else {
		open($pack_fh, '<', $pack_file) or error("can't open pack file: $!");
	}

	print(STDERR "processing pack file: $pack_file\n");

	$file_base =~ s/^.*\///;
	$file_base =~ s/\..*$//;

	push(@file_bases, $file_base);

	while (my $line = <$pack_fh>) {
		chomp($line);

		if ($line =~ /^from\.id/) {
			next();
		}

		my ($node1_id, $node1_orientation, $node2_id, $node2_orientation, $cov) = split(/\t/, $line);

		if ((defined($node1_id) && defined($node1_orientation) && defined($node2_id) && defined($node2_orientation) && defined($cov)) == 0) {
			error("invalid format found in pack edge table: $pack_file\n\t$line");
		}

		if ($use_model == 1) {
			if ($cov > 0) {
				$cov_histo{$cov}++;
			}
		}

		if (! exists($gfa_nodes{$node1_id}) || ! exists($gfa_nodes{$node2_id})) {
			next();
		}

		$pack_covs{$file_base}{$node1_id}{$node2_id} = $cov;
	}

	close($pack_fh);


	if ($use_model == 1) {
		my $hist_dir = "$model_dir/cov_hist";

		if (! -e $hist_dir) {
			mkdir($hist_dir) or error("$!");
		}

		my $cov_histo_file = "$hist_dir/$file_base.cov.histo";

		open(COV, '>', $cov_histo_file) or error("can't open cov histo file: $!");

		foreach my $cov (sort { $a <=> $b } keys %cov_histo) {
			print(COV "$cov\t$cov_histo{$cov}\n");
		}

		close(COV);


		my $gs_dir = "$model_dir/genomescope";

		if (! -e $gs_dir) {
			mkdir($gs_dir) or error("$!");
		}

		my $gs_base = "$gs_dir/$file_base";

		my $gs_cmd = "$gs_path -i $cov_histo_file -o $gs_dir -k $kmer_size -p $ploidy -n $file_base 1> $gs_base.gs.stdout 2> $gs_base.gs.stderr";

		system($gs_cmd);
	}

	return(0);
}


sub parse_gfa_var_file {
	open(GFA_VARS, '<', $gfa_var_file) or error("can't open gfa var file: $!");

	while (my $line = <GFA_VARS>) {
		chomp($line);

		if ($line =~ /^#/) {
			next();
		}

		my ($chr, $pos, $id, $ref, $alts, $qual, $filter, $info, $format, @gts) = split(/\t/, $line);

		if (! defined($id)) {
			next();
		}

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
		my @alleles = ($ref, @alts);
		my $head_node = $pos;
		my $rec_id = "$chr\t$pos\t$id";

		$rec_id_allele_counts{$rec_id} = $#alleles + 1;
		$gfa_nodes{$head_node}++;

		foreach my $gt (0..$#alleles) {
			my $node = $alleles[$gt];

			$gfa_nodes{$node}++;
			$rec_gt_nodes{$rec_id}{$gt} = $node;
			$rec_node_gts{$rec_id}{$node} = $gt;
		}
	}

	close(GFA_VARS);

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

	print("use option -h or --help to display help menu\n");

	exit(0);
}


sub parse_args {
	if ($#ARGV == -1) {
		pod2usage(-exitval => 0, -verbose => 1);
	}

	GetOptions ('v|var=s' => \$gfa_var_file,
				'p|pack=s{,}' => \@pack_files,
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

	if (defined($pack_list_file)) {
		open(PACKLIST, '<', $pack_list_file) or error("can't open pack list file: $!");

		foreach my $line (<PACKLIST>) {
			chomp($line);

			push(@pack_files, $line);
		}

		close(PACKLIST);
	}

	if (! @pack_files) {
		arg_error('at least 1 pack table file required');
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

 -v --var       gfa variants file (required)

 -p --pack      vg pack edge table(s)
                  ex: -p sample.1.pack.edge.table -p sample.2.pack.edge.table.gz

 --packlist     text file containing list of pack edge tables
                  (1 file per line)

1 or more vg (https://github.com/vgteam/vg) pack edge tables may be
specified using -p/--pack and/or --packlist. Pack edge tables are
generated using the `vg pack -D` command. Pack edge tables may be
uncompressed or compressed using either gzip or bzip2. (.gz or .bz2
file extension)

 --ploidy       1 (haploid) or 2 (diploid) currently supported
                  default: 1

 --rm_inv_head  remove variants with inverted head node
                  default: disabled
                 
Inversion variants can result in a negative sign prefix in the 'POS'
field head node id. Conventionally, this field represents the position
in the reference genome and negative values may cause issues with tools
that use vcfs. To remove such variants from vcf ouput, use --rm_inv_head.
Note, the reciprocal variant of the inversion should be called regardless,
so the variant information is still retained for most practical purposes. 

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

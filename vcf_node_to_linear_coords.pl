#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $vcf_file;
my $coords_file;
my @genos;
my $chr_delim = '\.';
my $chr_prefix;
my $help;

my %coords = ();
my @placeholder_seqs = ();

$placeholder_seqs[0] = 'A' x 10;
$placeholder_seqs[1] = 'C' x 10;
$placeholder_seqs[2] = 'G' x 10;
$placeholder_seqs[3] = 'T' x 10;

parse_args();
parse_coords();
parse_vcf();

exit(0);


sub parse_coords {
	my %coords_genos = ();
	my %coords_chrs = ();
	my $print_warning = 1;
	my $fh;

	if ($coords_file =~ /\.gz$/) {
		open($fh, '-|', "gzip -dc $coords_file") or error("can't read $coords_file: $!");
	}

	elsif ($coords_file =~ /\.bz2$/) {
		open($fh, '-|', "bzip2 -dc $coords_file") or error("can't read $coords_file: $!");
	}

	else {
		open($fh, '<', $coords_file) or error("can't open pack table file: $!");
	}

	while (my $line = <$fh>) {
		if ($line =~ /^path_name/) {
			next();
		}

		chomp($line);

		my ($path_name, $start, $stop, $node, $orientation) = split(/\t/, $line);
		my ($path_geno, $path_chr) = split(/$chr_delim/, $path_name);

		if (defined($chr_prefix)) {
			$path_chr = "${chr_prefix}${path_chr}";
		}

		if (! defined($path_chr)) {
			if ($print_warning == 1) {
				print(STDERR "warning: can't determine path chr from $path_name, using chr0\n");
				print(STDERR "\tif $coords_file contains chromosomes in the path names, use -d/--delim option to configure path name parsing\n");

				$print_warning = 0;
			}

			$path_chr = 'chr0';
		}

		$coords_chrs{$path_chr}++;
		$coords_genos{$path_geno}++;
		$coords{$node}{$path_geno} = "$path_chr:$start:$stop:$orientation";
	}

	close($fh);


	print(STDERR "$coords_file genotypes processed: ", join(', ', sort(keys %coords_genos)), "\n");
	print(STDERR "$coords_file chromosomes processd: ", join(', ', sort(keys %coords_chrs)), "\n");

	foreach my $geno (@genos) {
		if (! exists($coords_genos{$geno})) {
			error("genotype: $geno not found in $coords_file");
		}
	}

	return(0);
}


sub parse_vcf {
	open(VCF, '<', $vcf_file) or error("can't read $vcf_file: $!");

	while (my $line = <VCF>) {
		if ($line =~ /^#/) {
			print(STDOUT "$line");

			next();
		}

		chomp($line);

		my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split(/\t/, $line);
		my $head_node = $pos;
		my $var_nodes = $id;
		my $var_type = '';

		if ($ref eq '-' || $alt eq '-') {
			$var_type = 'indel';
		}

		elsif ($ref =~ s/^complex\.//) {
			$var_type = 'complex';
		}

		elsif ($ref =~ s/^long\.//) {
			$var_type = 'long';
		}

		elsif ($ref =~ s/^dense\.//) {
			$var_type = 'dense';
		}

		elsif ($ref =~ s/^multipath\.//i) {
			$var_type = 'multipath';
		}

		elsif (length($ref) == 1 && $ref =~ /[ACGTN]/i && length($alt) == 1 && $alt =~ /[ACGTN]/i) {
			$var_type = 'snp';
		}

		else {
			# multi-nucleotide variant
			$var_type = 'mnv';
		}

		my $orig_ref = $ref;
		my $orig_alt = $alt;

		$info = "VAR_TYPE=$var_type;HEAD_NODE=$head_node;VAR_NODES=$var_nodes;REF=$ref;ALT=$alt";

		if ($var_type eq 'complex' || $var_type eq 'long' || $var_type eq 'dense' || $var_type eq 'multipath') {
			$ref = $placeholder_seqs[0];

			my @alt_nodes = split(/,/, $alt);
			my @mod_alts = ();

			foreach my $index (1..($#alt_nodes + 1)) {
				push(@mod_alts, $placeholder_seqs[$index]);
			}

			$alt = join(',', @mod_alts);
		}


		my $found_linear_pos = 0;

		foreach my $path_geno (@genos) {
			if (exists($coords{$head_node}{$path_geno})) {
				my ($chr, $start, $stop, $orientation) = split(/:/, $coords{$head_node}{$path_geno});
				my $path_pos = $stop;

				if ($orientation eq '-') {
					$path_pos = $start;
				}

				$pos = $path_pos;

				# increment SNP pos to 1 bp beyond end of head node
				if ($var_type eq 'snp') {
					$pos++;
				}

				$id = "$path_geno:$chr:$pos";
				$found_linear_pos = 1;

				last();
			}

			#print(STDERR "coord not found for $path_geno: ", join("\t", $chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format), "\n");
		}

		if ($found_linear_pos == 1) {
			print(STDOUT join("\t", $chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples), "\n");
		}
	}

	close(VCF);

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

	GetOptions ('c|coords=s' => \$coords_file,
		'v|vcf=s' => \$vcf_file,
		'g|geno=s{,}' => \@genos,
		'd|delim=s' => \$chr_delim,
		'p|prefix=s' => \$chr_prefix,
		'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -verbose => 1);
	}

	if (! defined($coords_file)) {
		arg_error('coords file required');
	}

	if (! defined($vcf_file)) {
		arg_error('vcf file required');
	}

	if (scalar @genos == 0) {
		arg_error('at least one genotype is required');
	}

	return(0);
}


__END__

=head1 NAME

vcf_node_to_linear_coords.pl

=head1 SYNOPSIS

vcf_node_to_linear_coords.pl -c coords.file -v node.based.vcf -g genotype [options] > geno.coords.vcf

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -c --coords     coords file (required)
                   produced by gfa_nodes_to_linear_coords.pl
                   can be in gzip or bzip2 format

 -v --vcf        node-based vcf file (required)
                   produced by gfa_var_genotyper.pl

 -g --geno       genotype(s) used as reference for linear coordinates
                   one genotype required, can enter multiple genotypes

                   when multiple genotypes are provided, coordinates for the
                   first genotype will be output, when available, followed by
                   the second genotype, then third, etc...

                   ex: -g geno1
                   ex: -g geno1 -g geno2 ... -g genoX
                   

 -d --delim      pattern used to split genotype from chromosome in path names
                   default: '\.'

 -p --prefix     prefix prepended to chromosome names
                   useful if part of the chromosome name is used to split path
                   names and therefore discarded

                   ex: -d '\.chr' -p 'chr'
                   default: none

 -h --help       display help menu

=cut

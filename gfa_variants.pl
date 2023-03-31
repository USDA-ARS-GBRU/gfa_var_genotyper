#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $gfa_file;
my $chr_delim = '\.';
my $chr_prefix;
my $help;

my %paths = ();
my %seqs = ();
my %inv_head_nodes = ();
my @genos = ();
my $use_w_lines = 0;

parse_args();
parse_gfa_file();
call_vars();

exit(0);


sub call_vars {
	my $print_header = 1;

	foreach my $chr (sort keys %paths) {
		my %node_counts = ();
		my %links = ();
		my %geno_linked_nodes = ();
		my @ordered_nodes = ();

		foreach my $geno (sort keys %{$paths{$chr}}) {
			my ($rec_type, $path_name, $nodes, $overlaps);
			my @nodes = ();

			if ($use_w_lines == 0) {
				($rec_type, $path_name, $nodes, $overlaps) = split(/\t/, $paths{$chr}{$geno});
				@nodes = split(',', $nodes);
			}

			else {
				my ($rec_type, $sample, $hap_index, $seq_id, $seq_start, $seq_end, $walk) = split(/\t/, $paths{$chr}{$geno});
				my @walk_comps =  split(/(<|>)/,  $walk);

				$path_name = $seq_id;

				my $index = 1;

				while ($index < $#walk_comps) {
					my $orientation = '+';

					if ($walk_comps[$index] eq '<') {
						$orientation = '-';
					}

					my $node = "$walk_comps[$index + 1]"."$orientation";

					push(@nodes, $node);

					$index += 2;
				}
			}

			foreach my $node_index (0..$#nodes - 1) {
				my $node1 = $nodes[$node_index];
				my $node2 = $nodes[$node_index + 1];
				my $node1_orientation = '+';
				my $node2_orientation = '+';

				if ($node1 =~ s/([\-\+])$//) {
					$node1_orientation = $1;
				}

				if ($node2 =~ s/([\-\+])$//) {
					$node2_orientation = $1;
				}

				if ($node1_orientation eq '-') {
					if ($node2_orientation eq '-') {
						# double inversion, invert nodes and process as +/+ link
						my $tmp = $node1;
						$node1 = $node2;
						$node2 = $tmp;
					}

					elsif ($node2_orientation eq '+') {
						# negate node1 id to differentiate from +node1
						$node1 = "-$node1";

						$inv_head_nodes{$node1} = "path: $path_name\tnode1: $node1\tnode2: $node2";
					}
				}

				if (! exists($links{$node1}{$node2})) {
					$node_counts{$node1}++;
					$links{$node1}{$node2} = $node_counts{$node1} - 1;
					push(@ordered_nodes, $node1);
				}

				$geno_linked_nodes{$node1}{$geno} = $node2;
			}
		}


		if ($print_header == 1) {
			print(STDOUT "## fileformat=VCFvPanPipes\n");
			print(STDOUT "## https://github.com/USDA-ARS-GBRU/PanPipes\n");
			print(STDOUT "## input graph: $gfa_file\n");
			print(STDOUT '#', join("\t", 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', @genos), "\n");

			$print_header = 0;
		}


		# head_node -> var_node -> link_node
		my %proc_head_nodes = ();

		foreach my $head_node (@ordered_nodes) {
			if (exists($proc_head_nodes{$head_node})) {
				next();
			}

			$proc_head_nodes{$head_node}++;

			# monomorphic
			if ($node_counts{$head_node} <= 1) {
				next();
			}

			if (exists($inv_head_nodes{$head_node})) {
				print(STDERR "warning: inverted head node\t$inv_head_nodes{$head_node}\n");
			}

			my $ref_node;
			my $ref_seq;
			my @alt_nodes;
			my @alt_seqs;

			foreach my $var_node (sort { $a <=> $b } keys %{$links{$head_node}}) {
				if ($links{$head_node}{$var_node} == 0) {
					$ref_node = $var_node;
				}

				else {
					my $alt_index = $links{$head_node}{$var_node} - 1;
					$alt_nodes[$alt_index] = $var_node;
				}
			}


			my $alt_count = $#alt_nodes + 1;
			my $ref_link_count = keys %{$links{$ref_node}};
			my $alt_link_count = keys %{$links{$alt_nodes[0]}};

			if ($alt_count == 1 && $ref_link_count == 1 && $alt_link_count == 1) {
				my @ref_link_nodes = keys %{$links{$ref_node}};
				my $ref_link_node = $ref_link_nodes[0];
				my @alt_link_nodes = keys %{$links{$alt_nodes[0]}};
				my $alt_link_node = $alt_link_nodes[0];

				if ($ref_link_node eq $alt_nodes[0]) {
					# insertion as ref
					$ref_seq = $seqs{$ref_node};
					@alt_seqs = ('-');
				}

				elsif ($alt_link_node eq $ref_node) {
					# deletion as ref
					$ref_seq = ('-');
					@alt_seqs = ($seqs{$alt_nodes[0]});
				}

				elsif ($ref_link_node eq $alt_link_node) {
					# substitution
					$ref_seq = $seqs{$ref_node};
					@alt_seqs = ($seqs{$alt_nodes[0]});
				}

				else {
					$ref_seq = "long\.$ref_node";
					@alt_seqs = @alt_nodes;
				}
			}

			elsif ($alt_count > 1) {
				$ref_seq = "multipath\.$ref_node";
				@alt_seqs = @alt_nodes;
			}

			else {
				$ref_seq = "dense\.$ref_node";
				@alt_seqs = @alt_nodes;
			}


			my $vcf_rec = join("\t", $chr, $head_node, "$ref_node-$alt_nodes[0]", $ref_seq, join(',', @alt_seqs), '.', '.', '.', 'GT');

			foreach my $geno (@genos) {
				my $call = '.';

				if (exists($geno_linked_nodes{$head_node}{$geno})) {
					my $linked_node = $geno_linked_nodes{$head_node}{$geno};

					$call = $links{$head_node}{$linked_node};
				}

				$vcf_rec .= "\t$call";
			}

			print(STDOUT "$vcf_rec\n");
		}

		delete($paths{$chr});
	}

	return(0);
}


sub parse_gfa_file {
	my %genos = ();

	open(GFA, '<', $gfa_file) or error("can't open $gfa_file: $!");

	while (my $line = <GFA>) {
		chomp($line);

		if ($line =~ /^S\t/) {
			my ($rec_type, $seg, $seq) = split(/\t/, $line);

			$seqs{$seg} = $seq;
		}

		elsif ($line =~ /^P\t/) {
			my ($rec_type, $path_name, $nodes, $overlaps) = split(/\t/, $line);
			my ($path_geno, $path_chr) = split(/$chr_delim/, $path_name);

			if (defined($chr_prefix)) {
				$path_chr = "${chr_prefix}${path_chr}";
			}

			if (! defined($path_chr)) {
				print(STDERR "warning: can't determine path chr from $path_name, using chr0\n");
				print(STDERR "\tif $gfa_file contains chromosomes in the path names, use -d/--delim option to configure path name parsing\n");

				$path_chr = 'chr0';
			}

			$genos{$path_geno}++;
			$paths{$path_chr}{$path_geno} = $line;
		}

		elsif ($line =~ /^W\t/) {
			my ($rec_type, $sample, $hap_index, $seq_id, $seq_start, $seq_end, $walk) = split(/\t/, $line);

			if ($sample eq '_MINIGRAPH_') {
				next();
			}

			$genos{$sample}++;
			$paths{$seq_id}{$sample} = $line;
			$use_w_lines = 1;
		}
	}

	close(GFA);

	@genos = sort(keys %genos);

	print(STDERR "the following genotypes were found: ", join(', ', @genos), "\n");
	print(STDERR "the following chromosomes were found: ", join(', ', sort(keys %paths)), "\n");

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

	GetOptions ('g|gfa=s' => \$gfa_file,
				'd|delim=s' => \$chr_delim,
				'p|prefix=s' => \$chr_prefix,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -verbose => 1);
	}

	if (! defined($gfa_file)) {
		arg_error('gfa file required');
	}

	return(0);
}


__END__

=head1 NAME

gfa_variants.pl

=head1 SYNOPSIS

gfa_variants.pl -g graph.gfa [options] > graph.vcf

gfa_variants.pl writes node-based GFA graph variants to a VCF-like file

Variants are based on branch positions in the graph. gfa_variants.pl understands when a reverse-oriented variant is the same as a forward variant and so produces biologically appropriate variants within long inversions that are ortholgous despite their structural difference.  This behavior constrasts to available tools that we know of.  

Output is VCF-like, but instead of the 'POS' field being expressed in linear coordinates they are instead ordered graph node ids.  If you have used xmfa_tools.pl with sorting enabled (https://github.com/brianabernathy/xmfa_tools), the node order will reflect the primary sort reference that you specified (then secondary, etc.).  In other words, they are collinear with the genome sequence of that reference.  Though these node ids do not represent useful physical distances, they can be determined using gfa_nodes_to_linear_coords.pl and vcf_nodes_to_linear_coords.pl to convert graph-based vcfs to linear coordinates.

The gfa_variants.pl output POS field refers to the node id corresponding to the 'head node'.  'Allelic nodes', which serve as the REF and ALT fields, are the two possible branches extending from the head node.  If these two branches return to the same node over after a single allelic node, node sequences are used to determine the explicit base change. Otherwise, variants are encoded as 'long', 'dense' or 'multipath' with the relevant node id suffix.

Simple indels are encoded differently than conventional vcf format. Instead of giving the reference base and the insertion (or vice versus for deletions), gfa_variants.pl uses the '-' character and the indel sequence.  We find this more intuitive and in keeping with the variation graph concept, but this difference may break tools expecting a traditional vcf.

Inversion variants can result in a negative sign prefix in the 'POS' field head node id. Conventionally, this field represents the position in the reference genome and negative values may cause issues with tools that use vcfs.  Note, the reciprocal variant of the inversion should be called regardless, so the variant information is still retained for most practical purposes.

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -g --gfa        genome gfa file, vg-based (required)

 -d --delim      pattern used to split genotype from chromosome in path names
                   default: '\.'

 -p --prefix     prefix prepended to chromosome names
                   useful if part of the chromosome name is used to split path
                   names and therefore discarded
                   ex: -d '\.chr' -p 'chr'
                   default: none

 -h --help       display help menu

=cut

#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $gfa_file;
my $chr_delim_regex = '\.';
my $chr_delim_str = '.';
my $help;

parse_args();
parse_gfa_file();

exit(0);


sub parse_gfa_file {
	my %path_nodes = ();
	my %path_node_orientations = ();
	my %node_lens = ();

	open(GFA, '<', $gfa_file) or error("can't open $gfa_file: $!");

	print(STDOUT "path_name\tstart\tstop\tnode\torientation\n");

	while (my $line = <GFA>) {
		if ($line =~ /^P/) {
			chomp($line);

			my ($rec_type, $path_name, $nodes, $overlaps) = split(/\t/, $line);
			my @nodes = split(',', $nodes);

			foreach my $node (@nodes) {
				my $orientation = '+';

				if ($node =~ s/([\-\+])$//) {
					$orientation = $1;
				}

				$path_node_orientations{$path_name}{$node} = $orientation;
				$node_lens{$node} = 0;
			}

			push(@{$path_nodes{$path_name}}, @nodes);
		}

		elsif ($line =~ /^W\t/) {
			chomp($line);

			my ($rec_type, $sample, $hap_index, $seq_id, $seq_start, $seq_end, $walk) = split(/\t/, $line);

			if ($sample eq '_MINIGRAPH_') {
				next();
			}

			my @walk_comps =  split(/(<|>)/,  $walk);
			my $path_name = "${sample}${chr_delim_str}${seq_id}";

			my $index = 1;

			while ($index < $#walk_comps) {
				my $orientation = '+';

				if ($walk_comps[$index] eq '<') {
					$orientation = '-';
				}

				my $node = $walk_comps[$index + 1];

				$path_node_orientations{$path_name}{$node} = $orientation;
				$node_lens{$node} = 0;
				push(@{$path_nodes{$path_name}}, $node);

				$index += 2;
			}
		}
	}

	close(GFA);


	open(GFA, '<', $gfa_file) or error("can't open $gfa_file: $!");

	while (my $line = <GFA>) {
		if ($line =~ /^S/) {
			chomp($line);

			my ($rec_type, $node, $seq) = split(/\t/, $line);

			if (exists($node_lens{$node}) && defined($seq)) {
				$node_lens{$node} = length($seq);
			}
		}
	}

	close(GFA);


	foreach my $path_name (sort keys %path_nodes) {
		my $path_pos = 0;

		foreach my $index (0..$#{$path_nodes{$path_name}}) {
			my $node = $path_nodes{$path_name}[$index];

			if (! exists($node_lens{$node})) {
				print(STDERR "node: $node length not found\n");

				next();
			}

			my $len = $node_lens{$node};

			if ($len < 1) {
				print(STDERR "node: $node length: $len < 1\n");

				next();
			}

			if (! exists($path_node_orientations{$path_name}{$node})) {
				print(STDERR "path name: $path_name node: $node orientation not found\n");

				next();
			}

			my $orientation = $path_node_orientations{$path_name}{$node};
			my $start = $path_pos + 1;
			my $stop = $path_pos + $len;

			print(STDOUT "$path_name\t$start\t$stop\t$node\t$orientation\n");

			$path_pos = $stop;
		}
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

	GetOptions ('g|gfa=s' => \$gfa_file,
				'd|delim=s' => \$chr_delim_regex,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -verbose => 1);
	}

	if (! defined($gfa_file)) {
		arg_error('gfa file required');
	}

	$chr_delim_str = $chr_delim_regex;
	$chr_delim_str =~ s/\\//g;

	return(0);
}


__END__

=head1 NAME

gfa_nodes_to_linear_coords.pl

=head1 SYNOPSIS

gfa_nodes_to_linear_coords.pl [options] > gfa_node_coords.out

=head1 DESCRIPTION

gfa_nodes_to_linear_coords.pl uses an input GFA file to generate
linear coordinates for each reference and graph node

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -g --gfa       genome gfa file, vg-based (required)

 -d --delim     pattern used to split genotype from chromosome in path names
                  default: '\.'

 -h --help      display help menu

=cut

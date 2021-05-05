#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $pack_table_file;
my $enable_sort = 0;
my $help;

my %nodes = ();

parse_args();
parse_pack_table_file();

exit(0);


sub parse_pack_table_file {
	my $prev_node_id;
	my @covs = ();

	open(PACK, '<', $pack_table_file) or error("can't open pack table file: $!");

	print("node.id\tpos.covs\n");

	while (my $line = <PACK>) {
		chomp($line);

		if ($line =~ /^seq\.pos/) {
			next();
		}

		my ($seq_pos, $node_id, $node_offset, $cov) = split(/\t/, $line);

		if (defined($prev_node_id) && $node_id ne $prev_node_id) {
			if ($enable_sort == 1) {
				$nodes{$prev_node_id} = join(',', @covs);
			}

			else {
				print("$prev_node_id\t", join(',', @covs), "\n");
			}

			@covs = ();
		}

		$prev_node_id = $node_id;
		push(@covs, $cov);
	}

	close(PACK);

	if (defined($prev_node_id) && $#covs >= 0) {
		if ($enable_sort == 1) {
			$nodes{$prev_node_id} = join(',', @covs);
		}

		else {
			print("$prev_node_id\t", join(',', @covs), "\n");
		}
	}


	if ($enable_sort == 1) {
		foreach my $node_id (sort { $a <=> $b } keys %nodes) {
			print("$node_id\t$nodes{$node_id}\n");
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

	GetOptions ('p|pack=s' => \$pack_table_file,
				's|sort' => \$enable_sort,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -verbose => 1);
	}

	if (! defined($pack_table_file)) {
		arg_error('pack table file required');
	}

	return(0);
}


__END__

=head1 NAME

pack_table_to_seg_cov.pl

=head1 SYNOPSIS

pack_table_to_seg_cov.pl [options] > pack.seg.cov
pack_table_to_seg_cov.pl [options] | gzip > pack.seg.cov.gz

=head1 DESCRIPTION

pack_table_to_seg_cov.pl converts vg pack tables to compact segment cov files

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -p --pack    vg pack table file (required)

 -s --sort    sort by node (increases memory consumption and output file
              compressibility)
                default: disabled

 -h --help    display help menu

=cut

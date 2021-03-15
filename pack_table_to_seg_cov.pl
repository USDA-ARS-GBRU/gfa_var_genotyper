#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $pack_table_file;
my $help;

parse_args();
parse_pack_table_file();

exit(0);


sub parse_pack_table_file {
	my $prev_node_id;
	my $seg_len = 0;
	my $seg_cov = 0;

	open(PACK, '<', $pack_table_file) or error("can't open pack table file: $!");

	print("node.id\tlength\tcoverage\n");

	while (my $line = <PACK>) {
		chomp($line);

		if ($line =~ /^seq\.pos/) {
			next();
		}

		my ($seq_pos, $node_id, $node_offset, $cov) = split(/\t/, $line);

		if (defined($prev_node_id) && $node_id ne $prev_node_id) {
			if ($seg_len > 0) {
				my $avg_cov = int(($seg_cov / $seg_len) + 0.5);

				print("$prev_node_id\t$seg_len\t$avg_cov\n");
			}

			$seg_len = 0;
			$seg_cov = 0;
		}

		$seg_cov += $cov;
		$seg_len++;
		$prev_node_id = $node_id;
	}

	close(PACK);

	if (defined($prev_node_id) && $seg_len > 0) {
		my $avg_cov = int(($seg_cov / $seg_len) + 0.5);

		print("$prev_node_id\t$seg_len\t$avg_cov\n");
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

pack_table_to_seg_cov.pl converts vg pack tables to segment (avg) coverage files

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -p --pack    vg pack table file

 -h --help    display help menu

=cut

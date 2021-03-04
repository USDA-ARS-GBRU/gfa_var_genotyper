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
	my %seg_lens = ();
	my %seg_covs = ();

	open(PACK, '<', $pack_table_file) or error("can't open pack table file: $!");

	while (my $line = <PACK>) {
		chomp($line);

		if ($line =~ /^seq\.pos/) {
			next();
		}

		my ($seq_pos, $node_id, $node_offset, $cov) = split(/\t/, $line);

		if (exists($seg_lens{$node_id})) {
			if (($node_offset + 1) > $seg_lens{$node_id}) {
				$seg_lens{$node_id} = $node_offset + 1;
			}
		}

		else {
			$seg_lens{$node_id} = $node_offset + 1;
		}

		$seg_covs{$node_id} += $cov;
	}

	close(PACK);

	print("node_id\tcoverage\n");

	foreach my $node_id (sort { $a <=> $b } keys %seg_covs) {
		my $len = $seg_lens{$node_id};
		my $cov = $seg_covs{$node_id};
		my $avg_cov = int(($cov / $len) + 0.5);

		print("$node_id\t$avg_cov\n");
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

pack_table_to_seg_cov.pl [options]

=head1 DESCRIPTION

pack_table_to_seg_cov.pl converts vg pack tables to segment (avg) coverage files

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -p --pack      vg pack table file

 -h --help    display help menu

=cut

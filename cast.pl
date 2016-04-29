#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my $version = "0.1";
my %options = ();

getopts("a:b:th" , \%options);

my $qry = $options{a};                      # Query sequence file in FASTA format
my $ref = $options{b};                      # Reference genome sequence in FASTA format
my $dir = "./";
my $pam = "NGG";
my $target_length = 23;
my $pam_length = length $pam;
my $max_target_per_seq = 5;
my $target_file = "target.fa";

# Usage
&usage() if $options{h};

#-------------------------------------------------------------------------------
# Demo
#-------------------------------------------------------------------------------
if ($options{t}) {
	$qry = 'demo/cds.fa';
	$ref = 'demo/genome.fa';
}

#-------------------------------------------------------------------------------
# main
#-------------------------------------------------------------------------------
my $cds = load_fasta($qry);
my $genome = load_fasta($ref);
my %target = ();

cas_design($cds, $genome);


#-------------------------------------------------------------------------------
# Subroutines
#-------------------------------------------------------------------------------
sub cas_design {
	if ($pam =~ /^N/) {
		$pam = substr($pam, 1, $pam_length-1);
		$target_length++;
	}
	foreach my $seq_id (sort keys $cds) {
		find_pam($seq_id);
	}
	open (TARGET, ">$target_file") or die "Cannot open $target_file: $!\n";
	foreach my $tar_id (sort keys %target) {
		print TARGET ">$tar_id\n$target{$tar_id}\n";
	}
	close TARGET;
	run_cas9off($target_file, $ref);
}

sub run_cas9off {
	my $target_file = shift;
	my $ref = shift;
	system("$dir/sub/cas9off/offscan.pl -q $target_file -d $ref -a $dir/sub/cas9off/bin/seqmap");
}

sub find_pam {
	my $seq_id = shift;
	my $seq = $cds->{$seq_id};
	my $seq_len = length $seq;
	my $len = $target_length - $pam_length;
	my $i = 0;
	while ($seq =~ /(\S{$len}$pam)/) {
		$target{$seq_id . "." . $i} = $1;
		$i++;
		return 1 if $i > $max_target_per_seq;
	}
	return $i;
}

sub load_fasta {
	my $file = shift;
	my %fasta = ();
	open (FASTA, $file) or die "Cannot open fasta $file: $!\n";
	my $seq_id = undef;
	while (<FASTA>) {
		chomp;
		if (/^>(\S+)/) {
			$seq_id = $1;
		} else {
			s/\s+//g;
			$fasta{$seq_id} .= $_;
		}
	}
	close FASTA;
	return \%fasta;
}

sub usage {
	print <<USAGE;

cast $version

Usage:

    cast.pl <option>
    
    -h  Help
    -t  Test demo
    -a  Query sequence file, in FASTA format
    -b  Genome sequence file, in FASTA format

USAGE
	exit;
}

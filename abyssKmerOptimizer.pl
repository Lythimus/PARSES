#!/usr/bin/perl
use strict;
use Parallel::Iterator qw( iterate_as_array );
use File::Basename;

### author: Joseph Coco
# Determine the best size for the kmer to use for ABYSS, the DeNovo
# Assembler, and execute the assembly with that size for all data sets
# provided. Outputs several files containing the input file name prepended
# with kmer.contigs.fasta. Operates in parallel.
# 
# Optimizes ABYSS using the kmer which produces the best N50 score.
#
# NOTE: ABYSS must be within the executing user's PATH.
#
# ARGV0 = Regex to query files paths.
# ARGV1 = Minimum kmer to attempt.
# ARGV2 = Maximum kmer to attempt.
# ARGV3 = Data Type parameter. Either --illumina-quality or --standard-quality
#
###
my @reads = glob $ARGV[0];
my $sequence;
foreach my $seq (@reads){ # perform on each input file
	$sequence = $seq;
	# execute ABYSS with all kmer in parallel
	iterate_as_array( \&abyss, [$ARGV[1]..$ARGV[2]] );
#	foreach my $kmer ($ARGV[1]..$ARGV[2]){
#		print $kmer;
#	}
	# pull out best N50 score and make a copy of the file with best score
	my $dir = dirname($0);
	my @facBestResult = sort qx{ "$dir/fac.pl" $sequence.*.contigs.fa };
	$_ = $facBestResult[-1];
	split ('\t');
	chomp $_[1];
	`cp $_[1] $_[1].kmerOptimized.fa`;
}

sub abyss{
	my $kmer = shift;
	$kmer = $kmer+7;
	`ABYSS -k $kmer $ARGV[3] --out=$sequence.$kmer.kmer.contigs.fa $sequence --coverage-hist=$sequence.$kmer.kmer.contigs.coverage`;
}

#!/usr/bin/perl -w
###################################################################################################################
# This script takes files of mapped sqeuence reads in novo format as input.  For each input file, it produces a   #
# file as output with nonmapped (NM) reads and a counts file.                                                     #
###################################################################################################################

unless( $#ARGV >= 0 )
{
	print STDERR "usage: $0 mappedseqfile1 ...\n";
	exit;
}

for( $f = 0; $f <= $#ARGV; $f++ )
{
	# Open the input and output files
	open( INF, "<$ARGV[$f]" );
	open( NM, ">$ARGV[$f].NM.fastq" );
	open( NMA, ">$ARGV[$f].NM.fasta" );

	$tot = $uniq = $rep = $qc = $nm = $ql = 0;
	while( $line = <INF> )
	{
		# Skip any lines that aren't reads
		unless( $line =~ /^@/ )
		{ next; }

		++$tot;
		# Parse the line
		chomp($line);
		($head, $slr, $seq, $rdqual, $status, $alscore, $alqual, $chr, $pos, $strand) = split("\t", $line);
		$slr = $slr; $alscore = $alscore; $alqual = $alqual; $chr = $chr; $pos = $pos; $strand = $strand;

		if( $status eq "U" )
		{ ++$uniq; }
		elsif( $status eq "R" )
		{ ++$rep; }
		elsif( $status eq "QC" )
		{ ++$qc; }
		elsif( $status eq "NM" )
		{
			++$nm;
			print NM "$head\n$seq\n+\n$rdqual\n";
			$head =~ s/^@/>/;
			print NMA "$head\n$seq\n";
		}
		elsif( $status eq "QL" )
		{ ++$ql; }
		else
		{ print STDERR "Unknown status $status. Skipped this read.\n" }
	}

	close(INF);
	close(NM);
	close(NMA);

	open( COU, ">$ARGV[$f].counts" );
	print COU "Total Reads: $tot\n";
	print COU "Unique Maps: $uniq\n";
	print COU "Repeat Maps: $rep\n";
	print COU "No Mappings: $nm\n";
	print COU "Quality Con: $qc\n";
	print COU "Low Scoring: $ql\n";

	close(COU);
}


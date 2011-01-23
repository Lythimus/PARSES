#!/usr/bin/perl -w
###################################################################################################################
# This script takes fasta files of nonmapped reads and corresponding sam files of span reads from tophat and      #
# filters the span reads out of the nonmapped reads.  It prints a new fasta file with the span reads removed.     #
# There must be an even number of files specified on the command line and the respective paired files must match. #
# This is not checked by the script, the user should inspect the details file to ensure the files matched.        #
# For example: Filterspans.pl a.NM.fasta b.NM.fasta c.NM.fasta a.sam.spans b.sam.spans c.sam.spans                #
###################################################################################################################

# First check to make sure the command line parameters are correct; if not print usage to stderr and exit
unless (($#ARGV % 2) == 1)
{
	print STDERR "You must specify an even number of files on the command line.\n";
	print STDERR "usage: $0 file1.NM.fasta [file2.NM.fasta...] file1.sam.spans [file2.sam.spans...]\n";
	exit;
}

for( $fi = 0; $fi <= ($#ARGV / 2); $fi++ )
{
	$fastan = $ARGV[$fi];
	$spann = $ARGV[$fi + ($#ARGV / 2) + 1];

	unless( -r $fastan )
	{
		print STDERR "File $fastan is not readable.  Skipped this pair.\n";
		next;
	}
	unless( -r $spann )
	{
		print STDERR "File $spann is not readable.  Skipped this pair.\n";
		next;
	}
	unless( $fastan =~ /.NM.fasta$/ && $spann =~ /.sam.spans$/ )
	{
		print STDERR "Filename extensions are not in proper format.  Skipped this pair.\n";
		next;
	}

	print STDERR "Filtering spans in $spann from $fastan.\n";

	# Read the span reads into a hash
	%spanr = ();
	open( SPAN, "<$spann" );
	$spanreads = 0;
	while( $line = <SPAN> )
	{
		++$spanreads;
		($rname) = split " ", $line;
		$rname = ">" . $rname . "/1";
		if( defined $spanr{$rname} )
		{
			++$spanr{$rname};
			#print STDERR "FYI: Found a duplicated read name in the spans file: $rname\n";
		}
		else
		{ $spanr{$rname} = 1; }
	}
	close( SPAN );

	# Read the nonmapped reads from fasta and filter out those that are in the spanreads
	open( NM, "<$fastan" );
	open( NS, ">$fastan.nospans" );
	open( DET, ">$fastan.nospans.details" );
	print DET "Filtering spans in $spann from $fastan. Make sure these files match!\n";
	$nonmappedreads = 0;
	$filtered = 0;
	$kept = 0;
	$line = <NM>;
	chomp($line);
	$notdone = 1;
	do
	{
		$head = $line;

		# Read up the sequence
		$seq = <NM>;
		chomp($seq);
		do
		{
			if( $line = <NM> )
			{
				chomp( $line );
				unless( $line =~ /^>/ )
				{ $seq .= $line; }
			}
			else
			{
				$notdone = 0;
			}
		} until ( ($notdone == 0) || ($line =~ /^>/) );

		++$nonmappedreads;
		if( defined $spanr{$head} )
		{	# This read was one of the spans so filter it
			++$filtered;		
		}
		else
		{	# This read was not a span so keep it
			++$kept;
			print NS "$head\n$seq\n";
		}
	} while( $notdone );
	close( NM );
	close( NS );
	print DET "There were $spanreads reads in the span file.\n";
	print DET "There were $nonmappedreads reads in the nonmapped fasta file.\n";
	print DET "I filtered $filtered of these reads and kept $kept of these reads.\n";
	close( DET );
}

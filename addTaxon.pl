#!/usr/bin/perl -w
###################################################################################################################

# First collect the GI numbers from the blast hit file.  The taxonomy file is huge, so we only want to build the
# hash for those entries that appear in the blast hit file.
%GtoT = ();
open( BHF, "<$ARGV[1]" );
print STDERR "\nCollecting GI numbers from blast hit file...";
while( $hit = <BHF> )
{
	next if ($hit =~ /^#/);
	chomp($hit);
	($readid, $subjid) = split "\t", $hit;
	$readid = $readid;
	($gi, $gid) = split '\|', $subjid;
	$gi = $gi;
	$GtoT{$gid} = -1;
}
close(BHF);
print STDERR "done.\n";

# Now process the huge gi to taxid translation file, keeping translations for gis that we saw in the blast hit file
open( GTN, "<$ARGV[0]" );
print STDERR "\nBuilding gi to taxid translation hash...";
while( $line = <GTN> )
{
	chomp($line);
	($gin, $tin) = split " ", $line;
	if( defined $GtoT{$gin} )
	{
		if( $GtoT{$gin} == -1 )
		{
			$GtoT{$gin} = $tin;
		}
		elsif( $GtoT{$gin} != $tin )
		{ print STDERR "\nWARNING: GIN $gin has more than one TIN associated with it!\n"; }
	}
}
close(GTN);
print STDERR "done.\n";

# Now build the fasta hash
%fasta = ();
open( FASTA, "<$ARGV[2]" );
print STDERR "\nBuilding fasta hash...";
@mfa = <FASTA>;
close(FASTA);
for( $mi = 0; $mi <= $#mfa; $mi++ )
{
	$mine = $mfa[$mi];
	if( $mine =~ /^>/ )	# Check to see if it is a header line
	{
		chomp( $mine );
		$seq = $mfa[++$mi];
		chomp( $seq );
		while( ($mi < $#mfa) && ($mfa[$mi+1] =~ /^[-A-Za-z]/) )
                {       # Continue grabbing lines until the next non-sequence line
                        $tmine = $mfa[++$mi];
                        chomp($tmine);
                        $seq .= $tmine; # append this line to the end of the sequence
                }
                # Put this sequence into the hash
		$mine =~ s/>//;	# strip the > symbol off the sequence name
		$fasta{$mine} = $seq;
        }
}
print STDERR "done.\n";

open( TAXID, ">$ARGV[2].blast" );
open( BHF, "<$ARGV[1]" );
while( $hit = <BHF> )
{
	next if ($hit =~ /^#/);
	chomp($hit);
	($readid, $subjid) = split "\t", $hit;
	$readid = $readid;
	($gi, $gid) = split '\|', $subjid;
	$gi = $gi;
	print TAXID "$hit\t$GtoT{$gid}\n";
}
close(BHF);
foreach $key ( keys %fasta )
{
	print TAXID "$key	gi|000000|gb|AY000000.0|	0.00	0	0	0	0	0	0	0	0	0	-1\n";
}

close(TAXID);
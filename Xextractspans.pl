#!/usr/bin/perl -w

foreach $fn (@ARGV)
{
	unless( -r $fn )
	{
		print STDERR "\nFile $fn is not readable.  Skipping this one.\n";
		next;
	}

	open( INF, "<$fn" );
	open( OUF, ">$fn.spans" );
	open( NOR, ">$fn.mapped" );
	while( $line = <INF> )
	{
		if( $line =~ /^@/ )
		{ next; }
		@cols = split " ", $line;
		if( defined $cols[5] && length($cols[5]) > 3 )
		{ print OUF $line; }
		else
		{ print NOR $line; }
	}
	close(INF);
	close(OUF);
	close(NOR);
}


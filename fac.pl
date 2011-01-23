#!/usr/bin/perl
use strict;
use Getopt::Std qw'getopts';

my %opt;
getopts 'hHjt:', \%opt;
my $opt_threshold = defined $opt{'t'} ? $opt{'t'} : 100;
my $opt_filename = $opt{'H'} || (@ARGV > 1 && !$opt{'h'});
my $opt_jira = $opt{'j'};

sub eng($)
{
	my $x = shift;
	return $x if $x < 10000000;
	return substr($x / 1000000, 0, 5) . 'e6';
}

my ($short, $sum, $ambiguous, $any, $other);
my @x;

sub count($$)
{
	my $id = shift;
	my $seq = uc shift;
	#print $id, length $seq;
	#print "\n";
	my $x = $seq =~ tr/ACGT//;
	my $colourspace = $seq =~ tr/0123//;
	die unless $x == 0 || $colourspace == 0;
	$x = $colourspace if $x == 0;
	my $myambiguous = $seq =~ tr/KYSBWRDMHV//;
	my $myany = $seq =~ tr/N//;
	if ($x < $opt_threshold) {
		$short++;
		return;
	}
	$sum += $x;
	$ambiguous += $myambiguous;
	$any += $myany;
	$other += (length $seq) - $x - $myambiguous - $myany;
	push @x, $x;
}

sub fac($)
{
	my $path = shift;
	$short = $sum = 0;
	$ambiguous = $any = $other = 0;
	@x = ();

	my $id;
	my $seq;
	open IN, "<$path" or die "$path: $!\n";
	while (<IN>) {
		chomp;
		if (/^>/) {
			count $id, $seq if defined $id;
			$id = $_;
			$seq = '';
		} else {
			$seq .= $_;
			push @x, $seq;
		}
	}
	count $id, $seq if defined $id;
	close IN;

	my $n = @x;
	if ($n > 0) {
		my $mean = int $sum / $n;

		@x = sort { $a <=> $b } @x;
		my $min = $x[0];
		#my $q1 = $x[@x/4];
		my $q2 = $x[@x/2];
		#my $q3 = $x[@x*3/4];
		my $max = $x[-1];

		my ($nn50, $n50);
		my $n50sum = 0;
		while ($n50sum < $sum / 2) {
			$nn50++;
			$n50 = pop @x;
			$n50sum += $n50;
		}
		#my $np = int 100 * $n50sum / $sum;

		my $ntotal = eng $short + $n;
		my $sumeng = eng $sum;
print "$n50\t$path\n";
=pod
		print "$ntotal$," if $opt_threshold > 0;
		print $n, $nn50, $min, $q2, $mean, $n50, $max;
		#printf "$,%#.3g", $sum;
		print $, . eng($sum);
		print "$,ambig=$ambiguous" if $ambiguous > 0;
		print "$,any=$any" if $any > 0;
		print "$,other=$other" if $other > 0;
		print "$,$path" if $opt_filename;
		print "\n";
#=cut


		format Spaces =
@<<<<<<<@<<<<<<@<<<<<<@<<<<<<@<<<<<<@<<<<<<@<<<<<<@<<<<<<@<<<<<< ^<<<<<<<<<<<<<
$ntotal, $n, $nn50, $min, $q2, $mean, $n50, $max, $sumeng, $path
                                                                 ^<<<<<<<<<<<<<~~
$path
.
		format Pipes =
|@<<<<<<<|@<<<<<|@<<<<<|@<<<<<|@<<<<<|@<<<<<|@<<<<<|@<<<<<|@<<<<<<|@*|
$ntotal, $n, $nn50, $min, $q2, $mean, $n50, $max, $sumeng, $path
.
		$~ = $opt_jira ? 'Pipes' : 'Spaces';
		$^ = $opt_jira ? 'Pipes_TOP' : 'Spaces_TOP';
		$: = '/- ';
		write;
#=pod
		print "$,ambig=$ambiguous" if $ambiguous > 0;
		print "$,any=$any" if $any > 0;
		print "$,other=$other" if $other > 0;
		print "$,$path" if $opt_filename;
		print "\n";
=cut
	} else {
		print STDERR "warning: `$path' is empty\n";
	}
}

format Spaces_TOP =
n       n:@<<< n:N50  min    median   mean   N50    max    sum
$opt_threshold
.
format Pipes_TOP =
||n     ||n:@<<||n:N50||min  ||median  ||mean||N50  ||max  ||sum   ||
$opt_threshold
.

=pod
$, = "\t";
if ($opt_threshold > 0) {
	print 'n', "n:$opt_threshold$,";
} else {
	print "n$,";
}
print 'n:N50', 'min', 'median', 'mean', 'N50', 'max', "sum\n";
=cut

@ARGV = ('-') if @ARGV == 0;
fac $_ foreach @ARGV;

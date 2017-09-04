use strict;
use Switch;
use List::Util qw(sum);
use POSIX;
use Data::Dumper;
use List::MoreUtils 'any';
use Storable qw(dclone);

my $dataDrive = "/data/ohler/Dina/Drosophila/Research/data/";
my $resultDrive = "/data/ohler/Dina/Drosophila/Research/results/";

my @commands = @ARGV;
my $bodyFile = $commands[0];
my $TSSFile = $commands[1];
my $outFile = $commands[2];

my %body;
my %tss;

open IN, "<$bodyFile" or die "Can't open file:$bodyFile";
my $header1 = <IN>;
chomp($header1);
while (my $line = <IN>)
{
	chomp ($line);
	my ($gene, @rest) = split(/\s+/,$line);
	$body{$gene} = substr($line, length($gene)+1);
}
close IN;

open IN, "<$TSSFile" or die "Can't open file:$TSSFile";
my $header2 = <IN>;
chomp($header2);
$header2 = substr($header2,4);
while (my $line = <IN>)
{
	chomp ($line);
	my ($gene, @rest) = split(/\s+/,$line);
	$tss{$gene} = substr($line, length($gene)+1);
}
close IN;

open OUT, ">$outFile" or die "Can't open file:$outFile";
print OUT ("$header1\t$header2\n");
my %processed;
foreach my $g (keys %body)
{
	my $str = $body{$g};
	if (exists $tss{$g})
	{
		$str = $str. $tss{$g};
		$processed{$g}=1;
	}
	print OUT ("$g\t$str\n");
}

foreach my $p (keys %processed)
{
	delete $tss{$p};
}
foreach my $g (keys %tss)
{

	my $str = $tss{$g};
	print OUT ("$g\t$str\n");
}
close OUT;

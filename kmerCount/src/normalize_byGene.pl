use strict;
use List::Util qw(sum);
use POSIX;
use Data::Dumper;
use List::MoreUtils 'any';
use Storable qw(dclone);

my $dataDrive = "/data/ohler/Dina/Drosophila/Research/data/";
my $resultDrive = "/data/ohler/Dina/Drosophila/Research/results/";

my %geneCount;

#	my ($resultDirectory ,$seqFile , $kmerLength )= @_;
my $resultDirectory = $ARGV[0]; #"/data/ohler/Dina/Drosophila/Research/results/BaseLine_assign_DHS_to_closest_gene/ClustersDHS/" ;
my $seqFile =  $ARGV[1];
my $kmerLength =  $ARGV[2];

my $tabFile = $resultDirectory.$kmerLength."_mer/".$seqFile."_count.tab" ;



open IN, "<$tabFile" or die "Can not open file :$tabFile";
my %geneFeatures;
my %geneSeqLength;
my $kmerList = <IN>;
my %geneNoOfDHS;
while (my $line = <IN> )
{
	chomp($line);
	my @features = split(/\s+/,$line);
	my $dhs = $features[0];
	my @features2 = split(/\_/, $dhs);
	my $gene = $features2[0];
	my $dhsLength = $features2[@features2-1];
	#my $dhs_id_name = $stage."_".$no."_".$id."_".$idn;
	for(my $i=1; $i<@features;$i++)
	{	
		my ($index,$count) = split(/\:/, $features[$i]);
		$geneFeatures{$gene}{$index} +=$count;
		
	}
	#$dhsLength = $dhs{$dhs_id_name}{'end'} -$dhs{$dhs_id_name}{'start'}  +2;  #assuming that the data is in bed file format, to get the length end - start, but I add 1 extra bp before and after the actual length
	$geneSeqLength{$gene}+=$dhsLength;
	$geneNoOfDHS{$gene}++;
}

my $outputFile = $resultDirectory.$kmerLength."_mer/".$seqFile.".tab" ;
open OUT, ">$outputFile" or die "Can not open file :$outputFile";
print OUT ($kmerList);
foreach my $gene (keys %geneFeatures)
{
	print OUT ("$gene\t");
    foreach my $index (sort {$a <=> $b} keys %{ $geneFeatures{$gene} } ) 
    {  
    	
        my $normalizedCount = $geneFeatures{$gene}{$index} / ($geneSeqLength{$gene} - ($kmerLength*$geneNoOfDHS{$gene}) + $geneNoOfDHS{$gene});
		
		print OUT ("$index:$normalizedCount\t");
    }
    print OUT ("\n");
}
	

use strict;
use Switch;
use List::Util qw(sum);
use POSIX;
use Data::Dumper;
use List::MoreUtils 'any';
use Storable qw(dclone);

my $dataDrive = "/data/ohler/Dina/Drosophila/Research/data/";
my $resultDrive = "/data/ohler/Dina/Drosophila/Research/results/";

#my $dataDrive = "../../data/";
my $resultDrive = "/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/ClustersDHS/MC_5order/MC_ubiq_output/";

my %dhs;
for (my $i=12;$i<=39; $i++)
{
	#my $fasta = $resultDrive."MarkovChain_5order/".$i.".pos.fa";
	#my $bed = $resultDrive."MarkovChain_5order/".$i.".pos.bed";
	#my $fasta = $dataDrive."ClustersDHS/".$i."_pos.fa";
	#my $bed = $dataDrive."ClustersDHS/".$i."_pos.bed";
	readDHSs();	
	my $cluster =$i;
	my $fasta = $resultDrive.$cluster."_pos.fa";
	my $bed = $resultDrive.$cluster."_pos.bed";
	convertFastaToBed($fasta,$bed);
	
	my $sortedFile = $resultDrive.$cluster."_pos_sorted.bed";
	my $cmd_sort = "sort -k1,1 -k2,2n $bed >$sortedFile";
	system($cmd_sort);
	
	my $mergedFile = $resultDrive.$cluster."_pos_merged.bed";
	my $cmd_merge = "mergeBed -nms -i $sortedFile > $mergedFile";
	system($cmd_merge);
	
	my $faFileMerged = $resultDrive.$cluster."_pos_merged.fa";
	my $cmd_getFa= "twoBitToFa -bed=$mergedFile -noMask /data/ohler/Dina/Drosophila/Research/data/Genome/Drosophila_R5.2bit  $faFileMerged";
	system($cmd_getFa); 
	
	my $single =$resultDrive.$cluster."_pos_single.fa";
 		
	my $cmd = "perl /data/ohler/Dina/Drosophila/Research/code/scripts/convert_multi_fasta_to_single.pl $faFileMerged > $single";
 		system($cmd);
}

sub readDHSs 
{
	my $dhsFile = $dataDrive."DHS/DHS_all_R5_id.bed";
	open IN, "<$dhsFile" or die "Can't open file:$dhsFile";
	while (my $line = <IN>)
	{
		chomp ($line);
		my ($chr, $start, $end, $id) = split(/\s+/,$line);
		#$chr =~ s/\'//g;
 		#$end=~ s/\'//g;
 		#$start=~ s/\'//g;
 		#$id=~ s/\'//g;
 		
		$dhs{$id}{'chr'} = $chr;
		$dhs{$id}{'start'} = $start-1;
		$dhs{$id}{'end'} = $end+1;		
	} 
}
sub convertFastaToBed
{
	my ($fastaFile, $bedFile)=@_;
	
	open IN, "<$fastaFile" or die "Can't open file:$fastaFile";
 	open OUT, ">$bedFile" or die "Can't open file:$bedFile";
 	my $line;
 

 	while ($line = <IN>)
	{
		chomp($line);
		if($line =~ /^>/)
		{
			my ($gene, $stage,$no,$id,$idn,$time) = split(/\_/, $line);
			my $dhs_id = $stage."_".$no."_".$id."_".$idn;
			
			$gene=~ s/\>//g;
			my $chr  = $dhs{$dhs_id}{'chr'};
			my $start = $dhs{$dhs_id}{'start'};
			my $end = $dhs{$dhs_id}{'end'};	
			my $name = $dhs_id."_".$chr."_".$start."_".$end;
			print OUT ("$chr\t$start\t$end\t$gene\n"); 					
		}
	}
		
}

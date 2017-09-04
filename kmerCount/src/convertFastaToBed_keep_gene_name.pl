use strict;
use Switch;
use List::Util qw(sum);
use POSIX;
use Data::Dumper;
use List::MoreUtils 'any';
use Storable qw(dclone);

my $dataDrive = "/data/ohler/Dina/Drosophila/Research/data/";
my $resultDrive = "/data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_2order/Overlap_stage_2/";

my %dhs;
readDHSs();	

my @cutoffs = (5,8,10,12,14);
for (my $c = 0; $c < @cutoffs; $c++)
{
	for (my $i=11;$i<=39; $i++)
	{
	
		my $fasta = $resultDrive.$i."_pos_".$cutoffs[$c]."_single.fa";
		my $bed = $resultDrive."MergedDHS/".$i."_pos_".$cutoffs[$c]."_single.bed";	
		convertFastaToBed($fasta,$bed);
		splitBedMerge($bed);
		my $newFasta = $resultDrive."MergedDHS/".$i."_pos_".$cutoffs[$c]."_merged.fa";
		my $newFastaSingle = $resultDrive."MergedDHS/".$i."_pos_".$cutoffs[$c]."_single_merged.fa";
		convertBedToFasta($bed.".merged", $newFasta, $newFastaSingle);

		
	}
}

my $fasta = $resultDrive."ubiquitous_all_single.fa";
my $bed = $resultDrive."MergedDHS/ubiquitous_single.bed";		
convertFastaToBed($fasta,$bed);
splitBedMerge($bed);
my $newFasta = $resultDrive."MergedDHS/ubiquitous_all_merged.fa";
my $newFastaSingle = $resultDrive."MergedDHS/ubiquitous_all_single_merged.fa";
convertBedToFasta($bed.".merged", $newFasta, $newFastaSingle);

sub readDHSs 
{
	my $dhsFile = $dataDrive. "DHS_row_data/Jamm_peaks/allStages.filtered.peaks.idr.0.2.narrowPeaks";
	
 	open IN, "<$dhsFile" or die "Can not open file :$dhsFile";
 
 	while (my $line = <IN> )
 	{
 		chomp($line);
 		my ($chr, $start,$end, $id,@rest) = split(/\s+/, $line);
 		$dhs{$id}{'chr'} = $chr;
		$dhs{$id}{'start'} = $start;
		$dhs{$id}{'end'} = $end;		
 	}
 	close IN;
 
}


sub convertBedToFasta
{
	my ($bedFile, $fastaFile, $fastaSingleFile) = @_;
	my $cmd = "twoBitToFa -bed=$bedFile -noMask ".$dataDrive."/Genome/Drosophila_R5.2bit $fastaFile"; 
 	system($cmd);

 	my $cmd = "perl /data/ohler/Dina/Drosophila/Research/code/scripts/convert_multi_fasta_to_single.pl $fastaFile > $fastaSingleFile";
 	system($cmd);
	
}
sub convertFastaToBed
{
	my ($fastaFile, $bedFile) = @_;
	
	open IN, "<$fastaFile" or print "Can't open file:$fastaFile";
 	my $tempFile = $bedFile."temp";
 	open OUT, ">$tempFile" or print "Can't open file:$tempFile";
 	my $line;
 

 	while ($line = <IN>)
	{
		chomp($line);
		if($line =~ /^>/)
		{
			my ($gene, $chr, $start, $stageW,$stage,$time) = split(/\_/, $line);
			my $dhs_id = $chr."_".$start."_".$stageW."_".$stage;
			
			$gene=~ s/\>//g;
			my $chr  = $dhs{$dhs_id}{'chr'};
			my $start = $dhs{$dhs_id}{'start'};
			my $end = $dhs{$dhs_id}{'end'};	
			my $name = $dhs_id."_".$chr."_".$start."_".$end;
			#print ("$dhs_id\t");
		#	if ($time eq "initial")
			{
				print OUT ("$chr\t$start\t$end\t$gene\n");
			}
		}	
	}
	
	close OUT;
	my $cmd = "sort -k4,4 -k1,1 -k2,2n $tempFile > $bedFile";
	system($cmd);
	
	$cmd = "rm $tempFile";
	system($cmd);	
		
}

sub splitBedMerge
{
	my ($bedFile) = @_;
	my $cmd = "for gene in `cut -f 4 $bedFile | sort | uniq`; do grep -w \$gene $bedFile > ".$resultDrive."MergedDHS/split_results/\$gene.output.bed;  done";
	#print ($cmd."\n");
	system($cmd);
	sleep (0.1);
	$cmd = "for bed in ".$resultDrive."MergedDHS/split_results/*.output.bed; do   mergeBed -nms -i \$bed | cut -d';' -f 1 > \$bed.merged; done";
	#print ($cmd."\n");
	system($cmd);
	sleep (0.1);
	
	$cmd = "cat ".$resultDrive."MergedDHS/split_results/*.output.bed.merged > $bedFile.merged";
	#print ($cmd."\n");
	system ($cmd);

	$cmd = "rm -f ".$resultDrive."MergedDHS/split_results/*";
	system ($cmd);
	sleep (0.1);

}
sub getStrandInfoForGenes
{
	
	my $tssFile = $dataDrive. "Genome/FlyBase_Genes_map_uniq.tab";
	open IN, "<$tssFile" or die "Can not open file :$tssFile";
	my $line;
	my %gene_map;
	while ($line = <IN> )
	{
		chomp($line);
		my ($cgname,$fbname,$name2,$symbol) = split(/\s+/, $line);
	#	if (exists $genes{$gene})
		{
			$gene_map{$fbname} = $cgname;
 		}
	 }
	 
	 my $tssFile = $dataDrive. "Genome/fly_1st_exon.bed";
	open IN, "<$tssFile" or die "Can not open file :$tssFile";
	
	my %gene_strand;
	while ($line = <IN> )
	{
		chomp($line);
		my ($chr,$start,$end,$gene,$score,$strand) = split(/\s+/, $line);
		$gene_strand{$gene} = $strand;
	 }
	 
	my ($fastaFile, $bedFile)=@_;
	open IN, "<$fastaFile" or die "Can't open file:$fastaFile";
 	open OUT, ">$bedFile" or die "Can't open file:$bedFile";
 	my $line;
 

 	while ($line = <IN>)
	{
		chomp($line);
		#if($line =~ /^>/)
		{
			my ($d_chr,$d_start,$d_end,$dhs_id,$g_chr,$g_start,$g_end,$gene,$score,$strand) = split(/\s+/, $line);
			
			#my $gene_cgname = "";
			#if (exists $gene_map{$gene})
			{
			#	$gene_cgname = $gene_map{$gene};
			#	if (exists $gene_strand{$gene_cgname} )
			#	{
			#		my $strand = $gene_strand{$gene_cgname};
					
					my $chr  = $dhs{$dhs_id}{'chr'};
					my $start = $dhs{$dhs_id}{'start'};
					my $end = $dhs{$dhs_id}{'end'};	
					print OUT ("$chr\t$start\t$end\t$dhs_id\t$gene\t$strand\n");
			#	}
			#	else
			#	{
			#		print ("$gene does not exist in the strand\n");
				}
			}
			#else
			#{
		#		print ("$gene does not exist in the map\n");
		#	}
			
		#}	
	}
}


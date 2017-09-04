use strict;
use Switch;
use List::Util qw(sum);
use POSIX;
use Data::Dumper;
use List::MoreUtils 'any';
use Storable qw(dclone);

#CHANGE DHS IDS
# intersectBed -a allStages.filtered.peaks.idr.0.1.narrowPeaks -b ../../../data/Genome/Gene_TSSs_2bp.bed -u > idr.0.1/DHS_overlap_TSS_2bp.bed
# intersectBed -a allStages.filtered.peaks.idr.0.1.narrowPeaks -b ../../../data/Genome/Gene_TSSs_2bp.bed -v > idr.0.1/DHS_no_overlap_TSS_2bp.bed
# intersectBed -a allStages.filtered.peaks.idr.0.2.narrowPeaks -b ../../../data/Genome/Gene_TSSs_2bp.bed -u > idr.0.2/DHS_overlap_TSS_2bp.bed
# intersectBed -a allStages.filtered.peaks.idr.0.2.narrowPeaks -b ../../../data/Genome/Gene_TSSs_2bp.bed -v > idr.0.2/DHS_no_overlap_TSS_2bp.bed

#intersectBed -a DHS_no_overlap_TSS_2bp.bed -b ../../CRM/redfly_all_specified_modified.bed -wb  | sort -k4,4 -k14,14 --uniq > DHS_overlap_redfly_specified_genes.bed

#intersectBed -a DHS_no_overlap_TSS_2bp.bed -b ../../Stark_Data/enhancer_gene_assignments.bed -wb  |  sort -k4,4 -k15,15 --uniq > DHS_overlap_stark.bed

my $dataDrive = "/data/ohler/Dina/Drosophila/Research/data/";
my $resultDrive = "/data/ohler/Dina/Drosophila/Research/results_idr/";
my $codeVersion = "Initialization_RedFly_uniqGenes";

my %flymap;
mapGeneIDs();

$codeVersion = "Initialization_RedFly_uniqGenes";
splitDHSIntoClusters_specifiedGenes_redfly();
#concatenate_ubiquitous();
#getUnlabeled();  #needs to be changed
getmidpointDHS("/data/ohler/Dina/Drosophila/Research/data/DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed","/data/ohler/Dina/Drosophila/Research/data/DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp_midpoints.bed");

$codeVersion="Initialization_Stark_uniqGenes";
splitDHSIntoClusters();
#concatenate_ubiquitous();
#getUnlabeled(); 

##########initialize stark and redfly
$codeVersion = "Initialization_Stark_RedFly_uniqGenes";
concatenate_stark_redfly();
#concatenate_ubiquitous_stark_redfly();
#getUnlabeled_stark_redfly();

################TSS
#splitDHSIntoClusters_TSS();
#concatenate_ubiquitous_TSS();

############second stage MC
#These two functions are to get ubiquitous for all genes not only unique ones
#splitDHSIntoClusters_stark_ubiquitous_all_2nd_stage();
#splitDHSIntoClusters_redfly_specifiedGenes_ubiquitous_all_2nd_stage();
#getClosestPerClusterForUnspecifiedGenes_ubiquitous_all_2nd_stage(); --not used at all
#concatenate_ubiquitous_2nd_stage();

#getUnlabeled_2nd_stage();
################TSS
#splitDHSIntoClusters_TSS_all_genes();
#concatenate_ubiquitous_TSS_all_genes();

sub mapGeneIDs
{
	my $mapfile = $dataDrive."Genome/Gene_map_new_filtered_uniq.tab";
	open IN, "<$mapfile" or die "Can not open file :$mapfile";
	my $line = <IN> ;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($ensemble, $flybaseID, $chr, $start, $end, $strand) = split(/\s+/, $line);
 		if (exists $flymap{$ensemble} && ($flymap{$ensemble} ne $flybaseID))
 		{
 			print ("something is wrong \n");
 			
 		}
 		else
 		{
 			$flymap{$ensemble} =  $flybaseID;
 		}
 	}	
 	close IN;

}
sub getmidpointDHS
 {
 	#my $bedFile = "../../data/DHS/DHS_nooverlap_Redfly.bed";
 	#my $midFile = "../../data/DHS/DHS_nooverlap_Redfly_midpoints.bed";
 	my ($bedFile,$midFile) = @_;
 	open IN, "<$bedFile" or die "Can not open file :$bedFile";
 	open OUT, ">$midFile" or die "Can't open file : $midFile";
 	
 	my %dhs;
 	my $count = 0;
 	while (my $line = <IN> )
 	{
 			chomp($line);
 		my ($chr,$start,$end,$id) = split(/\s+/, $line);
 		my $mid = int(($start+$end) /2);
 		my $mid_end = $mid +1;
 		print OUT ("$chr\t$mid\t$mid_end\t$id\n");
 	}
 	close IN;
 	close OUT;		
 }
 
 
###initialize redfly
#intersectBed -a DHS_no_overlap_TSS_2bp.bed -b ../../CRM/redfly_all_specified_modified.bed -wb  | sort -k4,4 -k14,14 --uniq > DHS_overlap_redfly_specified_genes.bed
$codeVersion="Initialization_RedFly_uniqGenes";
sub splitDHSIntoClusters_specifiedGenes_redfly
{
 	my $DHSFile = $dataDrive. "DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";
	my %dhs;
	 
	open IN, "<$DHSFile" or die "Can't open file:$DHSFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id, @rest) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
 	close IN;
 	
	my $dhs_file = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_redfly_specified_genes.bed";
	open IN, "<$dhs_file" or die "Can not open file :$dhs_file";
	my $line ;
	my %gene_DHS;
	my %gene_start;
	my %gene_end;
	my %gene_chr;
	
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($chr,$dhs_start,$dhs_end,$dhs_id,$f1,$f2,$f3,$f4,$f5,$f6,$chr_r,$start_r,$end_r,$gene) = split(/\s+/, $line);
		push(@{$gene_DHS{$gene}},$dhs_id);
		push(@{$gene_start{$gene}},$dhs{$dhs_id}{'start'});
		push(@{$gene_end{$gene}},$dhs{$dhs_id}{'end'});
		push(@{$gene_chr{$gene}},$dhs{$dhs_id}{'chr'});
		
 	}
	close IN;
	
	
	
	my %geneCount;
	my %dhsCount;
	my %totalGeneCount;
	#my $all_DHS_Labeled = $resultDrive.$codeVersion."/ClustersDHS/redfly_dhs_specified_genes_used_in_initialization.bed";
	#open ALL, ">$all_DHS_Labeled" or die "Can not open file :$all_DHS_Labeled";
	
 	for (my $i=1; $i<=39; $i++)
 	{
 		my $clusterFile = $dataDrive."/geneClusters/$i.gene.cluster.uniq";
 		open IN, "<$clusterFile" or die "Can not open file :$clusterFile";
 		
 		my $clusterFile_fa = $resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize.bed";
 		open OUT, ">$clusterFile_fa" or die "Can not open file :$clusterFile_fa";
 		
 		while ($line = <IN> )
 		{
 			$totalGeneCount{$i}++;
 			chomp($line);
 			my ($cg_name,$gene) = split(/\s+/,$line);
 		
 				
 				if(exists($gene_DHS{$gene}))
	 			{
	 				
	 				for(my $j=0; $j<@{$gene_DHS{$gene}}; $j++)
	 				{
	 					my $s = $gene_start{$gene}[$j];
	 					my $e = $gene_end{$gene}[$j] ;
	 					print OUT ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene"."_".$gene_DHS{$gene}[$j]."_initial\n");
	 			#		print ALL ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene_DHS{$gene}[$j]\tcluster_$i\t$gene\n");
	 					$dhsCount{$i}{$gene_DHS{$gene}[$j]}=1;
	 					
	 				}
	 				$geneCount{$i}{$gene}++;
	 			}
 			
 		}
 		close OUT;
 		close IN;

		print ("cluster $i\n");
 		my $uniq_file =  $resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_uniq.bed";
 		my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq $clusterFile_fa > $uniq_file";
 		system($cmd);	

		my $fa_file = $resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_uniq.fa";
		my $cmd = "twoBitToFa -bed=$uniq_file -noMask ../../data/Genome/Drosophila_R5.2bit $fa_file"; 
 		system($cmd);
 		
 		my $single =$resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_uniq_single.fa";
 		my $cmd = "perl /data/ohler/Dina/Drosophila/Research/code/scripts/convert_multi_fasta_to_single.pl $fa_file > $single";
 		system($cmd);
 	
 	
 	}
 	##close ALL;
 	
 	print ("cluster#\t#genes_per_cluster\t#genes_per_cluster_with_DHS_overlap_redfly\t#DHSperCluster\n");
 	for (my $i=1;$i<=39; $i++)
	{
		my $c = keys %{$geneCount{$i}};
		my $dhs_c = keys %{$dhsCount{$i}};
		print ("cluster $i\t$totalGeneCount{$i}\t$c\t$dhs_c\n");
	}
	
		
 
}
###############Initialization Stark
#intersectBed -a DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed -b Stark_Data/enhancer_gene_assignments.bed -wb  |  sort -k4,4 -k14,14 --uniq > DHS_row_data/Jamm_peaks/DHS_overlap_stark.bed

sub splitDHSIntoClusters
{ 	
 	my $DHSFile = $dataDrive. "DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";
	my %dhs;
	 
	open IN, "<$DHSFile" or die "Can't open file:$DHSFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id, @rest) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
 	close IN;
 	
	my $dhs_file = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_stark.bed";
	open IN, "<$dhs_file" or die "Can not open file :$dhs_file";
	my $line ;
	my %gene_DHS;
	my %gene_start;
	my %gene_end;
	my %gene_chr;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($chr,$dhs_start,$dhs_end,$dhs_id,$f1,$f2,$f3,$f4,$f5,$f6,$chr_r,$start_r,$end_r,$enhancer_id,$gene,@rest) = split(/\s+/, $line);
		push(@{$gene_DHS{$gene}},$dhs_id);
		push(@{$gene_start{$gene}},$dhs{$dhs_id}{'start'});
		push(@{$gene_end{$gene}},$dhs{$dhs_id}{'end'});
		push(@{$gene_chr{$gene}},$dhs{$dhs_id}{'chr'});
 	}
	close IN;
	
	
	
	my %geneCount;
	my %dhsCount;
	my %totalGeneCount;
	#my $all_DHS_Labeled = $resultDrive.$codeVersion."/ClustersDHS/stark_specified_genes_used_in_initialization.bed";
	#open ALL, ">$all_DHS_Labeled" or die "Can not open file :$all_DHS_Labeled";
	
 	for (my $i=1; $i<=39; $i++)
 	{
 		my $clusterFile = $dataDrive."/geneClusters/$i.gene.cluster.uniq";
 		open IN, "<$clusterFile" or die "Can not open file :$clusterFile";
 		
 		my $clusterFile_fa = $resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize.bed";
 		open OUT, ">$clusterFile_fa" or die "Can not open file :$clusterFile_fa";
 		
 		while ($line = <IN> )
 		{
 			$totalGeneCount{$i}++;
 			chomp($line);
 			my ($cg_name,$gene) = split(/\s+/,$line);
 			
 			
 			
 				if(exists($gene_DHS{$gene}))
	 			{
	 				
	 				for(my $j=0; $j<@{$gene_DHS{$gene}}; $j++)
	 				{
	 					my $s = $gene_start{$gene}[$j];
	 					my $e = $gene_end{$gene}[$j] ;
	 					print OUT ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene"."_".$gene_DHS{$gene}[$j]."_initial\n");
	 	#				print ALL ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene_DHS{$gene}[$j]\tcluster_$i\t$gene\n");
	 					$dhsCount{$i}{$gene_DHS{$gene}[$j]}=1;
	 					
	 				}
	 				$geneCount{$i}{$gene}++;
	 			}
 		
 		}
 		close OUT;
 		close IN;
 	
 		print ("cluster $i\n");
 		my $uniq_file =  $resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_uniq.bed";
 		my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq $clusterFile_fa > $uniq_file";
 		system($cmd);	

		my $fa_file = $resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_uniq.fa";
		my $cmd = "twoBitToFa -bed=$uniq_file -noMask ../../data/Genome/Drosophila_R5.2bit $fa_file"; 
 		system($cmd);
 		
 		my $single =$resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_uniq_single.fa";
 		my $cmd = "perl /data/ohler/Dina/Drosophila/Research/code/scripts/convert_multi_fasta_to_single.pl $fa_file > $single";
 		system($cmd);
 	
 	}
 	close ALL;
 	
 	print ("cluster#\t#genes_per_cluster\t#genes_per_cluster_with_DHS_overlap_redfly\t#DHSperCluster\n");
 	for (my $i=1;$i<=39; $i++)
	{
		my $c = keys %{$geneCount{$i}};
		my $dhs_c = keys %{$dhsCount{$i}};
		print ("cluster $i\t$totalGeneCount{$i}\t$c\t$dhs_c\n");
	}
	
		
 
}

sub concatenate_ubiquitous
{
	my $output = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous.bed";
	my $cmd = "cat ";
	for (my $j=1; $j<= 10; $j++)
	{
		
		$cmd = $cmd. $resultDrive.$codeVersion."/ClustersDHS/".$j."_initialize_uniq.bed ";		
	}
	$cmd = $cmd . "> ".$output;
	print $cmd;
	print ("\n");
	system($cmd);
	
	my $ubiq_sorted = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous_uniq.bed";
	my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq $output  >$ubiq_sorted ";
	#print ($cmd);
	system($cmd);
		
		my $ubiq_fa = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous.fa";
		my $cmd = "twoBitToFa -bed=$ubiq_sorted -noMask ../../data/Genome/Drosophila_R5.2bit $ubiq_fa";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/ubiquitous_single.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $ubiq_fa > $single";
 		system($cmd);
	
}

sub getUnlabeled
{
		
	#get DHS used in initialization
	my %DHS_Initial;
	for (my $j=1; $j<= 39; $j++)
	{
		my $initializeFile =  $resultDrive.$codeVersion."/ClustersDHS/".$j."_initialize_uniq.bed";#_initialize_redfly_uniq.bed ";
		open IN, "<$initializeFile" or die "Can't open file:$initializeFile";	
		my $line;
		while ($line = <IN>)
		{
		
			chomp($line);
			my ($chr,$start,$end,$dhs_id, @rest) = split(/\s+/, $line);
			my ($gene,$chr,$id,$stageW,$stage,$state ) = split(/\_/,$dhs_id);
			my $theID = $chr."_".$id."_".$stageW."_".$stage;
			$DHS_Initial{$theID}++;
		}
		close IN;	
	}
	my $k = keys %DHS_Initial;
	print ("number of DHS in initialization = $k\n");
	my $all_DHS_overlap =  $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";
	my $unlabeledFile = $resultDrive.$codeVersion."/ClustersDHS/DHS_unlabeled.bed";
	open IN, "<$all_DHS_overlap" or die "Can't open file:$all_DHS_overlap";	
	open OUT, ">$unlabeledFile" or die "Can't open file:$unlabeledFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id, @rest) = split(/\s+/, $line);
		if(exists $DHS_Initial{$dhs_id})
		{
			#print ("DHS $dhs_id does  exist\n");
		}
		else{
			print OUT  ("$line\n");
		}
	}
	close IN;
	close OUT;
	
	my $midDHSFile =  $resultDrive. $codeVersion."/ClustersDHS/DHS_unlabeled_midpoints.bed";
	getmidpointDHS($unlabeledFile,$midDHSFile);
	
	my %dhs;
	 
	open IN, "<$unlabeledFile" or die "Can't open file:$unlabeledFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
	
	#Report all genes that are within 50000 bp upstream or downstream of CNVs.
	#bedtools window -a CNVs.bed -b genes.bed -w 10000
	for (my $i=1;$i<=39; $i++)
	{
		my $clusterFile = $dataDrive . "geneClusters/".$i."_cluster_TSS_uniq";
		my $output = $resultDrive. $codeVersion."/ClustersDHS/".$i."_unlabeled_window";
		my $cmd = "windowBed -a $clusterFile -b $midDHSFile -w 50000 > $output";
		print $cmd;
		system($cmd);
	}

	for (my $i=1;$i<=39; $i++)
	{
		my $prev_gene="";
		my $index = 1;
		
		my $window = $resultDrive. $codeVersion."/ClustersDHS/".$i."_unlabeled_window";
		my $bedIn = $resultDrive. $codeVersion."/ClustersDHS/".$i."_unlabeled_in.bed";
		open IN, "<$window" or die "Can't open file:$window";
		open OUTL, ">$bedIn" or die "Can't open file:$bedIn";
		my $line;
		
		while ($line = <IN>)
		{
			chomp($line);
			my ($chr,$start,$end,$genegc,$geneflybase,$strand,$chr_dhs,$start_dhs,$end_dhs,$dhs_id) = split(/\s+/, $line);
			my $s = $dhs{$dhs_id}{'start'};
			my $e = $dhs{$dhs_id}{'end'};
			print OUTL ("$chr_dhs\t$s\t$e\t$geneflybase"."_".$dhs_id."_unlabeled\n");
			$index++;
		}
		close OUT;
		
		
		my $unlabeledFileInUniq = $resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled_in_uniq.bed";
		my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq $bedIn  >$unlabeledFileInUniq ";
		#print ($cmd);
		system($cmd);
		
		my $unlabeled_out_file = $resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled.fa";
		my $cmd = "twoBitToFa -bed=$unlabeledFileInUniq -noMask ../../data/Genome/Drosophila_R5.2bit $unlabeled_out_file";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled_single.fa";
 		my $multi =$resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $multi > $single";
 		system($cmd);
	}	
	
}


########Initialization Stark and Redfly
sub concatenate_stark_redfly
{
	for (my $i=1; $i<=39; $i++)
	{
		my $output = $resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize.bed";
		my $cmd = "cat ".$resultDrive."Initialization_Stark_uniqGenes/ClustersDHS/".$i."_initialize_uniq.bed ".$resultDrive."Initialization_RedFly_uniqGenes/ClustersDHS/".$i."_initialize_uniq.bed  > ".$output;
		#my $cmd = "cat ".$resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_redfly.bed ".$resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_stark.bed > ".$output;
		
		print ("$cmd\n");
		system($cmd);
		
		
		my $ubiq_sorted = $resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_uniq.bed";
		my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq $output  >$ubiq_sorted ";
		system($cmd);
		
		my $ubiq_fa = $resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize.fa";
		my $cmd = "twoBitToFa -bed=$ubiq_sorted -noMask ../../data/Genome/Drosophila_R5.2bit $ubiq_fa";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_single.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $ubiq_fa > $single";
 		system($cmd);
		
	}
}

sub concatenate_ubiquitous_stark_redfly
{
	my $output = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous.bed";
	my $cmd = "cat ";
	for (my $j=1; $j<= 10; $j++)
	{
		
		$cmd = $cmd. $resultDrive.$codeVersion."/ClustersDHS/".$j."_initialize_uniq.bed ";
	}
	$cmd = $cmd . "> ".$output;
	print $cmd;
	print ("\n");
	system($cmd);
	
	my $ubiq_sorted = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous_uniq.bed";
	my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq $output  >$ubiq_sorted ";
	#print ($cmd);
	system($cmd);
		
		my $ubiq_fa = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous.fa";
		my $cmd = "twoBitToFa -bed=$ubiq_sorted -noMask ../../data/Genome/Drosophila_R5.2bit $ubiq_fa";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/ubiquitous_single.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $ubiq_fa > $single";
 		system($cmd);
	
}

sub getUnlabeled_stark_redfly
{

	my $unlabeledFile = "/data/ohler/Dina/Drosophila/Research/data/DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";
	my $midDHSFile =  $resultDrive. $codeVersion."/ClustersDHS/DHS_no_overlap_TSS_2bp_mid.bed";
	getmidpointDHS($unlabeledFile,$midDHSFile);
	
	my %dhs;
	 
	open IN, "<$unlabeledFile" or die "Can't open file:$unlabeledFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id, @rest) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
	
	#Report all genes that are within 50000 bp upstream or downstream of CNVs.
	#bedtools window -a CNVs.bed -b genes.bed -w 10000
	for (my $i=1;$i<=39; $i++)
	{
		my $clusterFile = $dataDrive . "geneClusters/".$i."_cluster_TSS_uniq";
		my $output = $resultDrive. $codeVersion."/ClustersDHS/".$i."_unlabeled_window";
		my $cmd = "windowBed -a $clusterFile -b $midDHSFile -w 50000 > $output";
		print $cmd;
		system($cmd);
	}

	for (my $i=1;$i<=39; $i++)
	{
		#get DHSs used in itialization for this cluster
		my %DHS_Initial;
		my $initializeFile =  $resultDrive.$codeVersion."/ClustersDHS/".$i."_initialize_uniq.bed";#_initialize_redfly_uniq.bed ";
		open IN, "<$initializeFile" or die "Can't open file:$initializeFile";	
		my $line;
		while ($line = <IN>)
		{
		
			chomp($line);
			my ($chr,$start,$end,$dhs_id, @rest) = split(/\s+/, $line);
			my ($gene,$chr,$id,$stageW,$stage,$state ) = split(/\_/,$dhs_id);
			my $theID = $chr."_".$id."_".$stageW."_".$stage;
			$DHS_Initial{$gene}{$theID}++;
			#print ("$gene\t");
		}
		close IN;	
		
		
		my $prev_gene="";
		my $index = 1;
		
		my $window = $resultDrive. $codeVersion."/ClustersDHS/".$i."_unlabeled_window";
		my $bedIn = $resultDrive. $codeVersion."/ClustersDHS/".$i."_unlabeled_in.bed";
		open IN, "<$window" or die "Can't open file:$window";
		open OUTL, ">$bedIn" or die "Can't open file:$bedIn";
		my $line;
		
		while ($line = <IN>)
		{
			chomp($line);
			my ($chr,$start,$end,$genegc,$geneflybase,$strand,$chr_dhs,$start_dhs,$end_dhs,$dhs_id) = split(/\s+/, $line);
			my $s = $dhs{$dhs_id}{'start'};
			my $e = $dhs{$dhs_id}{'end'};
			if(! exists ($DHS_Initial{$geneflybase}{$dhs_id}))
			{
				print OUTL ("$chr_dhs\t$s\t$e\t$geneflybase"."_".$dhs_id."_unlabeled\n");
				$index++;
			}
		}
		close OUT;
		
		
		my $unlabeledFileInUniq = $resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled_in_uniq.bed";
		my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq $bedIn  >$unlabeledFileInUniq ";
		#print ($cmd);
		system($cmd);
		
		my $unlabeled_out_file = $resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled.fa";
		my $cmd = "twoBitToFa -bed=$unlabeledFileInUniq -noMask ../../data/Genome/Drosophila_R5.2bit $unlabeled_out_file";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled_single.fa";
 		my $multi =$resultDrive.$codeVersion."/ClustersDHS/".$i."_unlabeled.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $multi > $single";
 		system($cmd);
	}	
	
}

#########################second stage MC
#This funciton is just for ubiquitous genes to get all genes for clusters 1-10 instead of uniq genes 
sub splitDHSIntoClusters_stark_ubiquitous_all_2nd_stage
{
	 	
 	my $DHSFile = $dataDrive. "DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";
	my %dhs;
	 
	open IN, "<$DHSFile" or die "Can't open file:$DHSFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
 	close IN;
 	
	my $dhs_file = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_stark.bed";
	open IN, "<$dhs_file" or die "Can not open file :$dhs_file";
	my $line ;
	my %gene_DHS;
	my %gene_start;
	my %gene_end;
	my %gene_chr;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($chr,$dhs_start,$dhs_end,$dhs_id,$f1,$f2,$f3,$f4,$f5,$f6,$chr_r,$start_r,$end_r,$enhancer_id,$gene) = split(/\s+/, $line);
		push(@{$gene_DHS{$gene}},$dhs_id);
		push(@{$gene_start{$gene}},$dhs{$dhs_id}{'start'});
		push(@{$gene_end{$gene}},$dhs{$dhs_id}{'end'});
		push(@{$gene_chr{$gene}},$dhs{$dhs_id}{'chr'});
 	}
	close IN;
	
	
	
	my %geneCount;
	my %dhsCount;
	my %totalGeneCount;
	#my $all_DHS_Labeled = $resultDrive.$codeVersion."/ClustersDHS/stark_specified_genes_used_in_initialization.bed";
	#open ALL, ">$all_DHS_Labeled" or die "Can not open file :$all_DHS_Labeled";
	
 	for (my $i=1; $i<=10; $i++)
 	{
 		my $clusterFile = $dataDrive."/geneClusters/$i.gene.cluster";
 		open IN, "<$clusterFile" or die "Can not open file :$clusterFile";
 		
 		my $clusterFile_fa = $resultDrive.$codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/".$i."_stark_all.bed";
 		open OUT, ">$clusterFile_fa" or die "Can not open file :$clusterFile_fa";
 		
 		while ($line = <IN> )
 		{
 			$totalGeneCount{$i}++;
 			chomp($line);
 			my $gene = $flymap{$line};
 				
 			if(exists($gene_DHS{$gene}))
	 		{
	 			
	 			for(my $j=0; $j<@{$gene_DHS{$gene}}; $j++)
	 			{
	 				my $s = $gene_start{$gene}[$j];
	 				my $e = $gene_end{$gene}[$j] ;
	 				print OUT ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene"."_".$gene_DHS{$gene}[$j]."_initial\n");
	 #				print ALL ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene_DHS{$gene}[$j]\tcluster_$i\t$gene\n");
	 				$dhsCount{$i}{$gene_DHS{$gene}[$j]}=1;
	 					
	 			}
	 			$geneCount{$i}{$gene}++;
	 		}
 			
 		}
 		close OUT;
 		close IN;
 	
 	
 	}
 	close ALL; 
}

sub splitDHSIntoClusters_redfly_specifiedGenes_ubiquitous_all_2nd_stage
{
	mapGeneIDs();
 	my $DHSFile = $dataDrive. "DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";
	my %dhs;
	 
	open IN, "<$DHSFile" or die "Can't open file:$DHSFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id,@rest) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
 	close IN;
 	
	my $dhs_file = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_redfly_specified_genes.bed";
	open IN, "<$dhs_file" or die "Can not open file :$dhs_file";
	my $line ;
	my %gene_DHS;
	my %gene_start;
	my %gene_end;
	my %gene_chr;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		 my ($chr,$dhs_start,$dhs_end,$dhs_id,$f1,$f2,$f3,$f4,$f5,$f6,$chr_r,$start_r,$end_r,$gene) = split(/\s+/, $line);
		push(@{$gene_DHS{$gene}},$dhs_id);
		push(@{$gene_start{$gene}},$dhs{$dhs_id}{'start'});
		push(@{$gene_end{$gene}},$dhs{$dhs_id}{'end'});
		push(@{$gene_chr{$gene}},$dhs{$dhs_id}{'chr'});
 	}
	close IN;
	
	
	
	my %geneCount;
	my %dhsCount;
	my %totalGeneCount;
	#my $all_DHS_Labeled = $resultDrive.$codeVersion."/ClustersDHS/redfly_dhs_specified_genes_used_in_initialization.bed";
	#open ALL, ">$all_DHS_Labeled" or die "Can not open file :$all_DHS_Labeled";
	
 	for (my $i=1; $i<=10; $i++)
 	{
 		my $clusterFile = $dataDrive."/geneClusters/$i.gene.cluster";
 		open IN, "<$clusterFile" or die "Can not open file :$clusterFile";
 		
 		my $clusterFile_fa = $resultDrive.$codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/".$i."_redfly_all.bed";
 		open OUT, ">$clusterFile_fa" or die "Can not open file :$clusterFile_fa";
 		
 		while ($line = <IN> )
 		{
 			$totalGeneCount{$i}++;
 			chomp($line);
 			
 			my $gene = $flymap{$line};
 		
 			if(exists($gene_DHS{$gene}))
	 		{
	 				
	 			for(my $j=0; $j<@{$gene_DHS{$gene}}; $j++)
	 			{
	 				my $s = $gene_start{$gene}[$j];
	 				my $e = $gene_end{$gene}[$j] ;
	 				print OUT ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene"."_".$gene_DHS{$gene}[$j]."_initial\n");
	 			#	print ALL ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene_DHS{$gene}[$j]\tcluster_$i\t$gene\n");
	 				$dhsCount{$i}{$gene_DHS{$gene}[$j]}=1;
	 					
	 			}
	 			$geneCount{$i}{$gene}++;
	 		}
	 			
 		}
 	
 		close OUT;
 		close IN;
 	
 	
 	
 	}
 	close ALL;
 	
 	print ("cluster#\t#genes_per_cluster\t#genes_per_cluster_with_DHS_overlap_redfly\t#DHSperCluster\n");
 	for (my $i=1;$i<=39; $i++)
	{
		my $c = keys %{$geneCount{$i}};
		my $dhs_c = keys %{$dhsCount{$i}};
		print ("cluster $i\t$totalGeneCount{$i}\t$c\t$dhs_c\n");
	}
	
}

sub getClosestPerClusterForUnspecifiedGenes_ubiquitous_all_2nd_stage
{
	#intersectBed -a DHS_no_overlap_TSS_10bp.bed -b ../../CRM/redfly_all_unspecified.bed -u  > DHS_overlap_redfly_unspecified_genes.bed
	#closestBed -a DHS_overlap_redfly_unspecified_genes_midpoints.bed -b ../../Genome/fly_TSS_10bp_uniq.bed -d -t all > DHS_overlap_redfly_unspecified_genes_closest.bed
	#read mapping data
	#flybase ID could be linked to more than one gene
	my $mapfile = $dataDrive."Genome/FlyBase_Genes_map.tab";

	my %flymap;
	open IN, "<$mapfile" or die "Can not open file :$mapfile";

	my $line = <IN> ;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($ensemble, $flybaseID, $annotationSymbo, $symbol) = split(/\s+/, $line);
 	
 
 		push (@{$flymap{$ensemble}}, $flybaseID);
 	
 	}	
 	close IN;
	
	
	
	my %geneClusters;
	my $line;
	for (my $i=1;$i<=10; $i++)
	{
		my $clusterFile = $dataDrive."geneClusters/".$i.".gene.cluster";
		open IN, "<$clusterFile" or die "Can't open file:$clusterFile";
		while ($line = <IN>)
		{
			chomp($line);
			#my @flybase = $flymap{$line};
			push (@{$geneClusters{$i}}, $line);
		}
	}
	close IN;
	
	#to get the actual start and end
	my %dhs;
	my $dhsFile = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_redfly_unspecified_genes.bed";
	 
	open IN, "<$dhsFile" or die "Can't open file:$dhsFile";
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
	
	my $closestFile = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_redfly_unspecified_genes_closest.bed";
	open IN, "<$closestFile" or die "Can't open file:$closestFile";
	#my %geneClusters_clone = dclone(%geneClusters);
	my %closestGenesDHS;
	my %GeneCount;
	my %DHSCount;
	my %dhs_gene;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$dhs_mid_start,$dhs_mid_end,$dhs_id,$chr_g,$gene_tss_s,$gene_tss_e,$gene,$score, $strand,$distance) = split(/\s+/, $line);
		
		for (my $i=1;$i<=10; $i++)
		{	
			#my $flybase = $flymap{$gene}[0];		
			my @output = grep {$gene eq $_} @{$geneClusters{$i}};
			if(@output >0)
			{
				$DHSCount{$i}{$dhs_id}=1;
				$GeneCount{$i}{$gene}++;
				push(@{$closestGenesDHS{$i}}, $dhs_id);
				$dhs_gene{$dhs_id}=$gene;
				#print ("$gene \t cluster $i\n");
			}
		}
		
	}
	
	for (my $i=1;$i<=10; $i++)
	{
		my $c = keys %{$GeneCount{$i}};
		my $s = keys %{$DHSCount{$i}};
		print ("cluster $i\t$c\t$s\n");
	}
	

	#my $allUNSPEC = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/redfly_dhs_unspecified_genes_used_in_initialization.bed";
	#open UNSPEC, ">$allUNSPEC" or die "Can't open file:$allUNSPEC";
	
	for (my $i=1;$i<=10; $i++)
	{
		my $index=1;
		my $outputFile = $resultDrive.$codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/".$i."_closest_dhs_unspecified_genes.bed";
		open OUT, ">$outputFile" or die "Can't open file:$outputFile";
		if(exists $closestGenesDHS{$i})
		{
			for(my $l=0; $l <@{$closestGenesDHS{$i}}; $l++)
			{
				my $id = $closestGenesDHS{$i}[$l];
				my $s = $dhs{$id}{'start'};
				my $e = $dhs{$id}{'end'};
				print OUT ("$dhs{$id}{'chr'}\t$s\t$e\t$flymap{$dhs_gene{$id}}[0]_".$id."_initial\n");
		#		print UNSPEC ("$dhs{$id}{'chr'}\t$dhs{$id}{'start'}\t$dhs{$id}{'end'}\t$id\tcluster_$i\t$flymap{$dhs_gene{$id}}[0]\n");
				$index++;
			}
		}
		close OUT;
	
		my $redFile = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_specified_redfly.bed";
		my $allFile = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_initialize_both.bed";
		my $starkFile = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_stark_all.bed";
		my $cmd = "cat $redFile $outputFile $starkFile> $allFile";
		system($cmd);
		
	}
	
		close UNSPEC;
		
		#my $cat = "cat ".$resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/redfly_dhs_specified_genes_used_in_initialization.bed ".$resultDrive.$codeVersion."/ClustersDHS/redfly_dhs_unspecified_genes_used_in_initialization.bed > ".$resultDrive.$codeVersion."/ClustersDHS/redfly_dhs_all_used_in_initialization.bed";
		#system($cat);
		
		
}

sub concatenate_ubiquitous_2nd_stage
{
	my $output = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/ubiquitous_all.bed";
	my $cmd = "cat ";
	for (my $j=1; $j<= 10; $j++)
	{
		
		$cmd = $cmd. $resultDrive.$codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/".$j."_stark_all.bed "; #   _initialize_both.bed ";	
		$cmd = $cmd. $resultDrive.$codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/".$j."_redfly_all.bed ";	
	}
	$cmd = $cmd . "> ".$output;
	print $cmd;
	print ("\n");
	system($cmd);
	
	my $ubiq_sorted = $resultDrive.$codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/ubiquitous_all_uniq.bed";
	my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq $output  >$ubiq_sorted ";
	#print ($cmd);
	system($cmd);
		
		my $ubiq_fa = $resultDrive.$codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/ubiquitous_all.fa";
		my $cmd = "twoBitToFa -bed=$ubiq_sorted -noMask ../../data/Genome/Drosophila_R5.2bit $ubiq_fa";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/ubiquitous_all_single.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $ubiq_fa > $single";
 		system($cmd);
	
}

sub getUnlabeled_2nd_stage
{
	my $all_DHS_overlap =  $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";
	my $midDHSFile =  $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp_midpoints.bed";
	
	my %dhs;
	 
	open IN, "<$all_DHS_overlap" or die "Can't open file:$all_DHS_overlap";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
	
	#Report all genes that are within 50000 bp upstream or downstream of CNVs.
	#bedtools window -a CNVs.bed -b genes.bed -w 10000
	for (my $i=11;$i<=39; $i++)
	{
		my $clusterFile = $dataDrive . "geneClusters/".$i."_cluster_TSS_non_uniq";
		my $output = $resultDrive. $codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/".$i."_unlabeled_window";
		my $cmd = "windowBed -a $clusterFile -b $midDHSFile -w 50000 > $output";
		print $cmd;
		system($cmd);
	}

	for (my $i=11;$i<=39; $i++)
	{
		#get DHSs used in stage 1 for this cluster
		my %DHS_Initial;
		my $initializeFile =  $resultDrive.$codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/".$i."_pos.fa";#_initialize_redfly_uniq.bed ";
		open IN, "<$initializeFile" or print "Can't open file:$initializeFile";	
		my $line;
		while ($line = <IN>)
		{
			chomp($line);
			if($line =~ "^\>")
			{
				my ($gene,$chr,$start,$stageW, $stage,$status) = split(/\_/, $line);
				my $theID = $chr."_".$start."_".$stageW."_".$stage;
				$DHS_Initial{$gene}{$theID}++;
			}
		}
		close IN;	
	
		
		my $prev_gene="";
		my $index = 1;
		
		my $window = $resultDrive. $codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/".$i."_unlabeled_window";
		my $bedIn = $resultDrive. $codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/".$i."_unlabeled_in.bed";
		open IN, "<$window" or die "Can't open file:$window";
		open OUTL, ">$bedIn" or die "Can't open file:$bedIn";
		my $line;
		
		while ($line = <IN>)
		{
			chomp($line);
			my ($chr,$start,$end,$gene,$flybase,$strand,$chr_dhs,$start_dhs,$end_dhs,$dhs_id) = split(/\s+/, $line);
			my $s = $dhs{$dhs_id}{'start'};
			my $e = $dhs{$dhs_id}{'end'};
			if(! exists ($DHS_Initial{$flybase}{$dhs_id}))
			{
				print OUTL ("$chr_dhs\t$s\t$e\t$flybase"."_".$dhs_id."_unlabeled2\n");
				$index++;
			}
		}
		close OUT;
		
		
		my $unlabeledFileInUniq = $resultDrive.$codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/".$i."_unlabeled_in_uniq.bed";
		my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq $bedIn  >$unlabeledFileInUniq ";
		print ($cmd);
		 system($cmd);
		
		my $unlabeled_out_file = $resultDrive.$codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/".$i."_unlabeled.fa";
		my $cmd = "twoBitToFa -bed=$unlabeledFileInUniq -noMask ../../data/Genome/Drosophila_R5.2bit $unlabeled_out_file";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/".$i."_unlabeled_single.fa";
 		my $multi =$resultDrive.$codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/".$i."_unlabeled.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $multi > $single";
 		system($cmd);
	}	
	
}


sub getUnlabeled_2nd_stage_old
{
	
	#get DHS used for uniq genes
	my %DHS_Initial;
	for (my $j=11; $j<= 39; $j++)
	{
		
		my $initializeFile =  $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$j."_pos.fa";
		open IN, "<$initializeFile" or print "Can't open file:$initializeFile";	
		my $line;
		while ($line = <IN>)
		{
			chomp($line);
			if($line =~ "^\>")
			{
			
				my ($gene,$chr,$start,$stageW, $stage,$status) = split(/\_/, $line);
				my $theID = $chr."_".$start."_".$stageW."_".$stage;
				$DHS_Initial{$theID}++;
			}
		}
		close IN;	
	}
	
	#read ubiquitous
	my $initializeFile =  $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/ubiquitous_all_single.fa";
	open IN, "<$initializeFile" or die "Can't open file:$initializeFile";	
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		if($line =~ "^\>")
		{
			my ($gene,$chr,$start,$stageW, $stage,$status) = split(/\_/, $line);
			my $theID = $chr."_".$start."_".$stageW."_".$stage;
			$DHS_Initial{$theID}++;
		}
	}
	close IN;	
	my $k = keys %DHS_Initial;
	print ("number of DHS after MC first stage = $k\n");
	my $all_DHS_overlap =  $dataDrive."DHS_row_data/Jamm_peaks/DHS_no_overlap_TSS_2bp.bed";
	my $unlabeledFile = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/DHS_unlabeled_second_stage.bed";
	open IN, "<$all_DHS_overlap" or die "Can't open file:$all_DHS_overlap";	
	open OUT, ">$unlabeledFile" or die "Can't open file:$unlabeledFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id, @rest) = split(/\s+/, $line);
		if(exists $DHS_Initial{$dhs_id})
		{
			#print ("DHS $dhs_id does  exist\n");
		}
		else{
			print OUT  ("$line\n");
		}
	}
	close IN;
	close OUT;
	
	my $midDHSFile =  $resultDrive. $codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/DHS_unlabeled_second_stage_midpoints.bed";
	getmidpointDHS($unlabeledFile,$midDHSFile);
	
	my %dhs;
	 
	open IN, "<$unlabeledFile" or die "Can't open file:$unlabeledFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
	
	#Report all genes that are within 50000 bp upstream or downstream of CNVs.
	#bedtools window -a CNVs.bed -b genes.bed -w 10000
	for (my $i=11;$i<=39; $i++)
	{
		my $clusterFile = $dataDrive . "geneClusters/".$i."_cluster_TSS_non_uniq";
		my $output = $resultDrive. $codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_unlabeled_window";
		my $cmd = "windowBed -a $clusterFile -b $midDHSFile -w 50000 > $output";
		print $cmd;
		system($cmd);
	}

	for (my $i=11;$i<=39; $i++)
	{
		my $prev_gene="";
		my $index = 1;
		
		my $window = $resultDrive. $codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_unlabeled_window";
		my $bedIn = $resultDrive. $codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_unlabeled_in.bed";
		open IN, "<$window" or die "Can't open file:$window";
		open OUTL, ">$bedIn" or die "Can't open file:$bedIn";
		my $line;
		
		while ($line = <IN>)
		{
			chomp($line);
			my ($chr,$start,$end,$gene,$flybase,$strand,$chr_dhs,$start_dhs,$end_dhs,$dhs_id) = split(/\s+/, $line);
			my $s = $dhs{$dhs_id}{'start'};
			my $e = $dhs{$dhs_id}{'end'};
			print OUTL ("$chr_dhs\t$s\t$e\t$flybase"."_".$dhs_id."_unlabeled2\n");
			$index++;
		}
		close OUT;
		
		
		my $unlabeledFileInUniq = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_unlabeled_in_uniq.bed";
		my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $bedIn  >$unlabeledFileInUniq ";
		#print ($cmd);
		system($cmd);
		
		my $unlabeled_out_file = $resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_unlabeled.fa";
		my $cmd = "twoBitToFa -bed=$unlabeledFileInUniq -noMask ../../data/Genome/Drosophila_R5.2bit $unlabeled_out_file";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_unlabeled_single.fa";
 		my $multi =$resultDrive.$codeVersion."/ClustersDHS/MC_3order/MC_ubiq_output/".$i."_unlabeled.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $multi > $single";
 		system($cmd);
	}	
	
}

##################TSS
sub splitDHSIntoClusters_TSS
{
 	my $DHSFile = $dataDrive. "DHS_row_data/Jamm_peaks/DHS_overlap_TSS_2bp.bed";
	my %dhs;
	 
	open IN, "<$DHSFile" or die "Can't open file:$DHSFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
 	close IN;
 	
	my $dhs_file = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_TSS_2bp_with_gene_names.bed";
	open IN, "<$dhs_file" or die "Can not open file :$dhs_file";
	my $line ;
	my %gene_DHS;
	my %gene_start;
	my %gene_end;
	my %gene_chr;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($chr,$dhs_start,$dhs_end,$dhs_id,$f1,$f2,$f3,$f4,$f5,$f6,$chr_r,$start_r,$end_r,$cg_gene, $gene, $strand) = split(/\s+/, $line);
		push(@{$gene_DHS{$gene}},$dhs_id);
		push(@{$gene_start{$gene}},$dhs{$dhs_id}{'start'});
		push(@{$gene_end{$gene}},$dhs{$dhs_id}{'end'});
		push(@{$gene_chr{$gene}},$dhs{$dhs_id}{'chr'});
 	}
	close IN;
	
	
	
	my %geneCount;
	my %dhsCount;
	my %totalGeneCount;
	
 	for (my $i=1; $i<=39; $i++)
 	{
 		my $clusterFile = $dataDrive."/geneClusters/$i.gene.cluster.uniq";
 		open IN, "<$clusterFile" or die "Can not open file :$clusterFile";
 		
 		my $clusterFile_fa = $resultDrive.$codeVersion."/ClustersDHS/".$i."_TSS.bed";
 		open OUT, ">$clusterFile_fa" or die "Can not open file :$clusterFile_fa";
 		
 		while ($line = <IN> )
 		{
 			$totalGeneCount{$i}++;
 			chomp($line);
 			my ($gene,$flybase) =  split(/\s+/, $line);
 			
 				if(exists($gene_DHS{$flybase}))
	 			{
	 				
	 				for(my $j=0; $j<@{$gene_DHS{$flybase}}; $j++)
	 				{
	 					my $s = $gene_start{$flybase}[$j];
	 					my $e = $gene_end{$flybase}[$j] ;
	 					print OUT ("$gene_chr{$flybase}[$j]\t$s\t$e\t$flybase"."_".$gene_DHS{$flybase}[$j]."_TSS\n");
#	 					print ALL ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene_DHS{$gene}[$j]\tcluster_$i\t$gene\n");
	 					$dhsCount{$i}{$gene_DHS{$flybase}[$j]}=1;
	 					
	 				}
	 				$geneCount{$i}{$flybase}++;
	 			
 			}
 		}
 		close OUT;
 		close IN;
 	
 		print ("cluster $i\n");
 		my $uniq_file =  $resultDrive.$codeVersion."/ClustersDHS/".$i."_TSS_uniq.bed";
 		my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq $clusterFile_fa > $uniq_file";
 		system($cmd);	

		my $fa_file = $resultDrive.$codeVersion."/ClustersDHS/".$i."_TSS_uniq.fa";
		my $cmd = "twoBitToFa -bed=$uniq_file -noMask ../../data/Genome/Drosophila_R5.2bit $fa_file"; 
 		system($cmd);
 		
 		my $single =$resultDrive.$codeVersion."/ClustersDHS/".$i."_TSS_uniq_single.fa";
 		my $cmd = "perl /data/ohler/Dina/Drosophila/Research/code/scripts/convert_multi_fasta_to_single.pl $fa_file > $single";
 		system($cmd);
 	
 	}
 
 	
 	print ("cluster#\t#genes_per_cluster\t#genes_per_cluster_with_DHS_overlap_redfly\t#DHSperCluster\n");
 	for (my $i=1;$i<=39; $i++)
	{
		my $c = keys %{$geneCount{$i}};
		my $dhs_c = keys %{$dhsCount{$i}};
		print ("cluster $i\t$totalGeneCount{$i}\t$c\t$dhs_c\n");
	}
	
}

sub concatenate_ubiquitous_TSS
{
	my $output = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous_TSS.bed";
	my $cmd = "cat ";
	for (my $j=1; $j<= 10; $j++)
	{
		
		$cmd = $cmd. $resultDrive.$codeVersion."/ClustersDHS/".$j."_TSS_uniq.bed ";  # _initialize_redfly_uniq.bed ";		
	}
	$cmd = $cmd . "> ".$output;
	print $cmd;
	print ("\n");
	system($cmd);
	
	my $ubiq_sorted = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous_TSS_uniq.bed";
	my $cmd = "sort -k1,1 -k2,2n -k3,3n -k4,4 --uniq $output  >$ubiq_sorted ";
	#print ($cmd);
	system($cmd);
		
		my $ubiq_fa = $resultDrive.$codeVersion."/ClustersDHS/ubiquitous_TSS.fa";
		my $cmd = "twoBitToFa -bed=$ubiq_sorted -noMask ../../data/Genome/Drosophila_R5.2bit $ubiq_fa";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/ubiquitous_TSS_single.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $ubiq_fa > $single";
 		system($cmd);
	
}
##################TSS stage 2
#intersectBed -a Embryo_all_stages.peaks -b ../../Genome/fly_TSS_10bp_uniq.bed -wb > DHS_overlap_TSS_10bp_with_gene_names.bed

sub splitDHSIntoClusters_TSS_all_genes
{
	mapGeneIDs();
 	my $DHSFile = $dataDrive. "DHS_row_data/Jamm_peaks/DHS_overlap_TSS_2bp.bed";
	my %dhs;
	 
	open IN, "<$DHSFile" or die "Can't open file:$DHSFile";
	my $line;
	while ($line = <IN>)
	{
		chomp($line);
		my ($chr,$start,$end,$dhs_id) = split(/\s+/, $line);
		$dhs{$dhs_id}{'start'} = $start;
		$dhs{$dhs_id}{'end'} = $end;
		$dhs{$dhs_id}{'chr'} = $chr;
		
	}
 	close IN;
 	
	my $dhs_file = $dataDrive."DHS_row_data/Jamm_peaks/DHS_overlap_TSS_2bp_with_gene_names.bed";
	open IN, "<$dhs_file" or die "Can not open file :$dhs_file";
	my $line ;
	my %gene_DHS;
	my %gene_start;
	my %gene_end;
	my %gene_chr;
 	while ($line = <IN> )
 	{
 		chomp($line);
 		my ($chr,$dhs_start,$dhs_end,$dhs_id,$f1,$f2,$f3,$f4,$f5,$f6,$chr_r,$start_r,$end_r,$cg_gene, $gene, $strand) = split(/\s+/, $line);
		push(@{$gene_DHS{$gene}},$dhs_id);
		push(@{$gene_start{$gene}},$dhs{$dhs_id}{'start'});
		push(@{$gene_end{$gene}},$dhs{$dhs_id}{'end'});
		push(@{$gene_chr{$gene}},$dhs{$dhs_id}{'chr'});
 	}
	close IN;
	
	
	
	my %geneCount;
	my %dhsCount;
	my %totalGeneCount;
	
 	for (my $i=1; $i<=39; $i++)
 	{
 		my $clusterFile = $dataDrive."/geneClusters/$i.gene.cluster";
 		open IN, "<$clusterFile" or die "Can not open file :$clusterFile";
 		
 		my $clusterFile_fa = $resultDrive.$codeVersion."/ClustersDHS/ModelSelection/ModelStage1/MC_multi_output/ModelStage2/MC_multi_output/".$i."_TSS_all.bed";
 		open OUT, ">$clusterFile_fa" or die "Can not open file :$clusterFile_fa";
 		
 		while ($line = <IN> )
 		{
 			$totalGeneCount{$i}++;
 			chomp($line);
 			
 			my $flybase = $flymap{$line};
 				if(exists($gene_DHS{$flybase}))
	 			{
	 				
	 				for(my $j=0; $j<@{$gene_DHS{$flybase}}; $j++)
	 				{
	 					my $s = $gene_start{$flybase}[$j];
	 					my $e = $gene_end{$flybase}[$j] ;
	 					print OUT ("$gene_chr{$flybase}[$j]\t$s\t$e\t$flybase"."_".$gene_DHS{$flybase}[$j]."_TSS\n");
#	 					print ALL ("$gene_chr{$gene}[$j]\t$s\t$e\t$gene_DHS{$gene}[$j]\tcluster_$i\t$gene\n");
	 					$dhsCount{$i}{$gene_DHS{$flybase}[$j]}=1;
	 					
	 				}
	 				$geneCount{$i}{$flybase}++;
	 			
 			}
 		}
 		close OUT;
 		close IN;
 	
 		print ("cluster $i\n");
 		my $uniq_file =  $resultDrive.$codeVersion."/ClustersDHS/ModelSelection/ModelStage1/MC_multi_output/ModelStage2/MC_multi_output/".$i."_TSS_all_uniq.bed";
 		my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $clusterFile_fa > $uniq_file";
 		system($cmd);	

		my $fa_file = $resultDrive.$codeVersion."/ClustersDHS/ModelSelection/ModelStage1/MC_multi_output/ModelStage2/MC_multi_output/".$i."_TSS_all_uniq.fa";
		my $cmd = "twoBitToFa -bed=$uniq_file -noMask ../../data/Genome/Drosophila_R5.2bit $fa_file"; 
 		system($cmd);
 		
 		my $single =$resultDrive.$codeVersion."/ClustersDHS/ModelSelection/ModelStage1/MC_multi_output/ModelStage2/MC_multi_output/".$i."_TSS_all_uniq_single.fa";
 		my $cmd = "perl /data/ohler/Dina/Drosophila/Research/code/scripts/convert_multi_fasta_to_single.pl $fa_file > $single";
 		system($cmd);
 	
 	}
 
 	
 	print ("cluster#\t#genes_per_cluster\t#genes_per_cluster_with_DHS_overlap_redfly\t#DHSperCluster\n");
 	for (my $i=1;$i<=39; $i++)
	{
		my $c = keys %{$geneCount{$i}};
		my $dhs_c = keys %{$dhsCount{$i}};
		print ("cluster $i\t$totalGeneCount{$i}\t$c\t$dhs_c\n");
	}
	
}

sub concatenate_ubiquitous_TSS_all_genes
{
	my $output = $resultDrive.$codeVersion."/ClustersDHS/MC_2order/MC_ubiq_output/ubiquitous_TSS_all.bed";
	my $cmd = "cat ";
	for (my $j=1; $j<= 10; $j++)
	{
		
		$cmd = $cmd. $resultDrive.$codeVersion."/ClustersDHS/ModelSelection/ModelStage1/MC_multi_output/ModelStage2/MC_multi_output/".$j."_TSS_all_uniq.bed ";  # _initialize_redfly_uniq.bed ";		
	}
	$cmd = $cmd . "> ".$output;
	print $cmd;
	print ("\n");
	system($cmd);
	
	my $ubiq_sorted = $resultDrive.$codeVersion."/ClustersDHS/ModelSelection/ModelStage1/MC_multi_output/ModelStage2/MC_multi_output/ubiquitous_TSS_all_uniq.bed";
	my $cmd = "sort -k1,1 -k2,2n -k3,3n --uniq $output  >$ubiq_sorted ";
	#print ($cmd);
	system($cmd);
		
		my $ubiq_fa = $resultDrive.$codeVersion."/ClustersDHS/ModelSelection/ModelStage1/MC_multi_output/ModelStage2/MC_multi_output/ubiquitous_TSS_all.fa";
		my $cmd = "twoBitToFa -bed=$ubiq_sorted -noMask ../../data/Genome/Drosophila_R5.2bit $ubiq_fa";
 		print ($cmd);
 		system($cmd);
		
		my $single =$resultDrive.$codeVersion."/ClustersDHS/ModelSelection/ModelStage1/MC_multi_output/ModelStage2/MC_multi_output/ubiquitous_TSS_all_single.fa";
 		my $cmd = "perl ../scripts/convert_multi_fasta_to_single.pl $ubiq_fa > $single";
 		system($cmd);
	
}

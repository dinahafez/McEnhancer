use strict;
use List::Util qw(sum);
use POSIX;
use Data::Dumper;
use List::MoreUtils 'any';
use Storable qw(dclone);

my $dataDrive = "";
my $resultDrive = "";

my $dhsLabeledFileName = $ARGV[0];
my $dhsUnLabeledFileName = $ARGV[1]; 
my $dhsBackgroundFileName = $ARGV[2];
my $dhsLabeledFileNameTemp = $ARGV[3]; 
my $order = $ARGV[4];  
my $max_itr = $ARGV[5];    

my $discount = $ARGV[6];
my $discount_n = $ARGV[7];
my $dif = $ARGV[8];

my %positiveDHS;
my %unlabeledDHS;
my @unlabeledDHSIDS=();
my @positiveDHSids=();
main();

sub main
{
	#calculate model parameters for labeled and background
        #positive MC
        my $positiveCountFileName = $resultDrive.$dhsLabeledFileNameTemp ;
        ngramCount ($dhsLabeledFileName, $positiveCountFileName, $order, "T");
        #negative MC 
        my $negativeCountFileName = $resultDrive.$dhsLabeledFileNameTemp.".neg";
       	ngramCount ($dhsBackgroundFileName,$negativeCountFileName, $order, "F");

	my $itr=1;
	#score unlabeled DHS
	my $positiveScores = scoreUnlabeledDHSFromFile ($dhsUnLabeledFileName, $positiveCountFileName, $order, "T");
	my $negativeScores = scoreUnlabeledDHSFromFile ($dhsUnLabeledFileName, $negativeCountFileName, $order, "F");
	my $prevNoOfUnlabeled = 0;
	my $noOfUnlabeled = 0;

	do{

		if($itr != 1)
		{
			#score unlabeled DHS
			my $unlabeledTempFile = $resultDrive.$dhsLabeledFileNameTemp.".unlabeled.temp";
			$positiveScores = scoreUnlabeledDHS ($unlabeledTempFile, $positiveCountFileName, $order);
			$negativeScores = scoreUnlabeledDHS ($unlabeledTempFile, $negativeCountFileName, $order);
			$prevNoOfUnlabeled = $noOfUnlabeled;
		}
		
		#update labeled and unlabeled 
		$noOfUnlabeled = updateLabeledAndUnlabeled($positiveScores, $negativeScores, $dhsLabeledFileName, $dhsUnLabeledFileName);
		#print Dumper @unlabeledDHSIDS;
		#reestimate positive DHSs
		$noOfUnlabeled = reestimatePositiveModel($dhsLabeledFileName, $positiveCountFileName,$negativeCountFileName,$order);
		
		print ("iteration number $itr\n");
		$itr++;
		
		print ("number of unlableled $noOfUnlabeled $prevNoOfUnlabeled\n");
		
	#}while ($itr <= $max_itr);
	}while (($prevNoOfUnlabeled != $noOfUnlabeled) && ($itr <= $max_itr));

	writePosFile();
	cleanFiles ();
}


sub ngramCount 
{
	
	my ($fastaFileName, $modelFileName, $order, $flag) =@_;

	open IN, "<$dataDrive$fastaFileName" or die "Can't open file:$dataDrive$fastaFileName";
	my $tempInputFile = $resultDrive.$modelFileName . ".temp";
	my $countFileName = $resultDrive.$modelFileName . ".count";
	open OUT, ">$tempInputFile" or die "Can't open file:$tempInputFile";
	my $line;
	my $id;
	while ($line = <IN>)
	{
		chomp($line);
		
		if ($line =~ /^>/) #save the dhs header
		{
			$id = $line;
		}
		else
		{
			#convert the sequence into a space separated characters
			my $cmd = "sed 's/\(.\)/\1 /g;s/ $//' $line";
			print OUT join(" ",split(//,$line));
			print OUT ("\n");
			if ($flag eq "T")
			{
				$positiveDHS{$id} = $line;
				push (@positiveDHSids,$id);
			}
		}
		
	}
	close (IN);
	close (OUT);

        #call ngram count from srilm
	my $cmd = "ngram-count -text ".$tempInputFile." -".$discount." ".$discount_n." -lm ".$countFileName."_1 -order 1 -no-sos -no-eos";
	print ("$cmd\n");      	
	system ($cmd);

        $cmd = "ngram-count -text ".$tempInputFile." -".$discount." ".$discount_n." -lm ".$countFileName."_2 -order 2  -no-sos -no-eos";
	print ("$cmd\n");      	
	system ($cmd);

	$cmd = "ngram-count -text ".$tempInputFile." -".$discount." ".$discount_n." -lm ".$countFileName."_3 -order 3 -no-sos -no-eos";
	print ("$cmd\n");      	
	system ($cmd);

	$cmd = "ngram-count -text ".$tempInputFile." -".$discount." ".$discount_n." -lm ".$countFileName."_4 -order 4 -no-sos -no-eos";
	print ("$cmd\n");      	
	system ($cmd);

	$cmd = "ngram-count -text ".$tempInputFile." -".$discount." ".$discount_n." -lm ".$countFileName."_5 -order 5 -no-sos -no-eos";
	print ("$cmd\n");      	
	system ($cmd);  
	              
}

sub scoreUnlabeledDHSFromFile 
{
	my ($dhsUnLabeledFileName, $modelFileName ,$order, $flag) = @_;
	my $tempFile = $resultDrive.$dhsLabeledFileNameTemp.".unlabeled.temp";
	my $countFileName = $resultDrive.$modelFileName.".count";

	if ($flag eq "T") #parse only one time
	{
		open IN, "<$dataDrive$dhsUnLabeledFileName" or die "Can't open file:$dataDrive$dhsUnLabeledFileName";
		open OUT, ">$tempFile" or die "Can't open file:$tempFile";
		my $line;
		my $id;
		while ($line = <IN>)
		{
			chomp($line);
			if ($line =~ /^>/) #if not the header for fasta
			{
				push (@unlabeledDHSIDS,$line);
				$id = $line;
			}		
			else		
			{
				#convert the sequence into a space separated characters
				my $cmd = "sed 's/\(.\)/\1 /g;s/ $//' $line";
				print OUT join(" ",split(//,$line));
				print OUT ("\n");
				$unlabeledDHS{$id} = $line;
			}
		
		}
		close (IN);
		close (OUT);
	}
	my $cmd = "";
	if ($order == 5)
	{
		$cmd = "ngram -lm ".$countFileName."_5  -lambda 0.2 -mix-lm ".$countFileName."_1 -mix-lm2 ".$countFileName."_2 -mix-lambda2 0.2 -mix-lm3 ".$countFileName."_3 -mix-lambda3 0.2 -mix-lm4 ".$countFileName."_4 -mix-lambda4 0.2 -bayes 0 -ppl ".$tempFile." -debug 1 -no-sos -no-eos | grep logprob ";
		print ("$cmd\n");  
		
	}  
	if ($order == 4)
	{
		$cmd = "ngram -lm ".$countFileName."_4  -lambda 0.3 -mix-lm ".$countFileName."_1 -mix-lm2 ".$countFileName."_2 -mix-lambda2 0.25 -mix-lm3 ".$countFileName."_3 -mix-lambda3 0.25 -bayes 0 -ppl ".$tempFile." -debug 1 -no-sos -no-eos | grep logprob ";
		print ("$cmd\n");  
		
	}  
	if ($order == 3)
	{
		$cmd = "ngram -lm ".$countFileName."_3  -lambda 0.3 -mix-lm ".$countFileName."_1 -mix-lm2 ".$countFileName."_2 -mix-lambda2 0.4 -bayes 0 -ppl ".$tempFile." -debug 1 -no-sos -no-eos | grep logprob ";
		print ("$cmd\n");  
		
	}  

	my @logprob = `$cmd`;
	return ( \@logprob);
}

sub scoreUnlabeledDHS 
{
	my ($testFileName, $modelFileName ,$order) = @_;
	my $countFileName = $resultDrive.$modelFileName.".count";
	
	my $cmd = "";
 	if ($order == 5)
	{
		$cmd = "ngram -lm ".$countFileName."_5  -lambda 0.2 -mix-lm ".$countFileName."_1 -mix-lm2 ".$countFileName."_2 -mix-lambda2 0.2 -mix-lm3 ".$countFileName."_3 -mix-lambda3 0.2 -mix-lm4 ".$countFileName."_4 -mix-lambda4 0.2 -bayes 0 -ppl ".$testFileName." -debug 1 -no-sos -no-eos | grep logprob ";
		print ("$cmd\n");  
		
	}  
	if ($order == 4)
	{
		$cmd = "ngram -lm ".$countFileName."_4  -lambda 0.3 -mix-lm ".$countFileName."_1 -mix-lm2 ".$countFileName."_2 -mix-lambda2 0.25 -mix-lm3 ".$countFileName."_3 -mix-lambda3 0.25 -bayes 0 -ppl ".$testFileName." -debug 1 -no-sos -no-eos | grep logprob ";
		print ("$cmd\n");  
		
	}  
	if ($order == 3)
	{
		$cmd = "ngram -lm ".$countFileName."_3  -lambda 0.3 -mix-lm ".$countFileName."_1 -mix-lm2 ".$countFileName."_2 -mix-lambda2 0.4 -bayes 0 -ppl ".$testFileName." -debug 1 -no-sos -no-eos | grep logprob ";
		print ("$cmd\n");  
		
	}  
	my @logprob = `$cmd`;
	return ( \@logprob);
}

sub updateLabeledAndUnlabeled
{
	my ( $positiveScores, $negativeScores, $posFastaFileName, $unlabeledFileName) = @_;

	my $tempInputFile = $resultDrive.$dhsLabeledFileNameTemp . ".temp";		#positive file temp
	open OUT, ">>$tempInputFile" or die "Can't open file:$tempInputFile";

	my $tempUnlabeledFile = $resultDrive.$dhsLabeledFileNameTemp . ".unlabeled.temp";	#unlabeled
	open OUTUNLABELED, ">$tempUnlabeledFile" or die "Can't open file:$tempUnlabeledFile";

	my @indexToDelete=();

	my $num = @unlabeledDHSIDS;
	#print Dumper @unlabeledDHSIDS;
	#print ("no of unlabeled DHSs $num\n");
	for (my $i=0; $i<@unlabeledDHSIDS; $i++ )
	{
		my @fields = split (/\s+/,$$positiveScores[$i]);
		my $posProb = $fields[3];
		my $posPpl = $fields[5];
 		
		@fields = split (/\s+/,$$negativeScores[$i]);
		my $negProb = $fields[3];
		my $negPpl = $fields[5];
		
		my $id = $unlabeledDHSIDS[$i];
		if ($posPpl =~ "undefined")
		{ $posPpl = 100000;}
		if ($negPpl =~ "undefined")
		{ $negPpl = 100000;}

		#if ($posPpl < $negPpl)   #positive DHS
		#{
			if (( $negPpl - $posPpl) > $dif) 
			{delete @unlabeledDHSIDS[$i];

			my $seq = $unlabeledDHS{$id};
			$positiveDHS{$id} = $seq;
			
			#print ("POSITIVE\n");
			#print ("$id\t$posPpl\t$negPpl\n");
			delete $unlabeledDHS{$id};

			print OUT join(" ",split(//,$positiveDHS{$id}));
			print OUT ("\n");

			push (@positiveDHSids,$id);
			
		}
		else
		{
			print OUTUNLABELED join(" ",split(//,$unlabeledDHS{$id}));
			print OUTUNLABELED ("\n");
			
		}
		
	}
	close OUT;
	close OUTUNLABELED;
	
	
	my @defined = grep {defined}  @unlabeledDHSIDS;
	@unlabeledDHSIDS = @defined;
	my $count = @unlabeledDHSIDS;
	return $count;
}

sub reestimatePositiveModel
{
	my ($posFastaFileName, $positiveModelFileName, $negativeModelFileName ,$order) = @_;
	my $tempInputFile = $resultDrive.$dhsLabeledFileNameTemp . ".temp";
	my $positiveCountFileName = $resultDrive.$positiveModelFileName . ".count";

       #call ngram count from srilm
	my $cmd = "ngram-count -text ".$tempInputFile." -".$discount." ".$discount_n." -lm ".$positiveCountFileName."_1 -order 1 -no-sos -no-eos";
	print ("$cmd\n");      	
	system ($cmd);

        $cmd = "ngram-count -text ".$tempInputFile." -".$discount." ".$discount_n." -lm ".$positiveCountFileName."_2 -order 2  -no-sos -no-eos";
	print ("$cmd\n");      	
	system ($cmd);

	$cmd = "ngram-count -text ".$tempInputFile." -".$discount." ".$discount_n." -lm ".$positiveCountFileName."_3 -order 3 -no-sos -no-eos";
	print ("$cmd\n");      	
	system ($cmd);

	$cmd = "ngram-count -text ".$tempInputFile." -".$discount." ".$discount_n." -lm ".$positiveCountFileName."_4 -order 4 -no-sos -no-eos";
	print ("$cmd\n");      	
	system ($cmd);

	$cmd = "ngram-count -text ".$tempInputFile." -".$discount." ".$discount_n." -lm ".$positiveCountFileName."_5 -order 5 -no-sos -no-eos";
	print ("$cmd\n");      	
	system ($cmd);  


	my $testFile = $resultDrive.$dhsLabeledFileNameTemp.".temp";
	my $positiveScores = scoreUnlabeledDHS ($testFile, $positiveModelFileName, $order);
	my $negativeScores = scoreUnlabeledDHS ($testFile, $negativeModelFileName, $order);


	my $posFileTemp = $resultDrive.$dhsLabeledFileNameTemp . ".temp";
	open OUT, ">$posFileTemp" or die "Can't open file:$posFileTemp";

	my $tempUnlabeledFile = $resultDrive.$dhsLabeledFileNameTemp . ".unlabeled.temp";
	open OUTUNLABELED, ">>$tempUnlabeledFile" or die "Can't open file:$tempUnlabeledFile";

	my @indexToDelete=();

	for (my $i=0; $i<@positiveDHSids; $i++ )
	{
		my @fields = split (/\s+/,$$positiveScores[$i]);
		my $posProb = $fields[3];
		my $posPpl = $fields[5];
 		
		@fields = split (/\s+/,$$negativeScores[$i]);
		my $negProb = $fields[3];
		my $negPpl = $fields[5];
		
		my $id = $positiveDHSids[$i];
		
		if ($posPpl < $negPpl)	#positive DHS
		{
			print OUT join(" ",split(//,$positiveDHS{$id}));
			print OUT ("\n");
		}
		else
		{
			
			delete @positiveDHSids[$i];
			
			my $seq = $positiveDHS{$id};
			$unlabeledDHS{$id} = $seq ;
			delete $positiveDHS{$id};

			print OUTUNLABELED join(" ",split(//,$unlabeledDHS{$id}));
			print OUTUNLABELED ("\n");

			push (@unlabeledDHSIDS,$id);
		}
		
	}
	close OUT;
	close OUTUNLABELED;

	my @fields = split (/\s+/,$$positiveScores[@positiveDHSids]);
	my $posProb = $fields[3];
	my $posPpl = $fields[5];
	my $totalPpl = $posPpl;


	my @defined = grep {defined}  @positiveDHSids;
	@positiveDHSids = @defined;

	my $count = @unlabeledDHSIDS;
	return $count;
		
}

sub writePosFile
{
	open OUT, ">$resultDrive$dhsLabeledFileNameTemp" or die "Can't open file:$resultDrive$dhsLabeledFileNameTemp";

	for (my $i=0; $i< @positiveDHSids; $i++)
	{
		print OUT ("$positiveDHSids[$i]\n");
		print OUT ("$positiveDHS{$positiveDHSids[$i]}\n");

	}
	close OUT;
}

sub cleanFiles
{
	`rm $resultDrive$dhsLabeledFileNameTemp".count_1";`; 
	`rm $resultDrive$dhsLabeledFileNameTemp".count_2";`; 
	`rm $resultDrive$dhsLabeledFileNameTemp".count_3";`; 
	`rm $resultDrive$dhsLabeledFileNameTemp".count_4";`; 
	`rm $resultDrive$dhsLabeledFileNameTemp".count_5";`; 
	`rm $resultDrive$dhsLabeledFileNameTemp".neg.count_1";`; 
	`rm $resultDrive$dhsLabeledFileNameTemp".neg.count_2";`;
	`rm $resultDrive$dhsLabeledFileNameTemp".neg.count_3";`;
	`rm $resultDrive$dhsLabeledFileNameTemp".neg.count_4";`;
	`rm $resultDrive$dhsLabeledFileNameTemp".neg.count_5";`; 
	`rm $resultDrive$dhsLabeledFileNameTemp".temp"`; 
	`rm $resultDrive$dhsLabeledFileNameTemp".neg.temp";`;
	`rm $resultDrive$dhsLabeledFileNameTemp".unlabeled.temp";`;
}


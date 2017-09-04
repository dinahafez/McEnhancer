!/bin/bash
#
#$ -S /bin/tcsh -cwd
#$ -o error2.out -j y
#$ -l h_vmem=10G
#$ -l mem_free=10G
#$ -e error2.log

dir=$1;		#output directory
order=$2;	#order of Markov chain 
itr=$3;		#max. number of iterations to keep estimating
cluster=$4;	#Cluster number to process

discount=$5;     #Type of discount used by srlim
discountn=$6;	 #value of the discount
dif=$7;		#Difference between the probability estimates of the positive and negative. DHSs with probability bellow this cutoff, will be put in a reject class. These are generally DHSs that are on the borders

if [ $# -ne 7 ]
  then
    echo "Please provide these arguments in this order:"
    echo "Directory where your files are stored. This will also be the output directory"
    echo "order of Markov chain  -- after model selection, I used 3" 
    echo "Max. number of iterations for the EM algorithm to keep running. Best is 3"
    echo "Cluster number to process. This is the positive class"
    echo "Type of discount to be used by srlim. Example wd"
    echo "Value of the discount -- 0.1"
    echo "Difference between the probability estimates of the positive and negative. DHSs with probability bellow this cutoff, will be put in a reject class. These are generally DHSs that are on the borders"

else
    #Markov chain
    ./pipeline_run_per_cluster.sh $dir $itr $cluster $order $discount $discountn $dif;
fi

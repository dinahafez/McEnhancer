# McEnhancer
Source code and results for McEnhancer: Predicting gene expression via semi-supervised assignment of enhancers to  target genes

 
This is the Readme file for running McEnhancer.

Data:

All DHSs and gene clusters data are saved in data folder. 
DHSs has row reads downloaded. Before processing, DHSs were mapped to dm3 genome. Then Jamm_v_1.0.6rev4 was used to call peaks on each embryo stage. The merge of all reproducible reads from all stages is in allStages.filtered.peaks.idr.0.2.narrowPeaks

DHSs were then  divided into two groups, those that overlap TSSs and those that are distal. 

CRMs coordinates reported in RedFly database along with their target genes are downloaded.

Data from Stark experiments, showing Vienna Tiles that showed enhancer expression and their matched gene targets is also available.


Names of genes in each gene cluster is located under data/geneClusters

-----------------------------------------------------------------------
Code:


Before running this script, you need to have your input files ready. All files should be in a fasta format, singe line. If your file is in a typical fasta format, where a line is up to 50bp, change it to a single line fasta. 
The sequence should be in a single line.
For example:

>DHS_ID_0001
CATTGGCATTGCCATCAAGACCCTCTGGACATAGGCACTGTTTCCTGTGATCGATCGTTTGACAGTGGGCATTGGTGCCGCAGGCAGTGGGATCTGCACATGGATCCACGCACTGCTGACCCACGCAGGACAATTCCGGCGGGCAGCCCTGATCGGACCGGCATCCAGGTACGCAGTTCAGGCCCAGACAGAGCTCTCCATGACCGCACTCATTGTCGTGGCGGCACAGGGGTTTGCAAACTCCCTGCTGGCATCGCTCGTTGGTCAAGCAGCCAGCATCGTCGGCGCATAGTGGTCGACATACCGACTCGAAGCAAGCCAATCCATTTCCACAGTCGCGATTCTCGCGGCATTCCAAGGGGGGAGAACGGACGCAGCCCACCTGGGGTGTGGGATTGGGCACCATACTCTCCAAACAACTGCAGCTAGCCCGGTGGTTGGACACCGAACAGGCGGCATTGGGACCACAGGGATTTTCCAGA

All files should be save in one directory. The following is a list with all the input files needed. 
1. File containing the DHSs for initialization. Those are the known examples, and we are quite certain that they regulate the genes. This should be names clusterNo_initialize_single.fa, like 15_initialize_single.fa

2. Files containing other DHSs that regulate other clusters. Those are the negative examples. Typically, those DHSs are known to regulate other clusters, and they are used for initialization for other clusters. So if you would like to learn the DHSs for cluster #15, and you would this against cluster 16 and 17. Then you would have two other files, one or 16 and one for 17. Their names should be 16_initialize_single.fa and 17_initialize_single.fa

3. File containing unlabeled DHSs, for which we would like to assign a label to. This would be named 15_unlabeled_single.fa

To run the program, write ./McEnhancer then supply these arguments in the following order
1. Directory where your files are stored. This will also be the output directory"
2. order of Markov chain  "-- after model selection, we decided to use 3, based on our dataset" 
3. Max. number of iterations for the EM algorithm to keep running. Best is 3"
4. Cluster number to process. This is the positive class --15 in our example"
5. Type of discount to be used by srlim. Example wd"
6. Value of the discount -- 0.1"
7. Difference between the probability estimates of the positive and negative. DHSs with probability bellow this cutoff, will be put in a reject class. These are generally DHSs that are on the borders"

Typically the code learns DHSs for one cluster (e.g cluster 15), against all other gene clusters from 1-18. You can definitely change this if you want.

After running this program, the output would be saved in a folder in the same directory called MC_15_output. Inside this folder, there will be a file for cluster 15 vs each of the other clusters. This shows the learnt DHSs for the positive class against each of the others. 

You should then find the overlapping DHSs, those that appeared more than a certain threshold.

You can then run the same script for stage2, but obviously adjust the input files

All the results for stage 1 are saved in results/ModelStage1.
Results for stage 2 are saved in results/ModelStage1/MC_multi_output/ModelStage2/

--------------------------------------------------------------------
Validation:

One way to validate the results is to run a classifier. 

Step 1:
To do this, first you will need to generate feature files for your positive and negative classes. Here, the feautre files are counts of 5-mers for all DHSs that are assigned to a specific gene, normalized by their lengths. 

To generate this feature files, you need to run "getKmer_normalize_byGene_Libsvm.py" located in "code/kmerCount/src/getKmer_normalize_byGene_Libsvm.py". This script expects these parameters in this order: 
1. path of input file
2. name of DHSs file (fasta single line, without the extension â€œ.faâ€)
3. output directory
4. no# kmers (we used 5)


Example:
python kmerCount/src/getKmer_normalize_byGene_Libsvm.py ../result/ModelStage1/ 15_initialize_single ../result/ModelStage1/ 5


Step 2:
For the classifier, we are basically using l1_logreg available at https://stanford.edu/~boyd/l1_logreg/index.html#download  
You can run the package on your generated files. 

Here we are providing a java program that does model selection, and cross validation. 
To use this Java program, you need to have l1_logreg installed, and the R library called â€˜ROCRâ€™
In this paper, for each gene cluster, we ran the classifier against each of the other clusters. The classifier will use the FILE.tab file generated in the previous step.
So assume your positive file is positive_single and your negative file is negative_file, you will call the classifier as in this example:

java -jar ClassifyKmerCounts_jar/ClassifyKmerCounts.jar Directory_where_files_are_saved positive_single negative_single /McEnhancer/code/




Amit's notebook-

White mold (Sclerotinia stem rot) is caused by a soilborne fungus, 
Sclerotinia sclerotiorum, one of the top ten devastating diseases of 
Soybean (Glycine max L.; Wrather et al., 2010). Sclerotinia 
sclerotiorum has a broad host range including over 400 host plants 
(Boland & Hall, 1994).challenging aspect to managing white mold is the 
diversity of the population of Sclerotinia sclerotiorum (Boland & Hall, 
1987; Tg et al., 2018). Despite being homothallic leading to what is 
considered a clonal population, there are chances of recombination leading 
to evolution of pathogen and subsequent diversity even within a field 
(Aldrich-Wolfe et al., 2015). The first step to find resistance to a 
complicated fungal population structure is to understand the genotypes of 
S. sclerotiorum present in the population. The objective of my study is 
to understand the population diversity of Sclerotinia sclerotiorum 
isolates collected from different regions of Midwest including Wisconsin, 
Ohio and Norrh Dakota. We have a total pf 100 isolates that we collcted 
from above mentioned states. For the trial, we are going to run it on 
three different isolates from wisconsin. For which I have extracted Dna 
and now have sent it for sequencing. The sequences will be done by short 
read illuma sequencing by UW-Biotech center. The data files we get after 
sequencing is fastq.gz. We have the refrence genome from NCBI database of 
sclerotinia sclerotiorum aligned to which the single nucleotide 
polymorphism will be checked.  After this, the next step would be to call of the varients using SNP Archer. Right now i am the process of figuring out 
how to download the SNPArcher. The NCBI sclerotinia refrence genome is 
38.9 MB with extension of .gna. 

Feburaray 20,2024 -I currently do not have any software downloaded as my computer says erroe"its new" so figuring out a way to download. I am going to work with assembling my genome to analyse the data. 

I do not have the alignmnets to work on now since I am working on the data. But i have got the gist of it. and perfomed this on the toy data.

Mrach 5, 2024

The packages that are required to install-
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)

 and loaded packages - library(ape)
library(adegenet)
library(phangorn)

loading the data set and Computing the genetic distances and then Get the NJ tree.
Before plotting, we can use the ladderize function which reorganizes the internal structure of the tree to get the ladderized effect when plotted
we can plot tree by using R command.
 ex-
 plot(tre, cex=.6)
title("A simple NJ tree")

and for parsimony menthod -
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)

and the load packages library(ape)
library(adegenet)
library(phangorn)

next is loading the sample data and convert it to phangorn object

 We need a starting tree for the search on tree space and compute the parsimony score of this tree 

example -
tre.ini <- nj(dist.dna(dna,model="raw"))
parsimony(tre.ini, dna2)

Now Search for the tree with maximum parsimony:

example-
> tre.pars <- optim.parsimony(tre.ini, dna2)
Final p-score 420 after  2 nni operations

Last- Plot tree

example- 
plot(tre.pars, cex=0.6)

Choose a maximum likelihood method that you like the best on your project 
dataset (it does not have to be RAxML or IQ-Tree, but you should read the 
paper and tutorial of whichever method you choose)


March 21, 2024

Although there are other software packages RAxML (Randomized Axelerated 
Maximum Likelihood) and IQ-TREE are both widely used software packages for phylogenetic analysis, especially useful 
in constructing phylogenetic trees based on genetic sequence data. When it 
comes to analyzing SNP (Single Nucleotide Polymorphism) variants of a 
particular organism, such as Sclerotinia sclerotiorum, the choice between 
RAxML and IQ-TREE depends on several factors including computational 
efficiency, model selection, user interface, and the specific goals of 
your phylogenetic analysis. RAxML has high computational efficacy, good 
performace abd better flexibility while it lacks in bulit model selection. 
while IQ-tree offers the above including built-in model selection tool, 
which can automatically determine the best-fitting model of evolution for 
your data. This is particularly useful for SNP data of Sclerotinia 
sclerotiorum as selecting an appropriate model can significantly impact 
the accuracy of the phylogenetic tree.
or analyzing SNP variants of Sclerotinia sclerotiorum would likely lean 
towards IQ-TREE due to its integrated model selection capability, which is 
crucial for SNP data analysis to ensure the accuracy of the phylogenetic 
tree. The automatic model selection simplifies the workflow and can lead 
to better-supported phylogenetic conclusions. Additionally, IQ-TREE's 
user-friendly approach and flexibility in handling various models make it 
an attractive option for researchers working on SNP data.

After SNP aclling and alignment of data. we choose IQ-Tree and then next 
step is to run it with command iqtree -s mydata.phy -m TEST. where Replace 
mydata.phy with the path to your data file. This command tells IQ-TREE 
to select the best-fit model automatically and then use that model to 
construct the phylogenetic tree.To assess the reliability of the branches 
in my phylogenetic tree, I might perform bootstrapping. Next step would be 
to interpret the result. 

April 11, 2024-  Performing Bayesian method in my data

*Downloaded  MrBayes from here-
brew reinstall mrbayes
sudo chown -R $(whoami) /usr/local/Cellar/open-mpi/4.1.0
brew reinstall mrbayes

* Once I have the alignment file (FASTA) for my data, I will have to convert it before running Mr. Bayes, since it only reads nexus files. 

* With the nexus file saved on my laptop, I will next create a mrbayes block in a separate txt file. In this file I will put the following items and save it: 

"begin mrbayes;
 set autoclose=yes;
 prset brlenspr=unconstrained:exp(10.0);
 prset shapepr=exp(1.0);
 prset tratiopr=beta(1.0,1.0);
 prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);
 lset nst=* rates=gamma ngammacat=4;
 mcmcp ngen=1,000,000 samplefreq=10 printfreq=100 nruns=1 nchains=3 savebrlens=yes;
 outgroup Anacystis_nidulans;
 mcmc;
 sumt;
end;" 

**** This is a file to insert the prior information to create the desired trees. For my data notice that I put the substitution model with a star  (not = *) because I will decide what model to use based on the one outputted by IQ Tree, once I know what model I will use I will do the proper * substitution to inform what model I will proceed using. Note that I also changed the generation number to 1,000,000 (ngen=1,000,000). For the rest of the lines, you noticed that in this case, I did not change the standard parameters because I don't have previous knowledge about my data to change the priors, so I will go ahead using the same. ****


* Next I have to "merge" my block file with my nex data file using the code shown below as an example: 
"cat algaemb.nex mbblock.txt > algaemb-mb.nex"

* Once merged, I can go ahead and run the command mb to run mrBayes on my "-mb.nex" file, like in the example below:

"mb algaemb-mb.nex"

* I can save the output in different formats, but for my purposes, I will choose to save the output of this proposed tree as .tre, so that I can proceed to alter its layout using software like R, for example, with R packages that allow me to manipulate phylogenetic trees.

April23,2024

In my work for the class project I am not going to do work based on coalescent methods, but it is good to learn about the software, procedure and codes for it
There are many computational tools and approaches available, depending on the specifics of your project and the nature of your data. 
For a project that requires modeling species and gene trees under a coalescent process, one widely used tool is BEAST (Bayesian Evolutionary Analsis Sampling Trees)
we learned about this in class.

BEAST is know for its excellent for analyzing complex evolutionary scenarios that involve sequence data from multiple genes.
 It has many advantages like in the estimation of phylogenies with explicit models of evolution, incorporating various molecular 
 clocks and coalescent processes, infering past-population dynamics and to deal with different types of sequence data and multiple
loci in a statistically robust Bayesian framework.

here’s a generic walkthrough of the steps involved-
Download and install BEAST 
 Format the sequence data in an XML file, which includes sequences, models of evolution, and settings for the MCMC analysis. This can be done using BEAUti (Bayesian Evolutionary Analysis Utility), which is part of the BEAST software suite.
 Configure the Coalescent Model: Within BEAUti:
fisrly Load the data, Set the substitution model, molecular clock model, and tree priors under the assumption of a coalescent process.
Choose from different coalescent models like constant size, exponential growth, or skyline plots based on your hypothesis about population dynamics
Once your model is specified, generate the XML input file for BEAST.
Next step is to run the beast and the cide is as-
beast -beagle -beagle_SSE -threads auto ConfigFile.xml
After BEAST finishes running, analyze the output using programs like Tracer (to examine MCMC trace files) and TreeAnnotator (to summarize the posterior sample of trees).

like for the lab we were to perfom this 
1. Prepared data in BEAUti and set up a coalescent model of constant population size.
2. Exported analysis settings to XML file named 'coalescent_toy_example.xml'.
3. Ran BEAST with the following command:
   beast -beagle -beagle_SSE -threads auto coalescent_toy_example.xml
4. Monitored run progress and convergence in Tracer.
5. Summarized posterior trees using TreeAnnotator with a burn-in of 10%.

I also read from online sourecs to write this information.
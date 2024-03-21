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

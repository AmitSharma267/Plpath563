Amit's notebook

workin with fungi sclerotinia sclerotiorum 
Currently I have 3 isolates from which i extracted DNA
now going to assemble genome.

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



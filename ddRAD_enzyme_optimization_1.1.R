##############################################
## T. Thurman                               ##
## Feb. 2016                                ##
## Using SimRAD to optimize enzyme choice   ##
## and size selection for ddRAD sequenceing ##
##############################################

###################################################################################
# HOW TO USE THIS SCRIPT
#
# This script requires 3 files
# 1) A reference genome for your organism, in FASTA format
#
# 2) A csv file of possible restriction enzymes, with their recognition sequences defined
# according to the convention in the SimRAD::insilico.digest help file. If you're using this script,
# you've probably got a copy of the "enzyme_list.csv" file. It works.
#
# 3) The ddRAD_opt_source.R file. It contains all the functions that run the analysis. If you're
# interested in figuring out the details, check it out.
#
# Besides that, the analysis requires the SimRAD R package and the dplyr package.
# Make sure you have them installed.
#
######################################################################################

########################
## Prep. for Analysis ##
########################

## Load the analysis functions into R
source("ddRAD_opt_source_1.1.R")

## Select genome file. File must be in FASTA format.
## Enter the filepath to the genome file below
genome.file <- "/Users/tthurman/Desktop/Gasterosteus_aculeatus.BROADS1.dna_sm.toplevel.fa"

## Select restriction enzyme file. File must be a .csv.
## Enter the filepath to the restriction enzyme file below
RE.file <- "enzyme_list.csv"

## List possible size-selection fragment sizes. Must be a vector of numerical values, in # of base pairs.
pos.frag.size <- c(350, 450)

## Specify the proportion of contigs you'd like to use for the simulations of RAD data.
## The more contigs you include, the more accurate the estimate will be, but the longer your simulations
## will take. Its best to analyze many enzymes at a low proportion (e.g., prop.contigs = 0.1) to find a
## few candidates for the best pairs, and then to run only one or a few enzymes with all contigs
## (prop.contigs = 1) to get a good estimate. If your genome is large enough (e.g., the sagrei genome),
## it will be too large to store as a single string in R, and you'll just have to use as high
## a proportion of contigs as possible.
prop.contigs <- 0.01


############################
## Analyze # of fragments ##
############################

## The "testEnzymes" function does all the analysis about # of fragments.
## By default, it will use the genome file, restriction enzyme file,
## list of possible fragments, and proportion of contigs that you specify above.
## It has 2 arguments you must define:

## 1) testAll. can be true or false (T/F). If T, will test all possible enzyme pairs in the list of
## restriction enzymes. If F, it will only test the enzyme pairs that you tell it to with the
## "selectedEnzymes" argument.

## 2) selectedEnzymes. Only need to specify this if testAll = F. This must be a list (the R class)
## of pairs of enzymes you wish to test. For example: list(c(3,4), c(5,6)) will test two sets
## of enzymes: the enzyme in row 3 with the enzyme in row 4, and the enzyme in row 5 with the
## enzyme in row 6. Be careful when counting rows: R does not import the header row, so the first
## enzyme is in row 1 of the R data frame, but in the second row when you look at the .csv in excel.

## Below is are two example function calls, one with enzymes specified, and one without

optimization.results <- testEnzymes(testAll = T)

optimization.results <- testEnzymes(testAll = F,
                                    selectedEnzymes = list(c(1,5), c(2,4), c(2,5)))


## It is best to write these results to a .csv so you don't lose them.
write.csv(optimization.results, file = "optimization_results.csv")


############################
## Put the # of fragments ##
## into a genomic context ##
############################

# Re-load the results from above, if necessary
optimization.results <- read.csv("optimization_results.csv")

## Now, we'll take the data on # of fragments and determine what that means for sequencing effort
## and genome coverage. These estimates are extrapolated to the genome level from the
## proportion of contigs sampled. They will be most accurate when prop.contigs = 1.

## The genomeStats function will add to the results dataframe two columns:

## 1) bp.seq.per.ind- The estimated # of base-pairs sequenced per individual (assuming 1x coverage).
## This is extrapolated from the prop.contigs #, and this number will be more accurate as prop.contigs
## approaches one

## 2) percent.genome- the estimated % of the genome contained in the sequenced portion of the
## estimated fragments

## The genomeStats function requires 4 arguments:

## 1) results- The results data frame generated above

## 2) paired.end- TRUE or FALSE- is sequencing paired-end or single-end?

## 3) seq.length- How many base pairs are being sequenced per read? Genome Quebec offers
## 50, 100, 150, or 250 on HiSeq and up to 300 on MiSeq (as of their website, checked Feb. 5 2016)

## 4) genome.size- The size of the genome, in base pairs. Helpful hint- gigabases = x10^9,
## megabases = 10^6

## An example function call, for single-end sequencing of 100 bp reads on the
## Anolis carolinensis genome (1.78 Gb) is shown below

optimization.results <- genomeStats(results = optimization.results,
                                    paired.end = F, seq.length = 100, genome.size = 1.78e9)

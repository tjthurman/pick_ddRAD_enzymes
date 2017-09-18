############################################## 
## T. Thurman                               ##
## Feb. 2016                                ##
## Using SimRAD to optimize enzyme choice   ##
## and size selection for ddRAD sequenceing ##
##############################################

## This is the source file, which contains the functions that actually do the analysis.
## Take a look, but be careful about changing things and breaking the script...

######################
## load libraries
library(SimRAD)
library(dplyr)



## testEnzymes is the top-level function.
## It checks some arguments, imports the genome and enzyme files,
## creates the list of enzyme pairs to be tested,
## and then runs a for loop on that list, calling
## the perEnzymePairSimulation function, which does the ddRAD simulation 
## for a pair of enzymes

testEnzymes <- function(genome = genome.file, enzymeFile = RE.file, testAll = F, selectedEnzymes = NULL,
                        propContigs = prop.contigs, possibleFragmentSizes = pos.frag.size) {
  ## First, test args, and define enzyme.pairs if the user has specified it.
  if (testAll == F) {
    if (is.null(selectedEnzymes) == T) {
      stop("If testAll = F, must provide a list of enzyme pairs!")
    }
    if (is.list(selectedEnzymes) == F) {
      stop("If testAll = F, the set of selected enzyme must be provided as a list of pairs!")
    }
    enzyme.pairs <- selectedEnzymes
  }
  
  # Load in the genome file, using only the specified proportion of contigs
  imported.genome <- ref.DNAseq(FASTA.file = genome.file, subselect.contigs = T, prop.contigs = propContigs)
  print("Genome imported")
  
  # load in the list of restriction enzymes
  re.list <- read.csv(enzymeFile, colClasses = c(rep("character", times = 3), NULL))
  
  # create vector of all possible comparisons, if user has specified they want to test all possible pairs
  if (testAll == T) {
    enzyme.pairs <- combn(x = 1:dim(re.list)[1], m = 2, simplify = F)
  }
  
  
  ## Then, set-up and run the for loop that will create the results
  total.results <- NULL
  for (pair in enzyme.pairs) {
    pair.results <- perEnzymePairSimulation(imported.genome, possibleFragmentSizes, enzymePair = as.vector(pair), RE.list = re.list)
    total.results <- rbind(total.results, pair.results)
  }
  rm(enzy1_3pSide, enzy1_5pSide, enzy1_name, enzy2_3pSide, enzy2_5pSide, enzy2_name, pos = 1)
  return(total.results)
}

## perEnzymePairSimulations calls other functions to estimate the # of fragments for a given enzyme pair
## across different size selection parameters
perEnzymePairSimulation <- function(imported.genome, possibleFragmentSizes, enzymePair, RE.list= re.list) {
  prepareEnzymeInfo(enzymePair, RE.list)
  x <- simulateFrags(imported.genome)
  y <- sizeSelectionDF(x, possibleFragmentSizes)
}

## prepareEnzymeInfo read the enzyme list csv for the cut recognition sites and names
## of the enzymes
prepareEnzymeInfo <- function(enzymePair, RE.list) {
  assign("enzy1_5pSide", value = as.character(RE.list[enzymePair[1], 2]), pos = 1)
  assign("enzy1_3pSide", value = as.character(RE.list[enzymePair[1], 3]), pos = 1)
  assign("enzy2_5pSide", value = as.character(RE.list[enzymePair[2], 2]), pos = 1)
  assign("enzy2_3pSide", value = as.character(RE.list[enzymePair[2], 3]), pos = 1)
  assign("enzy1_name", value = as.character(RE.list[enzymePair[1], 1]), pos = 1)
  assign("enzy2_name", value = as.character(RE.list[enzymePair[2], 1]), pos = 1)
}

## simulateFrags uses SimRAD functions to similate a double digest of the genome
## and to simulate the adaptor selection step in ddRAD. Returns a list of fragments
simulateFrags <- function(imported.genome) {
  # Simulate the digest
  digested.genome <- insilico.digest(imported.genome, enzy1_5pSide, enzy1_3pSide, enzy2_5pSide, enzy2_3pSide, verbose = F)
  # Simulate the adaptor selection step in ddRAD
  adaptor.fragments <- adapt.select(digested.genome, type = "AB+BA", enzy1_5pSide, enzy1_3pSide, enzy2_5pSide, enzy2_3pSide)
  print(c("# of ddRAD fragments = ", length(adaptor.fragments)))
  return(adaptor.fragments)
}

## sizeSelectionDF takes a list of fragments and a vector of possible fragment sizes as inputs,
## then creates a data frame with the results of how many fragments are expected per size selection parameters. 
## In reality, a lot of these combinations will be non-sensical (won't be doing size selection over more than 80-100 bp),
## to speed things up, may wish to modify the size selection list stuff so that you give it pairs of possible
## fragments.
sizeSelectionDF <- function(fragments, possibleFragmentSizes) {
  results.df <- matrix(nrow= length(possibleFragmentSizes)^2, ncol = 5)
  for (i1 in 1:length(possibleFragmentSizes)) {
    for (i2 in 1:length(possibleFragmentSizes)) {
      if (i1 < i2) {
        entry <- i1 * i2
        results.df[entry, 1] <- as.character(enzy1_name)
        results.df[entry, 2] <- as.character(enzy2_name)
        results.df[entry, 3] <- as.character(possibleFragmentSizes[i1])
        results.df[entry, 4] <- as.character(possibleFragmentSizes[i2])
        results.df[entry, 5] <- as.character(length(size.select(fragments, min.size = possibleFragmentSizes[i1], 
                                                   max.size = possibleFragmentSizes[i2], graph = F, verbose = F)))
      }
    }
  }
  results.df <- as.data.frame(results.df, stringsAsFactors = F)
  colnames(results.df) <-  c("enzyme.1", "enzyme.2","min.size.bp", "max.size.bp", "num.of.fragments")
  results.df$num.of.fragments <- as.numeric(results.df$num.of.fragments)
  results.df$min.size.bp <- as.numeric(results.df$min.size.bp)
  results.df$max.size.bp <- as.numeric(results.df$max.size.bp)
  results.df <- filter(results.df, !is.na(num.of.fragments))
  return(results.df)
}

############################
## Put the # of fragments ##
## into a genomic context ##
############################

## This function is is for analysis of # of fragments in the genomic context. 
## Uses dplyr to add two columns to the results dataframe generated by the above functions.

## 1) bp.seq.per.ind- The estimated # of base-pairs sequenced per individual (assuming 1x coverage).
## This is extrapolated from the prop.contigs #, and this number will be more accurate as prop.contigs
## approaches one

## 2) percent.genome- the estimated % of the genome contained in the sequenced portion of the 
## estimated fragments


genomeStats <- function(results, propContigs = prop.contigs, paired.end, seq.length, genome.size) {
  ## First, check some arguments
  if (paired.end == T) {
    multiplier <- 2
  } else if (paired.end == F) {
    multiplier <- 1
  } else {
    stop("paired end must be either T or F")
  }
  if (is.numeric(seq.length) == F) {
    stop("seq.length must be a numeric value")
  }
  if (genome.size < 1000) {
    warning("Check genome.size argument. Is the genome really that small?")
  }
  ## Check to see if the amount of sequenced DNA is bigger than the max fragment size
  ## if so, just use the max fragment size as the number of bp sequenced
  
  ## Problem- need to do this line by line with a for loop, can't use dplyr anymore, 
  ## as I have to check the value of max.size.bp for each line. 
  
  ## use a for loop for the bp.seq.per.ind value.
  for (row in 1:dim(results)[1]) {
    if (as.numeric(results$max.size.bp[row]) < seq.length*multiplier) {
      results$bp.seq.per.ind[row] <- as.numeric(results$num.of.fragments[row]/propContigs*results$max.size.bp[row])
    } else {
      results$bp.seq.per.ind[row] <- as.numeric(results$num.of.fragments[row]/propContigs*seq.length*multiplier)
    }
  }

  ## can still do the percent.genome with dplyr
  results <- mutate(results, percent.genome = round((bp.seq.per.ind/genome.size)*100, digits = 3))
  ## round the bp.seq.per.ind to something interpretable after the calculations are done
  results$bp.seq.per.ind <- as.numeric(format(results$bp.seq.per.ind, digits = 3))
  return(results)
}


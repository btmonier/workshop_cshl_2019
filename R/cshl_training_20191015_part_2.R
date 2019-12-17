#--------------------------------------------------------------------
# Script Name:   cshl_training_20191015_part_2.R
# Description:   CSHL Workshop 2019 - Part 2
# Author:        Brandon Monier and Guillaume Ramstein
# Created:       2019-10-08 at 14:15:16
# Last Modified: 2019-10-16 at 15:31:25
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to layout all possible
#    functions for the CSHL Workshop
#--------------------------------------------------------------------

# === Part 1 - PREAMBLE =============================================

### See `part_1` script for initial analyses...

## Load Java parameters ---
options(java.parameters = c("-Xmx10g")) # how much extra memory do you want?


## Load packages ---
library(dplyr)                # A grammar of data manipulation
library(foreach)              # Provides Foreach Looping Construct
library(magrittr)             # Pipe operator for R
library(readr)                # Read in data as tibble classes
library(rTASSEL)              # TASSEL in R!
library(SummarizedExperiment) # SummarizedExperiment container and methods
library(tibble)               # Simple data frames



# === Part 2 - rTASSEL DATA OBJECTS =================================

## Add family column to BLUEs data ---
phenoFamilyDF <- tasBLUE$BLUE %>%
    dplyr::mutate(
        family = gsub(
            pattern = "E.*",
            replacement = "",
            x = .$Taxa
        )
    ) %>%
    dplyr::mutate(family = factor(family)) %>%
    dplyr::select(Taxa, family, dplyr::everything())


## Convert to rTASSEL object ---
tasFamilyGenoPhenoBLUE <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj = tasGeno,
    phenoPathDFOrObj = phenoFamilyDF,
    taxaID = "Taxa",
    attributeTypes = c("factor", rep("data", 2))
)



# === Part 4 - CROSS VALIDATION =====================================

## Run and create cross validation object report ---
system.time(
    tasCV2Family <- rTASSEL::genomicPrediction(
        tasPhenoObj = tasFamilyGenoPhenoBLUE,
        kinship = tasKin,
        doCV = TRUE,
        kFolds = 5,
        nIter = 1
    )
)



# === Part 5 - LEAVE ONE FAMILY OUT CROSS VALIDATION ================

## (1) Create empty data frame object ---
LOFO <- data.frame()

t0 <- Sys.time()


## (2) Iterate through all families and leave one out (set to NA)
for (family in levels(phenoFamilyDF$family)) {

    message("LOFO: ", family)
    familyTaxa <- phenoFamilyDF %>%
        dplyr::filter(family == family) %>%
        dplyr::pull(Taxa) %>%
        as.factor()

    ## Converted iterated family to NA
    phenoDF_NA <- phenoFamilyDF
    phenoDF_NA[phenoDF_NA$family == family, sapply(phenoDF_NA, is.numeric)] <- NA

    ## Make rTASSEL phenotype object
    pheno_TASSEL <- rTASSEL::readPhenotypeFromDataFrame(
        phenotypeDF = phenoDF_NA[, c(1, 3, 4)],
        taxaID = "Taxa"
    )

    ## Genomic prediction
    GP <- rTASSEL::genomicPrediction(
        tasPhenoObj = pheno_TASSEL,
        kinship = tasKin,
        doCV = FALSE
    )

    ## Calculate accuracy
    out <- foreach::foreach(trait = unique(GP$Trait), .combine = rbind) %do% {

        GP_trait <- GP[GP$Trait == trait, ]

        r <- cor(
            phenoFamilyDF[match(familyTaxa, phenoFamilyDF$Taxa), trait],
            GP_trait[match(familyTaxa, GP_trait$Taxon), "Predicted"],
            use = "complete"
        )

        tibble::tibble(
            Trait    = trait,
            Family   = family,
            Accuracy = r
        )

    }

    LOFO <- rbind(LOFO, out)
    names(LOFO) <- c("Trait", "Family", "Accuracy")

}

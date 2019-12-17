#--------------------------------------------------------------------
# Script Name:   cshl_training_20191015_part_1.R
# Description:   CSHL Workshop 2019 - Part 1
# Author:        Brandon Monier and Guillaume Ramstein
# Created:       2019-10-08 at 14:15:16
# Last Modified: 2019-10-15 at 15:42:42
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to layout all possible
#    functions for the CSHL 2019 Workshop.
#--------------------------------------------------------------------

# === Part 1 - PREAMBLE =============================================

## Load Java parameters ---
options(java.parameters = c("-Xmx10g")) # how much extra memory do you want?


## Load packages ---
library(dplyr)                # A grammar of data manipulation
library(magrittr)             # Pipe operator for R
library(readr)                # Read in data as tibble classes
library(rTASSEL)              # TASSEL in R!
library(SummarizedExperiment) # SummarizedExperiment container and methods


## Start logging file ---
rTASSEL::startLogger(fullPath = "/home/bm646/Projects/cshl_training_2019")


## Get input variables ---

### Path to data directory
inputDir  <- "/home/bm646/Projects/cshl_training_2019/data/"

### Phenotype file
phenoFile <- paste0(inputDir, "phenotypic_data.csv")

### Genotype file
genoFile  <- paste0(inputDir, "AGPv4_NAM_subset.recode.vcf.gz")



# === Part 2 - rTASSEL DATA OBJECTS =================================

## Load and inspect genotype data ---

### Load data
tasGeno <- rTASSEL::readGenotypeTableFromPath(path = genoFile)

### Show data display on console
tasGeno

### Get data structure slot names
methods::slotNames(tasGeno)

### Inspect genotype slot - what does it show?
tasGeno@jGenotypeTable

### Get object size
object.size(tasGeno)


## Phenotype data ---

### Load `tibble`-class data frame
phenoDF <- readr::read_csv(file = phenoFile)

### Convert family column to factor
phenoDF <- phenoDF %>%
    dplyr::mutate(Location = factor(Location)) %>%
    dplyr::mutate(Year = factor(Year))

### Inspect data set - R
phenoDF

### Load data frame into rTASSEL data object
tasPheno <- rTASSEL::readPhenotypeFromDataFrame(
    phenotypeDF = phenoDF,
    taxaID = "Taxon",                                    # <- talk about this
    attributeTypes = c(rep("factor", 2), rep("data", 2)) # <- talk about this
)


## Create GenotypePhenotype object
tasGenoPheno <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj = tasGeno,
    phenoPathDFOrObj = tasPheno
)



# === Part 3 - ASSOCIATION and KINSHIP ANALYSES =====================

## BLUEs ---
tasBLUE <- rTASSEL::assocModelFitter(
    tasObj     = tasGenoPheno,
    formula    = . ~ Location + Year,
    fitMarkers = FALSE
)


## Inspect BLUEs ---
tasBLUE


## Convert BLUE output to rTASSEL data object ---
tasGenoPhenoBLUE <- rTASSEL::readGenotypePhenotype(
    genoPathOrObj    = tasGeno,
    phenoPathDFOrObj = tasBLUE$BLUE,
    taxaID           = "Taxa"
)


## Make TASSEL kinship object ---
tasKin <- rTASSEL::kinshipMatrix(tasObj = tasGenoPhenoBLUE)


## Create kinship matrix object (R) ---
## NOTE: relatively big object!
# tasKinR <- rTASSEL::kinshipToRMatrix(tasKin)


## GWAS - MLM ---
tasMLM <- rTASSEL::assocModelFitter(
    tasObj     = tasGenoPhenoBLUE,
    formula    = . ~ 1,
    fitMarkers = TRUE,
    kinship    = tasKin
)



# === Part 4 - VISUALIZE DATA =======================================

## Days to silking
manhattanDTS <- tasMLM$MLM_Stats %>%
    dplyr::filter(Trait == "DaysToSilk") %>%
    dplyr::filter(Marker != "None") %>%
    ggplot() +
    aes(x = Pos, y = -log10(p)) +
    geom_point(size = 0.7) +
    ggtitle(label = "Days to Silk") +
    xlab("SNP Position") +
    ylab(bquote(~-log[10]~ '('*italic(p)*'-value)')) +
    facet_grid(. ~ Chr, scales = "free_x") +
    theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )


## Ear weight
manhattanEW <- tasMLM$MLM_Stats %>%
    dplyr::filter(Trait == "EarWeight") %>%
    dplyr::filter(Marker != "None") %>%
    ggplot() +
    aes(x = Pos, y = -log10(p)) +
    geom_point(size = 0.7) +
    ggtitle(label = "Ear Weight") +
    xlab("SNP Position") +
    ylab(bquote(~-log[10]~ '('*italic(p)*'-value)')) +
    facet_grid(. ~ Chr, scales = "free_x") +
    theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

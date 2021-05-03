# global.R for variantVieweR_shiny_app

#----- Libraries -----#
source("./scripts/snakemake_helpers/snakemake_helpers.R")

requiredPackages <- c("shiny", "shinythemes", "shinyjs",
                      "DT", "data.table", "tidyverse")
for (package in requiredPackages) {
  TryInstall(package)
  
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    stop(cat("FATAL: Problem loading R package:", package, 
             "(global.R)\n\n"),
         call. = FALSE)
  }
}

requiredBioCPackages <- c("BiocManager", "Gviz")

for (package in requiredBioCPackages) {
  TryInstallBioconductor(package)
  
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    stop(cat("FATAL: Problem loading R package:", package, 
             "(global.R)\n\n"),
         call. = FALSE)
  }
}

#library("shiny")
#library("shinythemes")
#library("shinyjs")
#library("DT")
#library("Gviz")
#library("data.table")
#library("tidyverse")
#library("BiocManager")

#----- Global variables -----#

# Find results directory and dataSet.txt from config file
resultsPath <- readLines("VirusVariantViewR-RShiny-Application_config.txt",
                         warn = FALSE, n = 1)

dataSetFilePath <- paste0(resultsPath, "/dataSet.txt")

# Data set name from VirusVariantViewR_datasets.txt
dataSet <- readLines(dataSetFilePath, warn = FALSE)

previousMtxSize <- 0 # for determining when to re-paint the coverage plot

genome_size <- 7383 # For calculation/formatting avg. genome coverage

#----- Class Definitions -----#

# These classes are populated when the data set is selected

setClass("Sample",
         representation = representation(sample_name = "character",
                                         variant_list = "list"),
         prototype = prototype(sample_name = NA_character_,
                               variant_list = list()))

setClass("Variant",
         representation = representation(parent_sample = "character",
                                         position = "numeric",
                                         ref_allele = "character",
                                         alt_allele = "character",
                                         variant_df = "data.frame"))

#----- Function Definitions -----#

#----- Populating the Sample and Variant Objects -----#

# Returns a list of Sample objects, each of which contains a list of
#   Variant objects.
#dataSetSelect <- "orchard_CW3_pooled"
#testDataSet <- GenerateSampleData(dataSetSelect) # used to make dataSetSamples
#dataSetSamples <- as.character(testDataSet$Sample) # send to sampleVector

BuildSampleObjects <- function(dataSet, sampleVector) {
  
  # dataSet will come from input$dataSetSelect
  # sampleVector should automatically populate with all the samples available
  #   in the selected data set
  sampleClassList <- vector(mode = "list", length = length(sampleVector))
  
  # Build a Sample class object from the vector
  for (i in 1:length(sampleVector)) {
    currentSample <- sampleVector[i]
    # populate a Variant class object
    currentVCFFile <- GetVCF(dataSet, sample = currentSample)
    # Stop here if nothing is in the sample
    if (nrow(currentVCFFile) == 0) {
      warning(paste(currentSample, "has no variants detected"),
              call. = FALSE)
      newSample <- new("Sample",
                       sample_name = currentSample,
                       variant_list = list(0))
      sampleClassList[[i]] <- newSample
      names(sampleClassList)[i] <- currentSample
    } else {
      currentSampleVariantList <- vector(mode = "list", length = nrow(currentVCFFile))
      
      for (j in 1:nrow(currentVCFFile)) {
        positionString <- as.character(currentVCFFile[j, "Position"])
        referenceString <- as.character(currentVCFFile[j, "Reference"])
        alternativeString <- as.character(currentVCFFile[j, "Alternative"])
        variantID <- paste(positionString, referenceString, alternativeString,
                           sep = "_")
        variantReadableID <- paste(positionString, " ", referenceString, "/", alternativeString, sep = "")
        variantDF <- data.frame("Variation" = variantID,
                                "Variation_Readable" = variantReadableID,
                                "Reference_Allele" = referenceString,
                                "Reference_Codon" = currentVCFFile[j, "Reference Codon"],
                                "Reference_Protein" = currentVCFFile[j, "Reference Protein"],
                                "Mutant_Allele" = alternativeString,
                                "Mutant_Codon" = currentVCFFile[j, "Mutant Codon"],
                                "Mutant_Protein" = currentVCFFile[j, "Mutant Protein"],
                                "Mutation_Type" = currentVCFFile[j, "Mutation Type"])
        newVariant <- new("Variant",
                          parent_sample = currentSample,
                          position = as.numeric(positionString),
                          ref_allele = referenceString,
                          alt_allele = alternativeString,
                          variant_df = variantDF)
        currentSampleVariantList[[j]] <- newVariant
        names(currentSampleVariantList)[j] <- variantID
      }
      
      # Fill in the sample object
      newSample <- new("Sample",
                       sample_name = currentSample,
                       variant_list = currentSampleVariantList)
      # Add to sampleClassList
      sampleClassList[[i]] <- newSample
      names(sampleClassList)[i] <- currentSample
      #sampleClassList[[currentSample]] <- newSample
    }
    
  } 
  return(sampleClassList)
}

#----- Set up Sample Data Table -----#

# Alignment Count Data #
GenerateSampleData <- function(dataSet) {

  # Read in sample alignment count data, calculate percent MNV
  readCounts <- read.delim(paste0(resultsPath, "/", dataSet, "/alignment_counts.txt"),
                           colClasses = c("character", "numeric", "numeric", "numeric"))
  readCounts <- readCounts %>%
    mutate("percent_MNV" = round(100*((left_alignments + right_alignments)/total_reads), 
                                 digits = 2))
  
  # Genome Coverage Data #
  # Read in genome coverage count data, calculate avg. genome coverage
  rawFiles <- list.files(path = paste0(resultsPath, "/", dataSet, "/sample_data/genome_coverage/"), 
                         pattern = "*_coverage.txt")
  
  # column names determined from documentation at 
  # https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
  columnNames <- c("chromosome", "depth", "number_of_bases",
                   "chromosome_size", "fraction_of_bases")
  
  # Run calculation, store in genomeCovVec - will be combined with readCounts 
  genomeCovVec <- vector(mode = "integer", length = length(rawFiles))
  
  for (i in 1:length(rawFiles)) {
    
    currentFile <- rawFiles[i]
    currentFileName <- as.vector(strsplit(currentFile, "_coverage.txt")[[1]])
    
    currentData <- read.delim(paste0(resultsPath, "/", dataSet, 
                                     "/sample_data/genome_coverage/", 
                                     currentFile), 
                              header = FALSE,
                              col.names = columnNames)
    currentData <- currentData %>%
      dplyr::select(depth, number_of_bases) %>%
      dplyr::mutate(product = depth * number_of_bases)
    
    avg_cov <- round(sum(currentData$product)/genome_size,
                     digits = 0)
    
    genomeCovVec[i] <- avg_cov
    names(genomeCovVec)[i] <- currentFileName
    
  }
  
  # Try to Merge Metadata #
  # Add the average genome coverage numbers to the readCounts data
  sampleData <- readCounts %>%
    dplyr::mutate("avg_genome_cov" = genomeCovVec[sample]) %>%
    dplyr::select(sample, total_reads, percent_MNV, avg_genome_cov)
  names(sampleData) <- c("Sample", "Total Reads", "% MNV", "Average Coverage")
  sampleData$Sample <- as.character(sampleData$Sample)
  
  # Do we need to add metadata?
  #   Check for the existence of a metadata file.  If it exists, read-in,
  #   merge with the sample data and return.
  metadataFile <- paste0(resultsPath, "/", dataSet, "/", dataSet, "_metadata.txt")
  
  if (file.exists(metadataFile)) {
    
    metadata <- read.delim(metadataFile, check.names = FALSE,
                           colClasses = "character")
    
    sampleDataExtended <- tryCatch({ # in case the Sample column isn't correct 
      merge(metadata, sampleData, by = "Sample")},
      error = function(e) {
        # just return the un-merged sample data
        sampleData
      })
    return(sampleDataExtended)
    
  } else {
    return(sampleData)
  }

}


#----- Generate Sample Variant Table -----#

GetVCF <- function(dataSet, sample) {
  # Read in and format VCF file for selected sample
  vcfFile <- data.table::fread(file = paste0(resultsPath, "/", dataSet, 
                                             "/variants/annotated_variants/", 
                                             sample, "_variants_annotated.txt"),
                               header = TRUE, sep = "\t", 
                               stringsAsFactors = FALSE,
                               verbose = FALSE)
  
  # Get the primary alignment value from the current sample's alignment counts
  sampleAlignmentCounts <- GenerateSampleData(dataSet)
  samplePrimaryAlignments <- subset(sampleAlignmentCounts,
                                    Sample == sample)$`Primary Alignments`
  # Filter sampleAlignmentCounts by sample
  #samplePrimaryAlignments <- subset(GenerateSampleData,
                                    #sample == sample)$primary_alignments

  # basic columns to display, even if VCF file is empty:
  keepCols <- c("Reference Genome", "Position", "Reference",
                "Alternative", "Quality")

  # additional columns to display in cases where the VCF file is not empty:
  if (nrow(vcfFile) != 0) {
    keepCols <- c(keepCols, "Location", "Reference Codon", "Reference Protein", 
                  "Mutant Codon", "Mutant Protein", "Mutation Type", "Total Depth",
                  "Allelic Frequency (%)") 
  }

  vcfFileFormatted <- vcfFile %>%
    select(keepCols)
  
  comment(vcfFileFormatted) <- sample
  
  return(vcfFileFormatted)
  
}

#dataSetSelect <- "Combined_Data"
#testDataSet <- GenerateSampleData(dataSetSelect) # used to make dataSetSamples
#dataSetSamples <- as.character(testDataSet$Sample) # send to sampleVector
#vcfTest <- GetVCF(dataSet = dataSetSelect, sample = dataSetSamples[12])

#----- Generate Coverage Plot for Sample & Variants -----#

PlotCoverage <- function(dataSet, sample, positions = NULL, widths = 1) {

  # Get the coverage file
  # Colnames taken from 
  # https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
  bedgraphDT <- fread(paste0(resultsPath, "/", dataSet, "/alignment_files/", 
                             sample, "_sorted.bedGraph"),
                      col.names = c("chromosome", "start", "end", "value"))

  if (nrow(bedgraphDT) == 0) {
    stop("No coverage data to plot.", call. = FALSE)
  }
  
  # Generate the top axis track
  gtrack <- GenomeAxisTrack(fontsize = 20, fontcolor = "black", col = "black")
  # Generate the coverage track
  dtrack <- DataTrack(range = bedgraphDT, genome = "ModCR6", 
                      type = "histogram", name = " ",
                      background.title = "#2C3E50", col.histogram = "grey28",
                      # hex code matches flatly top bar, previously "slategrey"
                      fontsize = 18)
  
  # is positions null? 
  # if yes - plot tracks w/o highlights 
  # if no - plot tracks with highlights
  if (is.null(positions)) {
    coveragePlot <- plotTracks(list(gtrack, dtrack))
    #trackList <- list("gtrack" = gtrack, "dtrack" = dtrack)
    return(coveragePlot)
  } else {
    htrack <- HighlightTrack(trackList = list(dtrack),
                             start = positions,
                             width= widths,
                             inBackground = FALSE,
                             fill = "#FFE3E6b8")
    coveragePlot <- plotTracks(list(gtrack, htrack))
    #trackList <- list("gtrack" = gtrack, "htrack" = htrack)
    return(coveragePlot)
  }
}

#----- Build Common Variants Table (Exact Matching) -----#

# Need to handle condition where there's no data!
FindExactMatches <- function(sampleObjectList) {
  # Check if any samples have 0 variants and remove.  
  variantsAreNull <- vector(mode = "logical", length = length(sampleObjectList))
  for (i in 1:length(sampleObjectList)) {
    variantsAreNull[i]<- is.null(names(sampleObjectList[[i]]@variant_list))
  }
  
  selectedSampleObjectList <- sampleObjectList[!variantsAreNull]
  
  # Stop if all the samples selected have no variants in them.
  if (length(selectedSampleObjectList) == 0) {
    stop("No variants are detected in the selected samples.  Cannot identify common variants.",
         call. = FALSE)
  }
  
  #----- Exact Matching (Position & Variant) -----#
  
  # For all the variants between the selected samples, make a table showing which
  #   samples they are in.
  
  sampleVariantsDFList <- vector(mode = "list", length = length(selectedSampleObjectList))
  sampleVariantsPresenceList <- vector(mode = "list", length = length(selectedSampleObjectList))
  # pull extra info here and stick into currentSampleDF
  for (k in 1:length(selectedSampleObjectList)) {
    # Get current sample information
    currentSample <- selectedSampleObjectList[[k]]
    currentSampleName <- currentSample@sample_name
    
    # Build concatenated data frame from all the variants' variant_df slots
    #   get current Variants, put in list
    currentVariantList <- currentSample@variant_list
    #   pull out each variant's variant_df, put in list
    currentVariantDFList <- map(.x = seq(1:length(currentVariantList)),
                                .f = function(idx) currentVariantList[[idx]]@variant_df)
    #   reduce the list into one big df - (all variants' info for one sample only)
    combinedVariantDFs <- Reduce(f = function(df1, df2) {rbind(x = df1, y = df2)},
                                 x = currentVariantDFList)
    
    #   store the samples variants info in sampleVariantsDFList to be reduced later
    sampleVariantsDFList[[k]] <- combinedVariantDFs
    
    # Build variantPresenceDF that will be used for searching for common variants
    variantPresenceDF <- data.frame(names(currentSample@variant_list), 1)
    colnames(variantPresenceDF) <- c("variants", currentSampleName)
    sampleVariantsPresenceList[[k]] <- variantPresenceDF
  }
  
  # Merge the presence data frames together, convert NA's to 0's
  fullPresenceDF <- Reduce(f = function(df1, df2) {merge(x = df1, y = df2,
                                                         by = "variants", all = TRUE)},
                           x = sampleVariantsPresenceList)
  fullPresenceDF[is.na(fullPresenceDF)] <- 0
  
  # remove extra columns before attempting rowSums
  
  # What are the counts for each variant?
  fullPresenceDFModified <- column_to_rownames(fullPresenceDF, var = "variants")
  variantCounts <- rowSums(fullPresenceDFModified)
  
  # Which samples do you find each variant in?
  variantInSamples <- vector(mode = "character", length = nrow(fullPresenceDFModified))
  for (m in 1:nrow(fullPresenceDFModified)) {
    currentRow <- fullPresenceDFModified[m, ]
    currentVariantName <- rownames(currentRow)
    samplesVecString <- paste(colnames(currentRow)[currentRow != 0], collapse = ", ")
    variantInSamples[m] <- samplesVecString
    names(variantInSamples)[m] <- currentVariantName
  }
  
  # Change things to a DF merge
  variantCountsDF <- data.frame("Number_of_Samples" = variantCounts)
  variantInSamplesDF <- data.frame("Samples" = variantInSamples)
  # Make sure the dimensions of these match
  if (!(identical(dim(variantCountsDF), dim(variantInSamplesDF)))) {
    stop("Something is wrong with the variant/sample counting.",
         call. = FALSE)
  }
  # Otherwise, merge them.
  variantDataDF <- merge(variantCountsDF, variantInSamplesDF, by = "row.names")
  variantDataDF <- variantDataDF %>%
    filter(Number_of_Samples >= 2) %>%
    arrange(desc(Number_of_Samples)) %>%
    column_to_rownames(var = "Row.names") %>%
    rownames_to_column("Variation")
  
  # Put warning within app.R to deal with empty table.
  
  # Now pull extra information before displaying table - from Reducing information 
  #   from sampleVariantsDFList
  variantInformationDF <- Reduce(f = function(df1, df2) {merge(x = df1, y = df2, all = TRUE)},
                                 x = sampleVariantsDFList)
  
  # From variantInformationDF, append information to variantDataDF
  variantInformationFullDF <- merge(x = variantDataDF, variantInformationDF, 
                                    by = "Variation", all = FALSE)
  
  # Modify the table
  variantInformationFinal <- variantInformationFullDF %>%
    select(-c(Variation)) %>%
    select(Variation_Readable, everything()) %>%
    arrange(desc(Number_of_Samples))
  
  # Make more readable column names
  colnames(variantInformationFinal) <- c("Variation", "Number of Samples",
                                         "Samples", "Reference Allele", 
                                         "Reference Codon", "Reference Protein", 
                                         "Mutant Allele", "Mutant Codon", 
                                         "Mutant Protein", "Mutation Type")
  
  return(variantInformationFinal)
}




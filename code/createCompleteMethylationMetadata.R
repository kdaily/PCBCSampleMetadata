library(plyr)
library(dplyr)
library(reshape2)
library(synapseClient)
synapseLogin()

source("./addExtraColumns.R")

PROJECT_ID <- "syn1773109"
CELL_LINE_META_ID  <- 'syn2767694'
SAMPLE_META_ID = 'syn3131513'
ASSAY_META_ID = 'syn3105812'

## Get the cell line metadata
cellLineMetadataTable <- synGet(CELL_LINE_META_ID)
cellLineQuery <- paste("SELECT * FROM", cellLineMetadataTable$properties$id)
cellLineMetadata <- tbl_df(synTableQuery(cellLineQuery, loadResult=TRUE)@values)

# Get sample Metadata
sampleProcessMetadataTable <- synGet(SAMPLE_META_ID)
sampleProcessMetadataQuery <- synTableQuery(sprintf("SELECT * FROM %s", sampleProcessMetadataTable$properties$id))
sampleProcessMetadata  <- tbl_df(sampleProcessMetadataQuery@values)

sampleCols <- c("biologicalSampleName", "C4_Cell_Line_ID", 
                "Diffname_short", "SampleTypeIncomplete",
                "SampleType", "SampleTypeNotes",
                "Replicate", "ReplicateNotes",
                "ThawDate", "ThawDateNotes",
                "PassageAtThaw",  "PassageAtThawNotes",
                "DateDNAHarvested", "DateDNAHarvestedNotes", 
                "PassageAtDNAHarvest", "PassageAtDNAHarvestNotes",
                "DateDNAExtracted",  "DateDNAExtractedNotes",
                "DateDNASubmittedToCore", "Notes")

sampleProcessMetadata <- sampleProcessMetadata[, sampleCols]

sampleProcessMetadata <- filter(sampleProcessMetadata, tolower(DateDNASubmittedToCore) != "n/a")

## Get the assay metadata
assayMetadataTable <- synGet(ASSAY_META_ID)
assayQuery <- paste("SELECT * FROM", assayMetadataTable$properties$id)
assayMeta <- synTableQuery(assayQuery, loadResult=TRUE)
assayMetadata <- tbl_df(synTableQuery(assayQuery, loadResult=TRUE)@values)

tblAll <- left_join(assayMetadata, sampleProcessMetadata,
                    by=c("biologicalSampleName"))

tblAll <- left_join(tblAll, cellLineMetadata, by="C4_Cell_Line_ID")

tblAll <- tblAll %>% 
  extraColumns

tblAll %>% filter(is.na(C4_Cell_Line_ID)) %>% select(UID,biologicalSampleName,public)

# tc <- as.tableColumns(as.data.frame(tblAll))
# 
# schema <- TableSchema(name="Methylation Metadata", parent=PROJECT_ID, columns=tc$tableColumns)
# 
# tbl <- Table(tableSchema=schema,
#              values=as.data.frame(tblAll))
# 
# tbl <- synStore(tbl)

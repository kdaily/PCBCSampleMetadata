library(dplyr)
library(synapseClient)

source("./addExtraColumns.R")

synapseLogin()

PROJECT_ID <- "syn1773109"
CELL_LINE_META_ID <- 'syn2767694'
SAMPLE_META_ID = 'syn3131513'
ASSAY_META_ID = 'syn3105814'

## Get the cell line metadata
cellLineMetadataTable <- synGet(CELL_LINE_META_ID)
cellLineQuery <- paste("SELECT * FROM", cellLineMetadataTable$properties$id)
cellLineMetadata <- synTableQuery(cellLineQuery, loadResult=TRUE)@values

# Get sample Metadata
sampleProcessMetadataTable <- synGet(SAMPLE_META_ID)
sampleProcessMetadataQuery <- synTableQuery(sprintf("SELECT * FROM %s", sampleProcessMetadataTable$properties$id))
sampleProcessMetadata  <- sampleProcessMetadataQuery@values

sampleCols <- c("biologicalSampleName", "C4_Cell_Line_ID", 
                "Diffname_short", "SampleTypeIncomplete",
                "SampleType", "SampleTypeNotes",
                "Replicate", "ReplicateNotes",
                "ThawDate", "ThawDateNotes",
                "PassageAtThaw",  "PassageAtThawNotes",
                "DateRNAHarvested", "RNAHarvestedNotes", "PassageAtHarvest",
                "DatemiRNAextracted",  "DatemiRNAExtractedNotes",
                "DateRNASubmittedtoCore", "DateRNASubmittedToCoreNotes", "Notes")

sampleProcessMetadata <- sampleProcessMetadata[, sampleCols]

## Only keep samples that have a non-NA mRNA submitted to core date
sampleProcessMetadata <- filter(sampleProcessMetadata, tolower(DatemiRNAextracted) != "n/a")

## Get the assay metadata
assayMetadataTable <- synGet(ASSAY_META_ID)
assayQuery <- paste("SELECT * FROM", assayMetadataTable$properties$id)
assayMetadata <- synTableQuery(assayQuery, loadResult=TRUE)@values

tblAll <- left_join(assayMetadata, sampleProcessMetadata,
                    by=c("biologicalSampleName"))

tblAll <- left_join(tblAll, cellLineMetadata, by="C4_Cell_Line_ID")

tblAll <- tblAll %>% 
  extraColumns

tblAll %>% filter(is.na(C4_Cell_Line_ID)) %>% select(UID,biologicalSampleName,public,pass_qc)

# ## Clean up any duplicates because of problems in sample procesing table
# ## basically set all fields other than the sample info
# probCols <- setdiff(colnames(sampleProcessMetadata), 
#                     c("biologicalSampleName", "C4_Cell_Line_ID", "Diffname_short"))
# 
# problemSamples <- tblAll %>% count(UID) %>% filter(n>1)
# 
# fixSamples <- tblAll %>% 
#   filter(UID %in% problemSamples$UID) 
# 
# fixSamples[, probCols] <- NA
# fixSamples <- distinct(fixSamples)
# 
# tblAll <- tblAll %>% 
#   filter(!(UID %in% problemSamples$UID)) %>%
#   rbind(fixSamples)


# tc <- as.tableColumns(as.data.frame(tblAll))

# schema <- TableSchema(name="miRNA-Seq Metadata", parent=PROJECT_ID, columns=tc$tableColumns)
# 
# tbl <- Table(tableSchema=schema,
#              values=as.data.frame(tblAll))
# 
# tbl <- synStore(tbl)

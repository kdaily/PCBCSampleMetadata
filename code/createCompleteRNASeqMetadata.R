library(plyr)
library(dplyr)
library(reshape2)
library(synapseClient)
synapseLogin()

PROJECT_ID <- "syn1773109"
CELL_LINE_META_ID <- 'syn2767694'
SAMPLE_META_ID <- 'syn3131513'
ASSAY_META_ID <- 'syn3104413'

## Get the cell line metadata
cellLineMetadataTable <- synGet(CELL_LINE_META_ID)
cellLineQuery <- paste("SELECT * FROM", cellLineMetadataTable$properties$id)
cellLineMetadata <- tbl_df(synTableQuery(cellLineQuery, loadResult=TRUE)@values)

# Get sample Metadata
sampleProcessMetadataTable <- synGet(SAMPLE_META_ID)
sampleProcessMetadataQuery <- synTableQuery(sprintf("SELECT * FROM %s", sampleProcessMetadataTable$properties$id))
sampleProcessMetadata  <- tbl_df(sampleProcessMetadataQuery@values)

sampleProcessMetadata <- mutate(sampleProcessMetadata,
                                PassageAtThaw=as.numeric(PassageAtThaw),
                                PassageAtHarvest=as.numeric(PassageAtHarvest),
                                PassageAtDNAHarvest=as.numeric(PassageAtDNAHarvest))

## Get the assay metadata
assayMetadataTable <- synGet(ASSAY_META_ID)
assayQuery <- paste("SELECT * FROM", assayMetadataTable$properties$id)
assayMetadata <- tbl_df(synTableQuery(assayQuery, loadResult=TRUE)@values)
assayMetadata <- rename(assayMetadata, ReplicateIncomplete=replicate)

tblAll <- left_join(assayMetadata,
                    sampleProcessMetadata,
                    by=c("C4_Cell_Line_ID", "ReplicateIncomplete", "Diffname_short"))

tblAll <- left_join(tblAll, cellLineMetadata, by="C4_Cell_Line_ID")

tc <- as.tableColumns(as.data.frame(tblAll))

schema <- TableSchema(name="RNA-Seq Metadata", parent=PROJECT_ID, columns=tc$tableColumns)

tbl <- Table(tableSchema=schema,
             values=as.data.frame(tblAll))

# tbl <- synGet('syn3156503')

# act <- Activity(name="Table join",
#                 used=list(sampleProcessMetadataTable, cellLineMetadataTable, assayMetadataTable))

tbl <- synStore(tbl)

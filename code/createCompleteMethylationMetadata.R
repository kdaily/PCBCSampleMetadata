library(dplyr)
library(reshape2)
library(synapseClient)
synapseLogin()

PROJECT_ID <- "syn1773109"
CELL_LINE_META_ID  <-  'syn2767694'
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

## Rename replicate column to match assay data
sampleProcessMetadata <- rename(sampleProcessMetadata, BiologicalReplicate=Replicate)

## Get the assay metadata
assayMetadataTable <- synGet(ASSAY_META_ID)
assayQuery <- paste("SELECT * FROM", assayMetadataTable$properties$id)
assayMetadata <- tbl_df(synTableQuery(assayQuery, loadResult=TRUE)@values)

# In this case, there are two files for each methylation, so this needs to be
# accounted for by removing the columns and uniqifying
assayMetadata <- mutate(assayMetadata, File=NULL, Channel=NULL)
assayMetadata <- distinct(assayMetadata)

tblAll <- left_join(assayMetadata,
                    sampleProcessMetadata,
                    by=c("C4_Cell_Line_ID", "BiologicalReplicate", "Diffname_short"))

tblAll <- left_join(tblAll, cellLineMetadata, by="C4_Cell_Line_ID")

tc <- as.tableColumns(as.data.frame(tblAll))

schema <- TableSchema(name="Methylation Metadata", parent=PROJECT_ID, columns=tc$tableColumns)

tbl <- Table(tableSchema=schema,
             values=as.data.frame(tblAll))

tbl <- synStore(tbl)

tbl <- synStore(tbl,
                used=c(sampleProcessMetadataTable$properties$id, 
                       cellLineMetadataTable$properties$id,
                       assayMetadataTable$properties$id))



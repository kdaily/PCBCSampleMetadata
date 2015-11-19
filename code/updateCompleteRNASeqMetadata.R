library(plyr)
library(dplyr)
library(synapseClient)
synapseLogin()

tblCurrent <- synGet('syn3156503')
res <- synTableQuery(paste("SELECT * FROM", tblCurrent@properties$id))
bak <- res@values

source("./createCompleteRNASeqMetadata.R")

### CAREFUL! THIS WILL DELETE ALL ROWS.
synDeleteRows(res)

tblCurrent <- synGet('syn3156503')

tblNew <- Table(tableSchema=tblCurrent,
                values=as.data.frame(tblAll))

tblNew <- synStore(tblNew)

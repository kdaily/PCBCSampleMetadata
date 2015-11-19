library(plyr)
library(dplyr)
library(synapseClient)
synapseLogin()

tblCurrent <- synGet('syn3219876')
res <- synTableQuery(paste("SELECT * FROM", tblCurrent@properties$id))
bak <- res@values

source("./createCompleteMIRNASeqMetadata.R")

### CAREFUL! THIS WILL DELETE ALL ROWS.
synDeleteRows(res)
tblCurrent <- synGet('syn3219876')

tblNew <- Table(tableSchema=tblCurrent,
                values=as.data.frame(tblAll))

tblNew <- synStore(tblNew)

## Just in case, restore from backup of values
# tblNew <- Table(tableSchema=tblCurrent,
#                 values=bak)

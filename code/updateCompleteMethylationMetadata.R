library(plyr)
library(dplyr)
library(reshape2)
library(synapseClient)
synapseLogin()

tblCurrent <- synGet('syn3156828')
res <- synTableQuery(paste("SELECT * FROM", tblCurrent@properties$id))
bak <- res@values

source("./createCompleteMethylationMetadata.R")

### CAREFUL! THIS WILL DELETE ALL ROWS.
synDeleteRows(res)

tblCurrent <- synGet('syn3156828')

tblNew <- Table(tableSchema=tblCurrent,
                values=as.data.frame(tblAll))

tblNew <- synStore(tblNew)


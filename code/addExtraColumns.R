# Add extra columns to the metadata table

extraColumns <- function(x) {
  x %>% 
    dplyr::mutate(Cell_Type_of_Origin_Level2=ifelse(Cell_Type_of_Origin %in% c("CD34+ cells", "mononuclear"),
                                                    "blood",
                                                    as.character(Cell_Type_of_Origin)),
                  Reprogramming_Vector_Type_Level2=ifelse(Reprogramming_Vector_Type %in% c('lentivirus','retrovirus'),
                                                          "integrating",
                                                          as.character(Reprogramming_Vector_Type)),
                  Reprogramming_Vector_Type_Level2=ifelse(Reprogramming_Vector_Type_Level2 %in% c('plasmid','RNA','Sendai virus'),
                                                          "nonintegrating",
                                                          as.character(Reprogramming_Vector_Type_Level2)))
}

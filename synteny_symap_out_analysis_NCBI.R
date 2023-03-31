library(tidyr)
library(dplyr)
library(stringr)

#synteny_XT_XL_DR_MS_HS_df <- read.csv("synteny_HS_MS_XT_XL_DR_without_prod_ann.csv",header = TRUE)

#synteny_XT_XL_DR_MS_HS_df <- readxl::read_excel("Mod_synteny_XT_XL_DR_MS_HS_gped_data.xlsx")

#names(synteny_XT_XL_DR_MS_HS_ann_df) <- c("Description")

#synteny_XT_XL_DR_MS_HS_ann_df$Description[1]

#synteny_XT_XL_DR_MS_HS_d <- synteny_XT_XL_DR_MS_HS_df %>% separate(Description,into = c("XRow","Block","Block_Score",
 #       "NCBI_DR_Product","NCBI_HS_product","NCBI_MS_product","NCBI_XL_product","NCBI_XT_product","PgeneF","PgFSize",
  #      "HitIdx","Run_Size","NCBI_DR_Chr","NCBI_DR_Start","NCBI_DR_End","NCBI_DR_gene",
   #     "NCBI_HS_Chr","NCBI_HS_Start","NCBI_HS_End","NCBI_HS_gene",
    #    "NCBI_MS_Chr","NCBI_MS_Start","NCBI_MS_End","NCBI_MS_gene",
     #   "NCBI_XL_Chr","NCBI_XL_Start","NCBI_XL_End","NCBI_XL_gene",
      #  "NCBI_XT_Chr","NCBI_XT_Start","NCBI_XT_End","NCBI_XT_gene",
       # "NCBI_DR_ID","NCBI_XL_ID","NCBI_MS_ID","NCBI_HS_ID",
        #"NCBI_XL_All_Anno","NCBI_MS_All_Anno","NCBI_HS.All_Anno","NCBI_DR.All_Anno","NCBI_XT.ID","NCBI_XT.All_Anno"
        #),sep = ",")

#synteny_XT_XL_DR_MS_HS_df$X.Row. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.Row.)

##################### synteny details #########################

############### comparing outputs of old syntney file and new synteny file

synteny_XT_XL_DR_MS_HS_df <- read.csv("Updated_NCBI_Full_diff.csv",header = TRUE)

old_synteny <- read.csv("orient_NCBI_Full.csv",header = TRUE)

synteny_XT_XL_DR_MS_HS_df$X.Hit..Merge. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.Hit..Merge.)
synteny_XT_XL_DR_MS_HS_df$X.Hit..Merge. <- as.numeric(synteny_XT_XL_DR_MS_HS_df$X.Hit..Merge.)

max(synteny_XT_XL_DR_MS_HS_df$X.Hit..Merge.[!is.na(synteny_XT_XL_DR_MS_HS_df$X.Hit..Merge.) & synteny_XT_XL_DR_MS_HS_df$X.Hit..Merge. != "-"])

ab <- synteny_XT_XL_DR_MS_HS_df %>% filter(X.Hit..Merge. == 91)

unique(synteny_XT_XL_DR_MS_HS_df$X.Hit..Merge.)

synteny_XT_XL_DR_MS_HS_df$X.I_XT.Glen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Glen.)
synteny_XT_XL_DR_MS_HS_df$X.I_XT.Glen. <- gsub("-","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Glen.)
synteny_XT_XL_DR_MS_HS_df$X.I_XT.Glen. <- as.numeric(synteny_XT_XL_DR_MS_HS_df$X.I_XT.Glen.)
synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hlen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hlen.)
synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hlen. <- gsub("-","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hlen.)
synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hlen. <- as.numeric(synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hlen.)
min(synteny_XT_XL_DR_MS_HS_df$X.I_XT.Glen.[!is.na(synteny_XT_XL_DR_MS_HS_df$X.I_XT.Glen.)])
min(synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hlen.[!is.na(synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hlen.)])

synteny_XT_XL_DR_MS_HS_df$X.I_HS.Glen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Glen.)
synteny_XT_XL_DR_MS_HS_df$X.I_HS.Hlen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Hlen.)

filtered_data <- synteny_XT_XL_DR_MS_HS_df %>% filter(X.Block. != "\'-\'" &  (X.I_XT.ID.!= "\'-\'" | X.I_HS.ID.!= "\'-\'" | X.I_DR.ID.!= "\'-\'" | X.I_MS.ID.!= "\'-\'" | X.BI_L.ID.!= "\'-\'" | X.BI_S.ID.!= "\'-\'"))

write.csv(filtered_data,file = "filtered_data_synteny_raw.csv",row.names = FALSE,quote = FALSE)

####min of XT hit length is 36 and minimum of XT Gene length is 339

##### after the count analysis, restore the data to pass it into the function.

prepare_synteny_data <- function(){
  
synteny_XT_XL_DR_MS_HS_df <- read.csv("current_NCBI.csv",header = TRUE)
  
synteny_XT_XL_DR_MS_HS_df <- synteny_XT_XL_DR_MS_HS_df %>% filter(X.Block. != "\'-\'" &  (X.I_XT.ID.!= "\'-\'" | X.I_HS.ID.!= "\'-\'" | X.I_DR.ID.!= "\'-\'" | X.I_MS.ID.!= "\'-\'" | X.BI_L.ID.!= "\'-\'" | X.BI_S.ID.!= "\'-\'"))
  
synteny_XT_XL_DR_MS_HS_df$X.Hit.. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.Hit..)
#synteny_XT_XL_DR_MS_HS_df$X.Collinear. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.Collinear.)
#synteny_XT_XL_DR_MS_HS_df$X.Hit..Sim. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.Hit..Sim.)
#synteny_XT_XL_DR_MS_HS_df$X.Hit..Id. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.Hit..Id.)
#synteny_XT_XL_DR_MS_HS_df$X.Hit..Merge. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.Hit..Merge.)
synteny_XT_XL_DR_MS_HS_df$X.Block. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.Block.)
synteny_XT_XL_DR_MS_HS_df$X.Block.Score. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.Block.Score.)

######################################### DR #######################

synteny_XT_XL_DR_MS_HS_df$X.I_DR.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_DR.ID.)
synteny_XT_XL_DR_MS_HS_df$X.I_DR.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_DR.Chr.)
synteny_XT_XL_DR_MS_HS_df$X.I_DR.Gstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_DR.Gstart.)
synteny_XT_XL_DR_MS_HS_df$X.I_DR.Gend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_DR.Gend.)
synteny_XT_XL_DR_MS_HS_df$X.I_DR.Gst. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_DR.Gst.)
synteny_XT_XL_DR_MS_HS_df$X.I_DR.Gene.. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_DR.Gene..)
synteny_XT_XL_DR_MS_HS_df$X.I_DR.Hstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_DR.Hstart.)
synteny_XT_XL_DR_MS_HS_df$X.I_DR.Hend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_DR.Hend.)

############################ HS ############################

synteny_XT_XL_DR_MS_HS_df$X.I_HS.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Chr.)
synteny_XT_XL_DR_MS_HS_df$X.I_HS.Gstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Gstart.)
synteny_XT_XL_DR_MS_HS_df$X.I_HS.Gend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Gend.)
synteny_XT_XL_DR_MS_HS_df$X.I_HS.Gst. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Gst.)
synteny_XT_XL_DR_MS_HS_df$X.I_HS.Gene.. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Gene..)
synteny_XT_XL_DR_MS_HS_df$X.I_HS.Hstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Hstart.)
synteny_XT_XL_DR_MS_HS_df$X.I_HS.Hend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Hend.)
synteny_XT_XL_DR_MS_HS_df$X.I_HS.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.ID.)

############################# XL L ################################

synteny_XT_XL_DR_MS_HS_df$X.BI_L.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.Chr.)
synteny_XT_XL_DR_MS_HS_df$X.BI_L.Gstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.Gstart.)
synteny_XT_XL_DR_MS_HS_df$X.BI_L.Gend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.Gend.)
synteny_XT_XL_DR_MS_HS_df$X.BI_L.Gst. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.Gst.)
synteny_XT_XL_DR_MS_HS_df$X.BI_L.Gene.. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.Gene..)
synteny_XT_XL_DR_MS_HS_df$X.BI_L.Hstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.Hstart.)
synteny_XT_XL_DR_MS_HS_df$X.BI_L.Hend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.Hend.)
synteny_XT_XL_DR_MS_HS_df$X.BI_L.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.ID.)

############################## MS #################################

synteny_XT_XL_DR_MS_HS_df$X.I_MS.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.Chr.)
synteny_XT_XL_DR_MS_HS_df$X.I_MS.Gend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.Gend.)
synteny_XT_XL_DR_MS_HS_df$X.I_MS.Gstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.Gstart.)
synteny_XT_XL_DR_MS_HS_df$X.I_MS.Gst. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.Gst.)
synteny_XT_XL_DR_MS_HS_df$X.I_MS.Gene.. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.Gene..)
synteny_XT_XL_DR_MS_HS_df$X.I_MS.Hstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.Hstart.)
synteny_XT_XL_DR_MS_HS_df$X.I_MS.Hend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.Hend.)
synteny_XT_XL_DR_MS_HS_df$X.I_MS.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.ID.)

############################### XL S ##################################

synteny_XT_XL_DR_MS_HS_df$X.BI_S.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.Chr.)
synteny_XT_XL_DR_MS_HS_df$X.BI_S.Gstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.Gstart.)
synteny_XT_XL_DR_MS_HS_df$X.BI_S.Gend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.Gend.)
synteny_XT_XL_DR_MS_HS_df$X.BI_S.Gst. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.Gst.)
synteny_XT_XL_DR_MS_HS_df$X.BI_S.Gene.. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.Gene..)
synteny_XT_XL_DR_MS_HS_df$X.BI_S.Hstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.Hstart.)
synteny_XT_XL_DR_MS_HS_df$X.BI_S.Hend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.Hend.)
synteny_XT_XL_DR_MS_HS_df$X.BI_S.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.ID.)

############################# XT #############################

synteny_XT_XL_DR_MS_HS_df$X.I_XT.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Chr.)
synteny_XT_XL_DR_MS_HS_df$X.I_XT.Gstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Gstart.)
synteny_XT_XL_DR_MS_HS_df$X.I_XT.Gend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Gend.)
synteny_XT_XL_DR_MS_HS_df$X.I_XT.Gst. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Gst.)
synteny_XT_XL_DR_MS_HS_df$X.I_XT.Gene.. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Gene..)
synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hstart.)
synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hend.)
synteny_XT_XL_DR_MS_HS_df$X.I_XT.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.ID.)

############################### cleaning ID attribute and preparing grouping variable ################################

synteny_XT_XL_DR_MS_HS_df$X.I_XT.ID. <- sapply(synteny_XT_XL_DR_MS_HS_df$X.I_XT.ID., function(x)
  
  paste0(gsub(";(.)*","",str_split(x,pattern = "gene-")[[1]][2:length(str_split(x,pattern = "gene-")[[1]])]),collapse = ";")
  
)

synteny_XT_XL_DR_MS_HS_df_X.BI_S.ID. <- sapply(synteny_XT_XL_DR_MS_HS_df$X.BI_S.ID., function(x)
  
  paste0(gsub(";(.)*","",str_split(x,pattern = "gene-")[[1]][2:length(str_split(x,pattern = "gene-")[[1]])]),collapse = ";")
  
)

synteny_XT_XL_DR_MS_HS_df$X.BI_S.ID. <- synteny_XT_XL_DR_MS_HS_df_X.BI_S.ID.

synteny_XT_XL_DR_MS_HS_df_X.I_MS.ID. <- sapply(synteny_XT_XL_DR_MS_HS_df$X.I_MS.ID., function(x)
  
  paste0(gsub(";(.)*","",str_split(x,pattern = "gene-")[[1]][2:length(str_split(x,pattern = "gene-")[[1]])]),collapse = ";")
  
)

synteny_XT_XL_DR_MS_HS_df$X.I_MS.ID. <- synteny_XT_XL_DR_MS_HS_df_X.I_MS.ID.

synteny_XT_XL_DR_MS_HS_df_X.BI_L.ID. <- sapply(synteny_XT_XL_DR_MS_HS_df$X.BI_L.ID., function(x)
  
  paste0(gsub(";(.)*","",str_split(x,pattern = "gene-")[[1]][2:length(str_split(x,pattern = "gene-")[[1]])]),collapse = ";")
  
)

synteny_XT_XL_DR_MS_HS_df$X.BI_L.ID. <- synteny_XT_XL_DR_MS_HS_df_X.BI_L.ID.

synteny_XT_XL_DR_MS_HS_df_X.I_HS.ID. <- sapply(synteny_XT_XL_DR_MS_HS_df$X.I_HS.ID., function(x)
  
  paste0(gsub(";(.)*","",str_split(x,pattern = "gene-")[[1]][2:length(str_split(x,pattern = "gene-")[[1]])]),collapse = ";")
  
)

synteny_XT_XL_DR_MS_HS_df$X.I_HS.ID. <- synteny_XT_XL_DR_MS_HS_df_X.I_HS.ID.

synteny_XT_XL_DR_MS_HS_df_X.I_DR.ID. <- sapply(synteny_XT_XL_DR_MS_HS_df$X.I_DR.ID., function(x)
  
  paste0(gsub(";(.)*","",str_split(x,pattern = "gene-")[[1]][2:length(str_split(x,pattern = "gene-")[[1]])]),collapse = ";")
  
)

synteny_XT_XL_DR_MS_HS_df$X.I_DR.ID. <- synteny_XT_XL_DR_MS_HS_df_X.I_DR.ID.

#synteny_XT_XL_DR_MS_HS_df$X.I_XT.ID. <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.ID.)
#synteny_XT_XL_DR_MS_HS_df$X.BI_S.ID. <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.ID.)
#synteny_XT_XL_DR_MS_HS_df$X.I_MS.ID. <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.ID.)
#synteny_XT_XL_DR_MS_HS_df$X.BI_L.ID. <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.ID.)
#synteny_XT_XL_DR_MS_HS_df$X.I_HS.ID. <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.ID.)
#synteny_XT_XL_DR_MS_HS_df$X.I_DR.ID. <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$X.I_DR.ID.)

#synteny_XT_XL_DR_MS_HS_df$X.I_XT.ID. <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.ID.)
#synteny_XT_XL_DR_MS_HS_df$X.BI_S.ID. <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.ID.)
#synteny_XT_XL_DR_MS_HS_df$X.I_MS.ID. <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.ID.)
#synteny_XT_XL_DR_MS_HS_df$X.BI_L.ID. <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.ID.)
#synteny_XT_XL_DR_MS_HS_df$X.I_HS.ID. <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.ID.)
#synteny_XT_XL_DR_MS_HS_df$X.I_DR.ID. <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$X.I_DR.ID.)

synteny_XT_XL_DR_MS_HS_df <- synteny_XT_XL_DR_MS_HS_df %>% separate_rows(X.I_XT.ID.,sep = ";") %>% 
                 separate_rows(X.BI_S.ID.,sep = ";") %>% 
                     separate_rows(X.I_MS.ID.,sep = ";") %>% 
                       separate_rows(X.BI_L.ID.,sep = ";") %>% 
                         separate_rows(X.I_HS.ID.,sep = ";") %>% 
                          separate_rows(X.I_DR.ID.,sep = ";") 

#synteny_XT_XL_DR_MS_HS_df$ref_XT <- gsub("\\.(.)*","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.ID.)
synteny_XT_XL_DR_MS_HS_df$ref_XL_L <- gsub("\\.[L|S]","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.ID.)
synteny_XT_XL_DR_MS_HS_df$ref_XL_S <- gsub("\\.[L|S]","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.ID.)
#synteny_XT_XL_DR_MS_HS_df$ref_MS <- gsub("\\.(.)*","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.ID.)
#synteny_XT_XL_DR_MS_HS_df$ref_HS <- gsub("\\.(.)*","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.ID.)
#synteny_XT_XL_DR_MS_HS_df$ref_DR <- gsub("\\.(.)*","",synteny_XT_XL_DR_MS_HS_df$X.I_DR.ID.)

synteny_XT_XL_DR_MS_HS_df$ref_XL_L <- str_to_lower(synteny_XT_XL_DR_MS_HS_df$ref_XL_L)
synteny_XT_XL_DR_MS_HS_df$ref_XL_S <- str_to_lower(synteny_XT_XL_DR_MS_HS_df$ref_XL_S)
synteny_XT_XL_DR_MS_HS_df$ref_MS <- str_to_lower(synteny_XT_XL_DR_MS_HS_df$X.I_MS.ID.)
synteny_XT_XL_DR_MS_HS_df$ref_HS <- str_to_lower(synteny_XT_XL_DR_MS_HS_df$X.I_HS.ID.)
synteny_XT_XL_DR_MS_HS_df$ref_DR <- str_to_lower(synteny_XT_XL_DR_MS_HS_df$X.I_DR.ID.)
synteny_XT_XL_DR_MS_HS_df$ref_XT <- str_to_lower(synteny_XT_XL_DR_MS_HS_df$X.I_XT.ID.)

       ################## grouping them #########################

synteny_XT_XL_DR_MS_HS_gped_d <- synteny_XT_XL_DR_MS_HS_df %>% 
select(Synteny_Block_Location = X.Block.,
NCBI_DR_Chr = X.I_DR.Chr.,NCBI_DR_Gstart = X.I_DR.Gstart.,NCBI_DR_Gend = X.I_DR.Gend.,NCBI_DR_ID = X.I_DR.ID.,NCBI_DR_Str = X.I_DR.Gst.,DR_ref = ref_DR, 
NCBI_HS_Chr = X.I_HS.Chr.,NCBI_HS_Gstart = X.I_HS.Gstart.,NCBI_HS_Gend = X.I_HS.Gend.,NCBI_HS_ID = X.I_HS.ID.,NCBI_HS_Str = X.I_HS.Gst.,HS_ref = ref_HS, 
NCBI_MS_Chr = X.I_MS.Chr.,NCBI_MS_Gstart = X.I_MS.Gstart.,NCBI_MS_Gend = X.I_MS.Gend.,NCBI_MS_ID = X.I_MS.ID.,NCBI_MS_Str = X.I_MS.Gst.,MS_ref = ref_MS, 
NCBI_XL_L_Chr = X.BI_L.Chr.,NCBI_XL_L_Gstart = X.BI_L.Gstart.,NCBI_XL_L_Gend = X.BI_L.Gend.,NCBI_XL_L_ID = X.BI_L.ID.,NCBI_XL_L_Str = X.BI_L.Gst.,XL_L_ref = ref_XL_L, 
NCBI_XL_S_Chr = X.BI_S.Chr.,NCBI_XL_S_Gstart = X.BI_S.Gstart.,NCBI_XL_S_Gend = X.BI_S.Gend.,NCBI_XL_S_ID = X.BI_S.ID.,NCBI_XL_S_Str = X.BI_S.Gst.,XL_S_ref = ref_XL_S, 
NCBI_XT_Chr = X.I_XT.Chr.,NCBI_XT_Gstart = X.I_XT.Gstart.,NCBI_XT_Gend = X.I_XT.Gend.,NCBI_XT_ID = X.I_XT.ID.,NCBI_XT_Str = X.I_XT.Gst.,XT_ref = ref_XT) %>% distinct()

write.csv(synteny_XT_XL_DR_MS_HS_gped_d,file = "Mod_synteny_XT_XL_DR_MS_HS_gped_data.csv",row.names = F,quote = F)

return(synteny_XT_XL_DR_MS_HS_gped_d)

}

synteny_XT_XL_DR_MS_HS_df <- prepare_synteny_data()

##write.csv(ab,file = "pax8.csv",row.names = FALSE,quote = FALSE)

############### further parsing of grouping variable ############

#synteny_XT_XL_DR_MS_HS_df <- read.csv("Mod_synteny_XT_XL_DR_MS_HS_gped_data.csv",header = TRUE)

intermediate_grouping_data <- function(synteny_XT_XL_DR_MS_HS_df){

#pax6_data <- synteny_XT_XL_DR_MS_HS_df %>% filter(DR_ref == "pax6" | HS_ref == "pax6" | MS_ref == "pax6" |XT_ref == "pax6" | XL_L_ref == "pax6" | XL_S_ref == "pax6")

#write.csv(pax6_data,file = "pax6_data.csv",row.names = FALSE,quote = FALSE)

#synteny_XT_XL_DR_MS_HS_df$DR_ref <- gsub("\\.(.)*","",synteny_XT_XL_DR_MS_HS_df$DR_ref)
#synteny_XT_XL_DR_MS_HS_df$HS_ref <- gsub("\\.(.)*","",synteny_XT_XL_DR_MS_HS_df$HS_ref)
#synteny_XT_XL_DR_MS_HS_df$MS_ref <- gsub("\\.(.)*","",synteny_XT_XL_DR_MS_HS_df$MS_ref)
#synteny_XT_XL_DR_MS_HS_df$XL_L_ref <- gsub("\\.(.)*","",synteny_XT_XL_DR_MS_HS_df$XL_L_ref)
#synteny_XT_XL_DR_MS_HS_df$XL_S_ref <- gsub("\\.(.)*","",synteny_XT_XL_DR_MS_HS_df$XL_S_ref)

#synteny_XT_XL_DR_MS_HS_df$XT_reference <- synteny_XT_XL_DR_MS_HS_df$XT_ref

######### finding out each data file for data ##################
#### filtering out the genes of interest from each ortholog species and further
#### merging them into a data frame

######## removing NAs 
  
synteny_XT_XL_DR_MS_HS_df <- synteny_XT_XL_DR_MS_HS_df %>% filter(!NCBI_DR_ID == "NA" & NCBI_HS_ID != "NA" & NCBI_XL_L_ID != "NA" & NCBI_XL_S_ID != "NA" & NCBI_MS_ID != "NA" & NCBI_XT_ID != "NA")
  
Filtered_Gene_ls <- list()

unique_xb_gene_ids <- unique(synteny_XT_XL_DR_MS_HS_df$XT_ref[nchar(synteny_XT_XL_DR_MS_HS_df$XT_ref) > 0 & !is.na(synteny_XT_XL_DR_MS_HS_df$XT_ref) & !synteny_XT_XL_DR_MS_HS_df$XT_ref == '-' & !synteny_XT_XL_DR_MS_HS_df$XT_ref == "na"])  
  
synteny_XT_XL_DR_MS_HS_filtered_synteny <- synteny_XT_XL_DR_MS_HS_df %>% filter(Synteny_Block_Location != "-")

for(i in 1:length(unique_xb_gene_ids)){
  
  Filtered_Gene_ls[[i]] <- synteny_XT_XL_DR_MS_HS_filtered_synteny %>% filter(

    XT_ref == unique_xb_gene_ids[i] | XL_L_ref == unique_xb_gene_ids[i] | XL_S_ref == unique_xb_gene_ids[i] | DR_ref == unique_xb_gene_ids[i] | HS_ref == unique_xb_gene_ids[i] | MS_ref == unique_xb_gene_ids[i] 
  
)  %>% select(
  NCBI_DR_ID,
  NCBI_DR_Str,
  NCBI_MS_ID,
  NCBI_MS_Str,
  NCBI_HS_ID,
  NCBI_HS_Str,
  NCBI_XL_L_ID,
  NCBI_XL_L_Str,
  NCBI_XL_S_ID,
  NCBI_XL_S_Str,
  NCBI_XT_ID,
  NCBI_XT_Str) %>% distinct()
  
}

Mod_Parsed_Synteny_data <- do.call("rbind.data.frame",Filtered_Gene_ls)

#write.csv(Mod_Parsed_Synteny_data,file = "Intermediate_Mod_Parsed_Synteny_data.csv",row.names = FALSE,quote = FALSE)

Filtered_Gene_Mod_ls <- lapply(Filtered_Gene_ls, function(x)
  
  x %>% mutate(
    
    DR_Gene_Combined = paste0(NCBI_DR_ID,collapse = ";"),
    MS_Gene_Combined = paste0(NCBI_MS_ID,collapse = ";"),
    HS_Gene_Combined = paste0(NCBI_HS_ID,collapse = ";"),
    XL_L_Gene_Combined = paste0(NCBI_XL_L_ID,collapse = ";"),
    XL_S_Gene_Combined = paste0(NCBI_XL_S_ID,collapse = ";"),
    XT_Gene_Combined = paste0(NCBI_XT_ID,collapse = ";")
    
  ) %>% select(DR_Gene_Combined,DR_Str = NCBI_DR_Str,MS_Gene_Combined,MS_Str = NCBI_MS_Str,HS_Gene_Combined,HS_Str = NCBI_HS_Str,XL_L_Gene_Combined,XL_L_Str = NCBI_XL_L_Str,XL_S_Gene_Combined,XL_S_Str = NCBI_XL_S_Str,XT_Gene_Combined,XT_Str = NCBI_XT_Str) %>% distinct() 
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data <- do.call("rbind.data.frame",Filtered_Gene_Mod_ls)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$DR_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$DR_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][str_split(x,pattern = ";")[[1]] != "-"]),collapse = ";")
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$DR_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$DR_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][!is.na(str_split(x,pattern = ";")[[1]])]),collapse = ";")
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$HS_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$HS_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][str_split(x,pattern = ";")[[1]] != "-"]),collapse = ";")
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$HS_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$HS_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][!is.na(str_split(x,pattern = ";")[[1]])]),collapse = ";")
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$MS_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$MS_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][str_split(x,pattern = ";")[[1]] != "-"]),collapse = ";")
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$MS_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$MS_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][!is.na(str_split(x,pattern = ";")[[1]])]),collapse = ";")
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XL_L_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XL_L_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][str_split(x,pattern = ";")[[1]] != "-"]),collapse = ";")
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XL_L_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XL_L_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][!is.na(str_split(x,pattern = ";")[[1]])]),collapse = ";")
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XL_S_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XL_S_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][str_split(x,pattern = ";")[[1]] != "-"]),collapse = ";")
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XL_S_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XL_S_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][!is.na(str_split(x,pattern = ";")[[1]])]),collapse = ";")
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XT_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XT_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][str_split(x,pattern = ";")[[1]] != "-"]),collapse = ";")
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XT_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XT_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][!is.na(str_split(x,pattern = ";")[[1]])]),collapse = ";")
  
)

#Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d <- Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data %>%
 # separate_rows(XT_Gene_Combined, sep = ";") %>% distinct()

#Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$DR_Count <- sapply(Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$DR_Gene_Combined, function(x)
  
  
 # length(unique(str_split(x,pattern =  ";")[[1]]))
  
#) 

#Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$MS_Count <- sapply(Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$MS_Gene_Combined, function(x)
  
  
 # length(unique(str_split(x,pattern =  ";")[[1]]))
  
#)

#Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$HS_Count <- sapply(Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$HS_Gene_Combined, function(x)
  
  
 # length(unique(str_split(x,pattern =  ";")[[1]]))
  
#)

#Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$XL_L_Count <- sapply(Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$XL_L_Gene_Combined, function(x)
  
  
 # length(unique(str_split(x,pattern =  ";")[[1]]))
  
#)

#Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$XL_S_Count <- sapply(Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$XL_S_Gene_Combined, function(x)
  
  
 # length(unique(str_split(x,pattern =  ";")[[1]]))
  
#)

#Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$XT_Count <- sapply(Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$XT_Gene_Combined, function(x)
  
  
 # length(unique(str_split(x,pattern =  ";")[[1]])))

Filtered_Synteny_Genes <- Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data %>%
  select(XT = XT_Gene_Combined,HS = HS_Gene_Combined,XL_L = XL_L_Gene_Combined,XL_S = XL_S_Gene_Combined,MS = MS_Gene_Combined,DR = DR_Gene_Combined) %>% distinct()

Filtered_Synteny_Genes_D <- Filtered_Synteny_Genes %>%
  separate_rows(XT,sep = ";") %>%
  separate_rows(XL_L,sep = ";")  %>% 
  separate_rows(HS,sep = ";") %>% 
  separate_rows(MS,sep = ";")%>% 
  separate_rows(DR,sep = ";") %>% 
  separate_rows(XL_S,sep = ";") %>% distinct()


write.csv(Filtered_Synteny_Genes_D,file = "Filtered_Synteny_Genes_Data.csv",row.names = FALSE,quote = FALSE)


    return(Filtered_Synteny_Genes_D)


}

Filtered_Synteny_Genes_Data <- intermediate_grouping_data(synteny_XT_XL_DR_MS_HS_df)

#data_length <- nrow(Filtered_Synteny_Genes_Data)/2

#Filtered_Synteny_Genes_partitioned_D <- Filtered_Synteny_Genes_Data[1 : data_length,]
#Filtered_Synteny_Genes_partitioned_other_half <- Filtered_Synteny_Genes_Data[((data_length)+1) : nrow(Filtered_Synteny_Genes_Data),]


#Filtered_Synteny_Genes_partitioned_D <- Filtered_Synteny_Genes_partitioned_D %>%
 #separate_rows(XT,sep = ";") %>%
 #distinct()

#Filtered_Synteny_Genes_partitioned_D <- Filtered_Synteny_Genes_partitioned_D %>% 
 # separate_rows(HS,sep = ";") 

#Filtered_Synteny_Genes_partitioned_D <- Filtered_Synteny_Genes_partitioned_D %>% 
 # separate_rows(MS,sep = ";")

#Filtered_Synteny_Genes_partitioned_D <- Filtered_Synteny_Genes_partitioned_D %>% 
 # separate_rows(DR,sep = ";") %>% distinct()

#Filtered_Synteny_Genes_partitioned_D <- Filtered_Synteny_Genes_partitioned_D %>% 
 # separate_rows(DR,sep = ";") %>% distinct()

#Filtered_Synteny_Genes_partitioned_D <- Filtered_Synteny_Genes_partitioned_D %>% 
 # separate_rows(XL_L,sep = ";") %>% distinct()

#Filtered_Synteny_Genes_partitioned_D <- Filtered_Synteny_Genes_partitioned_D %>% 
 # separate_rows(XL_S,sep = ";") %>% distinct()


#ann_length <- round(nrow(Filtered_Synteny_Genes_Data) / 6,digits = 0)
#ann_length <- 6 ###choosing length for partitioning data
#ann_data <- 3033

#for (i in 2 : (ann_length - 1)){
  
 # filtered_ls[[i]] <- Filtered_Synteny_Genes_Data[(ann_data*i)+1 : ann_length*(i+1),]
#  sep_ls[[i]] <- filtered_ls[[i]] %>% 
 #   separate_rows(XT,sep = ";") %>%
  #  separate_rows(HS,sep = ";") %>% 
   # separate_rows(MS,sep = ";") %>% 
    #separate_rows(DR,sep = ";") %>%
    #separate_rows(XL_L,sep = ";") %>% 
    #separate_rows(XL_S,sep = ";") %>% distinct()

#}

#filtered_ls[[6]] <- Filtered_Synteny_Genes_Data[1 : 3033,] %>% 
 # separate_rows(XT,sep = ";") %>%
#  separate_rows(HS,sep = ";") %>% 
#  separate_rows(MS,sep = ";") %>% 
#  separate_rows(DR,sep = ";") %>%
 # separate_rows(XL_L,sep = ";") %>% 
#  separate_rows(XL_S,sep = ";") %>% distinct()

#final_parsed_d <- do.call("rbind.data.frame",filtered_ls)

#parsed_d <- rbind.data.frame(final_parsed_d,Filtered_Synteny_Genes_Data[18199:18202,])


#Filtered_Synteny_Genes_Data <- Filtered_Synteny_Genes_Data %>% 
 # separate_rows(XT,sep = ";") %>%
#  separate_rows(HS,sep = ";") %>% 
 # separate_rows(MS,sep = ";") %>% 
#  separate_rows(DR,sep = ";") %>%
#  separate_rows(XL_L,sep = ";") %>% 
#  separate_rows(XL_S,sep = ";") %>% distinct()

Filtered_Synteny_Genes_Data %>% filter(str_detect(XT,pattern = "LOC")) %>% select(XT) %>% distinct() %>% count()

#########4844 LOCs are present in total in XT 

Filtered_Synteny_Genes_Data %>% filter(str_detect(XL_L,pattern = "LOC")) %>% select(XL_L) %>% distinct() %>% count()

############### 4184 LOCs in Laevis L

Filtered_Synteny_Genes_Data %>% filter(str_detect(XL_S,pattern = "LOC")) %>% select(XL_S) %>% distinct() %>% count()

##############  2929 LOCs in Laevis S  

#abc <- ab %>% filter(XT == "pax8")

#abc1 <- ab %>% filter(str_detect(XT,pattern = "LOC"))

#ag <- abc1[1:"13000",]
  
############## adding Uniprot Accession, Entrez Gene IDs and XB IDs #######

parse_Gene_ID <- function(V10_Data,name){
  
  gene_data <- V10_Data %>% filter(type == "gene")
  
  gene_data$gene_id <- sapply(gene_data$attr,
                              
     function(x)  gsub("GeneID:","",str_extract(x,pattern = "GeneID:[0-9]+"))                             
                              
  )
  
  gene_data$gene_name <- sapply(gene_data$attr,
                                
    function(x)  gsub("gene=","",str_split(str_extract(x,pattern = "gene=.+;"),pattern = ";")[[1]][1])                             
                                
  )
  
  if(name == "Xenbase"){
    
    gene_data$XB_Gene_ID <- sapply(gene_data$attr,
                                   
      function(x)  gsub("Xenbase:","",str_extract(x,pattern = "Xenbase:XB-GENE-[0-9]+"))                             
                                   
    )
    
    gene_id_gene_name_mapping <- gene_data[,c("gene_id","gene_name","XB_Gene_ID")] %>% distinct()
    
    return(gene_id_gene_name_mapping)
    
  }
  
  else if(name == "DR" | name == "HS" | name == "MS" | name == "XL"){
    
    gene_id_gene_name_mapping <- gene_data[,c("gene_id","gene_name")] %>% distinct()
    
    return(gene_id_gene_name_mapping)
    
  }
}

#Filtered_Synteny_Genes_Data <- read.csv("Filtered_Synteny_Genes_Data.csv",header = TRUE) 

construct <- function(){

XT_Data <- read.delim("GCF_000004195.4_UCB_Xtro_10.0_genomic.gff",comment.char = "#",header = FALSE,sep = "\t") ###use the file created for Uniprot for trop
names(XT_Data) <- c("chrloc","source","type","start","end","score","strand","phase","attr")
parsed_XT_mapping <- parse_Gene_ID(XT_Data,"Xenbase")

Filtered_Syn_df <- Filtered_Synteny_Genes_Data %>% left_join(parsed_XT_mapping,by = c("XT" = "gene_name")) %>%
  select(DR,HS,MS,XL_L,XL_S,XT,XT_Gene_ID = gene_id,XB_ID = XB_Gene_ID)

#Filtered_Syn_df <- Filtered_d %>% left_join(parsed_XT_mapping,by = c("NCBI_XT_ID" = "gene_name")) %>%
# select(Synteny_Block_Location,NCBI_XT_Chr,NCBI_XT_ID,XT_Gene_ID = gene_id,XT_ID = XB_Gene_ID,NCBI_XL_Chr,NCBI_XL_ID,NCBI_HS_Chr,NCBI_HS_ID,NCBI_MS_Chr,NCBI_MS_ID,NCBI_DR_Chr,NCBI_DR_ID,XT_reference,DR_ref,HS_ref,MS_ref,XL_ref)

DR_Data <- read.delim("GCF_000002035.6_GRCz11_genomic.gff",comment.char = "#",header = FALSE,sep = "\t")
names(DR_Data) <- c("chrloc","source","type","start","end","score","strand","phase","attr")
parsed_DR_mapping <- parse_Gene_ID(DR_Data,"DR")

Filtered_Syn_d <- Filtered_Syn_df %>% left_join(parsed_DR_mapping,by = c("DR" = "gene_name")) %>%
  select(DR,DR_Gene_ID = gene_id,HS,MS,XL_L,XL_S,XT,XT_Gene_ID,XB_ID)
#select(Synteny_Block_Location,NCBI_XT_Chr,NCBI_XT_ID,XT_Gene_ID,XT_ID,NCBI_XL_Chr,NCBI_XL_ID,NCBI_HS_Chr,NCBI_HS_ID,NCBI_MS_Chr,NCBI_MS_ID,NCBI_DR_Chr,NCBI_DR_ID,DR_Gene_ID = gene_id,XT_reference,DR_ref,HS_ref,MS_ref,XL_ref)

HS_Data <- read.delim("GCF_000001405.40_GRCh38.p14_genomic.gff",comment.char = "#",header = FALSE,sep = "\t")
names(HS_Data) <- c("chrloc","source","type","start","end","score","strand","phase","attr")
parsed_HS_mapping <- parse_Gene_ID(HS_Data,"HS")

Filtered_Syn <- Filtered_Syn_d %>% left_join(parsed_HS_mapping,by = c("HS" = "gene_name")) %>%
  #select(Synteny_Block_Location,NCBI_XT_Chr,NCBI_XT_ID,XT_Gene_ID,XT_ID,NCBI_XL_Chr,NCBI_XL_ID,NCBI_HS_Chr,NCBI_HS_ID,HS_Gene_ID = gene_id,NCBI_MS_Chr,NCBI_MS_ID,NCBI_DR_Chr,NCBI_DR_ID,DR_Gene_ID,XT_reference,DR_ref,HS_ref,MS_ref,XL_ref)
  select(DR,DR_Gene_ID,HS,HS_Gene_ID = gene_id,MS,XL_L,XL_S,XT,XT_Gene_ID,XB_ID)

MS_Data <- read.delim("GCF_000001635.27_GRCm39_genomic.gff",comment.char = "#",header = FALSE,sep = "\t")
names(MS_Data) <- c("chrloc","source","type","start","end","score","strand","phase","attr")
parsed_MS_mapping <- parse_Gene_ID(MS_Data,"MS")

Filtered_Synteny <- Filtered_Syn %>% left_join(parsed_MS_mapping,by = c("MS" = "gene_name")) %>%
  #select(Synteny_Block_Location,NCBI_XT_Chr,NCBI_XT_ID,XT_Gene_ID,XT_ID,NCBI_XL_Chr,NCBI_XL_ID,NCBI_HS_Chr,NCBI_HS_ID,HS_Gene_ID,NCBI_MS_Chr,NCBI_MS_ID,MS_Gene_ID = gene_id,NCBI_DR_Chr,NCBI_DR_ID,DR_Gene_ID,XT_reference,DR_ref,HS_ref,MS_ref,XL_ref)
  select(DR,DR_Gene_ID,HS,HS_Gene_ID,MS,MS_Gene_ID = gene_id,XL_L,XL_S,XT,XT_Gene_ID,XB_ID)

XL_Data <- read.delim("GCF_017654675.1_Xenopus_laevis_v10.1_genomic.gff",comment.char = "#",header = FALSE,sep = "\t") ###use the file created for Uniprot for Laevis
names(XL_Data) <- c("chrloc","source","type","start","end","score","strand","phase","attr")
parsed_XL_mapping <- parse_Gene_ID(XL_Data,"Xenbase")

Filtered_Synteny_Final <- Filtered_Synteny %>% left_join(parsed_XL_mapping,by = c("XL_L" = "gene_name")) %>%
  #select(Synteny_Block_Location,NCBI_XT_Chr,NCBI_XT_ID,XT_Gene_ID,XT_ID,NCBI_XL_Chr,NCBI_XL_ID,XL_Gene_ID = gene_id,XL_ID = XB_Gene_ID,NCBI_HS_Chr,NCBI_HS_ID,HS_Gene_ID,NCBI_MS_Chr,NCBI_MS_ID,MS_Gene_ID,NCBI_DR_Chr,NCBI_DR_ID,DR_Gene_ID,XT_reference,DR_ref,HS_ref,MS_ref,XL_ref)
  select(DR,DR_Gene_ID,HS,HS_Gene_ID,MS,MS_Gene_ID,XL_L,XL_L_Gene_ID = gene_id,XL_L_ID = XB_Gene_ID,XL_S,XT,XT_Gene_ID,XB_ID)

Filtered_Synteny_Final <- Filtered_Synteny_Final %>% left_join(parsed_XL_mapping,by = c("XL_S" = "gene_name")) %>% 
  #select(Synteny_Block_Location,NCBI_XT_Chr,NCBI_XT_ID,XT_Gene_ID,XT_ID,NCBI_XL_Chr,NCBI_XL_ID,XL_Gene_ID = gene_id,XL_ID = XB_Gene_ID,NCBI_HS_Chr,NCBI_HS_ID,HS_Gene_ID,NCBI_MS_Chr,NCBI_MS_ID,MS_Gene_ID,NCBI_DR_Chr,NCBI_DR_ID,DR_Gene_ID,XT_reference,DR_ref,HS_ref,MS_ref,XL_ref)
  select(DR,DR_Gene_ID,HS,HS_Gene_ID,MS,MS_Gene_ID,XL_L,XL_L_Gene_ID,XL_L_ID,XL_S,XL_S_Gene_ID = gene_id,XL_S_ID = XB_Gene_ID,XT,XT_Gene_ID,XB_ID) %>% distinct()

write.csv(Filtered_Synteny_Final,file = "Filtered_Synteny_Finalized.csv",row.names = FALSE,quote = FALSE)

 return(Filtered_Synteny_Final)

}

Filtered_Synteny_Final <- construct()

Filtered_Synteny_Final$HS <- paste0("\'",Filtered_Synteny_Final$HS,"\'")

write.csv(Filtered_Synteny_Final,file = "Filtered_Synteny_Finalized.csv",row.names = FALSE)

################ adding Uniprot Data ##########################

#Uniprot_V10_XT_df <- read.delim("XT_v10_Proteome_Data.txt",sep = "\t",header = TRUE)

#Uniprot_V10_XT_df <- readxl::read_excel("XT_v10_Proteome_Data.xlsx")

#Uniprot_V10_XL_df <- read.delim("XL_v10_Proteome_Data.txt",sep = "\t",header = TRUE)

#names(Uniprot_V10_XT_df) <- c("Entry","Reviewed","Entry_Name","Protein_names","Gene_Names","Organism","Length","Ensembl","RefSeq","Xenbase","Xb_Gene_ID"
#)

#Uniprot_V10_XT <- Uniprot_V10_XT_df %>% select(Entry,Xb_Gene_ID)

#Uniprot_V10_XT <- Uniprot_V10_XT %>% separate_rows(Xb_Gene_ID, sep = ";") %>% filter(nchar(Xb_Gene_ID) > 0)

#Uniprot_V10_XL <- Uniprot_V10_XL %>% separate_rows(Xb_Gene_ID, sep = ";") %>% filter(nchar(Xb_Gene_ID) > 0)

#Uniprot_V10_XT$Xb_Gene_ID <- as.numeric(Uniprot_V10_XT$Xb_Gene_ID)

#Filtered_Synteny_Finalized$XT_Gene_ID <- as.numeric(Filtered_Synteny_Finalized$XT_Gene_ID)

#Uniprot_V10_XL$Xb_Gene_ID <- as.numeric(Uniprot_V10_XL$Xb_Gene_ID)

#parsed_XT_Filtered_Synteny_Data <- Filtered_Synteny_Finalized %>% left_join(Uniprot_V10_XT,by = c("XT_Gene_ID" = "Xb_Gene_ID")) %>%
  #select(Synteny_Block_Location,NCBI_XT_Chr,NCBI_XT_ID,XT_Gene_ID,XT_ID,NCBI_XL_Chr,NCBI_XL_ID,XL_Gene_ID,XL_ID,NCBI_HS_Chr,NCBI_HS_ID,HS_Gene_ID,NCBI_MS_Chr,NCBI_MS_ID,MS_Gene_ID,NCBI_DR_Chr,NCBI_DR_ID,DR_Gene_ID,Entry)
 # select(DR,DR_Gene_ID,HS,HS_Gene_ID,MS,MS_Gene_ID,XL_L,XL_L_Gene_ID,XL_L_ID,XL_S,XL_S_Gene_ID,XL_S_ID,XT,XT_Gene_ID,XB_ID,Entry) %>% distinct()

#parsed_XT_XL_Filtered_Synteny_Data <- Filtered_Synteny_Final %>% left_join(Uniprot_V10_XL,by = c("XT_Gene_ID" = "Xb_Gene_ID")) %>%
# select(DR,DR_Gene_ID,HS,HS_Gene_ID,MS,MS_Gene_ID,XL,XL_Gene_ID,XL_ID,XT,XT_Gene_ID,XB_ID,Entry,Protein_names)  

#write.csv(parsed_XT_Filtered_Synteny_Data,file = "Parsed_Filtered_Synteny_Genes_Data.csv",row.names = FALSE,quote = FALSE)

################################################################

#####creating a full additional file with synteny block details and pairwise alignment details of files

###############################################################


Synteny_Genes_Data <- read.csv("Mod_synteny_XT_XL_DR_MS_HS_gped_data.csv")

parse_full_chr_str_id_data <- function(Synteny_Genes_Data){

  Synteny_Genes_Data <- Synteny_Genes_Data %>% filter(Synteny_Block_Location != "-")

  ################## grouping them #########################
  
unique_xb_gene_ids <- unique(Synteny_Genes_Data$XT_ref[nchar(Synteny_Genes_Data$XT_ref) > 0 & !is.na(Synteny_Genes_Data$XT_ref) & !Synteny_Genes_Data$XT_ref == '-' & !Synteny_Genes_Data$XT_ref == "na"])  

Filtered_Gene_ls <- list()

for(i in 1:length(unique_xb_gene_ids)){
  
  Filtered_Gene_ls[[i]] <- Synteny_Genes_Data %>% filter( XT_ref == unique_xb_gene_ids[i] | XL_L_ref == unique_xb_gene_ids[i] | XL_S_ref == unique_xb_gene_ids[i] | DR_ref == unique_xb_gene_ids[i] | HS_ref == unique_xb_gene_ids[i] | MS_ref == unique_xb_gene_ids[i] 
) %>% distinct()

}

Filtered_d <- do.call("rbind.data.frame",Filtered_Gene_ls) %>% distinct()


#Filtered_d$NCBI_XT_ID <- sapply(Filtered_d$NCBI_XT_ID, function(x)
  
 # paste0(gsub(";(.)*","",str_split(x,pattern = "gene-")[[1]][2:length(str_split(x,pattern = "gene-")[[1]])]),collapse = ";")
  
#)

#Filtered_d$NCBI_XL_S_ID <- sapply(Filtered_d$NCBI_XL_S_ID, function(x)
  
 # paste0(gsub(";(.)*","",str_split(x,pattern = "gene-")[[1]][2:length(str_split(x,pattern = "gene-")[[1]])]),collapse = ";")
  
#)

#Filtered_d$NCBI_XL_L_ID <- sapply(Filtered_d$NCBI_XL_L_ID, function(x)
  
 # paste0(gsub(";(.)*","",str_split(x,pattern = "gene-")[[1]][2:length(str_split(x,pattern = "gene-")[[1]])]),collapse = ";")
  
#)

#Filtered_d$NCBI_MS_ID <- sapply(Filtered_d$NCBI_MS_ID, function(x)
  
 # paste0(gsub(";(.)*","",str_split(x,pattern = "gene-")[[1]][2:length(str_split(x,pattern = "gene-")[[1]])]),collapse = ";")
  
#)

#Filtered_d$NCBI_HS_ID <- sapply(Filtered_d$NCBI_HS_ID, function(x)
  
 # paste0(gsub(";(.)*","",str_split(x,pattern = "gene-")[[1]][2:length(str_split(x,pattern = "gene-")[[1]])]),collapse = ";")
  
#)

#Filtered_d$NCBI_DR_ID <- sapply(Filtered_d$NCBI_DR_ID, function(x)
  
 # paste0(gsub(";(.)*","",str_split(x,pattern = "gene-")[[1]][2:length(str_split(x,pattern = "gene-")[[1]])]),collapse = ";")
  
#)

write.csv(Filtered_d,file = "Parsed_Additional_Filtered_Synteny_Genes_Data.csv",row.names = FALSE,quote = FALSE)


#return(Filtered_d)

}

parse_full_chr_str_id_data(Synteny_Genes_Data)

#Filtered_d$NCBI_DR_ID <- gsub(";(.)*","",Filtered_d$NCBI_DR_ID)
#Filtered_d$NCBI_HS_ID <- gsub(";(.)*","",Filtered_d$NCBI_HS_ID)
#Filtered_d$NCBI_MS_ID <- gsub(";(.)*","",Filtered_d$NCBI_MS_ID)
#Filtered_d$NCBI_XL_L_ID <- gsub(";(.)*","",Filtered_d$NCBI_XL_L_ID)
#Filtered_d$NCBI_XL_S_ID <- gsub(";(.)*","",Filtered_d$NCBI_XL_S_ID)
#Filtered_d$NCBI_XT_ID <- gsub(";(.)*","",Filtered_d$NCBI_XT_ID)

#Filtered_d$NCBI_DR_ID <- gsub("gene-","",Filtered_d$NCBI_DR_ID)
#Filtered_d$NCBI_HS_ID <- gsub("gene-","",Filtered_d$NCBI_HS_ID)
#Filtered_d$NCBI_MS_ID <- gsub("gene-","",Filtered_d$NCBI_MS_ID)
#Filtered_d$NCBI_XL_L_ID <- gsub("gene-","",Filtered_d$NCBI_XL_L_ID)
#Filtered_d$NCBI_XL_S_ID <- gsub("gene-","",Filtered_d$NCBI_XL_S_ID)
#Filtered_d$NCBI_XT_ID <- gsub("gene-","",Filtered_d$NCBI_XT_ID)


####################################################################

######################### older version of symap synteny data ############

synteny_XT_XL_DR_MS_HS_df$X.Block. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.Block.)
synteny_XT_XL_DR_MS_HS_df$X.Block.Score. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.Block.Score.)
synteny_XT_XL_DR_MS_HS_df$X.PgeneF. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.PgeneF.)
synteny_XT_XL_DR_MS_HS_df$X.PgFSize. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.PgFSize.)
synteny_XT_XL_DR_MS_HS_df$X.HitIdx. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.HitIdx.)
synteny_XT_XL_DR_MS_HS_df$X.Run.Size. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.Run.Size.)

synteny_XT_XL_DR_MS_HS_df$X.NCBI_DR.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_DR.Chr.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_DR.Start. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_DR.Chr.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_DR.End. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_DR.End.)
#synteny_XT_XL_DR_MS_HS_df$X.NCBI_DR..Gene. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_DR..Gene.)

#synteny_XT_XL_DR_MS_HS_df$X.NCBI_HS..Gene. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_HS..Gene.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_HS.End. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_HS.End.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_HS.Start. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_HS.Start.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_HS.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_HS.Chr.)

synteny_XT_XL_DR_MS_HS_df$X.NCBI_MS.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_MS.Chr.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_MS.Start. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_MS.Start.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_MS.End. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_MS.End.)
#synteny_XT_XL_DR_MS_HS_df$X.NCBI_MS..Gene. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_MS..Gene.)

#synteny_XT_XL_DR_MS_HS_df$X.NCBI_XL..Gene. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_XL..Gene.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_XL.End. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_XL.End.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_XL.Start. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_XL.Start.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_XL.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_XL.Chr.)

#synteny_XT_XL_DR_MS_HS_df$X.NCBI_XT..Gene. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_XT..Gene.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_XT.End. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_XT.End.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_XT.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_XT.Chr.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_XT.Start. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_XT.Start.)

synteny_XT_XL_DR_MS_HS_df$X.NCBI_DR.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_DR.ID.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_HS.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_HS.ID.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_MS.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_MS.ID.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_XL.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_XL.ID.)
synteny_XT_XL_DR_MS_HS_df$X.NCBI_XT.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_XT.ID.)

synteny_XT_XL_DR_MS_HS_df$X.Run.Size. <- as.numeric(synteny_XT_XL_DR_MS_HS_df$X.Run.Size.)

############## ###################### ################# preparing grouping variable

synteny_XT_XL_DR_MS_HS_df$ref_DR <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_DR.ID.)
synteny_XT_XL_DR_MS_HS_df$ref_HS <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_HS.ID.)
synteny_XT_XL_DR_MS_HS_df$ref_MS <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_MS.ID.)
synteny_XT_XL_DR_MS_HS_df$ref_XL <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_XL.ID.)
synteny_XT_XL_DR_MS_HS_df$ref_XT <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$X.NCBI_XT.ID.)

synteny_XT_XL_DR_MS_HS_df$ref_XT <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$ref_XT)
synteny_XT_XL_DR_MS_HS_df$ref_XL <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$ref_XL)
synteny_XT_XL_DR_MS_HS_df$ref_MS <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$ref_MS)
synteny_XT_XL_DR_MS_HS_df$ref_HS <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$ref_HS)
synteny_XT_XL_DR_MS_HS_df$ref_DR <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$ref_DR)

synteny_XT_XL_DR_MS_HS_df$ref_XL <- gsub("\\.(.)*","",synteny_XT_XL_DR_MS_HS_df$ref_XL)

synteny_XT_XL_DR_MS_HS_df$ref_XL <- str_to_lower(synteny_XT_XL_DR_MS_HS_df$ref_XL)
synteny_XT_XL_DR_MS_HS_df$ref_MS <- str_to_lower(synteny_XT_XL_DR_MS_HS_df$ref_MS)
synteny_XT_XL_DR_MS_HS_df$ref_HS <- str_to_lower(synteny_XT_XL_DR_MS_HS_df$ref_HS)
synteny_XT_XL_DR_MS_HS_df$ref_DR <- str_to_lower(synteny_XT_XL_DR_MS_HS_df$ref_DR)

################## grouping genes ######################

synteny_XT_XL_DR_MS_HS_gped_d <- synteny_XT_XL_DR_MS_HS_df %>% group_by(ref_XT,ref_XL,ref_MS,ref_HS,ref_DR) %>% select(Synteny_Block_Location = X.Block.,
No_of_linkages_In_synteny_block = X.Block.Score.,
Run_Size_Length_of_Syn_Genes = X.Run.Size.,
NCBI_DR_Chr = X.NCBI_DR.Chr.,NCBI_DR_Start = X.NCBI_DR.Start.,NCBI_DR_End = X.NCBI_DR.End.,NCBI_DR_ID = X.NCBI_DR.ID.,DR_ref = ref_DR, 
NCBI_HS_Chr = X.NCBI_HS.Chr.,NCBI_HS_Start = X.NCBI_HS.Start.,NCBI_HS_End = X.NCBI_HS.End.,NCBI_HS_ID = X.NCBI_HS.ID.,HS_ref = ref_HS,
NCBI_MS_Chr = X.NCBI_MS.Chr.,NCBI_MS_Start = X.NCBI_MS.Start.,NCBI_MS_End = X.NCBI_MS.End.,NCBI_MS_ID = X.NCBI_MS.ID.,MS_ref = ref_MS,
NCBI_XL_Chr = X.NCBI_XL.Chr.,NCBI_XL_Start = X.NCBI_XL.Start.,NCBI_XL_End = X.NCBI_XL.End.,NCBI_XL_ID = X.NCBI_XL.ID.,XL_ref = ref_XL,
NCBI_XT_Chr = X.NCBI_XT.Chr.,NCBI_XT_Start = X.NCBI_XT.Start.,NCBI_XT_End = X.NCBI_XT.End.,NCBI_XT_ID = X.NCBI_XT.ID.,XT_reference = ref_XT)

write.csv(synteny_XT_XL_DR_MS_HS_gped_d,file = "Mod_synteny_XT_XL_DR_MS_HS_gped_data.csv",row.names = FALSE,quote = FALSE)

#############################################################

#### preparing the data requested 

#synteny_XT_XL_DR_MS_HS_df <- readxl::read_excel("Mod_synteny_XT_XL_DR_MS_HS_gped_data.xlsx")

synteny_XT_XL_DR_MS_HS_df <- read.csv("Mod_synteny_XT_XL_DR_MS_HS_gped_data.csv",header = TRUE)

pax6_data <- synteny_XT_XL_DR_MS_HS_df %>% filter(DR_ref == "pax6" | HS_ref == "pax6" | MS_ref == "pax6" |XT_reference == "pax6" | XL_ref == "pax6")

write.csv(pax6_data,file = "pax6_data.csv",row.names = FALSE,quote = FALSE)

synteny_XT_XL_DR_MS_HS_df$NCBI_DR_ID <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$NCBI_DR_ID)
synteny_XT_XL_DR_MS_HS_df$NCBI_HS_ID <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$NCBI_HS_ID)
synteny_XT_XL_DR_MS_HS_df$NCBI_MS_ID <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$NCBI_MS_ID)
synteny_XT_XL_DR_MS_HS_df$NCBI_XL_ID <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$NCBI_XL_ID)
synteny_XT_XL_DR_MS_HS_df$NCBI_XT_ID <- gsub(";(.)*","",synteny_XT_XL_DR_MS_HS_df$NCBI_XT_ID)

synteny_XT_XL_DR_MS_HS_df$NCBI_DR_ID <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$NCBI_DR_ID)
synteny_XT_XL_DR_MS_HS_df$NCBI_HS_ID <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$NCBI_HS_ID)
synteny_XT_XL_DR_MS_HS_df$NCBI_MS_ID <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$NCBI_MS_ID)
synteny_XT_XL_DR_MS_HS_df$NCBI_XL_ID <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$NCBI_XL_ID)
synteny_XT_XL_DR_MS_HS_df$NCBI_XT_ID <- gsub("gene-","",synteny_XT_XL_DR_MS_HS_df$NCBI_XT_ID)

synteny_XT_XL_DR_MS_HS_df$DR_ref <- gsub("\\.(.)*","",synteny_XT_XL_DR_MS_HS_df$DR_ref)
synteny_XT_XL_DR_MS_HS_df$HS_ref <- gsub("\\.(.)*","",synteny_XT_XL_DR_MS_HS_df$HS_ref)
synteny_XT_XL_DR_MS_HS_df$MS_ref <- gsub("\\.(.)*","",synteny_XT_XL_DR_MS_HS_df$MS_ref)
synteny_XT_XL_DR_MS_HS_df$XL_ref <- gsub("\\.(.)*","",synteny_XT_XL_DR_MS_HS_df$XL_ref)
synteny_XT_XL_DR_MS_HS_df$XT_reference <- gsub("\\.(.)*","",synteny_XT_XL_DR_MS_HS_df$XT_reference)

#synteny_XT_XL_DR_MS_HS_gped <- synteny_XT_XL_DR_MS_HS_df %>% filter(Synteny_Block_Location != "-") %>% group_by(across(c(DR_ref,HS_ref,MS_ref,XL_ref,XT_reference))) %>% summarise(
  
 # DR_Gene_Combined = paste0(unique(NCBI_DR_ID),collapse = ";"),
#  MS_Gene_Combined = paste0(unique(NCBI_MS_ID),collapse = ";"),
#  HS_Gene_Combined = paste0(unique(NCBI_HS_ID),collapse = ";"),
# XL_Gene_Combined = paste0(unique(NCBI_XL_ID),collapse = ";"),
#  XT_Gene_Combined = paste0(unique(NCBI_XT_ID),collapse = ";"),
#  .groups = "keep") 

#%>% select(XT_Gene_Combined,HS_Gene_Combined,XL_Gene_Combined,DR_Gene_Combined,MS_Gene_Combined) %>% distinct()

#synteny_data <- synteny_XT_XL_DR_MS_HS_gped %>% separate_rows(DR_Gene_Combined,sep = ";") %>% separate_rows(MS_Gene_Combined,sep = ";") %>% separate_rows(HS_Gene_Combined,sep = ";") %>% separate_rows(XT_Gene_Combined,sep = ";") %>% separate_rows(XL_Gene_Combined,sep = ";") %>% select(
  
#  XT = XT_Gene_Combined,
#  XL = XL_Gene_Combined,
#  HS = HS_Gene_Combined,
#  DR = DR_Gene_Combined,
#  MS = MS_Gene_Combined

#) %>% distinct()


####################################################################

######### finding out each data file for data ##################
#### filtering out the genes of interest from each ortholog species and further
#### merging them into a data frame

Filtered_Gene_ls <- list()

unique_xb_gene_ids <- unique(synteny_XT_XL_DR_MS_HS_df$XT_reference)

synteny_XT_XL_DR_MS_HS_filtered_synteny <- synteny_XT_XL_DR_MS_HS_df %>% filter(Synteny_Block_Location != "-")

for(i in 1:length(unique_xb_gene_ids)){
  
Filtered_Gene_ls[[i]] <- synteny_XT_XL_DR_MS_HS_filtered_synteny %>% filter(DR_ref == unique_xb_gene_ids[i] | 
                           HS_ref == unique_xb_gene_ids[i] | 
                           MS_ref == unique_xb_gene_ids[i] |
                           XT_reference == unique_xb_gene_ids[i] | 
                           XL_ref == unique_xb_gene_ids[i])  %>% select(DR_ref,NCBI_DR_Chr,NCBI_DR_ID,MS_ref,NCBI_MS_Chr,NCBI_MS_ID,HS_ref,NCBI_HS_Chr,NCBI_HS_ID,XL_ref,NCBI_XL_Chr,NCBI_XL_ID,XT_reference,NCBI_XT_Chr,NCBI_XT_ID) %>% distinct()
  
}

abc <- do.call("rbind",Filtered_Gene_ls)

Filtered_Gene_Mod_ls <- lapply(Filtered_Gene_ls, function(x)

  x %>% mutate(
  
  DR_Gene_Combined = paste0(NCBI_DR_ID,collapse = ";"),
  MS_Gene_Combined = paste0(NCBI_MS_ID,collapse = ";"),
  HS_Gene_Combined = paste0(NCBI_HS_ID,collapse = ";"),
  XL_Gene_Combined = paste0(NCBI_XL_ID,collapse = ";"),
  XT_Gene_Combined = paste0(NCBI_XT_ID,collapse = ";")

  ) %>% select(DR_Gene_Combined,DR_Chr = NCBI_DR_Chr,MS_Gene_Combined,MS_Chr = NCBI_MS_Chr,HS_Gene_Combined,HS_Chr = NCBI_HS_Chr,XL_Gene_Combined,XL_Chr = NCBI_XL_Chr,XT_Gene_Combined,XT_Chr = NCBI_XT_Chr) %>% distinct() 

)

#Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data <- rbind.data.frame(unlist(Filtered_Gene_Mod_ls))

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data <- do.call("rbind.data.frame",Filtered_Gene_Mod_ls)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$DR_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$DR_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][str_split(x,pattern = ";")[[1]] != "-"]),collapse = ";")

)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$HS_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$HS_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][str_split(x,pattern = ";")[[1]] != "-"]),collapse = ";")
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$MS_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$MS_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][str_split(x,pattern = ";")[[1]] != "-"]),collapse = ";")
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XL_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XL_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][str_split(x,pattern = ";")[[1]] != "-"]),collapse = ";")
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XT_Gene_Combined <- sapply(
  
  Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data$XT_Gene_Combined, function(x)
    
    paste0(unique(str_split(x,pattern = ";")[[1]][str_split(x,pattern = ";")[[1]] != "-"]),collapse = ";")
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d <- Mod_Parsed_Synteny_DR_HS_MS_XL_XT_data %>%
  separate_rows(XT_Gene_Combined, sep = ";") %>% distinct()

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$DR_Count <- sapply(Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$DR_Gene_Combined, function(x)
  
  
  length(unique(str_split(x,pattern =  ";")[[1]]))
  
) 

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$MS_Count <- sapply(Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$MS_Gene_Combined, function(x)
  
  
  length(unique(str_split(x,pattern =  ";")[[1]]))
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$HS_Count <- sapply(Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$HS_Gene_Combined, function(x)
  
  
  length(unique(str_split(x,pattern =  ";")[[1]]))
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$XL_Count <- sapply(Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$XL_Gene_Combined, function(x)
  
  
  length(unique(str_split(x,pattern =  ";")[[1]]))
  
)

Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$XT_Count <- sapply(Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d$XT_Gene_Combined, function(x)
  
  
  length(unique(str_split(x,pattern =  ";")[[1]])))

Filtered_Synteny_Genes <- Mod_Parsed_Synteny_DR_HS_MS_XL_XT_d %>% filter( !(DR_Count > 50 | MS_Count > 50 | HS_Count > 50 | XL_Count > 50 | XT_Count > 50)) %>%
  select(XT = XT_Gene_Combined,XT_Chr,HS = HS_Gene_Combined,HS_Chr,XL = XL_Gene_Combined,XL_Chr,MS = MS_Gene_Combined,MS_Chr,DR = DR_Gene_Combined,DR_Chr)

Filtered_Synteny_Genes_D <- Filtered_Synteny_Genes %>% 
  separate_rows(XL,sep = ";") %>% distinct()

Filtered_Synteny_Genes_D <- Filtered_Synteny_Genes_D %>% 
  separate_rows(HS,sep = ";") %>% distinct()

Filtered_Synteny_Genes_D <- Filtered_Synteny_Genes_D %>% 
  separate_rows(MS,sep = ";") %>% distinct()

Filtered_Synteny_Genes_Data <- Filtered_Synteny_Genes_D %>% 
  separate_rows(DR,sep = ";") %>% distinct()

#################################################################

##  Adding in Uniprot Accession,Entrez Gene IDs and XB Gene ID 

################################################################

parse_Gene_ID <- function(V10_Data,name){
  
  gene_data <- V10_Data %>% filter(type == "gene")
  
  gene_data$gene_id <- sapply(gene_data$attr,
                                    
   function(x)  gsub("GeneID:","",str_extract(x,pattern = "GeneID:[0-9]+"))                             
                                    
  )
  
  gene_data$gene_name <- sapply(gene_data$attr,
                                   
     function(x)  gsub("gene=","",str_split(str_extract(x,pattern = "gene=.+;"),pattern = ";")[[1]][1])                             
                                   
  )
  
  if(name == "Xenbase"){
    
    gene_data$XB_Gene_ID <- sapply(gene_data$attr,
                                  
      function(x)  gsub("Xenbase:","",str_extract(x,pattern = "Xenbase:XB-GENE-[0-9]+"))                             
                                  
    )
    
    gene_id_gene_name_mapping <- gene_data[,c("gene_id","gene_name","XB_Gene_ID")] %>% distinct()
    
    return(gene_id_gene_name_mapping)
    
  }
  
  else if(name == "DR" | name == "HS" | name == "MS" | name == "XL"){
  
  gene_id_gene_name_mapping <- gene_data[,c("gene_id","gene_name")] %>% distinct()
 
  return(gene_id_gene_name_mapping)
  
  }
}

Filtered_Synteny_Genes_Data <- read.csv("Filtered_Synteny_Genes_Data.csv",header = TRUE)

XT_Data <- read.delim("GCF_000004195.4_UCB_Xtro_10.0_genomic.gff",comment.char = "#",header = FALSE,sep = "\t")
names(XT_Data) <- c("chrloc","source","type","start","end","score","strand","phase","attr")
parsed_XT_mapping <- parse_Gene_ID(XT_Data,"Xenbase")

Filtered_Syn_df <- Filtered_Synteny_Genes_Data %>% left_join(parsed_XT_mapping,by = c("XT" = "gene_name")) %>%
select(DR,HS,MS,XL,XT,XT_Gene_ID = gene_id,XB_ID = XB_Gene_ID)

#Filtered_Syn_df <- Filtered_d %>% left_join(parsed_XT_mapping,by = c("NCBI_XT_ID" = "gene_name")) %>%
 # select(Synteny_Block_Location,NCBI_XT_Chr,NCBI_XT_ID,XT_Gene_ID = gene_id,XT_ID = XB_Gene_ID,NCBI_XL_Chr,NCBI_XL_ID,NCBI_HS_Chr,NCBI_HS_ID,NCBI_MS_Chr,NCBI_MS_ID,NCBI_DR_Chr,NCBI_DR_ID,XT_reference,DR_ref,HS_ref,MS_ref,XL_ref)

DR_Data <- read.delim("GCF_000002035.6_GRCz11_genomic.gff",comment.char = "#",header = FALSE,sep = "\t")
names(DR_Data) <- c("chrloc","source","type","start","end","score","strand","phase","attr")
parsed_DR_mapping <- parse_Gene_ID(DR_Data,"DR")

Filtered_Syn_d <- Filtered_Syn_df %>% left_join(parsed_DR_mapping,by = c("NCBI_DR_ID" = "gene_name")) %>%
  select(DR,DR_Gene_ID = gene_id,HS,MS,XL,XT,XT_Gene_ID,XB_ID)
  #select(Synteny_Block_Location,NCBI_XT_Chr,NCBI_XT_ID,XT_Gene_ID,XT_ID,NCBI_XL_Chr,NCBI_XL_ID,NCBI_HS_Chr,NCBI_HS_ID,NCBI_MS_Chr,NCBI_MS_ID,NCBI_DR_Chr,NCBI_DR_ID,DR_Gene_ID = gene_id,XT_reference,DR_ref,HS_ref,MS_ref,XL_ref)

HS_Data <- read.delim("GCF_000001405.40_GRCh38.p14_genomic.gff",comment.char = "#",header = FALSE,sep = "\t")
names(HS_Data) <- c("chrloc","source","type","start","end","score","strand","phase","attr")
parsed_HS_mapping <- parse_Gene_ID(HS_Data,"HS")

Filtered_Syn <- Filtered_Syn_d %>% left_join(parsed_HS_mapping,by = c("NCBI_HS_ID" = "gene_name")) %>%
  #select(Synteny_Block_Location,NCBI_XT_Chr,NCBI_XT_ID,XT_Gene_ID,XT_ID,NCBI_XL_Chr,NCBI_XL_ID,NCBI_HS_Chr,NCBI_HS_ID,HS_Gene_ID = gene_id,NCBI_MS_Chr,NCBI_MS_ID,NCBI_DR_Chr,NCBI_DR_ID,DR_Gene_ID,XT_reference,DR_ref,HS_ref,MS_ref,XL_ref)
  select(DR,DR_Gene_ID,HS,HS_Gene_ID = gene_id,MS,XL,XT,XT_Gene_ID,XB_ID)

MS_Data <- read.delim("GCF_000001635.27_GRCm39_genomic.gff",comment.char = "#",header = FALSE,sep = "\t")
names(MS_Data) <- c("chrloc","source","type","start","end","score","strand","phase","attr")
parsed_MS_mapping <- parse_Gene_ID(MS_Data,"MS")

Filtered_Synteny <- Filtered_Syn %>% left_join(parsed_MS_mapping,by = c("NCBI_MS_ID" = "gene_name")) %>%
  #select(Synteny_Block_Location,NCBI_XT_Chr,NCBI_XT_ID,XT_Gene_ID,XT_ID,NCBI_XL_Chr,NCBI_XL_ID,NCBI_HS_Chr,NCBI_HS_ID,HS_Gene_ID,NCBI_MS_Chr,NCBI_MS_ID,MS_Gene_ID = gene_id,NCBI_DR_Chr,NCBI_DR_ID,DR_Gene_ID,XT_reference,DR_ref,HS_ref,MS_ref,XL_ref)
  select(DR,DR_Gene_ID,HS,HS_Gene_ID,MS,MS_Gene_ID = gene_id,XL,XT,XT_Gene_ID,XB_ID)

XL_Data <- read.delim("GCF_017654675.1_Xenopus_laevis_v10.1_genomic.gff",comment.char = "#",header = FALSE,sep = "\t")
names(XL_Data) <- c("chrloc","source","type","start","end","score","strand","phase","attr")
parsed_XL_mapping <- parse_Gene_ID(XL_Data,"Xenbase")

Filtered_Synteny_Final <- Filtered_Synteny %>% left_join(parsed_XL_mapping,by = c("NCBI_XL_ID" = "gene_name")) %>%
  #select(Synteny_Block_Location,NCBI_XT_Chr,NCBI_XT_ID,XT_Gene_ID,XT_ID,NCBI_XL_Chr,NCBI_XL_ID,XL_Gene_ID = gene_id,XL_ID = XB_Gene_ID,NCBI_HS_Chr,NCBI_HS_ID,HS_Gene_ID,NCBI_MS_Chr,NCBI_MS_ID,MS_Gene_ID,NCBI_DR_Chr,NCBI_DR_ID,DR_Gene_ID,XT_reference,DR_ref,HS_ref,MS_ref,XL_ref)
  select(DR,DR_Gene_ID,HS,HS_Gene_ID,MS,MS_Gene_ID,XL,XL_Gene_ID = gene_id,XL_ID = XB_Gene_ID,XT,XT_Gene_ID,XB_ID)


##### after parsing gene id data and then adding in protein data #####

Uniprot_V10_XT_df <- read.delim("XT_v10_Proteome_Data.txt",sep = "\t",header = TRUE)

#Uniprot_V10_XL_df <- read.delim("XL_v10_Proteome_Data.txt",sep = "\t",header = TRUE)

names(Uniprot_V10_XT_df) <- c("Entry","Reviewed","Entry_Name","Protein_names","Gene_Names","Organism","Length","Ensembl","RefSeq","Xenbase","Xb_Gene_ID"
)

Uniprot_V10_XT <- Uniprot_V10_XT_df %>% select(Entry,Xb_Gene_ID)

Uniprot_V10_XT <- Uniprot_V10_XT %>% separate_rows(Xb_Gene_ID, sep = ";") %>% filter(nchar(Xb_Gene_ID) > 0)

#Uniprot_V10_XL <- Uniprot_V10_XL %>% separate_rows(Xb_Gene_ID, sep = ";") %>% filter(nchar(Xb_Gene_ID) > 0)

Uniprot_V10_XT$Xb_Gene_ID <- as.numeric(Uniprot_V10_XT$Xb_Gene_ID)

Filtered_Synteny_Final$XT_Gene_ID <- as.numeric(Filtered_Synteny_Final$XT_Gene_ID)

#Uniprot_V10_XL$Xb_Gene_ID <- as.numeric(Uniprot_V10_XL$Xb_Gene_ID)

parsed_XT_Filtered_Synteny_Data <- Filtered_Synteny_Final %>% left_join(Uniprot_V10_XT,by = c("XT_Gene_ID" = "Xb_Gene_ID")) %>%
  #select(Synteny_Block_Location,NCBI_XT_Chr,NCBI_XT_ID,XT_Gene_ID,XT_ID,NCBI_XL_Chr,NCBI_XL_ID,XL_Gene_ID,XL_ID,NCBI_HS_Chr,NCBI_HS_ID,HS_Gene_ID,NCBI_MS_Chr,NCBI_MS_ID,MS_Gene_ID,NCBI_DR_Chr,NCBI_DR_ID,DR_Gene_ID,Entry)
   select(DR,DR_Gene_ID,HS,HS_Gene_ID,MS,MS_Gene_ID,XL,XL_Gene_ID,XL_ID,XT,XT_Gene_ID,XB_ID,Entry)  
  
#parsed_XT_XL_Filtered_Synteny_Data <- Filtered_Synteny_Final %>% left_join(Uniprot_V10_XL,by = c("XT_Gene_ID" = "Xb_Gene_ID")) %>%
 # select(DR,DR_Gene_ID,HS,HS_Gene_ID,MS,MS_Gene_ID,XL,XL_Gene_ID,XL_ID,XT,XT_Gene_ID,XB_ID,Entry,Protein_names)  

write.csv(parsed_XT_Filtered_Synteny_Data,file = "Parsed_Filtered_Synteny_Genes_Data.csv",row.names = FALSE,quote = FALSE)


################################################################

#####creating a full additional file with synteny block details and pairwise alignment details of files

###############################################################

Synteny_Genes_Data <- read.csv("Mod_synteny_XT_XL_DR_MS_HS_gped_data.csv",header = TRUE)

Filtered_Synteny_Genes_Data <- Synteny_Genes_Data %>% filter(Synteny_Block_Location != "-")

unique_xb_ids <- unique(Filtered_Synteny_Genes_Data$XT_reference)

Filtered_ls <- list()

for(i in 1: length(unique_xb_ids)){

Filtered_ls[[i]] <- 
Filtered_Synteny_Genes_Data %>% filter(DR_ref == unique_xb_ids[i] | 
HS_ref == unique_xb_ids[i] | 
MS_ref == unique_xb_ids[i] |
XT_reference == unique_xb_ids[i] | 
XL_ref == unique_xb_ids[i]) 

}

Filtered_df <- do.call("rbind.data.frame",Filtered_ls)

Filtered_d <- Filtered_df %>% select(Synteny_Block_Location,NCBI_XT_Chr,NCBI_XT_ID,NCBI_XL_Chr,NCBI_XL_ID,NCBI_HS_Chr,NCBI_HS_ID,NCBI_MS_Chr,NCBI_MS_ID,NCBI_DR_Chr,NCBI_DR_ID,XT_reference,DR_ref,HS_ref,MS_ref,XL_ref)

Filtered_d$NCBI_DR_ID <- gsub(";(.)*","",Filtered_d$NCBI_DR_ID)
Filtered_d$NCBI_HS_ID <- gsub(";(.)*","",Filtered_d$NCBI_HS_ID)
Filtered_d$NCBI_MS_ID <- gsub(";(.)*","",Filtered_d$NCBI_MS_ID)
Filtered_d$NCBI_XL_ID <- gsub(";(.)*","",Filtered_d$NCBI_XL_ID)
Filtered_d$NCBI_XT_ID <- gsub(";(.)*","",Filtered_d$NCBI_XT_ID)

Filtered_d$NCBI_DR_ID <- gsub("gene-","",Filtered_d$NCBI_DR_ID)
Filtered_d$NCBI_HS_ID <- gsub("gene-","",Filtered_d$NCBI_HS_ID)
Filtered_d$NCBI_MS_ID <- gsub("gene-","",Filtered_d$NCBI_MS_ID)
Filtered_d$NCBI_XL_ID <- gsub("gene-","",Filtered_d$NCBI_XL_ID)
Filtered_d$NCBI_XT_ID <- gsub("gene-","",Filtered_d$NCBI_XT_ID)

write.csv(Filtered_d,file = "Parsed_Additional_Filtered_Synteny_Genes_Data.csv",row.names = FALSE,quote = FALSE)





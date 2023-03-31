
parse_full_chr_str_id_data <- function(){
  
  synteny_XT_XL_DR_MS_HS_df <- read.csv("current_NCBI.csv",header = TRUE)
  
  synteny_XT_XL_DR_MS_HS_df <- synteny_XT_XL_DR_MS_HS_df %>% filter(X.I_XT.ID.!= "\'-\'" | X.I_HS.ID.!= "\'-\'" | X.I_DR.ID.!= "\'-\'" | X.I_MS.ID.!= "\'-\'" | X.BI_L.ID.!= "\'-\'" | X.BI_S.ID.!= "\'-\'")
  
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
  synteny_XT_XL_DR_MS_HS_df$X.I_DR.Hstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_DR.Hstart.)
  synteny_XT_XL_DR_MS_HS_df$X.I_DR.Hend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_DR.Hend.)
  synteny_XT_XL_DR_MS_HS_df$X.I_DR.Glen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_DR.Glen.)
  synteny_XT_XL_DR_MS_HS_df$X.I_DR.Hlen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_DR.Hlen.)
  
  
  ############################ HS ############################
  
  synteny_XT_XL_DR_MS_HS_df$X.I_HS.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Chr.)
  synteny_XT_XL_DR_MS_HS_df$X.I_HS.Gstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Gstart.)
  synteny_XT_XL_DR_MS_HS_df$X.I_HS.Gend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Gend.)
  synteny_XT_XL_DR_MS_HS_df$X.I_HS.Gst. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Gst.)
  synteny_XT_XL_DR_MS_HS_df$X.I_HS.Hstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Hstart.)
  synteny_XT_XL_DR_MS_HS_df$X.I_HS.Hend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Hend.)
  synteny_XT_XL_DR_MS_HS_df$X.I_HS.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.ID.)
  synteny_XT_XL_DR_MS_HS_df$X.I_HS.Glen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Glen.)
  synteny_XT_XL_DR_MS_HS_df$X.I_HS.Hlen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_HS.Hlen.)
  
  
  ############################# XL L ################################
  
  synteny_XT_XL_DR_MS_HS_df$X.BI_L.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.Chr.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_L.Gstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.Gstart.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_L.Gend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.Gend.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_L.Gst. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.Gst.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_L.Hstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.Hstart.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_L.Hend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.Hend.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_L.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.ID.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_L.Glen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.Glen.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_L.Hlen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_L.Hlen.)
  
  
  ############################## MS #################################
  
  synteny_XT_XL_DR_MS_HS_df$X.I_MS.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.Chr.)
  synteny_XT_XL_DR_MS_HS_df$X.I_MS.Gend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.Gend.)
  synteny_XT_XL_DR_MS_HS_df$X.I_MS.Gstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.Gstart.)
  synteny_XT_XL_DR_MS_HS_df$X.I_MS.Gst. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.Gst.)
  synteny_XT_XL_DR_MS_HS_df$X.I_MS.Hstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.Hstart.)
  synteny_XT_XL_DR_MS_HS_df$X.I_MS.Hend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.Hend.)
  synteny_XT_XL_DR_MS_HS_df$X.I_MS.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.ID.)
  synteny_XT_XL_DR_MS_HS_df$X.I_MS.Glen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.Glen.)
  synteny_XT_XL_DR_MS_HS_df$X.I_MS.Hlen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_MS.Hlen.)
  
  
  ############################### XL S ##################################
  
  synteny_XT_XL_DR_MS_HS_df$X.BI_S.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.Chr.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_S.Gstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.Gstart.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_S.Gend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.Gend.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_S.Gst. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.Gst.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_S.Hstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.Hstart.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_S.Hend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.Hend.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_S.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.ID.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_S.Glen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.Glen.)
  synteny_XT_XL_DR_MS_HS_df$X.BI_S.Hlen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.BI_S.Hlen.)
  
  
  ############################# XT #############################
  
  synteny_XT_XL_DR_MS_HS_df$X.I_XT.Chr. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Chr.)
  synteny_XT_XL_DR_MS_HS_df$X.I_XT.Gstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Gstart.)
  synteny_XT_XL_DR_MS_HS_df$X.I_XT.Gend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Gend.)
  synteny_XT_XL_DR_MS_HS_df$X.I_XT.Gst. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Gst.)
  synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hstart. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hstart.)
  synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hend. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hend.)
  synteny_XT_XL_DR_MS_HS_df$X.I_XT.ID. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.ID.)
  synteny_XT_XL_DR_MS_HS_df$X.I_XT.Glen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Glen.)
  synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hlen. <- gsub("'","",synteny_XT_XL_DR_MS_HS_df$X.I_XT.Hlen.)
  
  
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
  
  #synteny_XT_XL_DR_MS_HS_df <- synteny_XT_XL_DR_MS_HS_df %>% separate_rows(X.I_XT.ID.,sep = ";") %>% 
  #  separate_rows(X.BI_S.ID.,sep = ";") %>% 
  #  separate_rows(X.I_MS.ID.,sep = ";") %>% 
  #  separate_rows(X.BI_L.ID.,sep = ";") %>% 
  #  separate_rows(X.I_HS.ID.,sep = ";") %>% 
  #  separate_rows(X.I_DR.ID.,sep = ";") 
  
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
  
  
  #Synteny_Genes_Data <- Synteny_Genes_Data %>% filter(Synteny_Block_Location != "-")
  
  ################## grouping them #########################
  
  unique_xb_gene_ids <- unique(synteny_XT_XL_DR_MS_HS_df$ref_XT[nchar(synteny_XT_XL_DR_MS_HS_df$ref_XT) > 0 & !is.na(synteny_XT_XL_DR_MS_HS_df$ref_XT) & !synteny_XT_XL_DR_MS_HS_df$ref_XT == '-' & !synteny_XT_XL_DR_MS_HS_df$ref_XT == "na"])  
  
  Filtered_Gene_ls <- list()
  
  for(i in 1:length(unique_xb_gene_ids)){
    
    Filtered_Gene_ls[[i]] <- synteny_XT_XL_DR_MS_HS_df %>% filter( ref_XT == unique_xb_gene_ids[i] | ref_XL_L == unique_xb_gene_ids[i] | ref_XL_S == unique_xb_gene_ids[i] | ref_DR == unique_xb_gene_ids[i] | ref_HS == unique_xb_gene_ids[i] | ref_MS == unique_xb_gene_ids[i] 
    ) %>% select(Synteny_Block_Location = X.Block.,Hit = X.Hit..,
    NCBI_DR_Chr = X.I_DR.Chr.,NCBI_DR_Gstart = X.I_DR.Gstart.,NCBI_DR_Gend = X.I_DR.Gend.,NCBI_DR_Glen = X.I_DR.Glen.,NCBI_DR_Hstart = X.I_DR.Hstart.,NCBI_DR_Hend = X.I_DR.Hend.,NCBI_DR_Hlen = X.I_DR.Hlen.,NCBI_DR_ID = X.I_DR.ID.,NCBI_DR_Str = X.I_DR.Gst.,ref_DR, 
    NCBI_HS_Chr = X.I_HS.Chr.,NCBI_HS_Gstart = X.I_HS.Gstart.,NCBI_HS_Gend = X.I_HS.Gend.,NCBI_HS_Glen = X.I_HS.Glen.,NCBI_HS_Hstart = X.I_HS.Hstart.,NCBI_HS_Hend = X.I_HS.Hend.,NCBI_HS_Hlen = X.I_HS.Hlen.,NCBI_HS_ID = X.I_HS.ID.,NCBI_HS_Str = X.I_HS.Gst.,ref_HS, 
    NCBI_MS_Chr = X.I_MS.Chr.,NCBI_MS_Gstart = X.I_MS.Gstart.,NCBI_MS_Gend = X.I_MS.Gend.,NCBI_MS_Glen = X.I_MS.Glen.,NCBI_MS_Hstart = X.I_MS.Hstart.,NCBI_MS_Hend = X.I_MS.Hend.,NCBI_MS_Hlen = X.I_MS.Hlen.,NCBI_MS_ID = X.I_MS.ID.,NCBI_MS_Str = X.I_MS.Gst.,ref_MS, 
    NCBI_XL_L_Chr = X.BI_L.Chr.,NCBI_XL_L_Gstart = X.BI_L.Gstart.,NCBI_XL_L_Gend = X.BI_L.Gend.,NCBI_XL_L_Glen = X.BI_L.Glen.,NCBI_XL_L_Hstart = X.BI_L.Hstart.,NCBI_XL_L_Hend = X.BI_L.Hend.,NCBI_XL_L_Hlen = X.BI_L.Hlen.,NCBI_XL_L_ID = X.BI_L.ID.,NCBI_XL_L_Str = X.BI_L.Gst.,ref_XL_L, 
    NCBI_XL_S_Chr = X.BI_S.Chr.,NCBI_XL_S_Gstart = X.BI_S.Gstart.,NCBI_XL_S_Gend = X.BI_S.Gend.,NCBI_XL_S_Glen = X.BI_S.Glen.,NCBI_XL_S_Hstart = X.BI_S.Hstart.,NCBI_XL_S_Hend = X.BI_S.Hend.,NCBI_XL_S_Hlen = X.BI_S.Hlen.,NCBI_XL_S_ID = X.BI_S.ID.,NCBI_XL_S_Str = X.BI_S.Gst.,ref_XL_S, 
    NCBI_XT_Chr = X.I_XT.Chr.,NCBI_XT_Gstart = X.I_XT.Gstart.,NCBI_XT_Gend = X.I_XT.Gend.,NCBI_XT_Glen = X.I_XT.Glen.,NCBI_XT_Hstart = X.I_XT.Hstart.,NCBI_XT_Hend = X.I_XT.Hend.,NCBI_XT_Hlen = X.I_XT.Hlen.,NCBI_XT_ID = X.I_XT.ID.,NCBI_XT_Str = X.I_XT.Gst.,ref_XT) %>% distinct()
    
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
  
  write.csv(Filtered_d,file = "Parsed_Additional_Filtered_Synteny_Genes_Data_Full.csv",row.names = FALSE,quote = FALSE)
  
  
  #return(Filtered_d)
  
}

parse_full_chr_str_id_data()


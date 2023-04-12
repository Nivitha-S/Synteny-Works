library(dplyr)
library(stringr)
library(tidyr)

######################################################################

Diopt_V10_Data <- read.csv("Mod_parsed_XT_HS_All_Diopt_Tool_Data.csv",header = TRUE)
 
Synteny_V10_Data <- read.csv("Filtered_Synteny_Finalized.csv",header = TRUE)

compare_synteny_diopt_data <- function(Synteny_V10_Data,Diopt_V10_Data){
  
parsed_syn_V10_D <- Synteny_V10_Data %>% select(XT,XL_L,XL_S,HS,MS,DR,XT_Gene_ID,XB_ID,XL_L_Gene_ID,XL_L_ID,XL_S_Gene_ID,XL_S_ID,HS_Gene_ID,MS_Gene_ID,DR_Gene_ID) %>% distinct()

#write.csv(parsed_syn_V10_D,file = "Parsed_Synteny_Data.csv",row.names = FALSE,quote = FALSE)

Synteny_V10_Data$HS <- gsub("\'",'',Synteny_V10_Data$HS)

Synteny_V10_Data$ref_synteny <- paste0(Synteny_V10_Data$XT,"_",Synteny_V10_Data$HS)

#parsed_syn_V10_D$HS <- gsub("\'",'',parsed_syn_V10_D$HS)

#parsed_syn_V10_D$ref_synteny <- paste0(parsed_syn_V10_D$XT,"_",parsed_syn_V10_D$HS)

Diopt_V10_Data$HS_Diopt <- gsub("\"",'',Diopt_V10_Data$HS)

Diopt_V10_Data$ref_diopt <- paste0(Diopt_V10_Data$XT,"_",Diopt_V10_Data$HS_Diopt)

#Compare_Data <- Synteny_V10_Data %>% full_join(Diopt_V10_Data,by = c("ref_synteny" = "ref_diopt")) %>% 
 # select(XT_synteny = XT.x,
  #       XT_Gene_ID_Synteny = XT_Gene_ID.x,
   #      XB_ID,
   #      HS_synteny = HS.x,
  #       HS_Gene_ID_Synteny = HS_Gene_ID.x,
  #       XL_L,
  #       XL_L_Gene_ID,
  #       XL_L_ID,
  #       XL_S,
  #       XL_S_Gene_ID,
  #       XL_S_ID,
  #       MS,
  #       MS_Gene_ID,
  #       DR,
  #       DR_Gene_ID,
  #       XT_Diopt = XT.y,
  #       HS_Diopt,
  #       ortho,
  #       swiftortho,
   #      phylome,
   #      proteinortho,
  #       Inparanoid,
  #       Sonic,
  #       FastOrtho,
  #       XT_Gene_ID_Diopt = XT_Gene_ID.y,
   #      HS_Gene_ID_Diopt = HS_Gene_ID.y,ref_synteny) %>% distinct()


Compare_Data <- parsed_syn_V10_D %>% full_join(Diopt_V10_Data,by = c("ref_synteny" = "ref_diopt")) %>% 
  select(XT_synteny = XT.x,
         XT_Gene_ID_Synteny = XT_Gene_ID.x,
         XB_ID,
         HS_synteny = HS.x,
         HS_Gene_ID_Synteny = HS_Gene_ID.x,
         XL_L,
         XL_L_Gene_ID,
         XL_L_ID,
         XL_S,
         XL_S_Gene_ID,
         XL_S_ID,
         MS,
         MS_Gene_ID,
         DR,
         DR_Gene_ID,
         XT_Diopt = XT.y,
         HS_Diopt = HS.y,
         ortho,
         swiftortho,
         phylome,
         proteinortho,
         Inparanoid,
         Sonic,
         FastOrtho,
         XT_Gene_ID_Diopt = XT_Gene_ID.y,
         HS_Gene_ID_Diopt = HS_Gene_ID.y,ref_synteny) %>% distinct()

############# creating a boolean value to equate the mappings and a score #######

#Compare_D <- Compare_Data %>% mutate(
  
 # ref_synteny = paste0(XT_synteny,HS_synteny,collapse = "_"),
  
#  ref_diopt =  paste0(XT_Diopt,HS_Diopt,collapse = "_")
  
#)

#Compare_Data$ref_synteny <- paste0(Compare_Data$XT_synteny,"_",Compare_Data$HS_synteny)

#Compare_Data$HS_Diopt <- gsub('"','',Compare_Data$HS_Diopt)

#Compare_Data$ref_diopt <- paste0(Compare_Data$XT_Diopt,"_",Compare_Data$HS_Diopt)

Compare_Data <- Compare_Data %>% mutate(
  
  XT_HS_Presence = ifelse(is.na(XT_synteny) | is.na(HS_synteny) | is.na(XT_Diopt) | is.na(HS_Diopt) | nchar(XT_synteny) == 0 | nchar(HS_synteny) == 0 | nchar(XT_Diopt) == 0 | nchar(HS_Diopt) == 0
,FALSE, TRUE)
  
) 


#%>% 
 # select(XT_synteny,
  ##      XL,
    #     MS,
     #    DR,
      #   DR_Gene_ID,
      #   XT_Gene_ID_Synteny,
      #   XB_ID,
       #  XL_Gene_ID,
      #   XL_ID,
      #   HS_Gene_ID_Synteny,
      #   MS_Gene_ID,
      #   XT_Diopt,
      #   HS_Diopt,
      #   ortho,
      #   swiftortho,
       #  phylome,
      #   proteinortho,
       #  Inparanoid,
      #   Sonic,
      #   FastOrtho,
       #  XT_Gene_ID_Diopt,
        # HS_Gene_ID_Diopt,
        # XT_HS_Presence,ref_synteny) %>% distinct()

#Final_Compared_Data <- Compare_Data %>% 
 # select(XT_synteny,
  #       HS_synteny,
   #      HS_Gene_ID_Synteny,
    #     XL,
     #    MS,
      #   DR,
       #  DR_Gene_ID,
        # XT_Gene_ID_Synteny,
      #   XB_ID,
      #   XL_Gene_ID,
      #   XL_ID,
      #   MS_Gene_ID,
      #   XT_Diopt,
      #   HS_Diopt,
      #   ortho,
      #   swiftortho,
       #  phylome,
      #   proteinortho,
      #   Inparanoid,
      #   Sonic,
      #   FastOrtho,
      #   XT_Gene_ID_Diopt,
      #   HS_Gene_ID_Diopt,
      #   XT_HS_Presence) %>% distinct()

#####################################################################

weightage_score_diopt <- list()
weightage_score_synteny <- list()
XT_HS_Presence_Synteny <- list()


for(i in 1: nrow(Compare_Data)){
  
     n <-  Compare_Data[i,c("ortho","swiftortho","phylome","proteinortho","Inparanoid","Sonic","FastOrtho")] == 1 & nchar(Compare_Data[i,c("ortho","swiftortho","phylome","proteinortho","Inparanoid","Sonic","FastOrtho")]) > 0 & !is.na(Compare_Data[i,c("ortho","swiftortho","phylome","proteinortho","Inparanoid","Sonic","FastOrtho")])

     n_Num_d <- length(n[n == TRUE])
     
     n_XT_HS <- nchar(Compare_Data[i,c("XT_synteny","HS_synteny")]) > 0 & !is.na(Compare_Data[i,c("XT_synteny","HS_synteny")])
     
     n_Num_XT_HS <- length(n_XT_HS[n_XT_HS == TRUE])
     
     n_syn <- nchar(Compare_Data[i,c("XT_synteny","HS_synteny","XL_L","XL_S","MS","DR")]) > 0 & !is.na(Compare_Data[i,c("XT_synteny","HS_synteny","XL_L","XL_S","MS","DR")])
       
     n_syn_Num <- length(n_syn[n_syn == TRUE])
     
    weightage_score_diopt[[i]] <- ifelse(is.numeric(n_Num_d) & n_Num_d != 0, ((n_Num_d*0.143) - 0.001) ,0) 
    weightage_score_synteny[[i]] <- ifelse(is.numeric(n_syn_Num) & n_syn_Num != 0, (n_syn_Num*0.17)-0.02, 0) 
    XT_HS_Presence_Synteny[[i]] <- ifelse(n_Num_XT_HS == 2,TRUE,FALSE)
    
}


weightage_score_diopt <- unlist(weightage_score_diopt)
weightage_score_syn <- unlist(weightage_score_synteny)
XT_HS_Presence_Synteny <- unlist(XT_HS_Presence_Synteny)

Compare_Data$weightage_score_diopt <- weightage_score_diopt
Compare_Data$weightage_score_syn <- weightage_score_syn
Compare_Data$XT_HS_Presence_Synteny <- XT_HS_Presence_Synteny

Compare_Data <- Compare_Data %>% mutate(
  
  Syn_Value = ifelse(XT_HS_Presence == TRUE,1,0)

)

#Compare_Data$HS_Diopt <- sapply(Compare_Data$HS_Diopt, function(x)
  
 # paste0('"',x,'"',collapse = "")
  
#)

write.csv(Compare_Data,file = "Compare_Data.csv",row.names = FALSE,quote = F)  


#return(Compare_Data)


}

compare_synteny_diopt_data(Synteny_V10_Data,Diopt_V10_Data)

Compare_Data <- read.csv("Compare_Data.csv")

prepare_full_data  <- function(Compare_Data){

Final_Compared_D <- Compare_Data %>% 
  select(XT_synteny,
         HS_synteny,
         XL_L,
         XL_S,
         MS,
         DR,
         DR_Gene_ID,
         XT_Gene_ID_Synteny,
         XB_ID,
         XL_L_Gene_ID,
         XL_L_ID,
         XL_S_Gene_ID,
         XL_S_ID,
         HS_Gene_ID_Synteny,
         MS_Gene_ID,
         XT_Diopt,
         HS_Diopt,
         ortho,
         swiftortho,
         phylome,
         proteinortho,
         Inparanoid,
         Sonic,
         FastOrtho,
         XT_Gene_ID_Diopt,
         HS_Gene_ID_Diopt,
         XT_HS_Presence,weightage_score_diopt,XT_HS_Presence_Synteny,weightage_score_syn) %>% distinct()

#filtered_Both_Compare_Df <- Compare_Data %>% select(ref_synteny,weightage_score_diopt,weightage_score_syn,XT_HS_Presence_Synteny,Syn_Value)

#filtered_Diopt_Compare_Df <- Compare_Data %>% filter(weightage_score_diopt != 0) %>% select(ref_synteny,weightage_score_diopt,ortho,
 #                                                  swiftortho,
  #                                                 phylome,
  #                                                 proteinortho,
  #                                                 Inparanoid,
  #                                                 Sonic,
   #                                                FastOrtho) %>% distinct()

#Score_Data <- Compare_Data %>% select(XT_synteny,
 #                                     HS_synteny,
  #                                    XL_L,
   #                                   XL_S,
    #                                  MS,
     #                                 DR,
      #                                DR_Gene_ID,
       #                               XT_Gene_ID_Synteny,
        #                              XB_ID,
         #                             XL_L_Gene_ID,
          #                            XL_L_ID,
           #                           XL_S_Gene_ID,
            #                          XL_S_ID,                                      
             #                         HS_Gene_ID_Synteny,
              #                        MS_Gene_ID,
               #                       XT_Diopt,
                #                      HS_Diopt,
                 #                     XT_Gene_ID_Diopt,
                  #                    HS_Gene_ID_Diopt,
                   #                   XT_HS_Presence,
                    #                  weightage_score_diopt,
                     #                 XT_HS_Presence_Synteny,weightage_score_syn) %>% distinct()

Final_Compared_D$HS_synteny <-  paste0("\'",Final_Compared_D$HS_synteny,"\'")

Final_Compared_D$HS_Diopt <- paste0("\'",Final_Compared_D$HS_Diopt,"\'")

write.csv(Final_Compared_D,file = "Comparison_Synteny_Diopt_Data_Full_Final.csv",row.names = FALSE)  

Final_Compared_Da <- Compare_Data %>% filter(str_detect(XT_synteny,"LOC") | str_detect(XL_L,"LOC") | str_detect(XL_S,"LOC") | str_detect(XT_Diopt,"LOC") | str_detect(XT_synteny,"provisional") | str_detect(XL_L,"provisional") | str_detect(XL_S,"provisional") | str_detect(XT_Diopt,"provisional") | str_detect(XT_synteny,"XB[0-9]+") | str_detect(XL_L,"XB[0-9]+") | str_detect(XL_S,"XB[0-9]+") | str_detect(XT_Diopt,"XB[0-9]+")) %>% select(-Syn_Value,-ref_synteny)

Final_Compared_Da$HS_synteny <- paste0("\'",Final_Compared_Da$HS_synteny,"\'")

Final_Compared_Da$HS_Diopt <- paste0("\'",Final_Compared_Da$HS_Diopt,"\'")

write.csv(Final_Compared_Da,file = "Comparison_Synteny_Diopt_Data_LOC_Final.csv",row.names = FALSE)  

#Score_Data$HS_synteny <- paste0("\'",Score_Data$HS_synteny,"\'")

#Score_Data$HS_Diopt <- paste0("\'",Score_Data$HS_synteny,"\'")

#Compare_Data$HS_synteny <- paste0("\'",Compare_Data$HS_synteny,"\'")

#Compare_Data$HS_Diopt <- paste0("\'",Compare_Data$HS_synteny,"\'")

#write.csv(Score_Data,file = "Score_Data.csv",row.names = FALSE,quote = FALSE)  

#write.csv(filtered_Diopt_Compare_Df,file = "filtered_Diopt_Data.csv",row.names = FALSE)  

#write.csv(filtered_Both_Compare_Df,file = "filtered_both_Value_Synteny_Diopt_Data.csv",row.names = FALSE,quote = FALSE)  

#write.csv(Final_Compared_Data,file = "Comparison_Synteny_Diopt_Data_Full_For_LOC.csv",row.names = FALSE)  

#return(list(filtered_Both_Compare_Df,Compare_Data))


}

prepare_full_data(Compare_Data)

#compare_data <- Full_compare_data[[2]]

#filtered_Both_Compare_Df <- Full_compare_data[[1]]

#compare_d$HS_synteny <- gsub("\\'","",compare_d$HS_synteny)
#compare_d$HS_Diopt <- gsub("\\'","",compare_d$HS_Diopt)

#count_file <- nrow(compare_df)/10

#compare_df <- compare_d %>% mutate(rank = dense_rank(XT_synteny))

#unique_XT_syn <- compare_d %>% filter(str_detect(XT_synteny,"LOC") | str_detect(XT_synteny,"provisional") | str_detect(XT_synteny,"XB[0-9]+") & XT_HS_Presence == FALSE ) %>% select(XT_synteny) %>% distinct() %>% count()
#unique_XT_diopt <- compare_d %>% filter(str_detect(XT_Diopt,"LOC") | str_detect(XT_Diopt,"provisional") | str_detect(XT_Diopt,"XB[0-9]+") & XT_HS_Presence == FALSE ) %>% select(XT_Diopt) %>% distinct() %>% count()
#overlapping_XT_syn_count <- compare_d %>% filter(str_detect(XT_synteny,"LOC") | str_detect(XT_synteny,"provisional") | str_detect(XT_synteny,"XB[0-9]+") & XT_HS_Presence == TRUE) %>% select(XT_synteny) %>% distinct() %>% count()

data_optimization <- function(data){
  
  library(dplyr)
  library(stringr)
  library(writexl)
  
compare_d <- read.csv("Comparison_Synteny_Diopt_Data_Full_Final.csv",header = TRUE)
  
if(data == "XT"){  
  
  #unique_XT_syn_d <- compare_d %>% filter(str_detect(XT_synteny,"LOC") | str_detect(XT_synteny,"provisional") | str_detect(XT_synteny,"XB[0-9]+") ) %>% select(XT_synteny) %>% distinct()
  #unique_XT_diopt_d <- compare_d %>% filter(str_detect(XT_Diopt,"LOC") | str_detect(XT_Diopt,"provisional") | str_detect(XT_Diopt,"XB[0-9]+") ) %>% select(XT_Diopt) %>% distinct()

  total_data <- data.frame(data = c(compare_d$XT_synteny,compare_d$XT_Diopt)) %>% distinct()
  
  unique_XT <- unique(total_data$data)
  
  xb_grouped_list <- list()
  
  for (i in 1 : length(unique_XT)) {
    
    a <- compare_d %>% filter(XT_synteny == unique_XT[i] | XT_Diopt == unique_XT[i]) 
    
    #%>% arrange(desc(weightage_score_synteny)) 
    
    a$xl_l_ref_syn <- str_to_lower(a$XL_L)

    a$xl_s_ref_syn <- str_to_lower(a$XL_S)
    
    a$hs_ref_syn <- str_to_lower(a$HS_synteny)
    
    a$dr_ref <- str_to_lower(a$DR)
    
    a$ms_ref <- str_to_lower(a$MS)
    
    #a$xt_ref_diopt <- str_to_lower(a$XT_Diopt)
    
    a$hs_ref_diopt <- str_to_lower(a$HS_Diopt)
    
    #unique_DR_id <- unique(a$DR)
    
    unique_HS_syn <- unique(a$hs_ref_syn)[!is.na(unique(a$hs_ref_syn)) | nchar(unique(a$hs_ref_syn)) != 0]
    
    unique_syn_XL <- c(unique(a$xl_l_ref_syn)[!is.na(unique(a$xl_l_ref_syn)) | nchar(unique(a$xl_l_ref_syn)) != 0],unique(a$xl_s_ref_syn)[!is.na(unique(a$xl_s_ref_syn)) & nchar(unique(a$xl_s_ref_syn)) != 0])
                       
    unique_syn_DR <- unique(a$dr_ref)[!is.na(unique(a$dr_ref)) | nchar(unique(a$dr_ref)) != 0]
    
    unique_syn_MS <- unique(a$ms_ref)[!is.na(unique(a$ms_ref)) | nchar(unique(a$ms_ref)) != 0]
    
    #unique_MS_id <- unique(a$MS)
    
    #unique_XT_syn <- unique(a$XT_synteny)
    
    #unique_XT_diopt <- unique(a$XT_Diopt)
    
    #unique_HS_diopt <- unique(a$hs_ref_diopt)
    
  if(nrow(a) > 3){
    
    filtered_data <- list()
    
    filtered_data_other <- list()
    
    if(!is.na(unique_HS_syn) & nchar(unique_HS_syn) > 0){
    
    for(j in 1: length(unique_HS_syn)){
      
      #inter_data <- unique_HS_syn[j] %in% a$hs_ref_syn | a$dr_ref == unique_HS_syn[j] | a$ms_ref == unique_HS_syn[j] | a$xt_ref_syn == unique_HS_syn[j] | a$xt_ref_diopt %in% unique_HS_diopt | a$hs_ref_diopt == unique_HS_diopt
      
      inter_data <- a %>% filter(hs_ref_syn == unique_HS_syn[j]) %>% mutate(
        
        hs_presence = ifelse(unique_HS_syn[j] %in% hs_ref_syn,1,0),
        dr_presence = ifelse(unique_HS_syn[j] %in% dr_ref,1,0),
        ms_presence = ifelse(unique_HS_syn[j] %in% ms_ref,1,0),
        xl_l_presence = ifelse(unique_HS_syn[j] == xl_l_ref_syn,1,0),
        xl_s_presence = ifelse(unique_HS_syn[j] == xl_l_ref_syn,1,0),
        total_present_count = hs_presence + dr_presence + ms_presence + xl_l_presence + xl_s_presence
        
      )
      
      inter_data <- inter_data %>% arrange(desc(total_present_count))
      
      filtered_data[[length(filtered_data)+1]] <- inter_data[1,]
      
    }
      
      if(a$XT_HS_Presence %in% FALSE){
        
      filtered_data_other[[length(filtered_data_other)+1]] <- a %>% filter(XT_HS_Presence == FALSE)    
        
      }
      
    filtered_other <- do.call("rbind.data.frame",filtered_data_other)  
      
    filtered_inter_d <- do.call("rbind.data.frame",filtered_data)
      
    filtered_inter_d <- filtered_inter_d %>% arrange(desc(weightage_score_syn))
    
    if(nrow(filtered_inter_d) > 3 & a$XT_HS_Presence %in% FALSE ){
    
    filtered_inter_data <- filtered_inter_d[1:3,] %>% select(-hs_presence,-ms_presence,-dr_presence,-xl_l_presence,-xl_s_presence,-total_present_count)  
    
    filtered_complete_data <- rbind.data.frame(filtered_inter_data,filtered_other) %>% distinct()  
    
    }
    
    else if(nrow(filtered_inter_d) <= 3 & a$XT_HS_Presence %in% FALSE){
      
      filtered_inter_data <- filtered_inter_d %>% select(-hs_presence,-ms_presence,-dr_presence,-xl_l_presence,-xl_s_presence,-total_present_count)   
      
      filtered_complete_data <- rbind.data.frame(filtered_inter_data,filtered_other) %>% distinct()  
      
    }
    
}
    
    #n_syn <- nchar(a[i,c("XT_synteny","HS_synteny","XL_L","XL_S","MS","DR")]) > 0 & !is.na(a[i,c("XT_synteny","HS_synteny","XL_L","XL_S","MS","DR")])
    
    #n_syn_Num <- length(n_syn[n_syn == TRUE])
    
    xb_grouped_list[[length(xb_grouped_list) + 1]] <- filtered_complete_data
    
    #write_xlsx(a,paste0("LOCXT/",unique_XT[i],"_file.xlsx"),col_names = TRUE)
    
  }
    
 else if(nrow(a) <= 3){
   
   xb_grouped_list[[length(xb_grouped_list) + 1]] <- a
   
 }
    
  }
  
  filtered_entire_data <- do.call("rbind.data.frame",xb_grouped_list)
  
  write.csv(filtered_entire_data,"Comparison_Synteny_Diopt_D_Final_XT.csv",row.names = F,col.names = T,quote = F)
  
  return(filtered_entire_data)
  
 #else if(nrow(a) > 3 & a$XT_HS_Presence %in% FALSE){
   
   if(!is.na(a$XL_L) & nchar(a$XL_L) != 0 | !is.na(a$XL_S) & nchar(a$XL_S) != 0){
      
     filtered_data <- list()
     
     for(j in 1: length(unique_syn_XL)){
       
       inter_data <- a %>% filter(xl_l_ref_syn == unique_syn_XL[j] | xl_s_ref_syn == unique_syn_XL[j] ) %>% mutate(
         
         hs_presence = ifelse(unique_syn_XL[j] %in% hs_ref_syn,1,0),
         dr_presence = ifelse(unique_syn_XL[j] %in% dr_ref,1,0),
         ms_presence = ifelse(unique_syn_XL[j] %in% ms_ref,1,0),
         xl_l_presence = ifelse(unique_syn_XL[j] == xl_l_ref_syn,1,0),
         xl_s_presence = ifelse(unique_syn_XL[j] == xl_l_ref_syn,1,0),
         total_present_count = hs_presence + dr_presence + ms_presence + xl_l_presence + xl_s_presence
         
       )
       
       inter_data <- inter_data %>% arrange(desc(total_present_count))
       
       filtered_data[[length(filtered_data)+1]] <- inter_data
       
     }
     
    filtered_inter_d <- do.call("rbind.data.frame",filtered_data)
     
    filtered_inter_d <- filtered_inter_d %>% arrange(desc(weightage_score_syn))
     
    filtered_inter_data <- filtered_inter_d[1:3,] %>% select(-hs_presence,-ms_presence,-dr_presence,-xl_l_presence,-xl_s_presence-total_present_count)  
     
     
   }
   
   else if(!is.na(a$DR) & nchar(a$DR) != 0){
     
     filtered_data <- list()
     
     for(j in 1: length(unique_syn_DR)){
       
       inter_data <- a %>% filter(dr_ref == unique_syn_DR[j] ) %>% mutate(
         
         hs_presence = ifelse(unique_syn_DR[j] %in% hs_ref_syn,1,0),
         dr_presence = ifelse(unique_syn_DR[j] %in% dr_ref,1,0),
         ms_presence = ifelse(unique_syn_DR[j] %in% ms_ref,1,0),
         xl_l_presence = ifelse(unique_syn_DR[j] == xl_l_ref_syn,1,0),
         xl_s_presence = ifelse(unique_syn_DR[j] == xl_l_ref_syn,1,0),
         total_present_count = dr_presence + ms_presence + xl_l_presence + xl_s_presence
         
       )
       
       inter_data <- inter_data %>% arrange(desc(total_present_count))
       
       filtered_data[[length(filtered_data)+1]] <- inter_data[1:3,]
       
     }
     
     filtered_inter_d <- do.call("rbind.data.frame",filtered_data)
     
     filtered_inter_d <- filtered_inter_d %>% arrange(desc(weightage_score_syn))
     
     filtered_inter_data <- filtered_inter_d[1:3,] %>% select(-hs_presence,-ms_presence,-dr_presence,-xl_l_presence,-xl_s_presence-total_present_count)  
     
   }
   
   else if(nchar(a$MS) == 0 & !is.na(a$MS)){
     
     filtered_data <- list()
     
     for(j in 1: length(unique_syn_MS)){
       
       inter_data <- a %>% filter(ms_ref == unique_syn_MS[j] ) %>% mutate(
         
        # hs_presence = ifelse(unique_syn_MS[j] %in% hs_ref_syn,1,0),
         dr_presence = ifelse(unique_syn_MS[j] %in% dr_ref,1,0),
         ms_presence = ifelse(unique_syn_MS[j] %in% ms_ref,1,0),
         xl_l_presence = ifelse(unique_syn_MS[j] == xl_l_ref_syn,1,0),
         xl_s_presence = ifelse(unique_syn_MS[j] == xl_l_ref_syn,1,0),
         total_present_count = dr_presence + ms_presence + xl_l_presence + xl_s_presence
         
       )
       
       inter_data <- inter_data %>% arrange(desc(total_present_count))
       
       filtered_data[[length(filtered_data)+1]] <- inter_data
       
     }
     
     filtered_inter_d <- do.call("rbind.data.frame",filtered_data)
     
     filtered_inter_d <- filtered_inter_d %>% arrange(desc(weightage_score_syn))
     
     filtered_inter_data <- filtered_inter_d[1:3,] %>% select(-hs_presence,-ms_presence,-dr_presence,-xl_l_presence,-xl_s_presence-total_present_count)  
     

   }
   
   else{
     
     
     filtered_inter_data <- a 
     
   }
   
   
   xb_grouped_list[[length(xb_grouped_list) + 1]] <- filtered_inter_data
   
   
  #}    

#else if(nrow(a) <= 3 & a$XT_HS_Presence %in% FALSE){
   
   xb_grouped_list[[length(xb_grouped_list) + 1]] <- a
   
   
 
   
#}
    
#}
  
  filtered_entire_data <- do.call("rbind.data.frame",xb_grouped_list)
  
  filtered_entire_data <- filtered_entire_data %>% select(XT_synteny,HS_synteny,XT_Gene_ID_Synteny,HS_Gene_ID_Synteny) %>% distinct()
  
  return(filtered_entire_data)
  
  write_xlsx(filtered_entire_data,file = "Comparison_Synteny_Diopt_Data_LOC_Final_XT.xlsx",col_names = T)
  
 #}

#if(data = "XL"){

unique_Xl_l_syn_d <- compare_d %>% filter(str_detect(XL_L,"LOC") | str_detect(XL_L,"provisional") | str_detect(XL_L,"XB[0-9]+")) %>% select(XL_L) %>% distinct()
unique_Xl_s_syn_d <- compare_d %>% filter(str_detect(XL_S,"LOC") | str_detect(XL_S,"provisional") | str_detect(XL_S,"XB[0-9]+")) %>% select(XL_S) %>% distinct()

total_data <- data.frame(data = c(unique_Xl_l_syn_d$XL_L,unique_Xl_s_syn_d$XL_S)) %>% distinct()

unique_XL_total <- unique(total_data$data)

unique_XL_L <- unique(unique_Xl_l_syn_d$XL_L)
unique_XL_S <- unique(unique_Xl_s_syn_d$XL_S)

xb_grouped_list <- list()

if(a$XT_HS_Presence == TRUE & nrow(a) > 3){
  
  for (i in 1 : length(unique_XL_total)) {
    
    a <- compare_d %>% filter(XL_L == unique_XL_total[i] | XL_S == unique_XL_total[i]) 
    
    a$xt_ref_syn <- str_to_lower(a$XT_synteny)
    
    a$hs_ref_syn <- str_to_lower(a$HS_synteny)
    
    a$dr_ref <- str_to_lower(a$DR)
    
    a$ms_ref <- str_to_lower(a$MS)
    
    a$xt_ref_diopt <- str_to_lower(a$XT_Diopt)
    
    a$hs_ref_diopt <- str_to_lower(a$HS_Diopt)
    
    #unique_DR_id <- unique(a$DR)
    
    #unique_HS_syn <- unique(a$hs_ref_syn)
    
    unique_HS_syn <- unique(a$hs_ref_syn)[!is.na(unique(a$hs_ref_syn)) & nchar(unique(a$hs_ref_syn)) != 0]
    
    unique_syn_XT <- unique(a$xt_ref_syn)[!is.na(unique(a$xt_ref_syn)) & nchar(unique(a$xt_ref_syn)) != 0]
    
    unique_syn_DR <- unique(a$dr_ref)[!is.na(unique(a$dr_ref)) & nchar(unique(a$dr_ref))]
    
    unique_syn_MS <- unique(a$ms_ref)[!is.na(unique(a$ms_ref)) & nchar(unique(a$ms_ref))]
    
    
    #unique_MS_id <- unique(a$MS)
    
    unique_XT_syn <- unique(a$XT_synteny)
    
    unique_XT_diopt <- unique(a$XT_Diopt)
    
    unique_HS_diopt <- unique(a$hs_ref_diopt)
    
    filtered_data <- list()
    
    if(!is.na(unique_HS_syn) & nchar(unique_HS_syn) > 0){
      
      for(j in 1: length(unique_HS_syn)){
        
        #inter_data <- unique_HS_syn[j] %in% a$hs_ref_syn | a$dr_ref == unique_HS_syn[j] | a$ms_ref == unique_HS_syn[j] | a$xt_ref_syn == unique_HS_syn[j] | a$xt_ref_diopt %in% unique_HS_diopt | a$hs_ref_diopt == unique_HS_diopt
        
        inter_data <- a %>% filter(hs_ref_syn == unique_HS_syn[j]) %>% mutate(
          
          hs_presence = ifelse(unique_HS_syn[j] %in% hs_ref_syn,1,0),
          dr_presence = ifelse(unique_HS_syn[j] %in% dr_ref,1,0),
          ms_presence = ifelse(unique_HS_syn[j] %in% ms_ref,1,0),
          #xl_l_presence = ifelse(unique_HS_syn[j] == xl_l_ref_syn,1,0),
          xt_ref_syn_presence = ifelse(unique_HS_syn[j] == xt_ref_syn,1,0),
          total_present_count = hs_presence + dr_presence + ms_presence + xl_l_presence + xl_s_presence
          
        )
        
        inter_data <- inter_data %>% arrange(desc(total_present_count))
        
        filtered_data[[j]] <- inter_data
        
      }
      
      filtered_inter_d <- do.call("rbind.data.frame",filtered_data)
      
      filtered_inter_d <- filtered_inter_d %>% arrange(desc(weightage_score_synteny))
      
      filtered_inter_data <- filtered_inter_d[1:3,] %>% select(-hs_presence,-ms_presence,-dr_presence,-xt_ref_syn_presence,-total_present_count)  
      
    }
    
    
    xb_grouped_list[[length(xb_grouped_list) + 1]] <- do.call("rbind.data.frame",filtered_inter_data)
    
    
  }
  
}  

else if(a$XT_HS_Presence == TRUE & nrow(a) <= 3 ){
  
  xb_grouped_list[[length(xb_grouped_list) + 1]] <- a
  
}

else if(a$XT_HS_Presence == FALSE & nrow(a) > 3){
  
  if(!is.na(a$XT_synteny) & nchar(a$XT_synteny) != 0){
    
    for(j in 1: length(unique_syn_XT)){
      
      inter_data <- a %>% filter(xt_ref_syn == unique_syn_XT[j] ) %>% mutate(
        
        #hs_presence = ifelse(unique_syn_XL[j] %in% hs_ref_syn,1,0),
        dr_presence = ifelse(unique_syn_XT[j] %in% dr_ref,1,0),
        ms_presence = ifelse(unique_syn_XT[j] %in% ms_ref,1,0),
        #xl_l_presence = ifelse(unique_syn_XT[j] == xl_l_ref_syn,1,0),
        xt_presence = ifelse(unique_syn_XT[j] == xt_ref_syn,1,0),
        total_present_count =  dr_presence + ms_presence + xt_presence
        
      )
      
      inter_data <- inter_data %>% arrange(desc(total_present_count))
      
      filtered_data[[length(filtered_data)]] <- inter_data
      
    }
    
    filtered_inter_d <- do.call("rbind.data.frame",filtered_data)
    
    filtered_inter_d <- filtered_inter_d %>% arrange(desc(weightage_score_synteny))
    
    filtered_inter_data <- filtered_inter_d[1:3,] %>% select(-ms_presence,-dr_presence,-xt_presence,-total_present_count)  
    
    
  }
  
  else if(!is.na(a$DR) & nchar(a$DR) != 0){
    
    for(j in 1: length(unique_syn_DR)){
      
      inter_data <- a %>% filter(dr_ref == unique_syn_DR[j] ) %>% mutate(
        
        #hs_presence = ifelse(unique_syn_DR[j] %in% hs_ref_syn,1,0),
        dr_presence = ifelse(unique_syn_DR[j] %in% dr_ref,1,0),
        ms_presence = ifelse(unique_syn_DR[j] %in% ms_ref,1,0),
        #xl_l_presence = ifelse(unique_syn_DR[j] == xl_l_ref_syn,1,0),
        xt_presence = ifelse(unique_syn_DR[j] == xt_ref_syn,1,0),
        total_present_count = dr_presence + ms_presence + xt_presence
        
      )
      
      inter_data <- inter_data %>% arrange(desc(total_present_count))
      
      filtered_data[[j]] <- inter_data
      
    }
    
    filtered_inter_d <- do.call("rbind.data.frame",filtered_data)
    
    filtered_inter_d <- filtered_inter_d %>% arrange(desc(weightage_score_synteny))
    
    filtered_inter_data <- filtered_inter_d[1:3,] %>% select(-ms_presence,-dr_presence,-xt_presence,-total_present_count)  
    
  }
  
  else if(nchar(a$MS) == 0 & !is.na(a$MS)){
    
    for(j in 1: length(unique_syn_MS)){
      
      inter_data <- a %>% filter(ms_ref == unique_syn_MS[j] ) %>% mutate(
        
        # hs_presence = ifelse(unique_syn_MS[j] %in% hs_ref_syn,1,0),
        dr_presence = ifelse(unique_syn_MS[j] %in% dr_ref,1,0),
        ms_presence = ifelse(unique_syn_MS[j] %in% ms_ref,1,0),
        xt_presence = ifelse(unique_syn_MS[j] == xt_ref_syn,1,0),
        #xl_s_presence = ifelse(unique_syn_MS[j] == xl_l_ref_syn,1,0),
        total_present_count = dr_presence + ms_presence + xt_presence 
        
      )
      
      inter_data <- inter_data %>% arrange(desc(total_present_count))
      
      filtered_data[[j]] <- inter_data
      
    }
    
    filtered_inter_d <- do.call("rbind.data.frame",filtered_data)
    
    filtered_inter_d <- filtered_inter_d %>% arrange(desc(weightage_score_synteny))
    
    filtered_inter_data <- filtered_inter_d[1:3,] %>% select(-hs_presence,-ms_presence,-dr_presence,-xt_presence,-total_present_count)  
    
    
  }
  
  else{
    
    
    filtered_inter_data <- a 
    
  }
  
  
  xb_grouped_list[[length(xb_grouped_list) + 1]] <- filtered_inter_data
  
  
}    

else if(a$XT_HS_Presence == FALSE & nrow(a) <= 3){
  
  xb_grouped_list[[length(xb_grouped_list) + 1]] <- a
  
  
}


filtered_entire_data <- do.call("rbind.data.frame",xb_grouped_list)

#write_xlsx(a,paste0("LOCXL/",unique_XT[i],"_file.xlsx"),col_names = TRUE)

write_xlsx(filtered_entire_data,file = "Comparison_Synteny_Diopt_Data_LOC_Final_XL.csv",row.names = FALSE,quote = FALSE)

#}  

#}



#}



}

}

filtered_entire_data <- data_optimization(data = "XT")

############### LOC count analysis #############

compare_df %>% filter(str_detect(XT_Diopt,"LOC") | 
                     str_detect(XT_Diopt,"provisional") | 
                      str_detect(XT_Diopt,"XB[0-9]+")) %>% 
                    select(XT_Diopt) %>% distinct() %>% count()


compare_df %>% filter(str_detect(XT_synteny,"LOC") | 
str_detect(XT_synteny,"provisional") | 
str_detect(XT_synteny,"XB[0-9]+")) %>% 
  select(XT_synteny) %>% distinct() %>% count()

compare_df %>% filter(str_detect(XL_L,"LOC") | 
                      str_detect(XL_L,"provisional") | 
                      str_detect(XL_L,"XB[0-9]+")) %>% 
                      select(XL_L) %>% distinct() %>% count() 

compare_d %>% filter(str_detect(XL_S,"LOC") | 
                     str_detect(XL_S,"provisional") |
                      str_detect(XL_S,"XB[0-9]+")) %>% 
                      select(XL_S) %>% distinct() %>% count()  


#########################################################

library(plotly)

#filtered_Diopt_Compare <- as.matrix(filtered_Both_Compare_Df[1:3,])

#abc <- filtered_Both_Compare_Df[1:3,c("ortho","Syn_Value","weightage_s")]

#m <- as.matrix(filtered_Both_Compare_Df[1:5,c("ortho","swiftortho","phylome","proteinortho","Inparanoid","Sonic","FastOrtho","weightage_s","Syn_Value")],nrow = 5, ncol = 9)

abc <- Compare_Data[1:10000,]

m <- as.matrix(abc[,c("weightage_score_diopt","weightage_score_syn","Syn_Value")],ncol = 3)

plot_ly(x = c("weightage_score_diopt","weightage_score_syn","Syn_Value"),y = abc$ref_synteny,z = m,type = "heatmap",colors = colorRamp(c("red","green","blue"))
)

#plot_ly(x = unique(filtered_Both_Compare_Df$ref_synteny)[1:5],y = c("ortho","swiftortho","phylome","portho","inp","Sonic","Fortho","Weightage","SynValue"),z = m,type = "heatmap",colors = colorRamp(c("red","green")))  

###########################################################

#### XT and HS mapping and count

overlapping_count_XT_HS <- length(unique(Compare_Data[Compare_Data$XT_HS_Presence == TRUE | str_detect(Compare_Data$XT_synteny,"LOC") | str_detect(Compare_Data$XT_synteny,"provisional") | str_detect(Compare_Data$XT_synteny,"XB[0-9]+"),]$XT_synteny))
     
total_count_synteny_XT <- length(unique(Synteny_V10_Data$XT)) 
#total_count_synteny_XL <- length(unique(Synteny_V10_Data$XL_L)) + length(unique(Synteny_V10_Data$XL_S)) 
total_count_diopt_XT <-   length(unique(Diopt_V10_Data$XT))

Unique_Synteny_Count <- total_count_synteny_XT - overlapping_count_XT_HS
Unique_Diopt_Count <- total_count_diopt_XT - overlapping_count_XT_HS

diopt_xt_loc <- Diopt_V10_Data %>% filter(str_detect(XT,"LOC") | 
                              str_detect(XT,"provisional") | 
                              str_detect(XT,"XB[0-9]+")) %>% select(XT) %>% distinct() %>% count()

syn_loc_xt <- Synteny_V10_Data %>% filter(str_detect(XT,"LOC") | 
                            str_detect(XT,"provisional") | 
                              str_detect(XT,"XB[0-9]+")) %>% select(XT) %>% distinct() %>% count()

syn_loc_xl_l <- Synteny_V10_Data %>% filter(str_detect(XL_L,"LOC") | 
                            str_detect(XL_L,"provisional") | 
                            str_detect(XL_L,"XB[0-9]+")) %>% select(XL_L) %>% distinct() %>% count() 

syn_loc_xl_s <- Synteny_V10_Data %>% filter(str_detect(XL_S,"LOC") | 
                            str_detect(XL_S,"provisional")
                             | str_detect(XL_S,"XB[0-9]+")) %>% select(XL_S) %>% distinct() %>% count()  

XT_HS_Synteny_Diopt_Data <- data.frame(description = c("overlapping_count_XT_HS","total_count_synteny","total_count_diopt","Unique_Synteny_Count","Unique_Diopt_Count","LOC_Count_Syn_XT","LOC_Count_Syn_XL_L","LOC_Count_Syn_XL_S","LOC_Count_Diopt_XT"),
                                       Value = c(overlapping_count_XT_HS,total_count_synteny_XT,total_count_diopt_XT,Unique_Synteny_Count,Unique_Diopt_Count,syn_loc_xt$n,syn_loc_xl_l$n,syn_loc_xl_s$n,diopt_xt_loc$n))

write.csv(XT_HS_Synteny_Diopt_Data,file = "XT_XL_HS_Synteny_Diopt_Data.csv",col.names = TRUE,quote = FALSE)  
  

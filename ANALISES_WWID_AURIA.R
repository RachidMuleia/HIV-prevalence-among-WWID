library(rio)
library(RDS)
library(gamlss)
library(tidyverse)
library(car)
library(ResourceSelection) # this is for goodness of fit hosmer and lemershow
library(rms) # this allows me to check for the calibration plot
library(Hmisc)
library(tidyverse)
library(sandwich)
library(lmtest)
setwd("/Users/rachidmuleia/Dropbox/INS/PID/ARTIGO_AURIA_SSR")
path
list.files(getwd())

# ---------------------- importacao de dados -x-------------------------------------

wid_df <- read.csv("DADOS_PID.csv", header = TRUE,na.strings = "")
wid_df <- wid_df |>
 filter(SEX_CAT == "2_FEMENINO")
View(wid_df)


wid_df <- wid_df |>
  mutate(
    #response_ssr_use = case_when(
      #RHCNNOW == "Sim" | CSCTTI %in% c("Nos_ultimos_3_mese", "4_6_meses") | PGMCOV1 == "Sim" | SINTOMAS_ITS == "2_NAO"  ~ "1_SIM"   ,
      #RHCNNOW == "Nao" & !(CSCTTI %in% c("Nos_ultimos_3_mese", "4_6_meses")) & PGMCOV1 == "Nao" & SINTOMAS_ITS == "1_SIM"  ~ "2_NAO"
      
      
      
      

    
    IDADE_SEX_CAT = case_when(
      IDENT43_LIMF1_AGE_MA < 15 ~ "1_<15",
      IDENT43_LIMF1_AGE_MA >= 15 ~ "2_>=15"
    ),
    
    SEXUAL_PATNERS = case_when(
      IDENT43_MSEX2_A <=1 ~ "0-1",
      IDENT43_MSEX2_A >1 & IDENT43_MSEX2_A < 3 ~ "2",
      IDENT43_MSEX2_A > 2 ~ "3+" 
    ),
    
    SELF_REPORTED_RISK = case_when(
      CSCTRS_CAT == "2_POSITIVO" ~"1_POSITIVE" ,
      VCTRISK1_CAT == "3_SEM_RESPOSTA" ~ "2_NAO_SABE",
      VCTRISK1_CAT == "1_SEM_RISCO_POUCO" ~ "3_NO_LESS_RISK",
      VCTRISK1_CAT == "2_RISCO_MODERADO_ALTO" ~ "4_MODER_HI"
    ),
    
    EDU_CAT1 = dplyr::recode(EDU_CAT1, "1_SEM_ESCOLA" = "2_PRIMARIO", "2_PRIMARIO" = "2_PRIMARIO", "3_SECUNDARIO_SUP" = "3_SECUNDARIO_SUP"),
    DERELIG_CAT = dplyr::recode(DERELIG_CAT, "1_CATOLICA" = "1_CRISTIAN", "2_PROTESTANTE" = "1_CRISTIAN", "3_MUCULMANA" = "3_MUCULMANA", "4_SEM_REL_ANIMISTA" = "4_SEM_REL_ANIMISTA")
    
  )



#  --------------- calcular AUDIT C --------------------------------------

calculate_audit_c <- function(freq_drink, typical_drinks, binge_drinking) {
  
  # Assign scores based on responses
  freq_drink_scores <- c("Nenhuma" = 0, "Mensal_ou_menos" = 1, "2_4_vezes_por_mes" = 2,
                         "2_3_vezes_por_semana" = 3, "4_ou_mais_vezes_por_semana" = 4,
                         "Semanal" = 3, "Diariamente" = 4, "Nao_sabe_ou_nao_lembra" = NA,
                         "Recusou_se_a_responder" = NA)
  
  typical_drinks_scores <- c("1_ou_2" = 0, "3_ou_4" = 1, "5_ou_6" = 2,
                             "6_a_9" = 3, "10+" = 4, "Nao_sabe_ou_nao_lembra" = NA,
                             "Recusou_se_a_responder" = NA)
  
  binge_drinking_scores <- c("Nunca" = 0, "Menos_de_uma_vez_por_mes" = 1,
                             "Pelo_menos_uma_vez_por_mes" = 2, 
                             "Pelo_menos_po_semana" = 3, 
                             "Diari_quase_todos_os_dias" = 4,
                             "Nao_sabe_nao_lembra" = NA,
                             "Recusou_se_a_responder" = NA)
  
  # Calculate the score by summing the assigned values
  score <- freq_drink_scores[freq_drink] + typical_drinks_scores[typical_drinks] + binge_drinking_scores[binge_drinking]
  
  return(score)
}


wid_df$audit_c_score <- apply(wid_df, 1, function(row) {
  calculate_audit_c(row["ALFRQ"], row["ALDAY"], row["ALBNGE"])
})


wid_df <- wid_df |>
  mutate( SEX_DROGAS = LASTREL7_CAT,
          PARTILHA_SERINGA = case_when(IDREL1 == "Sim" | IDSHARE == "Sim" ~ "1_SIM",
                                       IDSHARE == "Nao" & IDREL1 == "Nao" ~ "2_NAO",
                                       IDREL1 == "Nao" & IDSHARE == "Sim" ~ "1_SIM",
                                       IDREL1 == "Sim" & IDSHARE == "Sim" ~ "1_SIM",
                                       IDREL1 == "Nao" & is.na(IDSHARE) ~ "2_NAO"
          ),
          
          RESULTADO_ULTIMO_TESTE = case_when(
            CSCTRS == "Positivo" ~ "1_POSITIVO",
            CSCTRS == "Negativo" ~ "2_NEGATIVO",
            CSCTEV == "Nao" ~ "3_NUNCA_TESTOU"
            
            
          ),
          
          ACTUAL_GRAVIDA = case_when(
            RHTRYPRG == "Nao" ~ "2_NAO",
            RHTRYPRG == "Sim" ~ "1_SIM"
          ),
          
          HIST_GRAVIDEZ_PARTO = case_when(
            RHEVRPRG == "Sim" | NOTA26_RHPRCB == "Sim" ~ "1_SIM",
            RHEVRPRG == "Nao" ~ "2_NAO"
          ),
          
          HIST_ABORTO = case_when(
            RHPRGABO == "Sim" ~ "1_SIM",
            RHPRGABO == "Nao" | RHEVRPRG == "Nao" ~ "2_NAO"
          ),
          
          PLANO_GRAVIDEZ = case_when(
            RHPRENUM1 == "Sim" ~ "1_SIM",
            RHPRENUM1 == "Nao" ~ "2_NAO",
            RHTRYPRG  == "Nao" ~ "3_NUNCA_GRAVIDA"
          ),
          
          PROCUROU_MEDICO = case_when(
            STPHARM_CAT == "1_SIM" ~ "1_SIM",
            STPHARM_CAT == "2_NAO" ~ "2_NAO",
            SINTOMAS_ITS == "2_NAO" ~ "3_NUNCA_REPORTOU_SINTOMAS"
          ),
          
          BARREIRA_SAUDE = case_when(
            HEALTH3_CAT == "1_SIM" | HEALTH5_CAT == "1_SIM" ~ "1_SIM",
            HEALTH3_CAT == "2_NAO" & HEALTH5_CAT == "2_NAO" ~ "2_NAO",
            is.na(HEALTH3_CAT) ~ "2_NAO"
          ),
          
          UNPROTECTED_SEX = dplyr::recode(LASTREL6, "Nao" = "Sim", "Sim" = "Nao"),
          
          HIV_TEST_SIX_MONTH = case_when(
            CSCTTI %in% c("Nos_ultimos_3_mese", "4_6_meses") ~ "1_SIM",
            !(CSCTTI %in% c("Nos_ultimos_3_mese", "4_6_meses")) ~ "2_NAO"
          ),
          
          HEALTH1_CAT = case_when(
            HEALTH1 == "Sim" ~ "1_SIM",
            HEALTH1 == "Nao" ~ "2_NAO"
          ),
          ODNAX_CAT = dplyr::recode( ODNAX,
                                     "Sim" = "1_SIM", 
                                     "Nao" = "2_NAO",
                                     "Nao_sabe_ou_nao_se_lembra" = NULL,
                                     "Recusou_se_a_responder" = NULL
          ),
          TRATRDN_CAT = dplyr::recode( 
            TRATRDN,
            "Sim" = "1_SIM",
            "Nao"  = "2_NAO",
            "Recusou_se_a_responder" = NULL
            
            
          ),
          DEACT = factor(DEACT, levels = c("Sim", "Nao")),
          
          FREQ_INJEC_DROGAS = case_when(
            DRUGS6_1 == "Nao" ~ "3_NAOINJECTOU",
            ID6_FRQ %in% c("1_3_vezes_por_dia", "5+_vezes_por_dia") ~ "1_DIARIAMENTE",
            ID6_FRQ  %in% c("1_4_vezes_por_mes", "2_7_vezes_por_semana") ~ "2_ANUALMENSALSEMANAL"
          ),
          STABNF_CAT = case_when( # teve corrimento anormal ou teve ulcera  nos 12 meses antes do inquerito
            STABNF == 'Sim' | STULCM == 'Sim' | STAULC == 'SIM' ~ '1_SIM',
            STABNF == 'Nao' & STULCM == 'Nao' & STAULC == 'Nao' ~ '2_NAO',
            .default = NA
          ),
          
          
          
          STDIAG_CAT = case_when(  # alguem lhe informou que tinha ITS
            STDIAG == 'Sim' ~ '1_SIM',
            STDIAG == 'Nao' ~ '2_NAO',
            .default = NA
          ),
          
          SINTOMAS_ITS = case_when(
            STABNF == 'Sim' | STULCM == 'Sim' | STAULC == 'Sim'  ~ '1_SIM',
            STABNF == 'Nao' & STULCM == 'Nao' & STAULC == 'Nao'  ~ '2_NAO',
            STABNF == 'Nao' & is.na(STULCM) & STAULC == 'Nao' ~ '2_NAO',
            is.na(STABNF) & STULCM == 'Nao' & STAULC == 'Nao' ~ '2_NAO',
            STABNF == 'Nao' & STULCM == 'Nao' & is.na(STAULC) ~ '2_NAO',
            
            .default = NA
          ),
          
          SINTOMAS_ITS_INFO = case_when( # teve sintomas ou alguem lhe informou que tem uma ITS
            SINTOMAS_ITS == "1_SIM" | STDIAG == "Sim" ~ "1_SIM",
            SINTOMAS_ITS == "2_NAO" | STDIAG == "Nao" ~ "2_SIM",
            .default = NA
          ),
          
          
          STPHARM_CAT = case_when( # procurou aconselhamento medico quando teve sintomas de ITS
            STPHARM == 'Sim' & SINTOMAS_ITS == '1_SIM'~ '1_SIM',
            STPHARM == 'Nao' & SINTOMAS_ITS =='1_SIM' ~ '2_NAO',
            
            .default = NA
          ),
          STGPHRS_CAT = case_when( # Ja foi agredido fisicamente por ser pid
            STGPHRS == 'Sim_fisicamente' | STGPHRS == "Ambos" ~ '1_SIM',
            STGPHRS == 'Sim_varbalmente' | STGPHRS == 'Nao' ~ '2_NAO',
            
            .default = NA
          ),
          
          STGPFSX_CAT = case_when( # ja foi violado sexualmente por ser PID
            STGPFSX == 'Sim' ~ '1_SIM',
            STGPFSX == 'Nao' ~ '2_NAO',
            .default = NA
          ),
          
          STGXCLP_CAT = case_when( # teve alguma experiencia de discriminacao
            STGXCLP %in% c('Muitas_vezes', '1_vez', 'Poucas_vezes') ~ '1_SIM',
            STGXCLP == 'Nunca' ~ '2_NAO',
            .default = NA
          ),
          
          VIOLENCIA_CAT = case_when(
            STGPHRS_CAT == "1_SIM" | STGPFSX_CAT == "1_SIM"  ~ "1_SIM",
            STGPHRS_CAT == "2_NAO" & STGPFSX_CAT == "2_NAO"  ~ "2_NAO",
            .default = NA
          ),
          
          MSEX2_ACAT = case_when( # total de parceiros sexuais homens nos ultimos 12 meses
            IDENT43_MSEX2_A == 0 | LIMFVAG_M == 'Nao' ~ '0',
            IDENT43_MSEX2_A == 1 ~ '1',
            IDENT43_MSEX2_A > 1 ~ '2+'
          ),
          
          MSEX3_ACAT = case_when( # total de parceiros sexuais permanentes homens 
            IDENT43_MSEX3_A == 0 | LIMFVAG_M == 'Nao' ~ '0',
            IDENT43_MSEX3_A == 1 ~ '1',
            IDENT43_MSEX3_A > 1 ~ '2+'
          ),
          
          MSEX5_CAT = case_when( # Deu algum dinheiro em troca de sexo
            MSEX5 == 'Numero_de_parceiros' ~ '1_SIM',
            MSEX5 == 'Nenhum' | LIMFVAG_M == 'Nao' ~ '2_NAO',
            .default = NA
            
          ),
          
          MSEX6_CAT = case_when( # Recebeu dinheiro em troca de sexo
            MSEX6 == "Numero_de_parceiros" ~ "1_SIM",
            MSEX6 == "Nenhum" | LIMFVAG_M == "Nao" ~ "2_NAO",
            .default = NA
            
          ),
          
          LASTREL7_CAT = case_when( # Consumo de drogas injectáveis, drogas não injectáveis ou álcool a última vez que teve sexo
            LASTREL7 %in% c("Sim_alcool", "Sim_ambas", "Sim_drogas_nao_injetaveis") ~ "1_SIM",
            LASTREL7 == "Nao" | LIMFVAG_H == "Nao" | LIMFVAG_M == "Nao" | IDENT43_MSEX2_A == 0 | IDENT42_LIMFPART_A == 0  ~ "2_NAO",
            .default = NA
            
          ),
          
          IDNSTE_CAT1 = case_when(
            IDNSTE == "Sim" ~ "1_SIM",
            IDNSTE == "Nao" ~ "2_NAO",
            DRUGS6_1 == "Nao" ~ "3_NAO_INJECTOU "
          ),
        
          GRAVIDEZ = case_when(
            RHEVRPRG == "Sim" | RHTRYPRG == "Sim"~ "1_SIM",
            RHEVRPRG == "Nao" & RHTRYPRG == "Nao" ~ "2_NAO",
            .default = NA
          ),
          
          ACTUAL_GRAVIDA = case_when(
            RHTRYPRG == "Nao" ~ "2_NAO",
            RHTRYPRG == "Sim" ~ "1_SIM"
          ),
          ORIENTA_SEX1 = case_when(
            ORIENTA_SEX == "1_HETEROSSEXUAL" ~ "1_HETEROSSEXUAL",
            ORIENTA_SEX %in% c("2_BISSEXUAL", "3_HOMOSEX_OUTRO") ~ "2_BISEX_HOMOSEX",
            .default = NA
          ),
          MSEX7_CAT = case_when(
            MSEX7 == "Sim" ~ "1_SIM",
            MSEX7 == "Nao" | LIMFVAG_M == "Nao" ~ "2_NAO",
            .default = NA
          ),
          
          ALCOHOL1_PT = case_when(
            ALFRQ == "Nenhuma" ~ 0,
            ALFRQ == "Mensal_ou_menos" ~ 1,
            ALFRQ == "2_4_vezes_por_mes" ~ 2,
            ALFRQ == c("2_3_vezes_por_semana","Semanal") ~ 3,
            ALFRQ == c("4_ou_mais_vezes_por_semana", "Diariamente" ) ~ 4,
            .default = NA
            
          ),
          
          ALCOHOL2_PT = case_when(
            ALDAY == "1_ou_2"  ~ 0,
            ALDAY == "3_ou_4" ~ 1,
            ALDAY == "5_ou_6" ~ 3,
            ALDAY == "6_a_9" ~ 3,
            ALDAY == "10+" ~ 4,
            .default = NA
          ),
          
          ALCOHOL3_PT = case_when(
            ALBNGE == "Nunca" ~ 0,
            ALBNGE == "Menos_de_uma_vez_por_mes" ~ 1,
            ALBNGE == "Pelo_menos_uma_vez_por_mes" ~ 2,
            ALBNGE == "Pelo_menos_uma_vez_po_semana" ~ 3, 
            ALBNGE == "Diari_ou_quase_todos_os_dias" ~ 4,
            .default = NA
          )
          
  ) |> rowwise()|>
  mutate(ALCOHOL_PTSOMA = sum(c_across(ALCOHOL1_PT:ALCOHOL3_PT),na.rm = TRUE)) |>
  mutate(
    AUDIT_SCORE = case_when(
      ALCOHOL_PTSOMA < 3 ~ "2_NAO_ABUSIVO",
      ALCOHOL_PTSOMA >= 3  ~ "1_ABUSIVO",
      ALFRQ1 ==  "Nao" ~ "2_NAO_ABUSIVO",
      .default = NA
    ),
    
    IDNSTE1_CAT = case_when(
      IDNSTE1 %in% c("Facil", "Muito_facil" ) ~ "1_FACIL_MUITO_FACIL",
      IDNSTE1 %in% c("Dificil", "Muito_dificil" ) ~ "2_DIFICIL_MUITO_DIFICIL"
    ),
    
    ANOS_INJECT = ELAGEL_B - NOTA13_IDOLD_A,
    
    ANOS_INJECT_CAT = case_when(
      ANOS_INJECT == 0 ~ '<1_ANO',
      between(ANOS_INJECT, 1,3) ~ "2_1-3_ANOS",
      between(ANOS_INJECT, 4,6) ~ "3_4-6_ANOS",
      ANOS_INJECT >= 7 ~ "4_7+",
      .default = NA
      
    ),
    
    CONDOM1 = dplyr::recode(CONDOM1, "Sim" = "1_SIM" ,"Nao" = "2_NAO","Nao_sabe_ou_nao_se_lembra" = NULL)
    
)



# CONVERT THE DATA TO RDS FORMAT 

recruiter_id <- function(data, site_list, var_site) {
  data_city <- lapply(site_list,
                      \(x) {
                        data[data[[var_site]] == x, ]
                      })
  
  recruit_id <- lapply(data_city, function(x) {
    rid.from.coupons(
      data = x ,
      subject.coupon = 'coupon',
      coupon.variables = c("c1", "c2", "c3", "c4", "c5"),
      subject.id = 'CODPARTICIPANTE_dem'
    )
    
  })
  
  for (i in 1:length(data_city)) {
    data_city[[i]]$recruiter.id <- recruit_id[[i]]
  }
  
  all_df <- do.call(rbind, data_city)
  return(all_df)
  
}



# FUNCTION TO CONVERT ALL THE DATA SITE TO RDSAT FORMAT

wid_df2 <- recruiter_id(wid_df, site_list =  sort(unique(wid_df$SITE_CITY)), var_site = 'SITE_CITY') 
#View(pid_df2)


# this function converts data to rds format 
rdsat_format <- function(data, site_list, network_size, recruiter_ide,var_site,...){
  data_city <- lapply(
    site_list, 
    \(x){
      data[data[[var_site]] == x,]
    }
  )
  
  data_rdsat <- lapply(data_city, \(x){
    as.rds.data.frame(x,
                      recruiter.id =  recruiter_ide,
                      network.size = network_size,
                      ...)
  })
  
  return(data_rdsat)
}

wid_df2 <- as.data.frame(wid_df2)
for(i in 1: length(var_wid)){
  wid_df2[, var_wid[i]] <- factor( wid_df2[, var_wid[i]])
}

wid_df2$RESUL_HIV_PREV <-factor(wid_df2$RESUL_HIV_PREV)

data_list_rds <- rdsat_format(data = wid_df2,
                              site_list = sort(unique(wid_df2$SITE_CITY)),
                              network_size = 'ELINJMU',
                              recruiter_ide = 'recruiter.id',
                              var_site =  'SITE_CITY',
                              max.coupon = 5,
                              id = 'CODPARTICIPANTE_dem') 





var_wid <- c("AGE_CAT",  # this is the age
             "EDU_CAT1", # this is education level
             "DEMARSTA_CAT", # this is marital statis
             "DERELIG_CAT", # this is religion
             "IDADE_SEX_CAT", # Age at first sexual intercourse
             #"DEACT",
             "MSEX2_ACAT" # total number of sexual partners
             ,"MSEX6_CAT" # received money in exchange of sex
             ,"LASTREL6_CAT", # used condom with the last person had sex
             "CONDOM1" # ever used condom
             , "ORIENTA_SEX1", # sexual orientation
             "MSEX7_CAT", # ever had anal sex
             "AUDIT_SCORE", # abusive alcohol consumption
             "LASTREL7_CAT", # used alcohol  during the last time you had sex
             "PARTILHA_SERINGA" , # ever shared niddles
             "IDNSTE_CAT1", # used a new niddle in the last month
             "IDADE_DROGA_INJ" , # age that started using drugs
             "DRUG6_CAT", # droga injectavel com maior frequencia/ preferencia 
             "ANOS_INJECT_CAT", # for how long you inject drugs 
             "ID6_FRQ_CAT" , #  how often you inject drug in the last 30 days
             "IDNSTE1_CAT",  # access to niddles and seringes 
             "ODNAX_CAT",  # ever heard of nalaxon
             "TRATRDN_CAT", # received harm reduction services in the last 12 months
             "PEEREDU1_CAT",  # access to health services peer educator in the last 12 months
             "GRAVIDEZ",#  have you ever been pregnant
             "SINTOMAS_ITS", # symptons of STI in the last 12 months
             "STGPFSX_CAT",  # sexual violence
             "STGXCLP_CAT" # physical violence 
             )
             
    
# make all the variables factor 


# -------------------------- ANALISES DESCITIVAS -----------------------------------------------


descriptive_kp<-function(data, variables){
  #tab_list<-list()
  tab_list<-lapply(variables, function(x){t(t(table(data[,x])))})
  prop_list<-lapply(tab_list, prop.table)
  
  variavel_nam<-list()
  for(i in 1:length(tab_list)){
    variavel_nam[[i]]<-rep(variables[i], length(tab_list[[i]]))
  }
  data_table<-do.call(rbind, tab_list)
  data_prop<-do.call(rbind, prop_list)
  variavel_nam1<-unlist(variavel_nam)
  data_var<-data.frame(var_nam=variavel_nam1,Names=rownames(data_table), EST = paste(data_table[,1], '(', round(data_prop[,1]*100,1), ')'))
  return(data_var)
  
}


descriptive_analysis <- descriptive_kp(data = wid_df, variables = var_wid)

export_path <- "/Users/rachidmuleia/Dropbox/INS/PID/PAPER WWID_AURIA/TABELAS ANALISE"
export(descriptive_analysis, file = paste(export_path, "descriptive.xlsx",sep = "/"))




# compute the median and the quantiles 
quantile(wid_df$ELAGEL_B, na.rm = TRUE)

#-------------------------- analise bivariada tabulacao cruzada----------------------------------



perform_test <- function(table) {
  # Check if the table is a matrix or data frame
  if (!is.matrix(table) && !is.data.frame(table)) {
    stop("Input must be a table (matrix or data frame).")
  }
  
  # Compute the chi-squared test
  chisq_result <- chisq.test(table)
  
  # Check if any expected frequency is less than 5
  if (any(chisq_result$expected < 5)) {
    # Perform Fisher's exact test
    fisher_result <- fisher.test(table)
    return(list(test = "Fisher's Exact Test", result = fisher_result))
  } else {
    return(list(test = "Chi-Squared Test", result = chisq_result))
  }
}


cross_tabulation <- function(data, var_x, var_response, margin, digits = 1){
  tab_cross <- lapply(var_x, \(x){
    table(data[, x], data[, var_response])
  })
  
  chisq_cross <- lapply(tab_cross,  perform_test)
  
  prop_cross <- lapply(tab_cross, prop.table, margin = margin)
  
  table_ls <- list()
  for(i in 1: length(var_x)){
    col1<-paste(tab_cross [[i]][,1], " (",  round(prop_cross[[i]][,1]*100, digits = digits), ")",sep = "")
    
    col2<-paste( tab_cross [[i]][,2], " (",  round(prop_cross[[i]][,2]*100, digits =digits), ")",sep = "")
    
    table_ls[[i]] <- data.frame(var = rep(var_x[i],length(rownames(tab_cross[[i]]))),
                                category = rownames(tab_cross[[i]]),
                                col1 = col1, col2 = col2, 
                                p_value = rep(round(chisq_cross[[i]]$result$p.value,3),length(rownames(tab_cross[[i]])) ))
  }
  
  
  table_resul <- do.call(rbind, table_ls)
  names(table_resul)[c(3,4)] <- colnames(tab_cross[[1]])
  return(table_resul)
  
  
  
}



bivariate_table_HIV <- cross_tabulation(data = wid_df, var_x = var_wid, var_response = "RESUL_HIV_PREV", margin = 1) 

export(bivariate_table_HIV, file = paste(export_path, "tabela_bivariada_HIV.xlsx",sep = "/"))

xx<-cross_tabulation(data = wid_df, var_x = var_wid, var_response = "response_ssr", margin = 1)


# HCV by selected variables 
bivariate_table_HCV <- cross_tabulation(data = wid_df, var_x = var_wid, var_response = "RESUL_HCV_CAT", margin = 1) 
export(bivariate_table_HCV, file = paste(export_path, "tabela_bivariada_HCV.xlsx",sep = "/"))



# HBV by selected variables 
bivariate_table_HBV <- cross_tabulation(data = wid_df, var_x = var_wid, var_response = "RESUL_HBV_CAT", margin = 1) 
export(bivariate_table_HBV, file = paste(export_path, "tabela_bivariada_HBV.xlsx",sep = "/"))

#------------------------------------------ RDS estimates --------------------------------------------
tabulation_RDS <- function(data_list, var_outcome, dec_b, dec_rds, dec, pop_size, var_site = NULL, ...){
  
  # --- 1. Standardization: Define All Categories ---
  # Determine the complete set of unique categories across all data lists.
  # This prevents the "differing number of rows" error.
  all_categories <- unique(unlist(lapply(data_list, \(x) na.omit(unique(x[[var_outcome]])))))
  all_categories <- sort(all_categories)
  
  # --- 2. Brute/Unweighted Tabulation ---
  tab <- lapply(data_list, \(x){
    # Force the outcome column to a factor with all_categories as levels to standardize length.
    x_factor <- factor(x[[var_outcome]], levels = all_categories)
    table(x_factor)
  })
  
  tab_prop <- lapply(tab, \(x){
    round(prop.table(x) * 100, dec_b)
  })
  
  bruto_N <- lapply(tab, as.numeric)
  bruto_perc <- lapply(tab_prop, as.numeric)
  
  # --- 3. RDS Estimates and CI Calculation (Extraction) ---
  rds_estimates <- mapply(\(x,y){
    # Apply factor levels to the RDS input data for consistency
    x[[var_outcome]] <- factor(x[[var_outcome]], levels = all_categories) 
    
    RDS.bootstrap.intervals(x, outcome.variable = var_outcome, N = y, ...)
  }, data_list, pop_size, SIMPLIFY = FALSE)
  
  # Extract RAW RDS Estimate (Proportion) and SE for IVWA calculation
  raw_estimates <- lapply(rds_estimates, \(x){ x$estimate })
  raw_se <- lapply(rds_estimates, \(x){ attr(x$interval, 'bsresult')$se_estimate })
  
  # Extract RDS Percentage (Estimate)
  estimates_perc <- lapply(rds_estimates, \(x){
    round(x$estimate * 100, dec_rds)
  })
  
  # Extract Standard Error (SE) - Site-specific
  se <- lapply(rds_estimates, \(x){
    round(attr(x$interval, 'bsresult')$se_estimate * 100, dec)
  })
  
  # Extract CI Lower Bound - Site-specific
  ci_lower <- lapply(rds_estimates, \(x){
    round((x$estimate - 1.96 * attr(x$interval, 'bsresult')$se_estimate) * 100, dec)
  })
  
  # Extract CI Upper Bound - Site-specific
  ci_upper <- lapply(rds_estimates, \(x){
    round((x$estimate + 1.96 * attr(x$interval, 'bsresult')$se_estimate) * 100, dec)
  })
  
  # --- 4. INVERSE VARIANCE WEIGHTED AVERAGE (IVWA) CALCULATION ---
  # Ensure matrices are aligned by category for rowSums
  estimates_matrix <- do.call(cbind, raw_estimates)
  se_matrix <- do.call(cbind, raw_se)
  
  # Calculate Weights (W = 1 / SE^2)
  weights_matrix <- 1 / (se_matrix^2)
  
  # Calculate Pooled Prevalence and SE
  sum_P_times_W <- rowSums(estimates_matrix * weights_matrix)
  sum_W <- rowSums(weights_matrix)
  
  pooled_estimate <- sum_P_times_W / sum_W
  pooled_se <- sqrt(1 / sum_W)
  
  # Calculate Pooled CI Bounds
  pooled_ci_lower <- pooled_estimate - 1.96 * pooled_se
  pooled_ci_upper <- pooled_estimate + 1.96 * pooled_se
  
  # Format the Pooled Column (Pooled Estimate (CI Lower - CI Upper))
  pooled_prevalence_column <- paste0(
    round(pooled_estimate * 100, dec_rds), 
    ' (',
    round(pooled_ci_lower * 100, dec), 
    '-',
    round(pooled_ci_upper * 100, dec),
    ')'
  )
  
  # --- 5. Restructuring the Output (Site-specific results) ---
  all_site_data <- mapply(\(n, p, rds_est, rds_se, ci_l, ci_u){
    data.frame(
      N_BRUTO = n,
      PCT_BRUTO = p,
      RDS_PCT = rds_est,
      SE_RDS = rds_se,
      CI_LOWER = ci_l,
      CI_UPPER = ci_u
    )
  }, bruto_N, bruto_perc, estimates_perc, se, ci_lower, ci_upper, SIMPLIFY = FALSE)
  
  # *** FIX: Define all_prov_est here before it is used later ***
  all_prov_est <- do.call(cbind, all_site_data)
  
  # --- 6. Column Name Generation ---
  if (is.null(var_site)) {
    var_site <- paste0("Site", seq_along(data_list))
  }
  
  colnames_rds <- c()
  for (site in var_site) {
    colnames_rds <- c(colnames_rds, 
                      paste0(c('N_BRUTO', 'PCT_BRUTO', 'RDS_PCT', 'SE_RDS', 'CI_LOWER', 'CI_UPPER'), '_', site))
  }
  
  colnames(all_prov_est) <- colnames_rds
  
  # --- 7. Final Data Frame Assembly ---
  # Use the standardized 'all_categories' for the rows of the final output table.
  all_rds_est <- data.frame(
    VAR_NAME = rep(var_outcome, length(all_categories)),
    CAT = all_categories,
    all_prov_est,
    POOLED_RDS_IVWA = pooled_prevalence_column, # ADDED POOLED COLUMN
    stringsAsFactors = FALSE
  )
  
  return(all_rds_est)
}




RDS_many <- function(variables, data_list, ...){
  tab_all<- lapply(variables, \(x){
    tabulation_RDS(data_list = data_list, var_outcome = x, ...)
  })
  
  data_all <- do.call(rbind, tab_all)
  return(data_all)
}





info_dem <- RDS_many(variables = var_wid[6] , 
                     data_list = data_list_rds,dec_b=1,
                     dec_rds =1,
                     dec=2,
                     weight.type = "Gile's SS",
                     uncertainty = 'Gile',
                     confidence.level =0.95,
                     number.of.bootstrap.samples=15000, 
                     to.factor =TRUE, 
                     pop_size = pop_size, var_site = sort(unique(wid_df2$SITE_CITY)),
                     subset = RESUL_HIV_PREV == "1_POSITIVO"
)
lista_tabelas <-list()
 for(i in 1:length(data_list_rds)){
   for(j in 1:length(var_wid)){
     print(table(data_list_rds[[i]][,var_wid[j]], data_list_rds[[i]][, "RESUL_HIV_PREV"]))
   }
 }


RDS.bootstrap.intervals(data_list_rds[[5]],
                        outcome.variable = "RESUL_HIV_PREV",
                        weight.type = "Gile's SS",
                        uncertainty = 'Gile', confidence.level =0.95,
                        number.of.bootstrap.samples=1000,
                        subset = data_list_rds[[5]][, "AGE_CAT"] == "2_18-24",
                        N = pop_size[5])


RDS.bootstrap.intervals(data_list_rds[[5]],
                        outcome.variable = "RESUL_HIV_PREV",
                        weight.type = "Gile's SS",
                        uncertainty = 'Gile', confidence.level =0.95,
                        number.of.bootstrap.samples=1000,
                        subset = data_list_rds[[5]][, "AGE_CAT"] == "1_16-17",
                        N = pop_size[5])

prev_demo <- rds_prov_all(data_list = data_list_rds[5], var_vector = var_wid[1], var_outcome = 'RESUL_HIV_PREV',
                          pop_size = pop_size,dec=1, dec.i=1,site = 'SITE_CITY')


# A conceptual example of where the error might be in rds_prov_all:

rds_prov_all <- function(data_list, var_vector, var_outcome, pop_size, ..., site){
  all_results <- lapply(data_list, function(data) {
    # 1. Calculate prevalence for one province
    df <- rds_prev_all(data = data, vec_var = var_vector, var_outcome = var_outcome, pop_size = pop_size, ...)
    # 2. Add the site/city column
    df$SITE_CITY <- site # <<< Error likely happens if 'site' is a single string and df is NULL/empty
    return(df)
  })
  
  # 3. Combine all results (this is where the warnings come from)
  combined_df <- do.call(rbind, all_results) 
  
  # 4. Reorder columns (this is where the ERROR comes from)
  final_df <- combined_df %>% 
    relocate(SITE_CITY) # If combined_df is not a data frame, this fails
  
  return(final_df)
}
# Modify the rds_prev function as follows:
rds_prev <- function(data, x_var, var_outcome, population_size, dec, dec.i, ...){
  prevalence <- lapply(extract_lab(data = data, variable = x_var),\(x){
    # ... (RDS.bootstrap.intervals call remains the same)
  })
  
  # -----------------------------------------------------------
  # ADD CHECK for NULL/skipped result HERE
  # -----------------------------------------------------------
  
  p_value <- sapply(prevalence, \(x){
    if (is.null(x)) { 
      return(NA)  # Return NA if the prevalence calculation was skipped
    } else {
      # Only compare proportions if the result is valid
      return(RDS.compare.proportions(prevalence[[1]], x, M = 10000)[1,1])
    }
  })
  
  prev_estimate <- sapply(prevalence, \(x){
    if (is.null(x)) { 
      return(NA) # Return NA if the result is NULL
    } else {
      return(round(x$estimate[1]*100,dec))
    }
  })
  
  se_estimate <- sapply(prevalence, \(x){
    if (is.null(x)) { 
      return(0) # Return 0 SE for a 0% or 100% prevalence case
    } else {
      return(attr(x, 'bsresult')$se_estimate[1]*100)
    }
  })
  
  # ... (Rest of the function follows, but will now use NA/0 for skipped categories)
  # The 'prev_cat' calculation will also need to handle the NA values from 'prev_estimate' gracefully.
}



#------------------------------- logistic regression----------------------------------------------

# primeiro correr regressao bivariada 
wid_df <- as.data.frame(wid_df)
for(i in 1: length(var_wid)){
    wid_df2[, var_wid[i]] <- factor( wid_df[, var_wid[i]])
  }
 wid_df2$RESUL_HIV_PREV <-factor(wid_df2$RESUL_HIV_PREV)


logit_model <- list()

wid_df <- wid_df |>
  mutate(response =  case_when(
    RESUL_HIV_PREV == "2_NEGATIVO"~ 0,
    RESUL_HIV_PREV == "1_POSITIVO" ~ 1
  ),
  EDU_CAT1 = factor(EDU_CAT1),
  EDU_CAT1 = relevel(EDU_CAT1, ref = "3_SECUNDARIO_SUP"),
  
  
  
  )


var_wid1 <- c("SITE_CITY","AGE_CAT_INSIDA", "EDU_CAT1", "DEMARSTA_CAT", "DEACT",
             "IDADE_SEX_CAT", "SEXUAL_PATNERS","MSEX6_CAT","LASTREL6_CAT",
             "MSEX7_CAT","AUDIT_SCORE", "PARTILHA_SERINGA" ,"IDNSTE_CAT1",
             "IDADE_DROGA_INJ" ,"DRUG6_CAT","ANOS_INJECT_CAT" ,"ID6_FRQ_CAT" , "IDNSTE1_CAT", "ODNAX_CAT", "TRATRDN_CAT",
             "PEEREDU1_CAT", "GRAVIDEZ", "SINTOMAS_ITS" ,"STGPFSX_CAT", "STGXCLP_CAT")


robust_poisson_model <- list()
for (i in 1:(length(var_wid))){
  formula<-paste("response", var_wid[i], sep=" ~ ")
  #model_null[[i]] <- gamlss(response ~ 1, family = BI, data=na.omit(wid_df[, c('response', var_wid1[i])]))
  robust_poisson_model[[i]] <- coeftest(glm(as.formula(formula), family = poisson(link = "log"), data=na.omit(wid_df[, c('response', var_wid[i])])),
                               vcov = sandwich)
  
}

# extract OR from a glm object
extract_OR_glm <-function(x,dec_Or=2,dec=2){
  coeff<-x[,1]
  lab<-attr(x, 'dimnames')[[1]]
  li<-coeff-1.959964*x[,2]
  ls<-coeff+1.959964*x[,2]
  conf<-paste( paste( '(',round(exp(li),dec),sep=''),paste(round(exp(ls),dec),')',sep=''),sep='-')
  OR_conf<-paste(round(exp(coeff),dec_Or),conf, sep=' ' )
  data_OR<-data.frame(variabl_cat=lab, parameter=OR_conf, p_value=round(x[,4],3))
  return(data_OR)
  
}



extract_OR(model1_robust)

extract_OR<-function(x,dec_Or=1,dec=3){
  coeff<-na.omit(coef(x))
  lab<-labels(na.omit(coef(x)))
  li<-coeff-1.95*vcov(x, type="se")
  ls<-coeff+1.95*vcov(x, type="se")
  conf<-paste( paste( '(',round(exp(li),dec),sep=''),paste(round(exp(ls),dec),')',sep=''),sep='-')
  OR_conf<-paste(round(exp(coeff),dec_Or),conf, sep=' ' )
  data_OR<-data.frame(variabl_cat=lab, parameter=OR_conf, p_value=round(invisible(summary(x))[,4],3))
  return(data_OR)
  
}


extracted_OR <- lapply(
  robust_poisson_model, extract_OR_glm, dec_Or=1, dec=1)

Table_OR <- do.call(rbind, extracted_OR)



export(Table_OR, file = paste(export_path, "regressao_bivaariada_HIV_robust.xlsx",sep = "/"))






logit_model <- list()
model_null <- list()

for (i in 1:(length(var_wid1))){
  formula<-paste("response", var_wid1[i], sep=" ~ ")
  model_null[[i]] <- glm(response ~ 1, family = binomial, data=na.omit(wid_df[, c('response', var_wid1[i])]))
  logit_model[[i]] <- glm(as.formula(formula), family = binomial, data=na.omit(wid_df[, c('response', var_wid1[i])]))
  
}

anova_llrtest <-list()
for(i in 1:length(model_null)){
  anova_llrtest[[i]] <- anova(model_null[[i]], logit_model[[i]], test = "Chisq")
}



model_null <- gamlss(response ~ 1, family = BI, data=na.omit(wid_df[, c('response', var_wid1[i])]))

# function to extract OR 




# run the multivariate regression model 
var_wid2 <-  c("AGE_CAT", "SITE_CITY", "EDU_CAT1", "DEMARSTA_CAT", "DEACT",
                         "IDADE_SEX_CAT",  "LASTREL6_CAT",
                        
                          "ANOS_INJECT_CAT" ,"ID6_FRQ_CAT" , 
                         "PEEREDU1_CAT", "GRAVIDEZ", "SINTOMAS_ITS" , "STGXCLP_CAT")

"IDNSTE_CAT1"
"ANOS_INJECT_CAT"

"ID6_FRQ_CAT"

var_wid2 <-  c("AGE_CAT", "EDU_CAT1","SITE_CITY", "DEMARSTA_CAT", "DEACT",
               "IDADE_SEX_CAT",
               "ANOS_INJECT_CAT",
               "PEEREDU1_CAT", "GRAVIDEZ", "SINTOMAS_ITS" , "STGXCLP_CAT")

#"IDADE_DROGA_INJ"
eq <- as.formula(paste("response ~",paste( var_wid2, collapse = "+")))
#mult_model <- gamlss(eq, family = BI, data = na.omit(wid_df[, c("response", var_wid2)]))
mult_model <- glm(eq, family = poisson, data = na.omit(wid_df[, c("response", var_wid2)]))

mult_model_robust <- coeftest(mult_model, vcov = sandwich)

car::vif(mult_model)

summary(mult_model_robust)



# goodness of fit
null_model <- glm(response ~1 , family = poisson, data = na.omit(wid_df[, c("response", var_wid2)]))
null_model <- coeftest(null_model, vcov = sandwich)

anova(null_model, mult_model)

# check the goodness of fit using the hosmer and lemershow
HL <- ResourceSelection::hoslem.test(mult_model$y,
                                     fitted(mult_model))


HL


# check as well for calibration plot 
rms::val.prob(fitted(mult_model), mult_model$y)


extracted_OR_mult <- extract_OR_glm(mult_model_robust, dec_Or=1, dec=1)


export(extracted_OR_mult, file = paste(export_path, "regressao_multivariada_HIV_robust.xlsx",sep = "/"))

colSums(is.na(wid_df[, c("response", var_wid2)])) # check variables with large number of missings 

# These tables (network1, network2) should have been created using correlate() %>% stretch() %>% na.omit
# This function overlays two correlation tables. 
# This function also filters based on r value
# This function is very quick and dirty and not well documented.
# It is provided to the LMU Metagenomics course to streamline analysis.
# Niklas Schandry, 2020

overlay_corrnet <- function(network1 = NULL, 
                            network2 = NULL,
                            treat1 = "treat1",
                            treat2 = "treat2",
                            cutoff = 0.5 ) {

  
  full_net <- full_join(network1 %>% filter(abs(r) > cutoff),
                        network2 %>% filter(abs(r) > cutoff),
                        by = c("x", "y"), 
                        suffix =c("_1", "_2")) # This is much easier than tryin to use the treat1 and treat2 vars here, also the vars are not helpful at this stage.
  
  full_net %<>%
    mutate(corr_1 = case_when(r_1 < 0 ~ "Negative",
                              r_1 > 0 ~ "Positive",
                              TRUE ~ "NA"),
           corr_2 = case_when(r_2 < 0 ~ "Negative",
                              r_2 > 0 ~ "Positive",
                              TRUE ~ "NA") )
  
  full_net %<>%  mutate(
    edge = case_when(
      corr_1 == corr_2 ~ "No Change",
      ((corr_1 == "Negative" & corr_2 == "Positive") | (corr_1 == "Positive" && corr_2 == "Negative")) ~ "Sign change",
      (corr_2 == "NA") & (corr_1 != "NA") ~ paste(treat1, "Only"),
      (corr_2 != "NA") & (corr_1 == "NA") ~ paste(treat2, "Only"),
      TRUE ~ "check me"
    )) 
  
  full_net <- full_net %>%
    dplyr::select(x,y, edge) %>% 
    graph_from_data_frame() %>% 
    fortify
  return(full_net)
}
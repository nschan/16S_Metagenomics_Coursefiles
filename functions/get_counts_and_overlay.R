# This function creates two OTU level correlation networks, one for each treatment
# This function then overlays two correlation tables. 
# This function also filters based on r value and minimum counts.
# This function is very quick and dirty and not well documented.
# It is provided to the LMU Metagenomics course to streamline analysis.
# Niklas Schandry, 2020

get_counts_and_overlay <- function(phyl, # Phyloseq object
                                   treatvar = "sample_treatment", # Column that contains the treatments
                                   treat1, # First treatment
                                   treat2, # Second treatment
                                   OTUcol = "Genus", # Column that defines OTU
                                   mincounts = 99, # Minimum total counts per treatment to retain OTU for correlations
                                   mincorr = 0.5 # Correlation cutoff
) {
  counts <- phyl %>% 
    phyloseq::phyloseq_to_deseq2(reformulate(treatvar)) %>% 
    DESeq2::DESeq() %>% 
    DESeq2::counts() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("Amplicon")
  
  ## Add genus information
  counts %<>%  dplyr::left_join(phyl %>%    
                                  phyloseq::tax_table() %>%
                                  as.data.frame() %>%
                                  tibble::rownames_to_column("Amplicon") %>%
                                  dplyr::select(Amplicon, eval(OTUcol)), by = "Amplicon")  %>%
    # Drop amplicon column
    dplyr::select(-Amplicon) %>%
    # Use genus as rownames
    column_to_rownames(eval(OTUcol))
  
  
  counts %<>% 
    t() %>% # t() transposes. Transposition is switching rows are columns (sort of like rotating the table 90°, but not really)
    as.data.frame() %>% 
    tibble::rownames_to_column("Sample") %>% 
    dplyr::mutate(Sample = str_replace(Sample, "\\.","_")) # This should have no effect, is a remnant from course development. Also does not break anything
  
  ## We change this to do simpler subsetting
  counts %<>% dplyr::mutate(Treatment =
                              case_when(stringr::str_detect(Sample, treat1) ~ paste(treat1), 
                                        # Changed this so the treatment is detected from the sample name.
                                        # Maybe sub-optimal in general, but fine for this course..
                                        stringr::str_detect(Sample, treat2) ~ paste(treat2),
                                        TRUE ~ "Other"))
  
  # Format network 1
  
  network1 <- counts %>%
    dplyr::filter(Treatment == treat1) %>% 
    dplyr::select(-Sample, -Treatment)
  
  network1 <- network1[, colSums(network1) > mincounts] 
  
  network1 <- network1 %>% 
    corrr::correlate() %>% 
    corrr::stretch() %>% 
    na.omit()
  
  # Format network2
  
  network2 <- counts %>%
    dplyr::filter(Treatment == treat2) %>% 
    dplyr::select(-Sample, -Treatment)
  
  network2 <- network2[, colSums(network2) > mincounts] 
  
  network2 <- network2 %>% 
    corrr::correlate() %>% 
    corrr::stretch() %>% 
    na.omit()
  
  # Pass to overlay function
  
  overlay_corrnet(network1, network2, treat1 = treat1, treat2 = treat2, cutoff = mincorr)
  
}
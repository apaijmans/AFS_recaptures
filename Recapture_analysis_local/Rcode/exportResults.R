
exportResults <- function(results, msat_alleles, file_name_out){
  
  # Name list items numerically (serve later as group id)
  names(results) <- c(1:length(results))
  
  resultsDf <- bind_rows(results, .id="groups")
  
  results_df <- left_join(resultsDf, msat_alleles, by="SampleID") %>%
    mutate(groups = as.numeric(groups))
  
  
  
  #~~ Get number of mismatches and names of the loci that do not match
  mismatching_alleles <- results_df %>%
    group_by(groups) %>%
    #filter(n() > 1) %>%
    summarise(across(Pv9.a:Mang36.b, ~as.integer(n_distinct(.) == 1))) %>% # not ignoring NA
    #summarise(across(Pv9.a:Mang36.b, ~if(any(is.na(.))) NA else as.integer(n_distinct(.) == 1))) %>% # ignoring NA
    mutate(across(Pv9.a:Mang36.b, as.numeric)) 
  
  mismatches <- apply(mismatching_alleles, 1, function(x) which(x==0)  )
  
  mismatch_info <- cbind(data.frame(groups=1:length(mismatches)),
                         data.frame(Mismatch=sapply(mismatches, length)),
                         data.frame(Mismatching_alleles=sapply(mismatches, function(x) paste0(names(x), collapse = ", "))))
  
  
  
  #~~ Merge this info to the pup, then merge that back it into results
  results_pup <- results_df[!duplicated(results_df$groups),] %>% # select first row in each group (which is the sample all the others are compared to)
    full_join(mismatch_info, by="groups") %>%
    select(groups, SampleID, Incommon, Shared_loci, Matching_loci, Mismatching_alleles, PlateNumber, PlateLocation, Matches, Pv9.a:Mang36.b)
  
  results_info <- results_df %>%
    select(groups, SampleID, Incommon, Shared_loci, Matching_loci, PlateNumber, PlateLocation, Matches, Pv9.a:Mang36.b) %>%
    left_join(results_pup) %>%
    select(groups, SampleID, Incommon, Shared_loci, Matching_loci, Mismatching_alleles, PlateNumber, PlateLocation, Matches, Pv9.a:Mang36.b)
  
 
  
  #~~ Prepping for excel file 
  
  #insert blank row after each group
  results_excel <- do.call(rbind, by(results_info, results_info$groups, rbind, ""))
  
  #remove last blank row as it is pointless
  results_excel <- results_excel[1:(nrow(results_excel)-1),]
  
  
  
  #~~ Get indices for colored cells 
  
  # index of each pup (always 1st entry after empty row (+1) (+1 to allow for header in excel ) and remove last one)
  names(mismatches) <- c(1, which(results_excel[,1]=="")+1)+1
  
  indices <- unlist(mismatches)
  names(indices) <- gsub("\\..*","", names(indices) )
  
  cell_index <- data.frame(row_index = names(indices), col_index = indices)
  
  # correct for extra cols
  cell_index$col_index <- cell_index$col_index + (which(names(results_info) %in% "Pv9.a")-2) # gives index of Pv9.a in dataframe, corresponds to 2 in cell_index
  
  
  
  #~~ write this into an excel file 
  
  library(openxlsx)
  
  wb <- createWorkbook()
  
  addWorksheet(wb, "recaptures")
  writeData(wb = wb, sheet = "recaptures", x = results_excel)
  color <- createStyle(fgFill = "#9BC2E6") #BLUE
  addStyle(wb = wb, sheet = "recaptures", style = color, rows = cell_index$row_index, cols = cell_index$col_index)
  saveWorkbook(wb = wb, file = paste0(file_name_out, ".xlsx"), overwrite = TRUE)
  
}


#################
#   Libraries   #
#################

library(tidyverse)



#################
#   Variables   #
#################

# What comparison do you want to make? 
# 1. "pup-female" 
# 2. "pup-male" 
# 3. "pup-pup" 
# 4. "female-female" 
# 5. "male-male"
comp <- 2

file_date <- "_r159-168_31_03_25"
# file_name_in <- paste0("all_msat_genotypes_raw.xlsx")
# file_path_in <- "C:/Uni/03-Demography/Data/"
# # file_name_out is generated in function askComparison
# file_path_out <- "C:/Uni/03-Demography/Data/Processed/recaptures/"

file_name_in <- paste0("all_msat_genotypes_raw_r159-168.xlsx")
file_path_in <- "/grp/moleecol/Anneke/Recaptures/data/raw/"
# file_name_out is generated in function askComparison
file_path_out <- "/grp/moleecol/Anneke/Recaptures/data/processed/"

loci_9_names <- c("Pv9", "Hg6.3", "Hg8.10", "Hg1.3", "M11a", "PvcA", "Lw10", "Aa4", "PvcE")
threshold_loci_9 <- 7
threshold_total <- 23

msat <- NULL
counter <- 0
row_counter <- 0



#################
#   Functions   #
#################

# Function to prompt what comparison one wants to make
askComparison <- function(x, file_date) {
  
  cfg <- list(
    checkOnlyPuppies = FALSE,
    checkOnlyFemales = FALSE,
    checkOnlyMales   = FALSE,
    compareWithPup   = FALSE,
    compareWithFemale= FALSE,
    compareWithMale  = FALSE
  )
  
  if (x == 1) {
    cfg$checkOnlyPuppies <- TRUE
    cfg$compareWithFemale <- TRUE
    cfg$file_name_out <- paste0("results_pup_female", file_date)
    
  } else if (x == 2) {
    cfg$checkOnlyPuppies <- TRUE
    cfg$compareWithMale <- TRUE
    cfg$file_name_out <- paste0("results_pup_male", file_date)
    
  } else if (x == 3) {
    cfg$checkOnlyPuppies <- TRUE
    cfg$compareWithPup <- TRUE
    cfg$file_name_out <- paste0("results_pup_pup", file_date)
    
  } else if (x == 4) {
    cfg$checkOnlyFemales <- TRUE
    cfg$compareWithFemale <- TRUE
    cfg$file_name_out <- paste0("results_female_female", file_date)
    
  } else if (x == 5) {
    cfg$checkOnlyMales <- TRUE
    cfg$compareWithMale <- TRUE
    cfg$file_name_out <- paste0("results_male_male", file_date)
    
  } else stop("Invalid comparison")
  
  cfg
}


preprocessing <- function(msat, loci_9_names) {
  list(
    #change range indices
    range_data = 7:ncol(msat),
    #change loci indices
    loci_9 = which(names(msat) %in% loci_9_names)
  )
}


# Functions to create two datasets that will be compared with each other, eg pups against males, males against males etc.
generateData <- function(msat, cfg, mode = c("data1","data2")) {
  mode <- match.arg(mode)
  
  if (mode == "data1") {
    if (cfg$checkOnlyPuppies) return(msat[!is.na(msat$Pup) & msat$Pup == 1, ])
    if (cfg$checkOnlyFemales) return(msat[!is.na(msat$Female) & msat$Female == 1, ])
    if (cfg$checkOnlyMales) return(msat[!is.na(msat$Male) & msat$Male == 1, ])
  } else {
    if (cfg$compareWithPup) return(msat[!is.na(msat$Pup) & msat$Pup == 1, ])
    if (cfg$compareWithFemale) return(msat[!is.na(msat$Female) & msat$Female == 1, ])
    if (cfg$compareWithMale) return(msat[!is.na(msat$Male) & msat$Male == 1, ])
  }
  
  stop("Invalid dataset selection")
}


# executed in function compareData
checkNecessary <- function(id1, id2, cfg, group_of) {
  
  # no comparison necessary if id of row1 == row2 (ie don't compare a sample with itself)
  if (id1 == id2) return(FALSE)
  
  # check if IDs are already in same entry in the list of results (not necessary when comparing pup with male/female)
  # ie else if (checkOnlyPuppies != T | compareWithPup == T) in the earlier version of the script
  if (!cfg$checkOnlyPuppies || cfg$compareWithPup) {
    if (exists(id1, envir = group_of) && exists(id2, envir = group_of)) { # checks if row1 and row2 are in results
      return(group_of[[id1]] != group_of[[id2]]) # check if they are in the same group, if yes, comparison is NOT necessary
    }
  }
  
  # Otherwise comparison necessary
  TRUE
}



# executed in function compareData	
compareEntry <- function(row1, row2, loci_9, range_data,
                         threshold_loci_9, threshold_total) {
  
  v1 <- as.character(row1)
  v2 <- as.character(row2)
  
  # --- loci 9 ---
  r1 <- v1[loci_9]
  r2 <- v2[loci_9]
  
  L9 <- length(loci_9)
  
  # True matches: both non-NA and equal
  matches9 <- sum(!is.na(r1) & !is.na(r2) & r1 == r2)
  
  # Any non-match (including NA) is defined as a mismatch
  nonmatches_9 <- L9 - matches9
  
  # Maximum mismatches allowed
  allowed_mismatches_9 <- L9 - threshold_loci_9 # 9-7=2 loci 
  
  # Shared loci: both non-NA (matches + non-NA mismatches)
  shared9  <- sum(!is.na(r1) & !is.na(r2)) # number of loci shared EXCL NA
  
  if (nonmatches_9 > allowed_mismatches_9)
    return(0) # score 0 given back to compareData   
  
  
  # --- remaining loci ---
  other <- setdiff(range_data, loci_9)
  
  r1o <- v1[other]
  r2o <- v2[other]
  
  both_o <- !is.na(r1o) & !is.na(r2o)
  matches_o <- sum(r1o[both_o] == r2o[both_o]) # number of matching loci
  shared_o  <- sum(both_o) # number of loci shared
  
  if ((shared_o - matches_o) > ((length(other)) - threshold_total)) # max number of mismatches: (39-9)-23=7 loci 
    return(0) # score 0 given back to compareData
  # actual mismatches only, NAs are ignored because they are not in the shared_o vector
  
  # --- score ---
  total_shared  <- shared9 + shared_o
  total_matches <- matches9 + matches_o
  
  c(
    round(total_matches / total_shared * 100, 2), # % of matches
    total_shared, # n shared loci (so both genotyped but not necessarily a match)
    total_matches # n matches
  )
}



#key function
compareData <- function(msat1, msat2, cfg, pp,
                        threshold_loci_9, threshold_total) {
  
  results_env <- new.env(hash = TRUE, parent = emptyenv())
  group_of <- new.env(parent = emptyenv())
  
  counter <- 0L
  row_counter <- 0L
  
  for (i in seq_len(nrow(msat1))) {
    
    row1 <- msat1[i, ]
    id1 <- row1[[1]]
    
    row_counter <- row_counter + 1L
    
    # Print number of comparisons every 100 comparisons to see progress
    if (row_counter %% 100 == 0)
      message(Sys.time(), " ", row_counter, " ", counter)
    
    # compare row1 with row2
    for (j in seq_len(nrow(msat2))) {
      
      row2 <- msat2[j, ]
      id2 <- row2[[1]]
      
      # Check if comparison is necessary using function checkNecessary
      if (!checkNecessary(id1, id2, cfg, group_of)) next
      
      # Compared rows using function compareEntry
      scores <- compareEntry(
        row1, row2,
        pp$loci_9, pp$range_data,
        threshold_loci_9, threshold_total
      )
      
      if (scores[1] > 0) {
        
        if (!exists(id1, envir = results_env)) {
          
          results_env[[id1]] <- rbind(
            data.frame(
              SampleID = id1,
              Incommon = "",
              Shared_loci = "",
              Matching_loci = "",
              stringsAsFactors = FALSE
            ),
            data.frame(
              SampleID = id2,
              Incommon = paste(scores[1], " %"),
              Shared_loci = scores[2],
              Matching_loci = scores[3],
              stringsAsFactors = FALSE
            )
          )
          
          # record group membership explicitly
          group_of[[id1]] <- id1
          group_of[[id2]] <- id1
          
        }
        
        results_env[[id1]] <- rbind(
          results_env[[id1]],
          data.frame(
            SampleID = id2,
            Incommon = paste(scores[1], "%"),
            Shared_loci = scores[2],
            Matching_loci = scores[3],
            stringsAsFactors = FALSE
          )
        )
        
        counter <- counter + 1L
      }
    }
  }
  
  as.list(results_env)
  
}



#####################
#   start program   #
#####################

#~~ Read in csv data
msat_alleles <- readxl::read_excel(paste0(file_path_in, file_name_in), guess_max = 18000) %>%
  #filter(Best_genotype == 1) %>% # keep only samples with best genotypes: remove all duplicate pups and females
  filter(Pup==1 | Female==1 | Male ==1) %>% 
  #select(1, 16:93)
  select(SampleID,Pup,Female,Male,PlateNumber,PlateLocation,Best_genotype,Pv9.a:Mang36.b)

# bind cols with alleles into cols with loci
msat <- data.frame(msat_alleles[c(1:7)], 
                   mapply(paste0, 
                          msat_alleles[-c(1:7)][seq(1, ncol(msat_alleles)-7, by = 2)], 
                          sep="/",
                          msat_alleles[-c(1:7)][seq(2, ncol(msat_alleles)-6, by = 2)]) )

msat[msat == "NA/NA"] <- NA


names(msat) <- sub("\\.a", "", names(msat))

#msat[msat == "AGP20193 (REPEAT) (REPEAT)"] <- "AGP20193 (REPEAT)"
msat <- msat %>%
  mutate(SampleID = sub("\\?", "NOT", SampleID)) %>%
  mutate(SampleID = sub(" \\(REPEAT\\)", "_REPEAT", SampleID)) %>%
  mutate(SampleID = sub(" \\(POSITIVE CONTROL\\)", "_POSITIVECONTROL", SampleID)) %>%
  mutate(SampleID = sub(" \\(", "_", SampleID)) %>%
  mutate(SampleID = sub("\\)", "", SampleID)) %>%
  mutate(SampleID = gsub(" ", "", SampleID))

# Fix 1 sample ID that causes trouble
msat[which(msat$SampleID=="M1"),1] <- "dummy_M1"

# AGP01127 AGP97069 were causing trouble, because they matched 7/9 loci but not anymore when done for 39 loci
# Removed the best_genotype for the old samples so only the ones with 39 loci got compared. Manually check the old samples


# #~~ Generate some test data
# # test ind matching at 9 loci but mismatching at 2 other loci, see if it gets into results
# test_ind <- msat[which(msat$SampleID=="AGP05338"),]
# test_ind$SampleID <- "test"
# test_ind$Pup <- ""
# test_ind$Male <- "1"
# test_ind$Agaz3 <- "208/208"
# test_ind$Mang44 <- "147/147"
# 
# 
# msat <- rbind(msat[1:3,],
#               msat[which(msat$SampleID=="AGP94227"),],
#               msat[which(msat$SampleID=="AGM04020"),],
#               msat[22:48,],
#               msat[which(msat$SampleID=="AGM05031"),],
#               msat[3010:3060,],
#               msat[which(msat$SampleID=="AGM06103"),],
#               msat[which(msat$SampleID=="AGP94023"),],
#               msat[10:14,],
#               msat[which(msat$SampleID=="AGM06036"),],
#               test_ind,
#               msat[523:525,],
#               msat[which(msat$SampleID=="AGM06036"),],
#               msat[which(msat$SampleID=="AGM14060"),],
#               msat[which(msat$SampleID=="AGP05338"),],
#               msat[8750:8785,])
#
# msat <- msat[1:500,]

#~~ which comparison to make?
cfg <- askComparison(comp, file_date)


#~~ Preprocessing
pp <- preprocessing(msat %>% select(-Best_genotype), loci_9_names)

# Test data:
#msat <- msat[c(10:20, 270:280),] 

#~~ Generate datasets to compare
msat1 <- generateData(msat %>% 
                        filter(PlateNumber == "159" | 
                                 PlateNumber == "160" | 
                                 PlateNumber == "161" | 
                                 PlateNumber == "162" | 
                                 PlateNumber == "163" | 
                                 PlateNumber == "164" | 
                                 PlateNumber == "165" | 
                                 PlateNumber == "166" | 
                                 PlateNumber == "167" | 
                                 PlateNumber == "168") %>% 
                        select(-Best_genotype), cfg, mode = "data1") # only samples from last added plates


msat2 <- generateData(msat %>% 
                        #filter(Best_genotype == 1) %>% 
                        select(-Best_genotype), cfg, mode = "data2") # compare to all male genotypes



# #~~ Generate datasets to compare
# msat1_old <- generateData1(msat %>% 
#                              filter(PlateNumber == "159" | 
#                                       PlateNumber == "160" | 
#                                       PlateNumber == "161" | 
#                                       PlateNumber == "162" | 
#                                       PlateNumber == "163" | 
#                                       PlateNumber == "164" | 
#                                       PlateNumber == "165" | 
#                                       PlateNumber == "166" | 
#                                       PlateNumber == "167" | 
#                                       PlateNumber == "168") %>% 
#                              select(-Best_genotype)) # only samples from last added plates
# 
# 
# msat2_old <- generateData2(msat %>% 
#                              filter(Best_genotype == 1) %>% 
#                              select(-Best_genotype)) # compare to best genotypes
#msat1 <- msat1[9380:9400,]

# row1 <- msat2[9260,]
# row2 <- msat2[9261,]

#~~ Compare Data
results <- compareData(msat1, msat2, cfg, pp, threshold_loci_9, threshold_total) 

# results <- NULL

#~~ Save data so we do not have to rerun this part whenever we want to change postprocessing
save(results, file=paste0(file_path_out, cfg$file_name_out, ".Rdata"))

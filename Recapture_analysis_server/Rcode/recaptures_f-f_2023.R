
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
comp <- 4

file_date <- "25_09_23_incl_repeats"
# file_name_in <- paste0("all_msat_genotypes_raw.xlsx")
# file_path_in <- "C:/Uni/03-Demography/Data/"
# # file_name_out is generated in function askComparison
# file_path_out <- "C:/Uni/03-Demography/Data/Processed/recaptures/"

file_name_in <- paste0("all_msat_genotypes_raw.xlsx")
file_path_in <- "/grp/moleecol/Anneke/Recaptures/data/raw/"
# file_name_out is generated in function askComparison
file_path_out <- "/grp/moleecol/Anneke/Recaptures/data/processed/"

loci_9_names <- c("Pv9", "Hg6.3", "Hg8.10", "Hg1.3", "M11a", "PvcA", "Lw10", "Aa4", "PvcE")
threshold_loci_9 <- 7
threshold_total <- 23

msat <- NULL
counter <- 0
row_counter <- 0
counter_empty <- 0
results <- NULL 



#################
#   Functions   #
#################

# Function to prompt what comparison one wants to make
askComparison <- function(x){
  
  #x <- as.numeric(readline('What comparison do you want to make? \n 1. "pup-female" \n 2. "pup-male" \n 3. "pup-pup" \n 4. "female-female" \n 5. "male-male" \n'))
  
  if(x==1){
    checkOnlyPuppies <<- T
    checkOnlyFemales <<- F
    checkOnlyMales <<- F
    
    compareWithPup <<- F
    compareWithFemale <<- T
    compareWithMale <<- F 
    
    file_name_out <<- paste0("results_pup_female", file_date)
    
  } else if(x==2){
    checkOnlyPuppies <<- T
    checkOnlyFemales <<- F
    checkOnlyMales <<- F
    
    compareWithPup <<- F
    compareWithFemale <<- F
    compareWithMale <<- T 
    
    file_name_out <<- paste0("results_pup_male", file_date)
    
  } else if(x==3){
    checkOnlyPuppies <<- T
    checkOnlyFemales <<- F
    checkOnlyMales <<- F
    
    compareWithPup <<- T
    compareWithFemale <<- F
    compareWithMale <<- F 
    
    file_name_out <<- paste0("results_pup_pup", file_date)
    
  } else if(x==4){
    checkOnlyPuppies <<- F
    checkOnlyFemales <<- T
    checkOnlyMales <<- F
    
    compareWithPup <<- F
    compareWithFemale <<- T
    compareWithMale <<- F 
    
    file_name_out <<- paste0("results_female_female", file_date)
    
  } else if(x==5){
    checkOnlyPuppies <<- F
    checkOnlyFemales <<- F
    checkOnlyMales <<- T
    
    compareWithPup <<- F
    compareWithFemale <<- F
    compareWithMale <<- T 
    
    file_name_out <<- paste0("results_male_male", file_date)
    
  } else {print("Invalid comparison")}
  
}


# check how much needs to be deducted in new code, -13 gives negative indices
preprocessing <- function(msat) {
  
  #change range indices
  
  range_data <<- c(7:ncol(msat)) # <<- makes it available outside the function, alternative: assign("range_data", "new", envir = .GlobalEnv)
  
  
  #change loci indices
  
  loci_9 <<- which(names(msat) %in% loci_9_names)
  
  
  #change thresholds
  
  #threshold_loci_9 <<- threshold_loci_9*2
  #threshold_total <<- threshold_total*2
}



# FUnctions to create two datasets that will be compared with each other, eg pups against males, males against males etc.
generateData1 <- function(msat){
  
  if(checkOnlyPuppies == T) {
    msat1 <- msat %>% filter(Pup==1)
    
  } else if(checkOnlyFemales == T) {
    msat1 <- msat %>% filter(Female==1)
    
  } else if(checkOnlyMales == T) {
    msat1 <- msat %>% filter(Male==1)
    
  } else{
    print("Problem with dataset1")
  }
  
  return(msat1)
}

generateData2 <- function(msat){
  
  if(compareWithPup == T) {
    msat2 <- msat %>% filter(Pup==1)
    
  } else if(compareWithFemale == T) {
    msat2 <- msat %>% filter(Female==1)
    
  } else if(compareWithMale == T) {
    msat2 <- msat %>% filter(Male==1)
    
  } else{
    print("Problem with dataset2")
  }
  
  return(msat2)
}



# executed in function compareData		
checkNecessary <- function(row1, row2) {
  
  #not necessary if id of row1 == row2
  if(row1[[1]]==row2[[1]]) {
    return(F)
    
    # check if IDs are already in same entry in the list (not necessary when comparing pup with male/female)
  } else if(checkOnlyPuppies != T | compareWithPup==T){
    
    if(any(grepl(paste0("\\b",row1[[1]],"\\b"), results)) && any(grepl(paste0("\\b",row2[[1]],"\\b"), results))){ # checks if both row1 and row2 are in results
      
      if(grep(paste0("\\b",row1[[1]],"\\b"), results) == grep(paste0("\\b",row2[[1]],"\\b"), results)){ # checks if row1 and row2 are in the same sublist
        return(F)
      } else {
        return(T)}
      
    } else {
      return(T)}
    
    #otherwise necessary comparison
  } else {
    return(T)}
}



# executed in function compareData	
compareEntry <- function(row1, row2){
  
  #comparison just for all 9 loci 
  match_counter_loci_9 <- 0 
  shared_counter_loci_9 <- 0 
  mismatch_counter_loci_9 <- length(loci_9) - threshold_loci_9 #max number of mismatches: 9-7=2 loci 
  
  #comparison for samples with more than 9 loci
  match_counter_total <- 0
  mismatch_counter_total <- (range_data[length(range_data)] - range_data[1] + 1) - length(loci_9) - threshold_total #max number of mismatches: (39-9)-23=7 loci 
  
  counter_total <- 0
  
  
  #compare each value of row1 and row2 except id	
  
  #start with loci_9 values
  for(i in loci_9){
    
    if(row1[i]==row2[i] && !is.na(row1[i]) && !is.na(row2[i])) { #same and not empty then match
      match_counter_loci_9 <- match_counter_loci_9 + 1
      
    } else{
      mismatch_counter_loci_9 <- mismatch_counter_loci_9 - 1 #possible mismatch (mismatch and NA both considered mismatch): subtracts if empty or not the same value
      
      if(row1[i]!=row2[i] && !is.na(row1[i]) && !is.na(row2[i])) { #not empty, but not the same
        shared_counter_loci_9 <- shared_counter_loci_9 + 1
      }
      
      if(mismatch_counter_loci_9 < 0){ #max number of mismatches reached
        return(0) #score 0 given back to compareData 
      }
    }
  }
  
  shared_counter_loci_9 <- shared_counter_loci_9+match_counter_loci_9
  
  #continue with other values
  for(i in (range_data[1]:range_data[length(range_data)])[!range_data[1]:range_data[length(range_data)] %in% loci_9]) {
    
    if(!is.na(row1[i]) && !is.na(row2[i])) { #continue only if values in both samples
      counter_total <- counter_total + 1
      
      if(row1[i]==row2[i]) { #same (and not empty) then match
        match_counter_total <- match_counter_total + 1
        
        
      } else{
        mismatch_counter_total <- mismatch_counter_total - 1 #possible mismatch (actual mismatches only, NA's are ignored): subtracts if not the same value
        
        if(mismatch_counter_total < 0) { #max number of mismatches reached
          return(0) #score 0 given back to compareData   
        }
      }
    }
  }
  
  scores <- c(round((match_counter_loci_9+match_counter_total)/(shared_counter_loci_9+counter_total)*100, 2), #% of matches
              shared_counter_loci_9+counter_total, # n shared loci (so both genotyped but not necessarily a match)
              match_counter_loci_9 + match_counter_total) # n matches
  
  return(scores) #calculated score given back to compareData (in percentage and 3 decimals)
}

#key function
compareData <- function(){
  
  #compare all rows
  for(i in 1:nrow(msat1)) {
    
    row1 <- msat1[i,]
    
    row_counter <- row_counter + 1
    
    #print number of comparisons every 1000 comparisons to see progress
    if (row_counter %% 500 == 0){
      print(paste(Sys.time(), row_counter, counter))
    }
    
    
    #compare row1 with row2
    for(j in 1:nrow(msat2)){
      
      row2 <- msat2[j,]
      
      #function checkNecessary
      if(checkNecessary(row1, row2)==T){ 
        
        #function compareEntry
        scores <- compareEntry(row1, row2)
        
        if(scores[1] > 0){
          
          if(any(grepl(paste0("\\b",row1[[1]],"\\b"), results))){
            
            results[[grep(paste0("\\b",row1[[1]],"\\b"), results)]] <<- rbind(results[[grep(paste0("\\b",row1[[1]],"\\b"), results)]], 
                                                                              data.frame(SampleID = row2[[1]], 
                                                                                         Incommon = paste(as.character(scores[1]), " %"), 
                                                                                         Shared_loci = scores[2], 
                                                                                         Matching_loci = scores[3]))
            
            
          } else{
            results[[row1[[1]]]] <<- rbind(data.frame(SampleID = as.character(row1[[1]]), 
                                                      Incommon = as.character(""), 
                                                      Shared_loci = "", 
                                                      Matching_loci = "",
                                                      stringsAsFactors=FALSE), 
                                           data.frame(SampleID = as.character(row2[[1]]), 
                                                      Incommon = paste(as.character(scores[1]), " %"), 
                                                      Shared_loci = scores[2], 
                                                      Matching_loci = scores[3], 
                                                      stringsAsFactors=FALSE))
          }
          
          #increase counter
          counter <-  counter+1
          
        }
      }
    }
  }
  return(results<<-results)
}



#####################
#   start program   #
#####################

#~~ Read in csv data
msat_alleles <- readxl::read_excel(paste0(file_path_in, file_name_in), guess_max = 16000) %>%
  filter(Pup==1 | Female==1 | Male ==1) %>% 
  #select(1, 16:93)
  select(SampleID,Pup,Female,Male,PlateNumber,PlateLocation,Pv9.a:Mang36.b)

# bind cols with alleles into cols with loci
msat <- data.frame(msat_alleles[c(1:6)], 
                   mapply(paste0, 
                          msat_alleles[-c(1:6)][seq(1, ncol(msat_alleles)-7, by = 2)], 
                          sep="/",
                          msat_alleles[-c(1:6)][seq(2, ncol(msat_alleles)-6, by = 2)]) )

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
askComparison(comp)


#~~ Preprocessing
preprocessing(msat)

# Test data:
#msat <- msat[c(10:20, 270:280),] 

#~~ Generate datasets to compare
msat1 <- generateData1(msat)
msat2 <- generateData2(msat)

#msat1 <- msat1[9380:9400,]

# row1 <- msat2[9260,]
# row2 <- msat2[9261,]

#~~ Compare Data
results <- compareData()
# results <- NULL

#~~ Save data so we do not have to rerun this part whenever we want to change postprocessing
save(results, file=paste0(file_path_out, file_name_out, ".Rdata"))

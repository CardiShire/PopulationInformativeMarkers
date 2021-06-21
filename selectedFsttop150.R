# The three parameters that can be adjusted for the file. 
# Commented out arguments are for the .R file to be used with the file.
args = commandArgs(TRUE)

minReads <- args[1]
svmC <- args[2]
specificsample <- args[3]

selectedTopFst150 <- read.csv("~/src/Fstiminator/FstCalculations/HandTrainingData/selectedTopFst150.csv")


# Tidyverse is a collection of packages used for making data science faster and easier. 
# e1071 is a collection of packages that can do different machine learning
# classification and regression.
# Citation: H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
# https://CRAN.R-project.org/package=e107
suppressPackageStartupMessages(library(tidyverse))
library(e1071)
suppressPackageStartupMessages(library(reshape2))

# Selecting filenames from the current working directory matching a particular
# pattern in the parentheses. Currently using dates to select files based on when
# the Fst files were produced.
filenames <- list.files(pattern="*2019-11-26.fst") 

# Uniquify
#
# Given two columns makes sure they are unique
# 
# @param a column from a data frame
# @param b column from a data frame
# 
# @return the unique values from \code{a} and \code{b} in the data frame
#
# @example 
# uniquify(df$col1, df$col2)
uniquify <- function(a, b) {
  c(a,
    b) %>% unique()
}

# Remove columns from data frame that only contains zeros
#
# @param a data frame 
#
# @return data frame with the same or less columns
# 
# @examle
# remove_zero_cols(letsgowide)
#
remove_zero_cols <- function(df) {
  rem_vec <- NULL
  
  for(i in 1:ncol(df)){
    this_sum <- summary(df[,i])
    zero_test <- length(which(this_sum == 0))
    if(zero_test == 6) {
      rem_vec[i] <- names(df)[i]
    }
  }
  features_to_remove <- rem_vec[!is.na(rem_vec)]
  rem_ind <- which(names(df) %in% features_to_remove)
  df <- df[,-rem_ind]
  return(df)
}

# Opening one file and determine combinations for SVM analysis
#
# @param a single string containg a file name
#
# @return a data frame containing information about individual and replicate
# 
# @example
# onevsall(filenames[[1]])
onevsall <- function(filenames) {
  
  #make filenames into a tibble
  tibble(
    Filename=filenames
  ) -> tib
  
  #extract relevant information about filenames, e.g. sample # and replicate #
  tib %>%
    extract(Filename, c("S1", "S2"),
            regex="(S[0-9]+).*(S[0-9]+)",
            remove=FALSE) %>% 
    extract(Filename, c("R1", "R2"),
            regex="(R[1-3]+).*(R[1-3]+)",
            remove=FALSE) -> tib
  
  #make unique sample combinations based on sample #
  uniquify(
    tib$S1, tib$S2
  ) -> samplenumbers
  
  tib %>%
    separate(Filename, c("M1", "M2"), sep="_", remove=FALSE, extra="drop"
    ) -> tib
  
  uniquify(
    tib$M1, tib$M2
  ) -> samplenames 
  
  #update dataframe with column names and select useful columns
  bind_rows(
    select(tib, Filename, Sample=M1, Individual=S2, Replicate = R2),
    select(tib, Filename, Sample=M2, Individual=S1, Replicate = R1)
  ) -> longtibby
  
  separate(
    longtibby, Sample, c("Ss", "Bs", "T1", "R"), 
    sep = "-", remove = FALSE
  )-> longtibby
  
  subset(
    longtibby, R != Replicate
  ) -> longtibby
  
  subset(
    longtibby, 
    select = c("Filename", "Sample", "Individual", "Replicate")
  ) -> longtibby
  
  #expand grid makes unique combinations of sample and two individuals based on sample replicate
  expand.grid(
    samplenames, samplenumbers, 
    samplenumbers, stringsAsFactors = FALSE
  ) -> comb_filenames
  
  #make sure the replicate # does not match
  subset(
    comb_filenames, Var2 != Var3
  ) -> comb_filenames
  
  comb_filenames %>% 
    filter(Var2 < Var3
    ) -> comb_filenames
  
  comb_filenames %>% 
    extract(Var1, "R1", 
            regex="(R[1-3])", remove=FALSE
    ) -> comb_filenames
  
  comb_filenames %>% 
    mutate(
      Rvar2 = ifelse(
        R1 == "R1", "R2", ifelse(
          R1 == "R2", "R3", "R1"))) %>% 
    mutate(Rvar3 = ifelse(
      R1 == "R1", "R3", ifelse(
        R1 == "R2", "R1", "R2"))
    ) -> comb_filenames
  
  #join table multiple times to make all possible combinations of one sample
  #versus two individuals with the correct replicate numbers
  left_join(
    comb_filenames, longtibby, 
    by = c("Var1" = "Sample", 
           "Var2" = "Individual", 
           "Rvar2" = "Replicate")) %>%
    left_join(longtibby, by = c("Var1" = "Sample", 
                                "Var2" = "Individual", 
                                "Rvar3" = "Replicate")) %>%
    left_join(longtibby, by = c("Var1" = "Sample", 
                                "Var3" = "Individual", 
                                "Rvar2" = "Replicate")) %>%
    left_join(longtibby, by = c("Var1" = "Sample", 
                                "Var3" = "Individual", 
                                "Rvar3" = "Replicate")
    ) -> allComb
  
  #update column names to reflect output matrix information and the four files
  #names ot be used in each matrix
  colnames(allComb) <- c("Sample", "Sample_Replicate", 
                         "Individual_1", "Individual_2", 
                         "Replicate_A", "Replicate_B", 
                         "Filename_1", "Filename_2", 
                         "Filename_3", "Filename_4")
  paste(
    allComb$Individual_1, "_", 
    allComb$Replicate_A, " vs ", 
    allComb$Individual_2, "_", 
    allComb$Replicate_B, sep = ""
  ) -> allComb$IndComparison
  
  return(allComb)
}

# Make a list of lists containing samples from pairs of individuals
#
# @param a list of four files
# @param a string containing sample of interest
#
# @return makes a list of list for all combinations of four files for 
# two individuals compared to each other 
# 
# @example
# matrixlist(filenames, "S001-Hp-T1-R1")
#
matrixlist <- function(filenames, selectedsample = specificsample){
  
  #Using function `onevsall` on filenames
  #Filenames contains a list of file names
  onevsall(
    filenames
  ) -> allComb
  
  #onesample <- unique(allComb[1])
  
  specificsample <- allComb %>% filter(Sample == selectedsample)
  
  #Selecting just the columns with the filenames that are needed to make a matrix
  select(
    specificsample, c(Filename_1, Filename_2, 
                      Filename_3, Filename_4)
  ) -> matrixfilenames
  
  #Making a list of lists, aka as a large list, of the filenames needed 
  #to makes all the possible matrix needed for ML
  lapply(
    as.list(1:dim(
      matrixfilenames)[1]), 
    function(x) matrixfilenames[x[1],]
  ) -> ml
  
  return(ml)
}
# read one file in the list of four filenames
#
# @param a string
#
# @return a data frame containing Marker, Nucleotide Position, Reads_A, Reads_B, 
# Fst, Individual_A, Replicate_A, Individual_B, Replicate_B
#
# @example
# read_one_file(filenames[[1]])
#
read_one_file <- function(file) { 
  
  suppressMessages(read_table2(
    "MarkersForTargetedPanel.txt"
  )) -> MarkersForTargetedPanel
  
  read.delim2(
    "MicrobeNames.txt"
  ) -> MicrobeNames
  
  cols( 
    X1 = col_character(), X2 = col_double(), 
    X3 = col_double(), X4 = col_double(), 
    X5 = col_double()
  ) -> col_types
  
  read_delim(
    file, skip=1, delim=" ", 
    col_names = FALSE, 
    col_types = col_types
  ) -> temporarydata
  
  colnames(temporarydata) <- c(
    "Marker", "NucleotidePosition", 
    "Reads_A", "Reads_B", "Fst")  
  
  temporarydata$Filename <- file 
  
  temporarydata$sum_squares <- log10(
    (temporarydata$Reads_A - temporarydata$Reads_B)^2)
  
  temporarydata$sum_reads <- (
    temporarydata$Reads_A + temporarydata$Reads_B)
  
  temporarydata$read_dif <- abs(
    temporarydata$Reads_A - temporarydata$Reads_B)
  
  MicrobeNames$Clade <- as.character(
    MicrobeNames$Clade)
  
  left_join(
    temporarydata, MarkersForTargetedPanel, 
    by = "Marker"
  ) -> temporarydata
  
  left_join(
    temporarydata, MicrobeNames, 
    by = "Clade"
  ) -> temporarydata
  
  separate(
    temporarydata, "Filename", 
    c("Individual_A", "Bodysite_A", 
      "Timepoint_A", "Replicate_A", 
      "Individual_B", "Bodysite_B", 
      "Timepoint_B", "Replicate_B", 
      "Year", "Month", "Day", "fst", 
      "gz"), sep = "([\\-\\_\\.])", 
    remove = FALSE
  ) -> temporarydata
  
  subset(
    temporarydata, 
    select=-c(Year, Month, Day, fst, gz)
  ) -> temporarydata
  
  paste(
    temporarydata$Individual_A, "_", 
    temporarydata$Replicate_A, " vs ", 
    temporarydata$Individual_B, "_", 
    temporarydata$Replicate_B, sep = ""
  ) -> temporarydata$IndComparison
  
  return(temporarydata) 
}

# makes a list of lists with the filenames
ml <- matrixlist(filenames)

# read_four_files takes four file names and select SNPs based on 150 
# SNP predetermined panel
#
# @param one list in the list of lists of filenames
# @param integer
#
# @return matrix with selected SNPs
#
# @example
# input read_four_files(ml[[2]])
read_four_files <- function(ml, readmin = minReads) {
  
  xyz <- lapply(ml, read_one_file)
  
  bind_rows(
    xyz, .id = "column_label"
  ) -> bound
  
  
  suppressWarnings(right_join(bound, selectedTopFst150, 
                              by = c("Marker", "NucleotidePosition"))
  ) -> boundnew
  
  
  boundnew %>% 
    filter(Reads_A >= as.integer(readmin) | 
             Reads_B >= as.integer(readmin)
    ) -> boundnew
  
  subset(
    boundnew, select = c(
      IndComparison, MarkerPosition, 
      Marker, Fst)
  ) -> subsetData
  
  # subsetData %>% 
  #   group_by(IndComparison, Marker) %>% 
  #   arrange(desc(Fst)) %>% 
  #   slice(1:topNum
  #   ) -> top10
  # 
  # top10 %>% 
  #   ungroup(
  #   ) -> top10
  
  
  expand(
    subsetData, MarkerPosition, IndComparison
  ) -> x
  
  bound %>%
    mutate(Marker1 = Marker,
           Nucleotide1 = NucleotidePosition
    ) -> bound
  
  bound %>%
    unite(MarkerPosition,
          Marker1:Nucleotide1, sep= "_"
    ) -> bound
  
  suppressWarnings(right_join(
    subsetData, x, by = c(
      "IndComparison", "MarkerPosition"))
  ) -> AddNA
  
  suppressWarnings(left_join(
    AddNA, bound, by = c(
      "IndComparison", "MarkerPosition"))
  ) -> fillNA
  
  fillNA %>%
    subset(select = c(
      "IndComparison", "MarkerPosition",
      "Fst.y")
    ) -> fillNA
  
  separate(
    fillNA, "IndComparison",
    c("x", "y", "z", "w", "v"),
    sep = "([\\_\\ ])", remove = FALSE
  ) -> fillNA
  
  fillNA$NoReplicate <- paste(fillNA$x, "_", fillNA$w)
  
  fillNA %>%
    subset(select = c(
      "IndComparison", "MarkerPosition", "NoReplicate",
      "Fst.y")
    ) -> fillNA
  
  fillNA[is.na(fillNA)] <- 0
  
  fillNA[fillNA<0] <- 0
  
  matrix <- fillNA
  
  return(matrix)
}

#matrixnum <- read_four_files(ml[[3]], nFst, minReads)

# svm_model classifies unknown data point to one of two classes
# selected SNPs for the four finals are entered in the matrix
#
# @param list of lists
# @param integer
#
# @return svm output with binary classification and predication confidence
#
# @example
# svm_model(ml, svmcost = 100)
svm_model <- function(ml, svmcost = svmC) {
  
  suppressMessages(as.data.frame(
    lapply(ml, read_four_files)
  )) -> matrix
  
  matrix <- as.data.frame(ml)
  suppressMessages(as.data.frame(
    lapply(ml[1], read_four_files)
  )) -> matrix
  
  matrix %>% 
    spread(key = MarkerPosition , value = Fst.y
    ) -> letsgowide
  #removed when we removed marker as groupoing for nFST  
  # remove_zero_cols(
  #   letsgowide) -> letsgowide
  
  samples <- letsgowide[1]
  
  addtoallcomb <- as.data.frame(matrix(samples$IndComparison, ncol = 2, byrow = TRUE, dimnames = list(NULL, c('IndComparison_1', 'IndComparison_2'))))
  
  
  letsgowide <- letsgowide[,-1]
  training <- letsgowide[1:4,]
  
  as.list(
    rep(0, ncol(training))
  ) -> distribution
  
  data.frame(
    rbind(training, distribution)
  ) -> unknownrowadded
  
  colnames(unknownrowadded)<-colnames(training)
  
  rownames(unknownrowadded) <- c()
  
  unknownrowadded$NoReplicate[5] = "unknown"
  
  training$NoReplicate = factor(training$NoReplicate)
  
  test <- unknownrowadded[5,]
  
  x <- subset(training, select = -NoReplicate)
  y <- training$NoReplicate
  
  svm(
    NoReplicate ~ ., data = training, 
    type = "C", kernel = "linear", 
    scale = FALSE, probability = TRUE, 
    cross = 0, cost = svmcost
  ) -> model
  
  print(model)
  summary(model)
  
  predict(
    model, x,
    probability = TRUE, 
  ) -> pred
  
  predict(
    model, test[,-1],
    probability = TRUE
  ) -> testpred
  
  table(
    pred, y
  ) -> predtable
  
  table(
    testpred, true = test[,1]
  ) -> testpredtable
  
  test_pred <- attr(testpred, "probabilities")
  
  # list(model,
  #   samples, pred, test_pred
  #   ) -> return_list
  predic <- test_pred
  
  # predi <- do.call(rbind, predic)
  predi <- as_tibble(predic)
  
  long <- gather(predi, key = comparison, value = percentage)
  long$binary <- if_else(long$percentage > 0.5, 1, 0)
  long <- cbind(long, addtoallcomb)
  
  # binarydf <- cbind(addtoallcomb, long)
  # binarydf$svmcost <- svmcost
  # return(long)
  
  return(long)
}


# combined all SVM output with useful information about parameters
#
# @param list of lists
#
# @return data frame
#
# @example
# cp(ml)
cp <- function(ml2){
  
  allthingsrandom <- seq_along(ml2) %>% map(~svm_model(ml2[.]))
  
  combinedpred <- suppressWarnings(map_dfr(allthingsrandom, ~.x))
  
  filenames <- seq_along(ml2) %>% ml2[.]
  combinedpred$nFst <- "selectedTopFst150"
  combinedpred$minReads <- minReads
  combinedpred$svmcost <- svmC
  combinedpred$sample <- specificsample
  combinedpred$n <- nrow(combinedpred)
  
  binarycount <- combinedpred %>% filter(binary == 1) %>% count(comparison)
  binarycount <- summarise(group_by(combinedpred, comparison,binary),count =n())
  out <- merge(combinedpred, binarycount)
  
  meanout <- out %>% group_by(comparison, binary) %>% summarise(mean_class=mean(percentage))
  meanout <- out %>% group_by(comparison, binary) %>% summarise(mean_class=mean(percentage))
  sdout <- out %>% group_by(comparison, binary) %>% summarise(sd_class=sd(percentage))
  statsout <- merge(meanout, sdout)
  outadd <- merge(out, statsout)
  
  return(outadd)
}

# save results to file
write.csv(cp(ml), row.names = FALSE)
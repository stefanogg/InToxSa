# Script to remove duplicate mutations
myarg1 <- commandArgs()[6]
# if (length(commandArgs()) == 7) myarg2 <- commandArgs()[7]
source("~/R_functions/concatenate_snippy_functions.R")
remove_duplicate_mutations(file = myarg1,
                           write_file = T)
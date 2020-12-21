# [2_runScript.r]

# ------------------------------------------- ABOUT THE DATA ------------------------------------------
# 	- about:
#		> ...
#
#	- variables:
#		> ...

# 	- NOTE: ...

##-create a function that pauses R script till key is pressed                                                                                                                                                                                                                                                                                                            
wait <- function()
{                                                                                                                                                                                                                                                                                                            
	line <- readline(prompt = "\n\n\nPress [enter] to continue...")                                                                                                                                                                                                                                                                  
} 

## ---------------- load packages -----------------------------------------------------------
require(tidyr);  require(dplyr);  require(lattice)
require(latticeExtra);  require(ggplot2);  require(ggpubr)
require(tictoc);  require(MASS);  require(car)
require(UsingR);  require(extrafont); require(matahari)

message("\n | 1. Run regression model...")
message(" --------------------------------------------------------------------")
# ------------------------------------------- PART-1 --------------------------------------------------
# [1-read-data_I.r]
dance_start(value = FALSE, contents = FALSE)
source("1_regScript.r")
dance_save("./college_major_analysis.rds")
options(scipen = 0)
wait()

message("\n | 2. Test lm() with manual procedure (CASE I)...")
message(" --------------------------------------------------------------------")
# ------------------------------------------- PART-2 --------------------------------------------------
print(r.fit1$coef)
message("\n\n")
print(my.fit1$coefficients)
wait()

message("\n | 3. Varify model validation (CASE I)...")
message(" --------------------------------------------------------------------")
# ------------------------------------------- PART-3 --------------------------------------------------
print(mv1)
message("\n\n")
print(VIF)
wait()

message("\n | 4. Test lm() with manual procedure (CASE II)...")
message(" --------------------------------------------------------------------")
# ------------------------------------------- PART-4 --------------------------------------------------
print(r.fit2$coef)
message("\n\n")
print(my.fit3$coefficients)
wait()

# ------------------------------------------- PART-5 --------------------------------------------------
message("\n\n\n --------------------------------------------------------------------")
message(" Task complete...")
message(" Removing R objects from memory and clearing .mustache directory...")
message(" --------------------------------------------------------------------")
#clear_stash()
rm(list = ls())
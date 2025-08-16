

# Define required packages
required_packages <- c('lme4', 'readr', 'qcc', 'dplyr')

# Install and load required packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)  # Corrected: install.packages function
    library(pkg, character.only = TRUE)  # Load the package after installation
  } else {
    library(pkg, character.only = TRUE)  # Load the package if it's already installed
  }
}


# Load the parallel package

library(parallel)
options(warn = -1)



library(lme4)
library(dplyr)
library(readr)
library(qcc)

source('/Users/mehranmoazeni/Documents/SIMULATION_REPEAT/second_synth.R')

source('/Users/mehranmoazeni/Documents/SIMULATION_REPEAT/shift_imposer.R')
source('/Users/mehranmoazeni/Documents/SIMULATION_REPEAT/evaluator.R')


# Detect the number of cores
num_cores <- detectCores()

# Print the number of cores
print(paste("Number of cores available:", num_cores))


#nonpws = read_csv('sim_data/NonPWS_synthetic_2.csv.csv')
PWS = read_csv('/Users/mehranmoazeni/Documents/SIMULATION_REPEAT/simulation_data/merged_synth.csv')
# Example in R
#colnames(PWS) <- trimws(colnames(PWS))  # Removes leading and trailing spaces from column names
# nonpws$Complication_type<-NULL
# colnames(nonpws)[1]<-'ID'
colnames(PWS)[2]<-'ID'

#colnames(nonpws)[5]<-'norm_speed'
colnames(PWS)[6]<-'norm_speed'
PWS <- PWS %>% mutate_all(~replace(., is.na(.), 0))

PWS<- PWS[,c('ID', 'norm_speed', 'HCT', 'Flow', 'Motor_power')]
#nonpws <- nonpws %>% mutate_all(~replace(., is.na(.), 0))

# Set seed for reproducibility
#set.seed(145686)



# print(length(unique(pws$ID)))
# print(length(unique(nonpws$ID)))

trim_data <- function(data, id_column) {
  data %>%
    group_by(!!sym(id_column)) %>% 
    mutate(row_number = row_number()) %>% 
    filter(row_number <= 360) %>%
    select(-row_number) %>%  # Remove the auxiliary row_number column
    ungroup()
}


S1<- c(10,60,120,180, 240)
S2<-c(4,8,14,30)
N1<- c(5,10,25,50)

Var<-c('Flow', 'Motor_power' )

L1<- c(7, 5)
L2<- c(3.5,4.5)
LAMBDA<-c(0.2, 0.8)

#position<- 332
position<- 28
SHIFT<- c(-0.1, -0.2, -0.3, -0.4, -0.5)
SHIT_MP<- c(-0.1, -0.2, -0.3, -0.4, -0.5)
TYPE<-c("mean_change", 'random',"slow_trend_nonl", "slow_trend_l")



i = 0
J=0
print(paste0("Total simulation is ", length(S1)*length(S2)*length(LAMBDA)*length(L1)*length(L2)*length(TYPE)*length(SHIFT)*length(N1)*length(Var)))  

# Initialize a variable to store the start time of the entire process
start_time_total <- Sys.time()


library(doParallel)
library(foreach)
library(qcc)
# Set up the parallel backend
num_cores <- 6
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Start the parallel loop
foreach(v = Var, .packages = c('dplyr')) %:%
  foreach(s1 = S1) %:%
  foreach(s2 = S2) %:%
  foreach(n1 = N1) %dopar% {
    
    start_time_iter <- Sys.time()
    
    filename <- paste0("/Users/mehranmoazeni/Library/CloudStorage/OneDrive-UniversiteitUtrecht/1-PROJECTS/SIMULATION/code/Simulation/cloudcode/residuals2/",
                       "df", v, "_", s1, "_", s2, "_", n1, ".csv")
    df <- read.csv(filename)
    library(qcc)
    library(dplyr)
    library(readr)
    for (shift in SHIFT) {
      for (type in TYPE) {
        for (lambda in LAMBDA) {
          
          df_shifted <- shift_imposer(df, position = position, shift_mag = shift, var = v, type = type, lambda = lambda)
          print("I have given a shift")
          
          for (l1 in L1) {
            for (l2 in L2) {
              
              result_filename <- paste0("/Users/mehranmoazeni/Library/CloudStorage/OneDrive-UniversiteitUtrecht/1-PROJECTS/SIMULATION/code/Simulation/cloudcode/results_april_b/",
                                        "df_results", v, "_", s1, "_", s2, "_", n1, "_", l1, "_", l2, "_", lambda, "_", type, "_", shift, ".csv")
              
              if (file.exists(result_filename)) {
                next
              }
              position = 28
              df_appended <- ooc_optimized(df_shifted, L1 = l1, L2 = l2, lambda = lambda, Variable_name = v, s1 = s1)
              df_appended<- df_appended%>%group_by(ID)%>%slice(304:360)%>%ungroup
              df_evaluate <- evaluator(df_appended, s1, position = position, var = v)
              
              write.csv(df_evaluate, file = result_filename, row.names = FALSE)
              print(paste0("Finished file: ", result_filename))
            }
          }
        }
      }
    }
    
    end_time_iter <- Sys.time()
    duration_iter <- end_time_iter - start_time_iter
    print(paste0("Time taken for iteration: ", duration_iter, " seconds"))
  }

# Stop the parallel backend
stopCluster(cl)



for (v in Var) {
  for (s1 in S1) {
    for (s2 in S2) {
      for (n1 in N1) {
          start_time_iter <- Sys.time()
          
       # Define the filename
        filename <- paste0("/Users/mehranmoazeni/Library/CloudStorage/OneDrive-UniversiteitUtrecht/1-PROJECTS/SIMULATION/code/Simulation/cloudcode/residuals2/", "df", v, "_", s1, "_", s2, "_", n1, ".csv")
          df <- read.csv(filename)
        

        for (shift in SHIFT) {
          for (type in TYPE) {
            for (lambda in LAMBDA) {
                
                
                
              df_shifted <- shift_imposer(df, position = position, shift_mag = shift, var = v, type = type, lambda  = lambda)
              print("I have given a shift")
                


              for (l1 in L1) {
                for (l2 in L2) {
                    
                    filename <- paste0("/Users/mehranmoazeni/Library/CloudStorage/OneDrive-UniversiteitUtrecht/1-PROJECTS/SIMULATION/code/Simulation/cloudcode/results_april/", "df_results", v, "_", s1, "_", s2, "_", n1, "_", l1, "_", l2, "_", lambda, "_", type, "_", shift, ".csv")

                 if (file.exists(filename)) {
                     i = i+1
                     next 
                     
                     }
                    #print(filename)

                    df_appended <- ooc_optimized(df_shifted, L1 = l1, L2 = l2, lambda = lambda, Variable_name = v, s1 = s1)
                    df_evaluate <- evaluator(df_appended, s1, position = position, var = v)
  
                    write.csv(df_evaluate, file = paste0("/Users/mehranmoazeni/Library/CloudStorage/OneDrive-UniversiteitUtrecht/1-PROJECTS/SIMULATION/code/Simulation/cloudcode/results_april/", "df_results", v, "_", s1, "_", s2, "_", n1, "_", l1, "_", l2, "_", lambda, "_", type, "_", shift, ".csv"), row.names = FALSE)
                    print(paste0("Loop ", i, " is finished"))
                    i <- i + 1
                  }
                }
            }
          }
        }

        # Record the end time of the current loop iteration
        end_time_iter <- Sys.time()

        # Calculate and print the duration of the current loop iteration
        duration_iter <- end_time_iter - start_time_iter
        print(paste0("Time taken for iteration ", i, ": ", duration_iter, " seconds"))
      }
    }
  }
}

print(J)

#Record the end time of the entire process
end_time_total <- Sys.time()

# # Calculate and print the total duration of the entire process
duration_total <- end_time_total - start_time_total

duration_total



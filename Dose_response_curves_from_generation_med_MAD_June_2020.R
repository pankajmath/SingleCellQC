library(readxl)
library(dtw)
library(ggplot2)
library(fpca)
library(plotly)
library(fdasrvf)
library(readr)
library(fda.usc)
library(emdist)
library (plyr)
library(reshape2)
library(roahd)
library(dplyr)

# packageurl <- "https://cran.r-project.org/src/contrib/Archive/roahd//roahd_1.4.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")

#### Defined Functions to be used ###
#####################################
get_normalization_parameters <- function(x, mean_std, trim_tails, percentile_trim = 4){
  if (trim_tails){
    u_lim <-quantile(x, probs = (100- percentile_trim)/100, na.rm = TRUE)
    l_lim <-quantile(x, probs = percentile_trim/100, na.rm = TRUE)
    trimmed_x <- subset(x, x > l_lim & x < u_lim)
    
    if (mean_std){
      mu <- mean(trimmed_x)
      std <- sd(trimmed_x)
    } else {
      # get the median and MAD of Veh
      mu <- median(trimmed_x)
      # cat(sheet, '  Median: ', mu, '\n')
      std <- mad(trimmed_x, constant=cts)
      # cat(sheet, '  MAD: ', std, '\n')
    }
    
  } else {
    if (mean_std){
      mu <- mean(x)
      std <- sd(x)
    } else {
      # get the median and MAD of Veh
      mu <- median(x)
      # cat(sheet, '  Median: ', mu, '\n')
      std <- mad(x, constant=cts)
      # cat(sheet, '  MAD: ', std, '\n')
    }
  }
  return(unlist(list(mu, std)))
}

##### Define Functions #######
## Warping function and inverse warping function
warping_function <- function(align){
  phi_x <- align$index1
  phi_y <- align$index2
  f = matrix(data = NA, nrow = 1, ncol = max(phi_x))
  for (idx in 1:max(phi_x)){
    f[idx] = round(median(phi_y[which(phi_x == idx)]))
  }
  return(f)
}

inv_warping_function <- function(align){
  phi_x <- align$index2
  phi_y <- align$index1
  f = matrix(data = NA, nrow = 1, ncol = max(phi_x))
  for (idx in 1:max(phi_x)){
    f[idx] = round(median(phi_y[which(phi_x == idx)]))
  }
  return(f)
}

# function to obtain the warped quantile data from all the
# sheets for a given compound
drug_quantile_data <- function(quantile_data, drug_name, sheet_names, reference){
  comb_drug <- list()
  for (sheet in sheet_names){
    all_drugs <- colnames(quantile_data[[sheet]])
    if (is.element(drug_name, all_drugs)){
      Veh <- quantile_data[[sheet]]$Veh
      drug <- quantile_data[[sheet]][[drug_name]]
      normal <- reference
      
      # n_Veh <- (Veh - mean(Veh))/sd(Veh)
      # n_veh_normal <- dtw(n_Veh, normal, k = TRUE)
      
      n_veh_normal <- dtw(Veh, normal, k = TRUE)
      
      # print(c(n_veh_normal$normalizedDistance, n_veh_normal$distance))
      
      drug_Veh <- dtw(drug, Veh, k = TRUE)
      
      # print(c(drug_Veh$normalizedDistance, drug_Veh$distance))
      
      phi<- warping_function(n_veh_normal)
      
      psi_drug <- warping_function(drug_Veh)
      
      t_Veh <- normal[phi]
      t_drug <- t_Veh[psi_drug]
      
      comb_drug[[sheet]] <- t_drug
    }
  }
  return(comb_drug)
  
}

#### EMD distance between two curves curve
emd_dist_between_curves <- function(grid, curve1, curve2){
  c1 = t(rbind(grid,curve1))
  c2 = t(rbind(grid,curve2))
  return(emd(c1, c2))
} 

#### EMD distance from Median curve
emd_dist_from_median_curve <- function(data, fmed_all){
  # dist_mc = vector("list", length = 5)
  dist_mc <- setNames(as.list(matrix(NA, nrow = 1, ncol = length(data))),
                      names(data))
  for (ii in names(data)){
    samp_curve <- t(rbind(grid,data[[ii]]))
    dist_mc[[ii]]= emd(samp_curve, fmed_all)
  }
  return(dist_mc)
} 

# function to obtain the warped quantile data from all the
# sheets for a given compound
quantile_data_per_sheet <- function(quantile_data, sheet_name, reference){
  comb_drug <- list()
  all_drugs <- colnames(quantile_data[[sheet_name]])
  Veh <- quantile_data[[1]]$Veh
  normal <- reference
  for (drug_name in all_drugs){
    drug <- quantile_data[[sheet_name]][[drug_name]]
    # n_Veh <- (Veh - mean(Veh))/sd(Veh)
    
    # n_veh_normal <- dtw(n_Veh, normal, k = TRUE)
    n_veh_normal <- dtw(Veh, normal, k = TRUE)
    
    drug_Veh <- dtw(drug, Veh, k = TRUE)
    
    phi<- warping_function(n_veh_normal)
    
    psi_drug <- warping_function(drug_Veh)
    
    t_Veh <- normal[phi]
    t_drug <- t_Veh[psi_drug]
    
    comb_drug[[drug_name]] <- t_drug
  }
  
  return(comb_drug)
  
}



# function to obtain the warped quantile data from one sheet at a time when each sheet has its own vehicle

quantile_data_per_sheet_with_veh <- function(quantile_data, sheet_name, reference){
  comb_drug <- list()
  all_drugs <- colnames(quantile_data[[sheet_name]])
  Veh <- quantile_data[[sheet_name]]$Veh
  normal <- reference
  for (drug_name in all_drugs){
    drug <- quantile_data[[sheet_name]][[drug_name]]
    # n_Veh <- (Veh - mean(Veh))/sd(Veh)
    
    # n_veh_normal <- dtw(n_Veh, normal, k = TRUE)
    n_veh_normal <- dtw(Veh, normal, k = TRUE)
    
    drug_Veh <- dtw(drug, Veh, k = TRUE)
    
    phi<- warping_function(n_veh_normal)
    
    psi_drug <- warping_function(drug_Veh)
    
    t_Veh <- normal[phi]
    t_drug <- t_Veh[psi_drug]
    
    comb_drug[[drug_name]] <- t_drug
  }
  
  return(comb_drug)
  
}

##### user defined parameters and variables######
############################

## obtain the quantiles of standard normal
# data_frame <- read.csv('Q://psingh//Research//Pankaj_singlecellQC//Fabio//Data_Fabio//ER_heterogeneity_Veh_qunatiles.csv', header = F)
reference <- read.csv('D://QC_paper//standard_normal_quantiles.csv', header = F)# standard normal quantiles
# reference <- data_frame[,10] # standard normal quantiles
# rm(data_frame)

## load all the data from the excel file
data_path <- 'D://QC_paper//'
# file_name <- 'AllExperiments_Data_Pheno_PS.xlsx'
file_name <- 'All_plus_New_ExperimentsPankaj_PS.xlsx'
# get all the sheet names
sheet_names <- excel_sheets(paste(data_path, file_name, sep = ''))
# remove the first sheet name from the list(if this sheet does not contain
# single cell data)
# sheet_names<- sheet_names[2:length(sheet_names)]

num_sheets <- length(sheet_names) # number of sheets

# create a list containing all the data identified with sheet name
comb_data <- list()
for (sheet in sheet_names){
  comb_data[[sheet]] <-read_excel(paste(data_path, file_name, sep = ''), sheet)
}
# correct the column name in Exp 11 from 4-OHT to 4OHT
#colnames(comb_data$EXP11_ER)[8]<- '4OHT'

# constant value for MAD calculation
cts = 1.4826 # default constant value is 1.4826
mean_std <- FALSE # if FALSE, median and MAD are used instead of mean and standard deviation
trim_tails <- TRUE # if TRUE, tails are trimmed for mean and standard deviation calculation
percentile_trim <- 4 # percentile to trim from both ends

#### obtain the normalized Vehicles and quantile data (from 1st to 99th quantiles)
normed_data <- list()
quantile_data <- list()
quant_list <- seq(0.01, 0.99, 0.01) # sequence of 1st to 99th quantile
grid <- seq(1,99,1)

for (sheet in sheet_names){
  df <- comb_data[[sheet]]
  # remove the empty cells or NAs from Veh
  x = na.omit(df$Veh)
  parm <- get_normalization_parameters(x, mean_std, trim_tails, percentile_trim)
  mu <- parm[1]
  std <- parm[2]
  
  # list of all the compounds in the current sheet
  d_names <- colnames(df)
  
  normed_df <- list()
  # placeholder dataframe of appropriate size
  quantile_df <- data.frame()[1:length(quant_list), ] 
  # normalize each compound with mean and std from Veh of that experiment
  # also obtain the quantile values for the raw data for each compound
  for (ii in d_names){
    y <- na.omit(df[[ii]]) # remove the empty cells or NAs
    z <- (y - mu)/std # z-normalized with respect to Vehicle
    normed_df[[ii]] <- z
    quantile_df[[ii]] <- as.vector(quantile(z, probs = quant_list))
  }
  # collect z-normalized values and quantile values (from raw data)
  normed_data[[sheet]] <- normed_df
  quantile_data[[sheet]] <- quantile_df
  rm(df, normed_df, quantile_df, mu, std, x, y, ii)
}



Veh_data <- drug_quantile_data(quantile_data, 'Veh', sheet_names, reference)

## get the pairwise EMD distances for all the vehicle curves

veh_mat = do.call(rbind, Veh_data)
exp_names <- rownames(veh_mat)

emd_mat = matrix(NA, nrow = length(exp_names), ncol = length(exp_names))
rownames(emd_mat)<- exp_names

for (ii in 1:length(exp_names)){
  
  sam1 = veh_mat[ii,]
  
  for (jj in 1:length(exp_names)){
    
    sam2 = veh_mat[jj,]
    emd_mat[ii,jj] = emd_dist_between_curves(grid, sam1, sam2)
    
  }
}

## cluster them and cut tree with three clusters
dend <- hclust(dist(emd_mat, method = "euclidean"), method = "complete")
groups <- cutree(dend, k = 3)

# get the vehicle warped quantiles for group 1
group_1_Veh = Veh_data[names(groups[groups ==1])]
veh_mat = do.call(rbind, group_1_Veh) 

#### Functional Median curve 
fD = fData( grid, veh_mat )
fmed <- median_fData(fD)
fmed_all<- t(rbind(seq(1,fmed$P),fmed$values))




#### Distance from functional median curve
Veh_emd = emd_dist_from_median_curve(Veh_data, fmed_all)




#### Dose response curves from new data set; each sheet in the excel sheet has a different compound
#### except for the first sheet (which contains the vehicle and other controls)

# new_df <- read_excel('./Raw Data/MCF7/EXP31_Plate1.xlsx')
setwd("D://QC_paper//2020_0605_ATSDR_DRC//")
fname <- '2020_0605_ATSDR_MCF7_ER_DRC_Plate1_Pankaj.xlsx'

sheets <- excel_sheets(fname)
new_df <- list()
for (sheet in sheets){
  new_df[[sheet]] <-read_excel(fname, sheet)
}


df <- new_df[[sheets[1]]]
# remove the empty cells or NAs from Veh
x = na.omit(df$Veh)
# get the nomalization parmeters (mu and std)
parm <- get_normalization_parameters(x, mean_std, trim_tails, percentile_trim)
mu <- parm[1]
std <- parm[2]
cat(mu, std, '\n')

# list of all the compounds in the current sheet
d_names <- colnames(df)

new_quantile_data= list()

for (sheet in sheets){
  df <- new_df[[sheet]]
  # list of all the compounds in the current sheet
  d_names <- colnames(df)
  
  normed_df <- list()
  quantile_df <- data.frame()[1:length(quant_list), ]
  # normalize each compound with mean and std from Veh of that experiment
  # also obtain the quantile values for the raw data for each compound
  for (ii in d_names){
    y <- na.omit(df[[ii]]) # remove the empty cells or NAs
    z <- (y - mu)/std
    normed_df[[ii]] <- z
    quantile_df[[ii]] <- as.vector(quantile(z, probs = quant_list))
  }
  # collect z-normalized values and quantile values (from raw data)
  # normed_data[[sheet]] <- normed_df
  new_quantile_data[[sheet]] <- quantile_df
  # rm(df, normed_df, quantile_df, mu, std, x, y, ii)
}



# create figure directory if it does not exist
if (trim_tails) {
  if (mean_std){
    fig_dir_root = file.path('.', 'mean_std_emd_trimmed_tails')
  } else{
    fig_dir_root = file.path('.', 'median_mad_emd_trimmed_tails')
  }
  
} else{
  if (mean_std){
    fig_dir_root = file.path('.', 'mean_std_emd')
  } else{
    fig_dir_root = file.path('.', 'median_mad_emd')
  }
}


if (!dir.exists(fig_dir_root)){
  dir.create(fig_dir_root)
} else {
  print("Dir already exists!")
  
}
# create subdirectory for each experiment if it does not exist
sub_dir = tools::file_path_sans_ext(basename(fname))
dest_path <- file.path(fig_dir_root, sub_dir)

if (!dir.exists(dest_path)){
  dir.create(dest_path)
} else {
  print("Dir already exists!")
}

#dest_path <- 'D://psingh//Research//Pankaj_singlecellQC//Plots//DRC_3/'

ndict <- c(1.00E-24, 1.00E-21,
           1.00E-18, 
           1e-15,  1e-12, 
           1e-09,  1e-06)
names(ndict) <-c('yM', 'zM', 
                 'aM',
                 'fM','pM', 
                 'nM', 'uM')

#### EMD distance from Median curve
emd_dist_from_median <- function(data, fmed_all){
  # dist_mc = vector("list", length = 5)
  dist_mc <- list()
  for (ii in names(data)){
    samp_curve <- t(rbind(grid,data[[ii]]))
    dist_mc[[ii]]= emd(samp_curve, fmed_all)
  }
  return(dist_mc)
} 

all_drcs = list()
for (sheet in sheets){
  test <- quantile_data_per_sheet(new_quantile_data, sheet, reference)
  save_image_filename <- file.path(dest_path, paste(sheet, '_DRC.png', sep =''))
  ab <-emd_dist_from_median(test, fmed_all)
  dose_emd <- ldply (ab, data.frame)
  colnames(dose_emd) = c('Dose', 'EMD_dist')
  dose_emd$dose <-sapply(dose_emd$Dose, 
                         function(x) parse_number(x)*ndict[unlist(strsplit(x, parse_number(x)))[2]])
  dose_emd$Cmpd <- sheet
  all_drcs[[sheet]] <- dose_emd
  q<-ggplot(dose_emd, aes(x = reorder(Dose, dose), EMD_dist))+
    geom_point() +
    labs(x = 'Dose', y = 'EMD distance', 
         title = sheet)+ ylim(0,1.9)
  # ggplotly(q)
  print(q)
  
  ggsave(filename = save_image_filename)
  print(sheet)
}
ddf<- data.frame(Reduce(rbind, all_drcs))
ddf
csv_name <-  paste(tools::file_path_sans_ext(basename(fname)), '_EMD_distances.csv', sep ="")
write.csv(ddf, file.path(fig_dir_root, csv_name) , row.names =ddf$Cmpd)


##############################
######################
######## well to well Variation case
###############
##################


fname <- './Raw Data/Exp_20_Veh_O_wells.xlsx'
# fname <- './All_plus_New_ExperimentsPankaj_PS.xlsx'
sheets <- excel_sheets(fname)
new_df <- list()
for (sheet in sheets){
  new_df[[sheet]] <-read_excel(fname, sheet)
}


new_quantile_data= list()

for (sheet in sheets){
  df <- new_df[[sheet]]
  # list of all the compounds in the current sheet
  d_names <- colnames(df)
  
  normed_df <- list()
  quantile_df <- data.frame()[1:length(quant_list), ]
  
  x = na.omit(df$Veh)
  # get the nomalization parmeters (mu and std)
  parm <- get_normalization_parameters(x, mean_std, trim_tails, percentile_trim)
  mu <- parm[1]
  std <- parm[2]
  
  # normalize each compound with mean and std from Veh of that experiment
  # also obtain the quantile values for the raw data for each compound
  for (ii in d_names){
    # ### well to well variation case
    # x = na.omit(df[[ii]])
    # # get the nomalization parmeters (mu and std)
    # parm <- get_normalization_parameters(x, mean_std, trim_tails, percentile_trim)
    # mu <- parm[1]
    # std <- parm[2]
    
    ## otherwise
    y <- na.omit(df[[ii]]) # remove the empty cells or NAs
    z <- (y - mu)/std
    normed_df[[ii]] <- z
    quantile_df[[ii]] <- as.vector(quantile(z, probs = quant_list))
  }
  # collect z-normalized values and quantile values (from raw data)
  # normed_data[[sheet]] <- normed_df
  new_quantile_data[[sheet]] <- quantile_df
  # rm(df, normed_df, quantile_df, mu, std, x, y, ii)
}



# create figure directory if it does not exist
if (trim_tails) {
  if (mean_std){
    fig_dir_root = file.path('.', 'mean_std_emd_trimmed_tails')
  } else{
    fig_dir_root = file.path('.', 'median_mad_emd_trimmed_tails')
  }
  
} else{
  if (mean_std){
    fig_dir_root = file.path('.', 'mean_std_emd')
  } else{
    fig_dir_root = file.path('.', 'median_mad_emd')
  }
}


if (!dir.exists(fig_dir_root)){
  dir.create(fig_dir_root)
} else {
  print("Dir already exists!")
  
}
# create subdirectory for each experiment if it does not exist
sub_dir = tools::file_path_sans_ext(basename(fname))
dest_path <- file.path(fig_dir_root, sub_dir)

if (!dir.exists(dest_path)){
  dir.create(dest_path)
} else {
  print("Dir already exists!")
}

#dest_path <- 'D://psingh//Research//Pankaj_singlecellQC//Plots//DRC_3/'

ndict <- c(1.00E-24, 1.00E-21,
           1.00E-18, 
           1e-15,  1e-12, 
           1e-09,  1e-06)
names(ndict) <-c('yM', 'zM', 
                 'aM',
                 'fM','pM', 
                 'nM', 'uM')

#### EMD distance from Median curve
emd_dist_from_median <- function(data, fmed_all){
  # dist_mc = vector("list", length = 5)
  dist_mc <- list()
  for (ii in names(data)){
    samp_curve <- t(rbind(grid,data[[ii]]))
    dist_mc[[ii]]= emd(samp_curve, fmed_all)
  }
  return(dist_mc)
} 

all_drcs = list()
for (sheet in sheets){
  
  # # when the first sheet has vehicle and control compounds, rest of the sheets have different doses/compounds
  # test <- quantile_data_per_sheet(new_quantile_data, sheet, reference)
  
  # when each sheet has its own vehicle
  test <- quantile_data_per_sheet_with_veh(new_quantile_data, sheet, reference)
  
  save_image_filename <- file.path(dest_path, paste(sheet, '_DRC.png', sep =''))
  ab <-emd_dist_from_median(test, fmed_all)
  dose_emd <- ldply (ab, data.frame)
  colnames(dose_emd) = c('Dose', 'EMD_dist')
  dose_emd$dose <-sapply(dose_emd$Dose, 
                         function(x) parse_number(x)*ndict[unlist(strsplit(x, parse_number(x)))[2]])
  dose_emd$Cmpd <- sheet
  all_drcs[[sheet]] <- dose_emd
  q<-ggplot(dose_emd, aes(x = reorder(Dose, dose), EMD_dist))+
    geom_point() +
    labs(x = 'Dose', y = 'EMD distance', 
         title = sheet)+ ylim(0,1.6)
  # ggplotly(q)
  print(q)
  
  ggsave(filename = save_image_filename)
  print(sheet)
}
ddf<- data.frame(Reduce(rbind, all_drcs))

csv_name <-  paste(tools::file_path_sans_ext(basename(fname)), '_EMD_distances.csv', sep ="")
write.csv(ddf, file.path(fig_dir_root, csv_name) , row.names =ddf$Cmpd )

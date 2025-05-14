#### Bog-Fen Redox Code


# Environment Set up: -----------------------------------------------------
library(tidyverse)
#Helpful Functions
#Pareto Scaling:
PS_helper <- function(x) {
  (x - mean(x)) / sqrt(sd(x, na.rm = T))
}	

pareto_scale <- function(x){
  mtb_scaled <- data.frame(apply(x, 2, PS_helper))
  return(mtb_scaled)
}

#Auto Scaling:
AS_helper <- function(x) {
  (x - mean(x)) / sd(x, na.rm = T)
} 

auto_scale <- function(x){
  mtb_scaled <- apply(x, 2, AS_helper) 
  return(mtb_scaled)
}

#Log Transformation Functions:
log_helper <- function(x, min.val) {
  log2((x + sqrt(x ^ 2 + min.val ^ 2)) / 2)
}

#Log Scaling:
log_transform <- function(x){
  x_nz <- x[ ,which(apply(x, 2, sum) != 0)]
  min.val <- min(abs(x_nz[x_nz!=0]))/10
  x_log_trans <- data.frame(apply(x_nz, 2, log_helper, min.val))
  return(x_log_trans)
}

#Coefficient of Variation
cv <- function(x){
  sd(x)/mean(x)
}

#Color palettes
pal_site <- c('#009900', '#3333FF')
names(pal_site) <- c('bog', 'fen') 

pal_sxd <- c('#2E603D', '#06B25C', '#402BCA', '#20A9FE')
names(pal_sxd) <- c('bog_ti', 'bog_tf', 'fen_ti', 'fen_tf')


# LC-MS Exploratory Analysis: ---------------------------------------------

## Data Load ---------------------------------------------------------------

### HILIC -------------------------------------------------------------------

#Load in duplicate features which need to be removed from the data
hilic_bad <- read_lines('Data/LCMS/HILIC/HILIC_Duplicate_Features.txt')
hilic_IS_frags <- read_lines('Data/LCMS/HILIC/IS_Remove_HILIC.txt')
hilic_badFormulas <- read_lines('Data/LCMS/HILIC/hilic_badFormula.txt')

#Make the data ready for analysis
hilic_raw <- read_csv('Data/LCMS/HILIC/HILIC_Data_wQCs.csv') %>%
  #Remove norm area from file names
  rename_with(~gsub('Norm. Area: ', '', .x)) %>%
  #Remove .raw from file names
  rename_with(~gsub('.raw', '', .x)) %>%
  #Remove the file number of file names
  rename_with(~gsub(" \\(F[[:digit:]]+\\)", '', .x)) %>%
  #Remove teh duplicate features:
  filter(!Compound_ID %in% hilic_bad) %>%
  filter(!Compound_ID %in% hilic_IS_frags) %>%
  filter(!Compound_ID %in% hilic_badFormulas) %>%
  #Make the compound IDs the row names
  column_to_rownames('Compound_ID') %>%
  #Transpose the matrix to make it so samples are rows and compounds are columns
  t()  %>%
  as.data.frame() %>%
  #Add FT at the front of every compound to make life easier
  rename_with(~paste0('FT_', .x)) 

testsamp <- sample(1:ncol(hilic_raw), 100)

#Grab QCs
hilic_qcs <- hilic_raw %>%
  rownames_to_column('sample') %>%
  filter(str_detect(sample, 'QC')) %>%
  filter(str_detect(sample, 'MS1')) %>%
  column_to_rownames('sample')

#Keep good CVs
hilic_qc_cv <- apply(hilic_qcs, 2, cv)
hilic_keep <- names(hilic_qc_cv[hilic_qc_cv < 0.3])

#Finalized data to keep
hilic_clean <- hilic_raw %>%
  log_transform() %>%
  pareto_scale() %>%
  dplyr::select(all_of(hilic_keep))


### RP ----------------------------------------------------------------------

#Load in duplicate features which need to be removed from the data
rp_bad <- read_lines('Data/LCMS/RP/RP_Duplicate_Features.txt')
rp_IS_frags <- read_lines('Data/LCMS/RP/IS_Remove_RP.txt')
rp_badFormulas <- read_lines('Data/LCMS/RP/rp_badFormula.txt')

#Load in the data
rp_raw <- read_csv('Data/LCMS/RP/RP_Data.csv') %>%
  #Remove norm area from file names
  rename_with(~gsub('Norm. Area: ', '', .x)) %>%
  #Remove .raw from file names
  rename_with(~gsub('.raw', '', .x)) %>%
  #Remove the file number of file names
  rename_with(~gsub(" \\(F[[:digit:]]+\\)", '', .x)) %>%
  #Remove duplicate features
  filter(!Compound_ID %in% rp_bad) %>%
  filter(!Compound_ID %in% rp_IS_frags) %>%
  filter(!Compound_ID %in% rp_badFormulas) %>%
  #Make the compound IDs the row names
  column_to_rownames('Compound_ID') %>%
  #Transpose the matrix to make it so samples are rows and compounds are columns
  t()  %>%
  as.data.frame() %>%
  rename_with(~paste0('FT_', .x))

#Check the distribution - Clearly bad
graphics.off()
#plot.new()
par(mfrow = c(1,1))
par(mar = c(1,1,1,1))
boxplot(rp_raw[testsamp])

#Grab QCs
rp_qcs <- rp_raw %>%
  rownames_to_column('sample') %>%
  filter(str_detect(sample, 'QC')) %>%
  filter(str_detect(sample, 'MS1')) %>%
  column_to_rownames('sample')

#Set up filtering criteria
rp_qc_cv <- apply(rp_qcs, 2, cv)
rp_keep <- names(rp_qc_cv[rp_qc_cv < 0.3])

#Normalize the data:
rp_clean <- rp_raw %>%
  log_transform() %>%
  pareto_scale() %>%
  dplyr::select(all_of(rp_keep))

#Much better
boxplot(rp_raw[testsamp])
boxplot(rp_clean[testsamp])

#Load in metadata
meta_hilic <- read_csv('Data/LCMS/HILIC/metadata_HILIC.csv') %>%
  #Fix the day type for fen
  mutate(Day = ifelse(Site == 'fen' & Day == 'D28', 'D35', Day))

meta_rp <- read_csv('Data/LCMS/RP/metadata_RP.csv') %>%
  #Fix the day type for fen
  mutate(Day = ifelse(Site == 'fen' & Day == 'D28', 'D35', Day))

#Combine metadata with metabolite data
rp_final <- rp_clean %>%
  rownames_to_column('Filename') %>%
  left_join(meta_rp) %>%
  filter(str_detect(Filename, 'D'))  %>%
  modify_at('Day', factor, levels = c('D0', 'D7', 'D14', 'D21', 'D28', 'D35','D42')) %>%
  modify_at('Site', factor, levels = c('fen', 'bog'))


hilic_final <- hilic_clean %>%
  rownames_to_column('Filename') %>%
  left_join(meta_hilic) %>%
  filter(str_detect(Filename, 'D')) %>%
  modify_at('Day', factor, levels = c('D0', 'D7', 'D14', 'D21', 'D28', 'D35','D42')) %>%
  modify_at('Site', factor, levels = c('fen', 'bog'))


## Exploratory PCAs --------------------------------------------------------
library(mixOmics)

### RP ----------------------------------------------------------------------

#Extract out just metabolite data to computea  PCA
rp_control <- rp_final %>%
  filter(Status == 'Alive') %>%
  filter(Day %in% c('D0', 'D7', 'D14', 'D28', 'D35')) %>%
  filter(Treatment == 'unamended')

X_rp_control <- rp_control %>%
  column_to_rownames('Filename') %>%
  dplyr::select(starts_with('FT_')) 

meta_rp_control <- rp_control %>%
  dplyr::select(-starts_with('FT_')) %>%
  modify_at('Day', factor, levels = c('D0', 'D7', 'D14', 'D28', 'D35'))

pca_rp_control <- pca(X_rp_control)
#Some seperation by Site but clearly something else at play here
rp_elip <- plotIndiv(pca_rp_control, group = rp_control$Site, ind.names = F, legend = T, ellipse = T)$df.ellipse
plotIndiv(pca_rp_control, group = rp_control$Day, ind.names = F, legend = T)
rp_pca_data <- plotIndiv(pca_rp_control)$df %>%
  cbind(meta_rp_control)

ggplot(rp_pca_data, aes(x, y, color = Day, shape = Site)) +
  geom_point(size = 6)  +
  xlab('PC1: 59%') +
  ylab('PC2: 12%') +
  geom_polygon(data = rp_elip, aes(x = Col1, y = Col2), inherit.aes = F, color = pal_site[2], fill = NA, size = 1) +
  geom_polygon(data = rp_elip, aes(x = Col3, y = Col4), inherit.aes = F, color = pal_site[1], fill = NA, size = 1) +
  ggtitle('RP') +
  cowplot::theme_cowplot()

### HILIC -------------------------------------------------------------------

hilic_control <- hilic_final %>%
  filter(Status == 'Alive') %>%
  filter(Day %in% c('D0', 'D7', 'D14', 'D28', 'D35')) %>%
  filter(Treatment == 'unamended')

X_hilic_control <- hilic_control %>%
  column_to_rownames('Filename') %>%
  dplyr::select(starts_with('FT_')) 

meta_hilic_control <- hilic_control %>%
  dplyr::select(-starts_with('FT_')) %>%
  modify_at('Day', factor, levels = c('D0', 'D7', 'D14', 'D28', 'D35'))

pca_hilic_control <- pca(X_hilic_control)
#Some seperation by Site but clearly something else at play here
hilic_elip <- plotIndiv(pca_hilic_control, group = hilic_control$Site, ind.names = F, legend = T, ellipse = T)$df.ellipse
plotIndiv(pca_hilic_control, group = hilic_control$Day, ind.names = F, legend = T)
hilic_pca_data <- plotIndiv(pca_hilic_control)$df %>%
  cbind(meta_hilic_control)

ggplot(hilic_pca_data, aes(x, y, color = Day, shape = Site)) +
  geom_point(size = 6)  +
  #scale_color_manual(values = day_pal) +
  xlab('PC1: 52%') +
  ylab('PC2: 17%') +
  geom_polygon(data = hilic_elip, aes(x = Col1, y = Col2), inherit.aes = F, color = pal_site[2], fill = NA, size = 1) +
  geom_polygon(data = hilic_elip, aes(x = Col3, y = Col4), inherit.aes = F, color = pal_site[1], fill = NA, size = 1) +
  ggtitle('HILIC') +
  cowplot::theme_cowplot() 


## Differential Abundance Testing ------------------------------------------

#Refactor to categorical to make it easier:
lm_data_rp <- rp_control %>%
  filter(Day %in% c('D0', 'D7', 'D14', 'D28', 'D35')) %>%
  mutate(t_cat = ifelse(Day == 'D0', 'T0', 
                        ifelse(Day == 'D7', 'T1',
                               ifelse(Day == 'D14', 'T2', 'T3')))) %>%
  modify_at('t_cat', factor, levels = c('T0', 'T1', 'T2', 'T3')) %>%
  pivot_longer(starts_with('FT_'), names_to = 'ft', values_to = 'area') %>%
  group_by(ft) %>%
  nest() 

lm_rp <- lm_data_rp %>%
  mutate(mod = purrr::map(data, function(x) lm(area ~ Site + t_cat + Site*t_cat, data = x)))

meta_rp2 <- meta_rp %>%
  mutate(t_cat = ifelse(Day == 'D0', 'T0', 
                        ifelse(Day == 'D7', 'T1',
                               ifelse(Day == 'D14', 'T2', 'T3')))) %>%
  modify_at('t_cat', factor, levels = c('T0', 'T1', 'T2', 'T3')) 

#Setup Contrasts:
#Bog vs Fen
T0_BvF <- c(0,1,0,0,0,0,0,0)
T1_BvF <- c(0,1,0,0,0,1,0,0)
T2_BvF <- c(0,1,0,0,0,0,1,0)
T3_BvF <- c(0,1,0,0,0,0,0,1)

#Bog - Time 
B_T0vT1 <- c(0,0,1,0,0,1,0,0)
B_T0vT2 <- c(0,0,0,1,0,0,1,0)
B_T0vT3 <- c(0,0,0,0,1,0,0,1)
B_T1vT2 <- c(0,0,-1,1,0,-1,1,0)
B_T1vT3 <- c(0,0,-1,0,1,-1,0,1)
B_T2vT3 <- c(0,0,0,-1,1,0,-1,1)

#Fen - Time 
F_T0vT1 <- c(0,0,1,0,0,0,0,0)
F_T0vT2 <- c(0,0,0,1,0,0,0,0)
F_T0vT3 <- c(0,0,0,0,1,0,0,0)
F_T1vT2 <- c(0,0,-1,1,0,0,0,0)
F_T1vT3 <- c(0,0,-1,0,1,0,0,0)
F_T2vT3 <- c(0,0,0,-1,1,0,0,0)

lm_rp_contra <- lm_rp %>%
  mutate(p_T0_BvF = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(T0_BvF, 1)))$test$pvalues[1])) %>%
  mutate(p_T1_BvF = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(T1_BvF, 1)))$test$pvalues[1])) %>%
  mutate(p_T2_BvF = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(T2_BvF, 1)))$test$pvalues[1])) %>%
  mutate(p_T3_BvF = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(T3_BvF, 1)))$test$pvalues[1])) %>%
  mutate(p_B_T0vT1 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(B_T0vT1, 1)))$test$pvalues[1])) %>%
  mutate(p_B_T0vT2 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(B_T0vT2, 1)))$test$pvalues[1])) %>%
  mutate(p_B_T0vT3 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(B_T0vT3, 1)))$test$pvalues[1])) %>%
  mutate(p_B_T1vT2 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(B_T1vT2, 1)))$test$pvalues[1])) %>%
  mutate(p_B_T1vT3 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(B_T1vT3, 1)))$test$pvalues[1])) %>%
  mutate(p_B_T2vT3 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(B_T2vT3, 1)))$test$pvalues[1])) %>%
  mutate(p_F_T0vT1 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(F_T0vT1, 1)))$test$pvalues[1])) %>%
  mutate(p_F_T0vT2 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(F_T0vT2, 1)))$test$pvalues[1])) %>%
  mutate(p_F_T0vT3 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(F_T0vT3, 1)))$test$pvalues[1])) %>%
  mutate(p_F_T1vT2 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(F_T1vT2, 1)))$test$pvalues[1])) %>%
  mutate(p_F_T1vT3 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(F_T1vT3, 1)))$test$pvalues[1])) %>%
  mutate(p_F_T2vT3 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(F_T2vT3, 1)))$test$pvalues[1])) 

lm_rp_contra %>%
  column_to_rownames('ft') %>%
  dplyr::select(starts_with('p_')) %>%
  apply(., 2, function(x) sum(x < 0.05))

lm_rp_contra_adj <- lm_rp_contra %>%
  pivot_longer(starts_with('p_'), names_to = 'test') %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(value, method = 'BH')) %>%
  mutate(new_test = paste0(gsub('p_', '', test), '_adj')) %>%
  dplyr::select(ft, new_test, p_adj) %>%
  pivot_wider(names_from = 'new_test', values_from = 'p_adj')

lm_rp_contra_adj %>%
  column_to_rownames('ft') %>%
  dplyr::select(ends_with('_adj')) %>%
  apply(., 2, function(x) sum(x < 0.05))

calc_log2fc <- function(x, factor1, factor2, factor1_value, factor2_groups){
  if(length(factor2_groups) != 2){
    stop('Must Be Exactly 2 Groups')
  }
  val_mat <- x[x[factor1] == factor1_value,]
  group1 <- val_mat[val_mat[factor2] == factor2_groups[1], ][['area']]
  group2 <- val_mat[val_mat[factor2] == factor2_groups[2], ][['area']]
  l2fc <- log2(mean(group2)/mean(group1))
  return(l2fc)
}

rp_fold <- rp_raw %>%
  rownames_to_column('Filename') %>%
  left_join(meta_rp2) %>%
  filter(Status == 'Alive') %>%
  filter(Treatment == 'unamended') %>%
  filter(t_cat %in% c('T0', 'T1', 'T2', 'T3')) %>%
  pivot_longer(cols = starts_with('FT_'), names_to = 'ft', values_to = 'area') %>%
  group_by(ft) %>%
  nest() %>%
  mutate(T0_l2fc = map_dbl(data, function(x) calc_log2fc(x, factor1 = 't_cat', factor2 = 'Site', factor1_value = 'T0', factor2_groups = c('fen', 'bog')))) 

lm_data_hilic <- hilic_control %>%
  filter(Day %in% c('D0', 'D7', 'D14', 'D28', 'D35')) %>%
  mutate(t_cat = ifelse(Day == 'D0', 'T0', 
                        ifelse(Day == 'D7', 'T1',
                               ifelse(Day == 'D14', 'T2', 'T3')))) %>%
  modify_at('t_cat', factor, levels = c('T0', 'T1', 'T2', 'T3')) %>%
  pivot_longer(starts_with('FT_'), names_to = 'ft', values_to = 'area') %>%
  group_by(ft) %>%
  nest()

lm_hilic <- lm_data_hilic %>%
  mutate(mod = purrr::map(data, function(x) lm(area ~ Site + t_cat + Site*t_cat, data = x)))

lm_hilic_contra <- lm_hilic %>%
  mutate(p_T0_BvF = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(T0_BvF, 1)))$test$pvalues[1])) %>%
  mutate(p_T1_BvF = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(T1_BvF, 1)))$test$pvalues[1])) %>%
  mutate(p_T2_BvF = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(T2_BvF, 1)))$test$pvalues[1])) %>%
  mutate(p_T3_BvF = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(T3_BvF, 1)))$test$pvalues[1])) %>%
  mutate(p_B_T0vT1 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(B_T0vT1, 1)))$test$pvalues[1])) %>%
  mutate(p_B_T0vT2 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(B_T0vT2, 1)))$test$pvalues[1])) %>%
  mutate(p_B_T0vT3 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(B_T0vT3, 1)))$test$pvalues[1])) %>%
  mutate(p_B_T1vT2 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(B_T1vT2, 1)))$test$pvalues[1])) %>%
  mutate(p_B_T1vT3 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(B_T1vT3, 1)))$test$pvalues[1])) %>%
  mutate(p_B_T2vT3 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(B_T2vT3, 1)))$test$pvalues[1])) %>%
  mutate(p_F_T0vT1 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(F_T0vT1, 1)))$test$pvalues[1])) %>%
  mutate(p_F_T0vT2 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(F_T0vT2, 1)))$test$pvalues[1])) %>%
  mutate(p_F_T0vT3 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(F_T0vT3, 1)))$test$pvalues[1])) %>%
  mutate(p_F_T1vT2 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(F_T1vT2, 1)))$test$pvalues[1])) %>%
  mutate(p_F_T1vT3 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(F_T1vT3, 1)))$test$pvalues[1])) %>%
  mutate(p_F_T2vT3 = map_dbl(mod, function(x) summary(multcomp::glht(x, linfct = matrix(F_T2vT3, 1)))$test$pvalues[1])) 

lm_hilic_contra %>%
  column_to_rownames('ft') %>%
  dplyr::select(starts_with('p_')) %>%
  apply(., 2, function(x) sum(x < 0.05))

lm_hilic_contra_adj <- lm_hilic_contra %>%
  pivot_longer(starts_with('p_'), names_to = 'test') %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(value, method = 'BH')) %>%
  mutate(new_test = paste0(gsub('p_', '', test), '_adj')) %>%
  dplyr::select(ft, new_test, p_adj) %>%
  pivot_wider(names_from = 'new_test', values_from = 'p_adj')

lm_hilic_contra_adj %>%
  column_to_rownames('ft') %>%
  dplyr::select(ends_with('_adj')) %>%
  apply(., 2, function(x) sum(x < 0.05))

meta_hilic2 <- meta_hilic %>%
  mutate(t_cat = ifelse(Day == 'D0', 'T0', 
                        ifelse(Day == 'D7', 'T1',
                               ifelse(Day == 'D14', 'T2', 'T3')))) %>%
  modify_at('t_cat', factor, levels = c('T0', 'T1', 'T2', 'T3')) 

hilic_fold <- hilic_raw %>%
  rownames_to_column('Filename') %>%
  left_join(meta_hilic2) %>%
  filter(Status == 'Alive') %>%
  filter(Treatment == 'unamended') %>%
  filter(t_cat %in% c('T0', 'T1', 'T2', 'T3')) %>%
  pivot_longer(cols = starts_with('FT_'), names_to = 'ft', values_to = 'area') %>%
  group_by(ft) %>%
  nest() %>%
  mutate(T0_l2fc = map_dbl(data, function(x) calc_log2fc(x, factor1 = 't_cat', factor2 = 'Site', factor1_value = 'T0', factor2_groups = c('fen', 'bog')))) 

canopus_rp <- read_csv('Data/LCMS/RP/Canopus_RP_Consensus.csv') %>%
  modify_at('progenesis_id', ~paste0('FT_', .x))

canopus_hilic <- read_csv('Data/LCMS/HILIC/Canopus_HILIC_Consensus.csv') %>%
  modify_at('progenesis_id', ~paste0('FT_', .x))

link_rp <- read_csv('Data/LCMS/RP/CompoundLink_RP.csv')  %>%
  modify_at('Compound_ID', ~paste0('FT_', .x)) %>%
  mutate(Name = ifelse(abs(`Annot. DeltaMass [ppm]`) > 5, NA, Name)) %>%
  mutate(annotation_type = ifelse(!is.na(Name) & `mzVault Search` == 'Full match', 'L1',
                                  ifelse(!is.na(Name) & `mzVault Search` != 'Full match' & `mzCloud Search` != 'No results', 'L2', '')))
link_hilic <- read_csv('Data/LCMS/HILIC/CompoundLink_HILIC_Fix.csv') %>%
  modify_at('Compound_ID', ~paste0('FT_', .x)) %>%
  mutate(Name = ifelse(abs(`Annot. DeltaMass [ppm]`) > 5, NA, Name)) %>%
  mutate(annotation_type = ifelse(!is.na(Name) & `Annot. Source: mzVault Search` == 'Full match', 'L1',
                                  ifelse(!is.na(Name) & `Annot. Source: mzVault Search` != 'Full match' & `Annot. Source: mzCloud Search` != 'No results', 'L2', '')))

## chemRich Analysis -------------------------------------------------------

source('ChemRichFunctions.R')

lm_rp_cr <- lm_rp_contra_adj %>%
  dplyr::select(ft, T0_BvF_adj) %>%
  left_join(., dplyr::select(rp_fold, -data)) %>%
  left_join(., dplyr::select(canopus_rp, -contains('Probability'), -id), by = c('ft' = 'progenesis_id')) %>%
  dplyr::select(ft, T0_l2fc, T0_BvF_adj, `ClassyFire#level 5`) %>%
  dplyr::rename('compound_name' = ft, 'effect_size' = T0_l2fc, 'pvalue' = T0_BvF_adj, 'set' = `ClassyFire#level 5`)

lm_hilic_cr <- lm_hilic_contra_adj %>%
  dplyr::select(ft, T0_BvF_adj) %>%
  left_join(., dplyr::select(hilic_fold, -data)) %>%
  left_join(., dplyr::select(canopus_hilic, -contains('Probability'), -id), by = c('ft' = 'progenesis_id')) %>%
  dplyr::select(ft, T0_l2fc, T0_BvF_adj, `ClassyFire#level 5`) %>%
  dplyr::rename('compound_name' = ft, 'effect_size' = T0_l2fc, 'pvalue' = T0_BvF_adj, 'set' = `ClassyFire#level 5`)

cr_input <- rbind(lm_rp_cr, lm_hilic_cr) %>%
  filter(!set %in% c('Amino acids and derivatives', 'Peptides')) %>%
  chemRICHFix()

plot_set <- c('Monosaccharides', 'Glycosyl compounds', 'Aminosaccharides', 'Branched fatty acids', 
              'Sugar acids and derivatives', 'Coumaric acids and derivatives', '7-hydroxycoumarins', 
              'Medium-chain fatty acids', 'Oligosaccharides', 'Ketones', 'Disaccharides', 'Benzoic acid esters', 
              'Substituted imidazoles', 'Hydroquinones', 'Flavanones', 'Benzoic acids', 'Chromones',
              'Phthalic acid and derivatives', 'Alkyl phosphates')


plotChemRich_custom(cr_input, colors = c('#3333FF', 'yellow3', '#009900'), sets = plot_set, binary = T) +
  cowplot::theme_cowplot() +
  scale_x_continuous(labels = NULL) +
  xlab(NULL) +
  ylab('—log(p value)') +
  theme(legend.position = c(0.75, 0.9),
        text = element_text(size = 20),
        axis.text.y = element_text(size = 18))


# Figure 1: ---------------------------------------------------------------

nosc_plot <- ggplot(nosc_data, aes(x = unique, y = NOSC, fill = unique)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.05, fill = 'grey') +
  scale_fill_manual(values = pal_site_cap, name = NULL) +
  cowplot::theme_cowplot() +
  xlab(NULL) +
  theme(legend.position = c(0.8, 0.9))

dif_plot <- ggplot(comp_plot, aes(y = unique, x = p_class, fill = Class_new)) +
  geom_col(color = 'black', linewidth = 1) +
  scale_fill_manual(values = comp_colors, name = 'Class') +
  cowplot::theme_cowplot() +
  xlab('Relative Abundance') +
  scale_x_continuous(labels = scales::label_percent()) +
  ylab(NULL) +
  theme(legend.position = 'bottom')

rp_pca_data <- mixOmics::plotIndiv(pca_rp_control)$df %>%
  cbind(meta_rp_control) %>%
  mutate(sxd = paste0(Site, '_', Day)) %>%
  modify_at('sxd', factor, levels = c('fen_D0', 'bog_D0', 'fen_D7', 'bog_D7', 'fen_D14', 'bog_D14', 'fen_D35', 'bog_D28'))

fen_day <- c('#01267A', '#024AF2', '#018ADF', '#5DC0FE')
names(fen_day) <- c('fen_D0', 'fen_D7', 'fen_D14', 'fen_D35')
bog_day <- c('#276501', '#0F8D01', '#3BDE00','#6FFF47')
names(bog_day) <- c('bog_D0', 'bog_D7', 'bog_D14', 'bog_D28')

pca_plot <- ggplot(rp_pca_data, aes(x, y, color = sxd)) +
  geom_polygon(data = rp_elip, aes(x = Col1, y = Col2), inherit.aes = F, color = pal_site[2], fill = NA, size = 1) +
  geom_polygon(data = rp_elip, aes(x = Col3, y = Col4), inherit.aes = F, color = pal_site[1], fill = NA, size = 1) +
  geom_point(size = 6)  +
  scale_color_manual(values = c(fen_day, bog_day), 
                     labels = c('Fen: D0', 'Bog: D0', 'Fen: D7', 'Bog: D7', 'Fen: D14', 'Bog: D14', 'Fen: D35', 'Bog: D28'),
                     name = 'Site Day') +
  xlab('PC1: 59%') +
  ylab('PC2: 12%') +
  cowplot::theme_cowplot() 
pca_plot

#Note: the order of X-axis is random so figure will turn out slightly different.
cr_plot <- plotChemRich_custom(cr_input, colors = c('#3333FF','yellow3','#009900'), sets = plot_set, binary = T) +
  cowplot::theme_cowplot()  +
  scale_x_continuous(labels = NULL) +
  scale_y_continuous(limits = c(0, 20)) +
  xlab(NULL) +
  ylab('—log(p value)') +
  theme(legend.position = c(0.75, 0.9))

classification <- tribble(
  ~Class, ~OC_low, ~OC_high, ~HC_low, ~HC_high,
  'Lipid', 0, 0.3, 1.5, 2.5,
  NULL, 0, 0.125, 0.8, 1.5,
  'Cond. HC', 0, 0.95, 0.2, 0.8,
  NULL , 0.3, 0.55, 1.5, 2.3,
  'Amino sugar', 0.55, 0.7, 1.5, 2.2,
  'Carbohydrate', 0.7, 1.5, 1.5, 2.5,
  'Lignin', 0.125, 0.65, 0.8, 1.5,
  'Tannin', 0.65, 1.1, 0.8, 1.5, 
) %>% 
  mutate(label_x = (OC_low + OC_high) / 2,
         label_y = HC_high - 0.1)

## Compound class rectangles (for plotting of Van Krevelen diagrams) 

class_rect <-  geom_rect(data = classification,
                         aes(xmin = OC_low,
                             xmax = OC_high,
                             ymin = HC_low,
                             ymax = HC_high),
                         color = 'black',
                         fill = NA,
                         linewidth = 1,
                         inherit.aes = FALSE, 
                         linetype = 'dashed')
rect_label <- geom_label(data = classification,
                         aes(x = label_x,
                             y = label_y,
                             label = Class),
                         inherit.aes = FALSE,
                         size = 3)

cont_bog <- ggplot(bog_vk, aes(x = OC, y = HC, color = unique)) + 
  stat_density_2d(geom = 'polygon', contour = TRUE,
                  aes(fill = after_stat(level)), colour = 'black') +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  labs(x = 'O:C',
       y = 'H:C') +
  cowplot::theme_cowplot() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2.55)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1.55)) +
  theme(legend.position = 'none') +
  class_rect  +
  rect_label 

cont_fen <- ggplot(fen_vk, aes(x = OC, y = HC, color = unique)) + 
  stat_density_2d(geom = 'polygon', contour = TRUE,
                  aes(fill = after_stat(level)), colour = 'black') +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(x = 'O:C',
       y = 'H:C') +
  cowplot::theme_cowplot() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2.55)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1.55)) +
  theme(legend.position = 'none') +
  theme(plot.margin = margin(unit(c(5.5, 12.5, 5.5, 0), 'pt'))) +
  class_rect  +
  rect_label 

cont_plot <- cowplot::plot_grid(cont_bog, cont_fen, nrow = 1)

p1 <- cowplot::plot_grid(dif_plot, cr_plot, ncol = 1, rel_heights = c(5,10), labels = c('A', 'D'))
p2 <- cowplot::plot_grid(nosc_plot, pca_plot, align = 'h', axis = 'tb', rel_widths = c(5,8), labels = c('B', 'C'))

p3 <- cowplot::plot_grid(p2, cont_plot, ncol = 1, labels = c('', 'E'), rel_heights = c(5,4))
cowplot::plot_grid(p1, p3, nrow = 1, rel_widths = c(5,5))  



# Metatranscriptomics Data Prep: ------------------------------------------

library(edgeR)
## Data Load In ------------------------------------------------------------

### Bog ---------------------------------------------------------------------

#Read in the raw count data for bog
bog_counts_raw <- read_delim('Data/MetaT/htseq_2304MAGs_UCD_bog_100M_97_REVSTRANDED_no0s.txt', delim = '\t', 
                             col_names = c("gene","STM_0716_S_M_E014","STM_0716_S_M_E015",
                                           "STM_0716_S_M_E016","STM_0716_S_M_E086","STM_0716_S_M_E087",
                                           "STM_0716_S_M_E088","STM_0716_S_M_E097","STM_0716_S_M_E098",
                                           "STM_0716_S_M_E099")) %>%
  #Remove ambiguous and unmapped genes:
  filter(!gene %in% c('__no_feature', '__ambiguous')) %>%
  column_to_rownames('gene') 

saveRDS(bog_counts_raw, 'Data/MetaT/bog_counts_raw.rds')

### Fen ---------------------------------------------------------------------

#Read in the raw count data for the fen:
fen_counts_raw <- read_delim('Data/metaT/htseq_2304MAGs_100M_97_REVSTRANDED_no0s.txt', delim = '\t', 
                             col_names = c("gene","STM_0716_E_M_E002","STM_0716_E_M_E003",
                                           "STM_0716_E_M_E004","STM_0716_E_M_E025","STM_0716_E_M_E027",
                                           "STM_0716_E_M_E029","STM_0716_E_M_E030","STM_0716_E_M_E031",
                                           "STM_0716_E_M_E033","STM_0716_E_M_E034","STM_0716_E_M_E035",
                                           "STM_0716_E_M_E050","STM_0716_E_M_E051","STM_0716_E_M_E052",
                                           "STM_0716_E_M_E054","STM_0716_E_M_E055","STM_0716_E_M_E056",
                                           "STM_0716_E_M_E058","STM_0716_E_M_E059","STM_0716_E_M_E060",
                                           "STM_0716_E_M_E062","STM_0716_E_M_E063","STM_0716_E_M_E064",
                                           "STM_0716_E_M_E066","STM_0716_E_M_E067","STM_0716_E_M_E068",
                                           "STM_0716_E_M_E070","STM_0716_E_M_E071","STM_0716_E_M_E072",
                                           "STM_0716_E_M_E121","STM_0716_E_M_E122","STM_0716_E_M_E123",
                                           "STM_0716_E_M_E125","STM_0716_E_M_E126","STM_0716_E_M_E127",
                                           "STM_0716_E_M_E129","STM_0716_E_M_E130","STM_0716_E_M_E131")) %>%
  #Remove ambiguous and unmapped genes:
  filter(!gene %in% c('__no_feature', '__ambiguous')) %>%
  column_to_rownames('gene')

### Metadata ---------------------------------------------------------------

metadata <- readRDS('Data/MetaT/MetaDataFinal.RDS') %>%
  mutate(tcat = ifelse(day == 'day0', 'ti', 'tf')) %>%
  modify_at('tcat', factor, levels = c('ti', 'tf')) 

metadata_fen <- metadata %>%
  filter(site == 'fen') 

metadata_bog <- metadata %>%
  filter(site == 'bog')

pal_site <- c('#009900', '#3333FF')
names(pal_site) <- c('bog', 'fen') 

## Filtering only on SAMPLE type -------------------------------------------

### Bog --------------------------------------------------------------------

bog_filter_group <- metadata_bog$day %in% c('day0', 'day28')
names(bog_filter_group) <- metadata_bog$sample

bog_small <- bog_counts_raw[,bog_filter_group]
bog_small_nz <- bog_small[apply(bog_small, 1, sum) != 0,] #370,113

### Fen --------------------------------------------------------------------

names_filter_fen <- metadata_fen$treatment == 'unamended' & metadata_fen$day %in% c('day0', 'day35')
fen_ua <- metadata_fen$sample[names_filter_fen]
fen_filter_group <- colnames(fen_counts_raw) %in% fen_ua

fen_small <- fen_counts_raw[,fen_filter_group]
fen_small_nz <- fen_small[apply(fen_small, 1, sum) != 0,] #1,249,206

## Calculate geTMMs --------------------------------------------------------

gene_lengths <- read_delim('Data/MetaT/gene_lengths.txt', col_names = c('gene', 'length'))

#Load in the MAGs:
gene_tax <- readRDS('Data/MetaT/MAG_Tax.RDS')

mag_tax <- gene_tax %>%
  dplyr::select(-gene, -GTDB, -MAGarb) %>%
  distinct()

#Load in annotations
dram <- read_delim('Data/MetaT/reactorEMERGE_annotations.tsv')  %>%
  dplyr::rename('gene' = `...1`) 

### Bog --------------------------------------------------------------------

bog_small_gl <- bog_small_nz %>%
  rownames_to_column('gene') %>%
  left_join(gene_lengths) %>%
  column_to_rownames('gene')

bog_small_norm <- sweep(bog_small_gl, 1, bog_small_gl$length, '/') %>%
  dplyr::select(-length)

bog_group_small <- rep('A', ncol(bog_small_norm))

dge_bog_small <- edgeR::DGEList(counts = bog_small_norm, group = bog_group_small)

dge_bog_norm_small <- edgeR::calcNormFactors(dge_bog_small)

bog_getmm <- as.data.frame(edgeR::cpm(dge_bog_norm_small))

### Fen --------------------------------------------------------------------

fen_small_gl <- fen_small_nz %>%
  rownames_to_column('gene') %>%
  left_join(gene_lengths) %>%
  column_to_rownames('gene')

fen_small_norm <- sweep(fen_small_gl, 1, fen_small_gl$length, '/') %>%
  dplyr::select(-length)

fen_group_small <- rep('A', ncol(fen_small_norm))

dge_fen_small <- edgeR::DGEList(counts = fen_small_norm, group = fen_group_small)

dge_fen_norm_small <- edgeR::calcNormFactors(dge_fen_small)

fen_getmm <- as.data.frame(edgeR::cpm(dge_fen_norm_small))

## Combine Bog and Fen -----------------------------------------------------

bf_getmm <- full_join(rownames_to_column(bog_getmm, 'gene'), rownames_to_column(fen_getmm, 'gene')) %>%
  pivot_longer(cols = starts_with('STM_'), names_to = 'sample', values_to = 'getmm') %>%
  replace_na(list('getmm' = 0)) %>%
  pivot_wider(names_from = 'sample', values_from = 'getmm')

bf_counts <- full_join(rownames_to_column(bog_small_nz, 'gene'), rownames_to_column(fen_small_nz, 'gene')) %>%
  pivot_longer(cols = starts_with('STM_'), names_to = 'sample', values_to = 'getmm') %>%
  replace_na(list('getmm' = 0)) %>%
  pivot_wider(names_from = 'sample', values_from = 'getmm')

## Time Course Fen ---------------------------------------------------------

tc_fen_filter <- metadata_fen$treatment == 'unamended'
fen_tc_ua <- metadata_fen$sample[tc_fen_filter]
fen_tc_group <- colnames(fen_counts_raw) %in% fen_tc_ua

fen_tc <- fen_counts_raw[,fen_tc_group]
fen_tc_nz <- fen_tc[apply(fen_tc, 1, sum) != 0,] #1,543,777

fen_tc_gl <- fen_tc_nz %>%
  rownames_to_column('gene') %>%
  left_join(gene_lengths) %>%
  column_to_rownames('gene')

fen_tc_norm <- sweep(fen_tc_gl, 1, fen_tc_gl$length, '/') %>%
  dplyr::select(-length)

fen_group_tc <- rep('A', ncol(fen_tc_norm))

dge_fen_tc <- edgeR::DGEList(counts = fen_tc_norm, group = fen_group_tc)

dge_fen_norm_tc <- edgeR::calcNormFactors(dge_fen_tc)

fen_getmm_tc <- as.data.frame(edgeR::cpm(dge_fen_norm_tc))

## Time Course Bog ---------------------------------------------------------

bog_tc_gl <- bog_counts_raw %>%
  rownames_to_column('gene') %>%
  left_join(gene_lengths) %>%
  column_to_rownames('gene')

bog_tc_norm <- sweep(bog_tc_gl, 1, bog_tc_gl$length, '/') %>%
  dplyr::select(-length)

bog_group_tc <- rep('A', ncol(bog_tc_norm))

dge_bog_tc <- edgeR::DGEList(counts = bog_tc_norm, group = bog_group_tc)

dge_bog_norm_tc <- edgeR::calcNormFactors(dge_bog_tc)

bog_getmm_tc <- as.data.frame(edgeR::cpm(dge_bog_norm_tc))

## Only to Bog Day 28 ------------------------------------------------------

bf_tc_getmm <- full_join(rownames_to_column(bog_getmm, 'gene'), rownames_to_column(fen_getmm_tc, 'gene')) %>%
  pivot_longer(cols = starts_with('STM_'), names_to = 'sample', values_to = 'getmm') %>%
  replace_na(list('getmm' = 0)) %>%
  pivot_wider(names_from = 'sample', values_from = 'getmm')

bf_tc_counts <- full_join(rownames_to_column(bog_small_nz, 'gene'), rownames_to_column(fen_tc_nz, 'gene')) %>%
  pivot_longer(cols = starts_with('STM_'), names_to = 'sample', values_to = 'getmm') %>%
  replace_na(list('getmm' = 0)) %>%
  pivot_wider(names_from = 'sample', values_from = 'getmm')

#Final geTMM data can be downloaded from Supplementary Data 5


# Supplemental Figure 1: Community Composition: ---------------------------
active_getmm <- bf_tc_getmm %>%
  left_join(gene_tax) %>%
  filter(!is.na(MAG)) %>%
  group_by(MAG)  %>%
  pivot_longer(cols = starts_with('STM_'), names_to = 'sample', values_to = 'getmm') %>%
  left_join(metadata) %>%
  group_by(site, phylum) %>%
  summarise(total_exp = sum(getmm)) %>%
  group_by(site) %>%
  mutate(rel_prop = total_exp/sum(total_exp)) %>%
  mutate(phylum_alt = ifelse(rel_prop < 0.01, 'Other', phylum))


active_n <- bf_tc_getmm %>%
  left_join(gene_tax) %>%
  filter(!is.na(MAG)) %>%
  pivot_longer(cols = starts_with('STM_'), names_to = 'sample', values_to = 'getmm') %>%
  left_join(metadata) %>%
  group_by(MAG, site) %>%
  summarise(total_exp = sum(getmm)) 

ind_mags <- active_n %>% 
  filter(total_exp != 0) %>%
  group_by(site) %>%
  nest() 

ggvenn::ggvenn(data = list(bog = ind_mags$data[[1]]$MAG, fen = ind_mags$data[[2]]$MAG))

active_mags <- active_n %>%
  ungroup() %>%
  filter(total_exp != 0) %>%
  left_join(mag_tax) %>%
  filter(!is.na(MAG)) %>%
  group_by(site, phylum) %>%
  summarise(n_mag = n()) %>%
  group_by(site) %>%
  mutate(rel_prop = n_mag/sum(n_mag)) %>%
  mutate(phylum_alt = ifelse(rel_prop < 0.01, 'Other', phylum)) 

active_u <- unique(c(active_mags$phylum_alt, active_n$phylum_alt))
active_u_s <- c(active_u[-3], active_u[3])


pal_long <- as.vector(paletteer::paletteer_d("ggthemes::Tableau_20"))
active_pal <- c(pal_long[-c(13,14)], pal_long[c(13,14)])
names(active_pal) <- active_u_s

active_getmm <- active_getmm %>%
  modify_at('phylum_alt', factor, levels = active_u_s)

rel_tran <- ggplot(active_getmm, aes(x = site, y = rel_prop, fill = phylum_alt)) +
  geom_col() +
  scale_fill_manual(values = active_pal,
                    name = 'Phylum') +
  cowplot::theme_cowplot() +
  scale_y_continuous(expand = c(0,0),
                     labels = scales::label_percent()) +
  scale_x_discrete(labels = c('Bog', 'Fen')) +
  ylab('Proportion of Expressed Transcripts') +
  xlab(NULL)

active_mags <- active_mags %>%
  modify_at('phylum_alt', factor, levels = active_u_s)

rel_mag <- ggplot(active_mags, aes(x = site, y = rel_prop, fill = phylum_alt)) +
  geom_col() +
  scale_fill_manual(values = active_pal,
                    name = 'Phylum') +
  cowplot::theme_cowplot() +
  scale_y_continuous(expand = c(0,0),
                     labels = scales::label_percent()) +
  scale_x_discrete(labels = c('Bog', 'Fen')) +
  ylab('Proportion of Active MAGs') +
  xlab(NULL)

cowplot::plot_grid(rel_mag, rel_tran, labels = c('A', 'B'))


# Within MAG Balancing: ---------------------------------------------------
mag_expression <- bf_getmm %>%
  left_join(dram) %>%
  left_join(gene_tax) %>%
  group_by(MAG) %>%
  nest() %>%
  filter(!is.na(MAG)) %>%
  mutate(expressed_genes = purrr::map(data, function(x){
    x %>%
      pivot_longer(starts_with('STM_'), names_to = 'sample', values_to = 'getmm') %>%
      group_by(gene) %>%
      filter(any(getmm > 0)) %>%
      dplyr::select(gene, sample, getmm, contains('_id'), contains('hit')) 
  })) %>%
  mutate(per_mag_expression = purrr::map(expressed_genes, function(x){
    x %>%
      dplyr::rename('label' = ko_id) %>%
      group_by(sample) %>%
      nest() %>%
      mutate(ED_exp = map_dbl(data, function(x) recurse_calc('K00036+(K01057,K07404)+K01690+K01625', data = x))) %>%
      mutate(EMP_exp = map_dbl(data, function(x) recurse_calc('(K01810,K06859,K13810,K15916)+(K00850,K16370,K21071,K00918)+(K01623,K01624,K11645,K16305,K16306)+K01803', data = x))) %>%
      mutate(AC_exp = map_dbl(data, function(x) recurse_calc('(K13788,K04020)+(K01512)', data = x))) %>%
      mutate(but_exp = map_dbl(data, function(x) recurse_calc('K00626+(K00074,K00022,K07516,K01825,K01782,K07514,K00023)+(K17865,K01692,K07515,K07511)+(K17829,K00209)+(K00634,K00929,K01034+K01035,K19709)', data = x))) %>%
      mutate(lac_exp = map_dbl(data, function(x) recurse_calc('K00101,K00016,K00102,K03777', data = x))) %>%
      mutate(prop_exp = map_dbl(data, function(x) recurse_calc('(K01958,K01959,K01960)+K00024+K01676+(K00239,K00244,K18209,K18556)+(K01899,K01902)+K01847+K05606+(K11264,K03416)', data = x))) %>%
      mutate(etOH_new_exp = map_dbl(data, function(x) recurse_calc('(K00132,K04072,K04073,K18366,K04021)+(K13953,K13954,K00001,K00121,K04072,K00002,K12957,K13979,K00114,K04022)', data = x))) %>%
      pivot_longer(cols = ends_with('_exp'), names_to = 'path', values_to = 'expression')
  })) 

mag_expression_unnested <- mag_expression %>%
  dplyr::select(MAG, per_mag_expression) %>%
  unnest('per_mag_expression')

mag_exp <- mag_expression_unnested %>%
  dplyr::select(-data) %>%
  pivot_wider(names_from = 'path', values_from = 'expression')

mag <- mag_exp %>%
  left_join(metadata) %>%
  pivot_longer(cols = ends_with('_exp'), names_to = 'path', values_to = 'expression') %>%
  #filter(expression != 0) %>%
  group_by(path, sample) %>%
  summarise(n_mag = length(unique(MAG))) %>%
  left_join(metadata)

mag_balance <- mag_exp %>%
  mutate(ox_full = ED_exp + EMP_exp) %>%
  mutate(ferm_full = etOH_new_exp + but_exp + lac_exp + prop_exp + AC_exp) %>%
  left_join(metadata) %>%
  mutate(sxd = paste0(site, '_', day)) %>%
  modify_at('sxd', factor, levels = c('bog_day0', 'fen_day0', 'fen_day7', 'fen_day14', 'fen_day21',
                                      'bog_day28', 'fen_day35')) %>%
  filter(ox_full != 0 & ferm_full != 0)

ggplot(mag_balance, aes(x = ox_full, y = ferm_full, color = sxd)) +
  geom_point() +
  facet_wrap(~sxd, scales = 'free') 

mag_balance %>%
  group_by(sxd) %>%
  nest() %>%
  mutate(r = purrr::map_dbl(data, function(x) cor(x$ox_full, x$ferm_full))) %>%
  mutate(rho = purrr::map_dbl(data, function(x) cor(x$ox_full, x$ferm_full, method = 'sp'))) 


# Metatranscriptomics Pathway Mapping: ------------------------------------

## Load in the manually curated genes: -----
pmo <- readxl::read_excel('Data/MetaT/yanni_gene_curation.xlsx')

pmo_small <- pmo %>%
  filter(gene %in% bf_tc_getmm$gene) %>%
  filter(annotation %in% 'PMO')

pmo_bad <- pmo %>%
  filter(gene %in% bf_tc_getmm$gene) %>%
  filter(str_detect(annotation, 'exclude'))

amo_small <- pmo %>%
  filter(annotation %in% 'AMO')

nar <- readxl::read_excel('Data/MetaT/yanni_gene_curation.xlsx', sheet = 2)

nar_small <- nar %>%
  filter(str_detect(recommendation, 'reduction'))

nxr_small <- nar %>%
  filter(str_detect(recommendation, 'oxidation'))

## Prep data for timecourse analysis: -----
#Load in getmm normalized data
bf_dram_tc <- bf_tc_getmm %>%
  #Join on annotations and taxonomy
  left_join(dram) %>%
  left_join(gene_tax) %>%
  #Trim data to only important information
  dplyr::select(gene, starts_with('STM_'), ko_id) %>%
  mutate(ko_id = strsplit(as.character(ko_id), ',')) %>%
  unnest(ko_id) %>%
  #Nest for pathway mapping
  mutate(ko_custom = ifelse(gene %in% amo_small$gene & ko_id == 'K10944', 'K10944_amo',
                            ifelse(gene %in% amo_small$gene & ko_id == 'K10945', 'K10945_amo', 
                                   ifelse(gene %in% pmo_small$gene & ko_id == 'K10944', 'K10944_pmo',
                                          ifelse(gene %in% pmo_small$gene & ko_id == 'K10945', 'K10945_pmo',
                                                 ifelse(gene %in% pmo_small$gene & ko_id == 'K10946', 'K10946_pmo',
                                                        ifelse(gene %in% nar_small$gene & ko_id == 'K00370', 'K00370_nar',
                                                               ifelse(gene %in% nxr_small$gene & ko_id == 'K00370', 'K00370_nxr', ko_id)))))))) %>%
  dplyr::select(starts_with('STM_'), ko_custom) %>%
  pivot_longer(starts_with('STM_'), names_to = 'sample', values_to = 'getmm') %>%
  rename('label' = ko_custom) %>%
  group_by(sample) %>%
  nest() 

## Pathway Mapping: -----
#The recurse_calc() function can't be shared however it follows the rules explained by KEGG: https://www.genome.jp/kegg/module.html
#Briefly, KOs inside parentheses are worked from the inside out
#A + sign designates AND and values for those KOs are averaged
#A , sign designates OR and values for those KOs are summed
#K1,K2+K3 is interpretted as K1 OR (K1 AND K2)

#Nitrogen Metabolism
nitrogen_metab <- bf_dram_tc %>%
  #Reductive Processes:
  mutate(DNR_exp = map_dbl(data, function(x) recurse_calc('(K00370_nar+K00371_nar+K00374_nar)+(K00362+K00363,K03385+K15876)', data = x))) %>%
  mutate(ANR_exp = map_dbl(data, function(x) recurse_calc('(K00367,K10534,K00372,K00360)+(K00366,K17877,K26139+K26138,K00361)', data = x))) %>%
  mutate(DeNit_exp = map_dbl(data, function(x) recurse_calc('(K00370_nar+K00371_nar+K00374_nar)+(K04561+K02305)+K00376', data = x))) %>%
  #Oxidative Processes:
  mutate(Nit_exp = map_dbl(data, function(x) recurse_calc('(K10944_amo+K10945_amo+K10946_amo)+K10535+(K00370_nxr+K00371_nxr)', data = x))) %>%
  mutate(Annamox_exp = map_dbl(data, function(x) recurse_calc('(K20932+K20933+K20934)+K20935', data = x))) %>%
  dplyr::select(-data) %>%
  pivot_longer(cols = ends_with('_exp'), names_to = 'path', values_to = 'exp') %>%
  left_join(metadata) %>%
  mutate(sxd = paste0(site, '_', tcat)) %>%
  modify_at('sxd', factor, levels = c('bog_ti', 'bog_tf', 'fen_ti', 'fen_tf')) 


#Methane oxidation:
methane_ox_metab <- bf_dram_tc %>%
  mutate(methox_exp = map_dbl(data, function(x) recurse_calc('((K10944_pmo+K10945_pmo+K10946_pmo),(K16157+K16158+K16159+K16160+K16161+K16162))+((K14028,K14029),K23995)', data = x))) %>%
  dplyr::select(-data) %>%
  pivot_longer(cols = ends_with('_exp'), names_to = 'path', values_to = 'exp') %>%
  left_join(metadata) %>%
  mutate(sxd = paste0(site, '_', tcat)) %>%
  modify_at('sxd', factor, levels = c('bog_ti', 'bog_tf', 'fen_ti', 'fen_tf')) 

#Sulfur metabolism
sulfur_metab <- bf_dram_tc %>%
  #Reduction
  mutate(ASR_exp = map_dbl(data, function(x) recurse_calc('(((K13811,K00958+K00860,K00955+K00957,K00956+K00957+K00860)+K00390),(K13811+K05907))+(K00380+K00381,K00392)', data = x))) %>%
  mutate(DSR_exp = map_dbl(data, function(x) recurse_calc('K00958+(K00394+K00395)+(K11180+K11181,K27196)+(K27187+K27188+K27189+K27190+K27191)', data = x))) %>%
  #Oxidation
  mutate(SOX_exp = map_dbl(data, function(x) recurse_calc('K17222+K17223+K17224,K17225,K22622+K17226+K17227', data = x))) %>%
  dplyr::select(-data) %>%
  pivot_longer(cols = ends_with('_exp'), names_to = 'path', values_to = 'exp') %>%
  left_join(metadata) %>%
  mutate(sxd = paste0(site, '_', tcat)) %>%
  modify_at('sxd', factor, levels = c('bog_ti', 'bog_tf', 'fen_ti', 'fen_tf')) 

#Aromatic compound metabolism:
aromatics <- bf_tc_getmm %>%
  left_join(dram) %>%
  dplyr::select(gene, starts_with('STM_'), ko_id, camper_id) %>%
  mutate(ko_id = strsplit(as.character(ko_id), ',')) %>%
  unnest(ko_id) %>%
  pivot_longer(starts_with('STM_'), names_to = 'sample', values_to = 'getmm') %>%
  rename('label' = ko_id) %>%
  group_by(sample) %>%
  nest()  %>%
  #Benzoyl-CoA Degrdation Pathway
  mutate(bzcoa_exp = map_dbl(data, function(x) recurse_calc('(K04112+K04113+K04114+K04115,K19515+K19516)+K07537+K07538+K07539', data = x))) %>%
  #Switch IDs to CAMPER
  nplyr::nest_rename(data, 'ko_id' = label) %>%
  nplyr::nest_rename(data, 'label' = camper_id) %>%
  #Phloroglucinol Degradation
  mutate(phloro_exp = map_dbl(data, function(x) recurse_calc('D00023+D00042+D00043+D00044', data = x)))  %>%
  #Caffeic Acid Respiration
  mutate(caf_exp = map_dbl(data, function(x) recurse_calc('D00018+D00017+(D00016+D00015+D00014)', data = x)))  %>%
  dplyr::select(-data) %>%
  pivot_longer(cols = ends_with('_exp'), names_to = 'path', values_to = 'exp') %>%
  left_join(metadata) %>%
  mutate(sxd = paste0(site, '_', tcat)) %>%
  modify_at('sxd', factor, levels = c('bog_ti', 'bog_tf', 'fen_ti', 'fen_tf')) 

#For methanogenesis since the CODH/ACS complex is mutually shared across acetogenic bacteria
#We first isolate only the Archaea and then go from there:
arc_methane <- bf_tc_getmm %>%
  left_join(dram) %>%
  left_join(gene_tax) %>%
  #Pull only Archaea
  filter(domain == 'Archaea') %>%
  dplyr::select(gene, starts_with('STM_'), ko_id) %>%
  mutate(ko_id = strsplit(as.character(ko_id), ',')) %>%
  unnest(ko_id) %>%
  dplyr::select(starts_with('STM_'), ko_id) %>%
  pivot_longer(starts_with('STM_'), names_to = 'sample', values_to = 'getmm') %>%
  rename('label' = ko_id) %>%
  group_by(sample) %>%
  nest() %>%
  #Acetoclastic Methanogenesis
  mutate(acetometh_exp = map_dbl(data, function(x) recurse_calc('((K00925+K00625),K01895)+(K00193+K00197+K00194)', data = x))) %>%
  #Hydrogenotrophic Methanogenesis
  mutate(hydrometh_exp = map_dbl(data, function(x) recurse_calc('(K00200+K00201+K00202+K00203,K11261+(K00205,K11260,K00204))+K00672+K01499+(K00319,K13942)+K00320', data = x))) %>%
  dplyr::select(-data) %>%
  pivot_longer(cols = ends_with('_exp'), names_to = 'path', values_to = 'exp') %>%
  left_join(metadata) %>%
  mutate(sxd = paste0(site, '_', tcat)) %>%
  modify_at('sxd', factor, levels = c('bog_ti', 'bog_tf', 'fen_ti', 'fen_tf'))

#Fermentation Metabolism:
ferm <- bf_dram_tc %>%
  mutate(etOH_exp = map_dbl(data, function(x) recurse_calc('(K00132,K04072,K04073,K18366,K04021)+(K13953,K13954,K00001,K00121,K04072,K00002,K12957,K13979,K00114,K04022)', data = x))) %>%
  mutate(AC_exp = map_dbl(data, function(x) recurse_calc('(K13788,K04020)+(K01512)', data = x))) %>%
  mutate(but_exp = map_dbl(data, function(x) recurse_calc('K00626+(K00074,K00022,K07516,K01825,K01782,K07514,K00023)+(K17865,K01692,K07515,K07511)+(K17829,K00209)+(K00634,K00929,K01034+K01035,K19709)', data = x))) %>%
  mutate(lac_exp = map_dbl(data, function(x) recurse_calc('K00101,K00016,K00102,K03777', data = x))) %>%
  mutate(prop_exp = map_dbl(data, function(x) recurse_calc('(K01958,K01959,K01960)+K00024+K01676+(K00239,K00244,K18209,K18556)+(K01899,K01902)+K01847+K05606+(K11264,K03416)', data = x))) %>%
  dplyr::select(-data) %>%
  pivot_longer(cols = ends_with('_exp'), names_to = 'path', values_to = 'exp') %>%
  left_join(metadata) %>%
  mutate(sxd = paste0(site, '_', tcat)) %>%
  modify_at('sxd', factor, levels = c('bog_ti', 'bog_tf', 'fen_ti', 'fen_tf'))

#Oxidative Metabolism
#Only upstream parts of pathways, not shared components
ox <- bf_dram_tc %>%
  mutate(ED_exp = map_dbl(data, function(x) recurse_calc('K00036+(K01057,K07404)+K01690+K01625', data = x))) %>%
  mutate(EMP_exp = map_dbl(data, function(x) recurse_calc('(K01810,K06859,K13810,K15916)+(K00850,K16370,K21071,K00918)+(K01623,K01624,K11645,K16305,K16306)+K01803', data = x))) %>%
  dplyr::select(-data) %>%
  pivot_longer(cols = ends_with('_exp'), names_to = 'path', values_to = 'exp') %>%
  left_join(metadata) %>%
  mutate(sxd = paste0(site, '_', tcat)) %>%
  modify_at('sxd', factor, levels = c('bog_ti', 'bog_tf', 'fen_ti', 'fen_tf'))

#Hemicellulase_genes
depoly <- bf_dram_tc %>%
  mutate(hemicell_xgca_exp = map_dbl(data, function(x) recurse_calc('K18651', data = x))) %>%
  dplyr::select(-data) %>%
  pivot_longer(cols = ends_with('_exp'), names_to = 'path', values_to = 'exp') %>%
  left_join(metadata) %>%
  mutate(sxd = paste0(site, '_', tcat)) %>%
  modify_at('sxd', factor, levels = c('bog_ti', 'bog_tf', 'fen_ti', 'fen_tf'))

#Now combine all together into a massive table
tea_tc <- arc_methane %>%
  rbind(methane_ox_metab) %>%
  rbind(nitrogen_metab) %>%
  rbind(sulfur_metab) %>%
  rbind(aromatics) %>%
  rbind(ferm) %>%
  rbind(ox) %>%
  rbind(depoly) %>%
  pivot_wider(names_from = 'path', values_from = 'exp') %>%
  #Summary_stats
  mutate(methane_all = acetometh_exp + hydrometh_exp) %>%
  mutate(NR_all = DNR_exp + ANR_exp + DeNit_exp) %>%
  mutate(SR_all = ASR_exp + DSR_exp)  %>%
  mutate(ferm_all = etOH_exp + AC_exp + but_exp + lac_exp + prop_exp) %>%
  mutate(methox_all = methox_exp) %>%
  mutate(car_all = caf_exp) %>%
  mutate(ox_all = EMP_exp + ED_exp) %>%
  mutate(hemi_all =  hemicell_xgca_exp) %>%
  mutate(sox_all =  SOX_exp) %>%
  mutate(nox_all =  Annamox_exp + Nit_exp) %>%
  modify_at('day', factor, levels = c('day0', 'day7', 'day14', 'day21', 'day28', 'day35')) %>%
  pivot_longer(cols = ends_with('_all'), names_to = 'tea')

# Figure 3: ---------------------------------------------------------------
## Fen ------
tea_tc_fen <- tea_tc %>%
  filter(site == 'fen')

tea_tc_fen_ranged <- tea_tc_fen %>%
  pivot_wider(names_from = 'tea') %>%
  dplyr::select(sample, ends_with('_all')) %>%
  column_to_rownames('sample') %>%
  auto_scale() %>%
  rownames_to_column('sample') %>%
  left_join(metadata) %>%
  modify_at('day', factor, levels = c('day0', 'day7', 'day14', 'day21', 'day28', 'day35')) %>%
  pivot_longer(cols = ends_with('_all'), names_to = 'tea')

ggplot(tea_tc_fen_ranged, aes(x = day, y = value, color = tea, group = tea, fill = tea)) +
  geom_smooth(linewidth = 2, se = F)

ggplot(tea_tc_fen, aes(x = day, y = value, color = tea, group = tea, fill = tea)) +
  geom_smooth(linewidth = 2, se = F)

tea_pal <- c('#EF476F', '#F8961E', '#2176FF', '#F9C80E', '#57C4E5', '#730071', '#08F7B7', '#073B4C', '#D1D646')
names(tea_pal) <- c('methane_all', 'ferm_all', 'NR_all', 'SR_all', 'nox_all', 'sox_all', 'ox_all', 'methox_all', 'hemi_all')

#A panel: Nitrate and Sulfate:
NS_actual_fen <- tea_tc_fen %>%
  filter(tea %in% c('NR_all', 'SR_all', 'sox_all', 'nox_all')) %>%
  modify_at('tea', factor, levels = c('NR_all', 'SR_all', 'nox_all', 'sox_all'))

ns_plot_fen_a <- ggplot(NS_actual_fen, aes(x = day, y = value, color = tea, group = tea, fill = tea)) +
  geom_smooth(linewidth = 2, se = T) +
  scale_color_manual(values = tea_pal, labels = c('Nitrate Reduction', 'Sulfate Reduction', 
                                                  'Nitrification', 'Sulfur Oxidation'),
                     name = 'Process') +
  scale_fill_manual(values = tea_pal, labels = c('Nitrate Reduction', 'Sulfate Reduction', 
                                                 'Nitrification', 'Sulfur Oxidation'),
                    name = 'Process') +
  scale_x_discrete(labels = NULL) +
  #scale_x_discrete(labels = c('Day 0', 'Day 7', 'Day 14', 'Day 21', 'Day 35')) + 
  cowplot::theme_minimal_grid()  +
  ylab('Expression (geTMM)') +
  xlab(NULL) +
  theme(legend.position = 'none')

NS_ranged_fen <- tea_tc_fen_ranged %>%
  filter(tea %in% c('NR_all', 'SR_all', 'sox_all', 'nox_all')) %>%
  modify_at('tea', factor, levels = c('NR_all', 'SR_all', 'nox_all', 'sox_all'))

ns_plot_fen_r <- ggplot(NS_ranged_fen, aes(x = day, y = value, color = tea, group = tea, fill = tea)) +
  geom_smooth(linewidth = 2, se = F) +
  scale_color_manual(values = tea_pal, labels = c('Nitrate Reduction', 'Sulfate Reduction', 
                                                  'Nitrification', 'Sulfur Oxidation'),
                     name = 'Process') +
  scale_fill_manual(values = tea_pal, labels = c('Nitrate Reduction', 'Sulfate Reduction', 
                                                 'Nitrification', 'Sulfur Oxidation'),
                    name = 'Process') +
  scale_x_discrete(labels = NULL) +
  #scale_x_discrete(labels = c('Day 0', 'Day 7', 'Day 14', 'Day 21', 'Day 35')) +
  cowplot::theme_minimal_grid()  +
  ylab('Autoscaled Expression') +
  xlab(NULL) 

fen_a_leg <- cowplot::get_legend(ns_plot_fen_r)
ns_plot_fen_r_nl <- ns_plot_fen_r +
  theme(legend.position = 'none')

#B panel: Carbon Processes:

C_actual_fen <- tea_tc_fen %>%
  filter(tea %in% c('ox_all', 'ferm_all', 'hemi_all', 'methane_all', 'methox_all')) %>%
  modify_at('tea', factor, levels = c('ox_all', 'ferm_all', 'hemi_all', 'methane_all', 'methox_all'))

c_plot_fen_a <- ggplot(C_actual_fen, aes(x = day, y = value, color = tea, group = tea, fill = tea)) +
  geom_smooth(linewidth = 2, se = T) +
  scale_color_manual(values = tea_pal) +
  scale_fill_manual(values = tea_pal) +
  scale_x_discrete(labels = c('Day 0', 'Day 7', 'Day 14', 'Day 21', 'Day 35')) +
  cowplot::theme_minimal_grid()  +
  ylab('Expression (geTMM)') +
  xlab(NULL) +
  theme(legend.position = 'none')

c_ranged_fen <- tea_tc_fen_ranged %>%
  filter(tea %in% c('ox_all', 'ferm_all', 'hemi_all', 'methane_all', 'methox_all')) %>%
  modify_at('tea', factor, levels = c('ox_all', 'ferm_all', 'hemi_all', 'methane_all', 'methox_all'))

c_plot_fen_r <- ggplot(c_ranged_fen, aes(x = day, y = value, color = tea, group = tea, fill = tea)) +
  geom_smooth(linewidth = 2, se = F) +
  scale_color_manual(values = tea_pal, labels = c('Oxidative \n Metabolism', 'Fermentation', 
                                                  'Hemicellulase', 'Methanogenesis', 'Methane Oxidation'),
                     name = 'Process') +
  scale_fill_manual(values = tea_pal, labels = c('Oxidative \n Metabolism', 'Fermentation', 
                                                 'Hemicellulase', 'Methanogenesis', 'Methane Oxidation'),
                    name = 'Process') +
  scale_x_discrete(labels = c('Day 0', 'Day 7', 'Day 14', 'Day 21', 'Day 35')) +
  cowplot::theme_minimal_grid()  +
  ylab('Autoscaled Expression') +
  xlab(NULL) 

b_fen_leg <- cowplot::get_legend(c_plot_fen_r)
c_plot_fen_r_nl <- c_plot_fen_r +
  theme(legend.position = 'none')


cowplot::plot_grid(ns_plot_fen_a, ns_plot_fen_r_nl, fen_a_leg, c_plot_fen_a, c_plot_fen_r_nl, b_fen_leg,
                   rel_widths = c(10, 10, 3),
                   align = 'hv', axis = 'tblr',
                   labels = c('A', '', '', 'B', '', ''))

## Bog ------
tea_tc_bog <- tea_tc %>%
  filter(site == 'bog')

tea_tc_bog_ranged <- tea_tc_bog %>%
  pivot_wider(names_from = 'tea') %>%
  dplyr::select(sample, ends_with('_all')) %>%
  column_to_rownames('sample') %>%
  auto_scale() %>%
  rownames_to_column('sample') %>%
  left_join(metadata) %>%
  modify_at('day', factor, levels = c('day0', 'day28')) %>%
  pivot_longer(cols = ends_with('_all'), names_to = 'tea')

ggplot(tea_tc_bog_ranged, aes(x = day, y = value, color = tea, group = tea, fill = tea)) +
  geom_smooth(linewidth = 2, se = F)

ggplot(tea_tc_bog, aes(x = day, y = value, color = tea, group = tea, fill = tea)) +
  geom_smooth(linewidth = 2, se = F)

#A panel: Nitrate and Sulfate:
NS_actual_bog <- tea_tc_bog %>%
  filter(tea %in% c('NR_all', 'SR_all', 'sox_all', 'nox_all')) %>%
  modify_at('tea', factor, levels = c('NR_all', 'SR_all', 'nox_all', 'sox_all'))

ns_plot_bog_a <- ggplot(NS_actual_bog, aes(x = day, y = value, color = tea, group = tea, fill = tea)) +
  geom_smooth(linewidth = 2, se = T) +
  scale_color_manual(values = tea_pal, labels = c('Nitrate Reduction', 'Sulfate Reduction', 
                                                  'Nitrification', 'Sulfur Oxidation'),
                     name = 'Process') +
  scale_fill_manual(values = tea_pal, labels = c('Nitrate Reduction', 'Sulfate Reduction', 
                                                 'Nitrification', 'Sulfur Oxidation'),
                    name = 'Process') +
  scale_x_discrete(labels = NULL) +
  #scale_x_discrete(labels = c('Day 0', 'Day 7', 'Day 14', 'Day 21', 'Day 35')) + 
  cowplot::theme_minimal_grid()  +
  ylab('Expression (geTMM)') +
  xlab(NULL) +
  theme(legend.position = 'none')

NS_ranged_bog <- tea_tc_bog_ranged %>%
  filter(tea %in% c('NR_all', 'SR_all', 'sox_all', 'nox_all')) %>%
  modify_at('tea', factor, levels = c('NR_all', 'SR_all', 'nox_all', 'sox_all'))

ns_plot_bog_r <- ggplot(NS_ranged_bog, aes(x = day, y = value, color = tea, group = tea, fill = tea)) +
  geom_smooth(linewidth = 2, se = F) +
  scale_color_manual(values = tea_pal, labels = c('Nitrate Reduction', 'Sulfate Reduction', 
                                                  'Nitrification', 'Sulfur Oxidation'),
                     name = 'Process') +
  scale_fill_manual(values = tea_pal, labels = c('Nitrate Reduction', 'Sulfate Reduction', 
                                                 'Nitrification', 'Sulfur Oxidation'),
                    name = 'Process') +
  scale_x_discrete(labels = NULL) +
  #scale_x_discrete(labels = c('Day 0', 'Day 7', 'Day 14', 'Day 21', 'Day 35')) +
  cowplot::theme_minimal_grid()  +
  ylab('Autoscaled Expression') +
  xlab(NULL) 

bog_a_leg <- cowplot::get_legend(ns_plot_bog_r)
ns_plot_bog_r_nl <- ns_plot_bog_r +
  theme(legend.position = 'none')

#B panel: Carbon Processes:

C_actual_bog <- tea_tc_bog %>%
  filter(tea %in% c('ox_all', 'ferm_all', 'hemi_all', 'methane_all', 'methox_all')) %>%
  modify_at('tea', factor, levels = c('ox_all', 'ferm_all', 'hemi_all', 'methane_all', 'methox_all'))

c_plot_bog_a <- ggplot(C_actual_bog, aes(x = day, y = value, color = tea, group = tea, fill = tea)) +
  geom_smooth(linewidth = 2, se = T) +
  scale_color_manual(values = tea_pal) +
  scale_fill_manual(values = tea_pal) +
  scale_x_discrete(labels = c('Day 0', 'Day 28')) +
  cowplot::theme_minimal_grid()  +
  ylab('Expression (geTMM)') +
  xlab(NULL) +
  theme(legend.position = 'none')

c_ranged_bog <- tea_tc_bog_ranged %>%
  filter(tea %in% c('ox_all', 'ferm_all', 'hemi_all', 'methane_all', 'methox_all')) %>%
  modify_at('tea', factor, levels = c('ox_all', 'ferm_all', 'hemi_all', 'methane_all', 'methox_all'))

c_plot_bog_r <- ggplot(c_ranged_bog, aes(x = day, y = value, color = tea, group = tea, fill = tea)) +
  geom_smooth(linewidth = 2, se = F) +
  scale_color_manual(values = tea_pal, labels = c('Oxidative \n Metabolism', 'Fermentation', 
                                                  'Hemicellulase', 'Methanogenesis', 'Methane Oxidation'),
                     name = 'Process') +
  scale_fill_manual(values = tea_pal, labels = c('Oxidative \n Metabolism', 'Fermentation', 
                                                 'Hemicellulase', 'Methanogenesis', 'Methane Oxidation'),
                    name = 'Process') +
  scale_x_discrete(labels = c('Day 0', 'Day 28')) +
  cowplot::theme_minimal_grid()  +
  ylab('Autoscaled Expression') +
  xlab(NULL) 

bog_b_leg <- cowplot::get_legend(c_plot_bog_r)
c_plot_bog_r_nl <- c_plot_bog_r +
  theme(legend.position = 'none')


cowplot::plot_grid(ns_plot_bog_a, ns_plot_bog_r_nl, bog_a_leg, c_plot_bog_a, c_plot_bog_r_nl, bog_b_leg,
                   rel_widths = c(10, 10, 3),
                   align = 'hv', axis = 'tblr',
                   labels = c('A', '', '', 'B', '', ''))


# Stress_genes -------------------------------------------------------------
box_data <- tea_tc %>%
  ungroup() %>%
  filter(day %in% c('day0', 'day28', 'day35')) %>%
  dplyr::select(sample, site, sxd, EMP_exp, ED_exp, etOH_exp, lac_exp, prop_exp, but_exp, AC_exp) %>%
  pivot_longer(cols = ends_with('_exp'), names_to = 'path') %>%
  modify_at('path', factor, levels = c('EMP_exp', 'ED_exp', 'etOH_exp', 'lac_exp', 'prop_exp', 'but_exp', 'AC_exp')) %>%
  distinct()

#Stress genes
bf_pfam <- bf_getmm %>%
  left_join(dram) %>%
  pivot_longer(cols = starts_with('STM_'), names_to = 'sample', values_to = 'getmm') %>%
  dplyr::select(sample, pfam_hits, getmm) %>%
  filter(!is.na(pfam_hits)) %>%
  group_by(sample) %>%
  nest()

expression <- bf_pfam %>%
  mutate(expression = purrr::map(data, function(x) data.frame(
    Thioredoxin_exp = mean(x$getmm[str_detect(x$pfam_hits, 'PF00085')]),
    GSHPx_exp = mean(x$getmm[str_detect(x$pfam_hits, 'PF00255')]),
    Catalase_exp = mean(x$getmm[str_detect(x$pfam_hits, 'PF00199')])
  ))) %>%
  dplyr::select(-data) %>%
  unnest('expression')


# Supplemental Figure 2 Code - Nitrate and Sulfate Transporters: ----------
#Nitrogen transport genes:
n_transport <- bf_dram_tc %>%
  mutate(amt_exp = map_dbl(data, function(x) recurse_calc('K03320', data = x))) %>%
  mutate(GDH_exp = map_dbl(data, function(x) recurse_calc('K15371', data = x))) %>%
  mutate(GS_exp = map_dbl(data, function(x) recurse_calc('K01915', data = x))) %>%
  mutate(GlnKB_exp = map_dbl(data, function(x) recurse_calc('K04752+K04751', data = x))) %>%
  mutate(DNR_u_exp = map_dbl(data, function(x) recurse_calc('(K00362+K00363,K03385+K15876)', data = x))) %>%
  mutate(ANR_exp = map_dbl(data, function(x) recurse_calc('(K00367,K10534,K00372,K00360)+(K00366,K17877,K26139+K26138,K00361)', data = x))) %>%
  mutate(DeNit_u_exp = map_dbl(data, function(x) recurse_calc('(K04561+K02305)+K00376', data = x))) %>%
  mutate(Nit_u_exp = map_dbl(data, function(x) recurse_calc('K10535', data = x))) %>%
  mutate(Annamox_u_exp = map_dbl(data, function(x) recurse_calc('(K20932+K20933+K20934)+K20935', data = x))) %>%
  dplyr::select(-data) %>%
  mutate(NR_exp = DNR_u_exp + ANR_exp + DeNit_u_exp) %>%
  mutate(nox_exp =  Annamox_u_exp + Nit_u_exp) %>%
  dplyr::select(-DNR_u_exp, -ANR_exp, -DeNit_u_exp, -Annamox_u_exp, -Nit_u_exp) %>%
  pivot_longer(cols = ends_with('_exp'), names_to = 'path', values_to = 'exp') %>%
  left_join(metadata) %>%
  modify_at('day', factor, levels = c('day0', 'day7', 'day14', 'day21', 'day28', 'day35')) 

n_transport_fen <- n_transport %>%
  filter(site == 'fen')

n_t_plot <- ggplot(n_transport_fen, aes(x = day, y = exp, color = path, group = path, fill = path)) +
  geom_smooth(linewidth = 2, se = F) +
  scale_x_discrete(labels = c('Day 0', 'Day 7', 'Day 14', 'Day 21', 'Day 35')) +
  cowplot::theme_minimal_grid()  +
  ylab('Expression (geTMM)') +
  xlab(NULL) 

#Sulfur transport genes
s_transport <- bf_dram_tc %>%
  #Sulfate Transport:
  mutate(cysP_exp = map_dbl(data, function(x) recurse_calc('K02048', data = x))) %>%
  #Sulfate Transport:
  mutate(sbp_exp = map_dbl(data, function(x) recurse_calc('K23163', data = x))) %>%
  #Sulfate Transport:
  mutate(SULP_exp = map_dbl(data, function(x) recurse_calc('K03321', data = x))) %>%
  #Sulfonate Transport:
  mutate(ssuA_exp = map_dbl(data, function(x) recurse_calc('K15553', data = x))) %>%
  #Sulfur cycle:
  mutate(ASR_exp = map_dbl(data, function(x) recurse_calc('(((K13811,K00958+K00860,K00955+K00957,K00956+K00957+K00860)+K00390),(K13811+K05907))+(K00380+K00381,K00392)', data = x))) %>%
  mutate(DSR_exp = map_dbl(data, function(x) recurse_calc('K00958+(K00394+K00395)+(K11180+K11181,K27196)+(K27187+K27188+K27189+K27190+K27191)', data = x))) %>%
  mutate(SOX_exp = map_dbl(data, function(x) recurse_calc('K17222+K17223+K17224,K17225,K22622+K17226+K17227', data = x))) %>%
  dplyr::select(-data) %>%
  mutate(SR_exp = ASR_exp + DSR_exp) %>%
  dplyr::select(-ASR_exp, -DSR_exp) %>%
  pivot_longer(cols = ends_with('_exp'), names_to = 'path', values_to = 'exp') %>%
  left_join(metadata) %>%
  modify_at('day', factor, levels = c('day0', 'day7', 'day14', 'day21', 'day28', 'day35')) 

s_transport_fen <- s_transport %>%
  filter(site == 'fen')

s_t_plot <- ggplot(s_transport_fen, aes(x = day, y = exp, color = path, group = path, fill = path)) +
  geom_smooth(linewidth = 2, se = F) +
  scale_x_discrete(labels = c('Day 0', 'Day 7', 'Day 14', 'Day 21', 'Day 35')) +
  cowplot::theme_minimal_grid()  +
  ylab('Expression (geTMM)') +
  xlab(NULL) 

cowplot::plot_grid(n_t_plot, s_t_plot, align = 'h', axis = 'tb', 
                   labels = c('A', 'B'))



# Supplemental Figure 3 Methanogenesis Code: ------------------------------
meth_fen <- arc_methane %>%
  filter(site == 'fen') %>%
  modify_at('day', factor, levels = c('day0', 'day7', 'day14', 'day21', 'day35')) 

meth_pal <- c('#1B9E77', '#D95F02')
names(meth_pal) <- c('acetometh_exp', 'hydrometh_exp')

arc_methane_both <- bf_tc_getmm %>%
  left_join(dram) %>%
  left_join(gene_tax) %>%
  #Pull only Archaea
  filter(domain == 'Archaea') %>%
  dplyr::select(gene, starts_with('STM_'), ko_id) %>%
  mutate(ko_id = strsplit(as.character(ko_id), ',')) %>%
  unnest(ko_id) %>%
  dplyr::select(starts_with('STM_'), ko_id) %>%
  pivot_longer(starts_with('STM_'), names_to = 'sample', values_to = 'getmm') %>%
  rename('label' = ko_id) %>%
  group_by(sample) %>%
  nest() %>%
  #Acetoclastic Methanogenesis
  mutate(acetometh_exp = map_dbl(data, function(x) recurse_calc('((K00925+K00625),K01895)+(K00193+K00197+K00194)', data = x))) %>%
  #Hydrogenotrophic Methanogenesis
  mutate(hydrometh_exp = map_dbl(data, function(x) recurse_calc('(K00200+K00201+K00202+K00203,K11261+(K00205,K11260,K00204))+K00672+K01499+(K00319,K13942)+K00320', data = x))) %>%
  dplyr::select(-data) %>%
  pivot_longer(cols = ends_with('_exp'), names_to = 'path', values_to = 'exp') %>%
  left_join(metadata) %>%
  mutate(sxd = paste0(site, '_', tcat)) %>%
  filter(day %in% c('day0', 'day28', 'day35')) %>%
  modify_at('sxd', factor, levels = c('bog_ti', 'bog_tf', 'fen_ti', 'fen_tf'))

pb <- ggplot(arc_methane_both, aes(x = path, y = exp, fill = sxd)) +
  geom_boxplot() +
  scale_fill_manual(values = pal_sxd,
                    name = NULL,
                    labels = c('Bog Initial',
                               'Bog Final',
                               'Fen Initial',
                               'Fen Final')) +
  cowplot::theme_minimal_grid() +
  xlab(NULL) +
  scale_x_discrete(labels = c('Acetoclastic \n Methanogenesis', 'Hydrogenotrophic \n Methanogenesis')) +
  theme(legend.position.inside = c(0.2,0.9),
        legend.position = 'inside',
        legend.background = element_rect(fill = 'white')) +
  ylab('Metatranscriptome Expression (geTMM)')

pl <- ggplot(meth_fen, aes(x = day, y = exp, color = path, group = path, fill = path)) +
  geom_smooth(linewidth = 2, se = T) +
  scale_color_manual(values = meth_pal, labels = c('Acetoclastic \n Methanogenesis', 'Hydrogenotrophic \n Methanogenesis'),
                     name = 'Pathway') +
  scale_fill_manual(values = meth_pal, labels = c('Acetoclastic \n Methanogenesis', 'Hydrogenotrophic \n Methanogenesis'),
                    name = 'Pathway') +
  scale_x_discrete(labels = c('Day 0', 'Day 7', 'Day 14', 'Day 21', 'Day 35')) +
  cowplot::theme_minimal_grid()  +
  ylab('Metatranscriptome Expression (geTMM)') +
  xlab(NULL)  +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.1, 0.8),
        legend.background = element_rect(fill = 'white'))

cowplot::plot_grid(pl, pb, labels = c('A', 'B'), nrow = 1)


# Figure 2:  --------------------------------------------------------------
#First pull the MAGs that are actually capable of doing metabolism
per_mag_tc <- bf_tc_counts %>%
  left_join(dram) %>%
  left_join(gene_tax) %>%
  group_by(MAG) %>%
  nest() %>%
  filter(!is.na(MAG)) %>%
  #Quantify expressed genes by summing up across each MAG
  mutate(expressed_genes = purrr::map(data, function(x){
    x %>%
      pivot_longer(starts_with('STM_'), names_to = 'sample', values_to = 'count') %>%
      group_by(gene) %>%
      filter(any(count > 0)) %>%
      dplyr::select(gene, sample, count) 
  })) %>%
  #Now how much is each gene expressed?
  mutate(gene_sums = purrr::map(expressed_genes, function(x){
    x %>%
      group_by(sample) %>%
      summarise(total_exp = sum(count))
  }))

mag_expression_final <- per_mag_tc %>%
  dplyr::select(MAG, gene_sums) %>%
  unnest('gene_sums') %>%
  left_join(metadata) %>%
  group_by(site, sample) %>%
  nest() %>%
  #Only MAGS with more than 5 counts are considered expressing
  mutate(active_mags = purrr::map(data, function(x){
    x %>%
      group_by(MAG) %>%
      summarise(total = sum(total_exp)) %>%
      filter(total > 15)
  }))

#Load in the database
db_carbon <- readxl::read_excel('Data/MetaG/metabolism_summary.xlsx',
                                sheet = 'carbon utilization') %>%
  dplyr::select(-gene_description, -module, -header, -subheader) %>%
  distinct() %>%
  column_to_rownames('gene_id') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('MAG')

ox_coded <- mag_expression_final %>%
  unnest('active_mags') %>%
  ungroup() %>%
  dplyr::select(-data, -total)  %>%
  left_join(db_carbon) %>%
  drop_na() %>%
  mutate(step1_EMP = ifelse(K01810 != 0 | K06859 != 0 | K13810 != 0 | K15916 != 0, 1, 0)) %>%
  mutate(step2_EMP = ifelse(K00850 != 0 | K16370 != 0 | K00918 != 0, 1, 0)) %>%
  mutate(step3_EMP = ifelse(K01623 != 0 | K01624 != 0 | K11645 != 0 | K16305 != 0 | K16305 != 0, 1, 0)) %>%
  mutate(step4_EMP = ifelse(K01803 != 0, 1, 0))  %>%
  mutate(step1_ED = ifelse(K00036 != 0, 1, 0))  %>%
  mutate(step2_ED = ifelse(K01057 != 0 | K07404 != 0, 1, 0))  %>%
  mutate(step3_ED = ifelse(K01690 != 0, 1, 0))  %>%
  mutate(step4_ED = ifelse(K01625 != 0, 1, 0)) %>%
  dplyr::select(sample, site, MAG, starts_with('step')) %>%
  group_by(sample, site, MAG) %>%
  reframe(EMP_total = step1_EMP + step2_EMP + step3_EMP + step4_EMP,
          ED_total = step1_ED + step2_ED + step3_ED + step4_ED) %>%
  mutate(EMP_coded = ifelse(EMP_total >= 3, 1, 0)) %>%
  mutate(ED_coded = ifelse(ED_total >= 3, 1, 0)) %>%
  group_by(sample, site) %>%
  reframe(per_EMP = sum(EMP_coded)/length(EMP_coded),
          per_ED = sum(ED_coded)/length(ED_coded))  %>%
  pivot_longer(cols = starts_with('per_'), names_to = 'path')

ox_summary <- ox_coded %>%
  group_by(site, path) %>%
  reframe(m = mean(value),
          sd = sd(value))

ggplot(ox_summary, aes(x = site, y = m, fill = site)) +
  geom_bar(stat = 'summary') +
  geom_errorbar(aes(ymin = m-sd, ymax = m+sd)) +
  facet_wrap(~path, nrow = 1) 

ox_coded %>%
  group_by(path) %>%
  nest() %>%
  mutate(t_test = map_dbl(data, function(x) wilcox.test(value ~ site, data = x)$p.value))

db_all <- read_delim('Data/MetaG/annotations.tsv')

db_ko <- db_all %>%
  dplyr::select(fasta, ko_id) %>%
  dplyr::rename('MAG' = fasta) %>%
  filter(!is.na(ko_id)) %>%
  group_by(MAG, ko_id) %>%
  summarise(n_copy = n()) %>%
  pivot_wider(names_from = 'ko_id', values_from = 'n_copy', values_fill = 0)

ferm_expressed <- mag_expression_final %>%
  unnest('active_mags') %>%
  dplyr::select(-data, -total)  %>%
  left_join(db_ko) %>%
  drop_na() 

ferm_coded <- ferm_expressed %>%
  #Ethanol - 2 Steps 
  mutate(step1_etOH = ifelse(K00132 != 0 | K04072 != 0 | K04073 != 0 | K04073 != 0 | K18366 != 0 | K04021 != 0, 1, 0)) %>%
  mutate(step2_etOH = ifelse(K13953 != 0 | K13954 != 0 | K00001 != 0 | K00121 != 0 | K04072 != 0 | K00002 != 0 | K12957 != 0 | K13979 != 0 | K00114 != 0 | K04022 != 0, 1, 0)) %>%
  #Acetate - 2 Steps
  mutate(step1_AC = ifelse(K13788 != 0, 1, 0)) %>%
  mutate(step2_AC = ifelse(K01512 != 0, 1, 0)) %>%
  #Butyrate - 5 Steps
  mutate(step1_but = ifelse(K00626 != 0, 1, 0)) %>%
  mutate(step2_but = ifelse(K00074 != 0 | K00022 != 0 | K07516 != 0 | K01825 != 0 | K01782 != 0 | K00023 != 0, 1, 0)) %>%
  mutate(step3_but = ifelse(K17865 != 0 | K01692 != 0 , 1, 0)) %>%
  mutate(step4_but = ifelse(K17829 != 0 | K00209 != 0 , 1, 0)) %>%
  mutate(step5a_but = ifelse(K00634 != 0 | K00929 != 0 | K19709 != 0, 1, 0)) %>%
  mutate(step5b_but = ifelse(K01034 != 0 & K01035 != 0, 1, 0)) %>%
  mutate(step5_but = ifelse(step5a_but != 0 | step5b_but != 0, 1, 0)) %>%
  #Lactate - 1 Step
  mutate(step1_lac = ifelse(K00101 != 0 | K00016 != 0 | K00102 != 0 | K03777 != 0, 1, 0)) %>%
  #Propionate - 8 Steps
  mutate(step1_prop = ifelse(K01958 != 0 | K01959 != 0 | K01960 != 0, 1, 0)) %>%
  mutate(step2_prop = ifelse(K00024 != 0 , 1, 0)) %>%
  mutate(step3_prop = ifelse(K01676 != 0 , 1, 0)) %>%
  mutate(step4_prop = ifelse(K00239 != 0 | K00244 != 0 | K18209 != 0, 1, 0)) %>%
  mutate(step5_prop = ifelse(K01902 != 0, 1, 0)) %>%
  mutate(step6_prop = ifelse(K01847 != 0, 1, 0)) %>%
  mutate(step7_prop = ifelse(K05606 != 0, 1, 0)) %>%
  mutate(step8_prop = ifelse(K11264 != 0 | K03416 != 0, 1, 0)) %>%
  group_by(sample, site, MAG) %>%
  reframe(etOH_total = step1_etOH + step2_etOH,
          AC_total = step1_AC + step2_AC,
          but_total = step1_but + step2_but + step3_but + step4_but + step5_but,
          lac_total = step1_lac,
          prop_total = step1_prop + step2_prop + step3_prop + step4_prop + step5_prop + step6_prop + step7_prop + step8_prop) %>%
  mutate(etOH_coded = ifelse(etOH_total == 2, 1, 0),
         AC_coded = ifelse(AC_total == 2, 1, 0),
         but_coded = ifelse(but_total >= 3, 1, 0),
         lac_coded = ifelse(lac_total == 1, 1, 0),
         prop_coded = ifelse(prop_total >= 5, 1, 0)) %>%
  group_by(sample, site) %>%
  reframe(per_etOH = sum(etOH_coded)/length(etOH_coded),
          per_AC = sum(AC_coded)/length(AC_coded),
          per_but = sum(but_coded)/length(but_coded),
          per_lac = sum(lac_coded)/length(lac_coded),
          per_prop = sum(prop_coded)/length(prop_coded)) %>%
  pivot_longer(cols = starts_with('per_'), names_to = 'path')

ferm_summary <- ferm_coded %>%
  group_by(site, path) %>%
  reframe(m = mean(value),
          sd = sd(value))

ggplot(ferm_summary, aes(x = site, y = m, fill = site)) +
  geom_col() +
  geom_errorbar(aes(ymin = m-sd, ymax = m+sd)) +
  facet_wrap(~path, nrow = 1) 

ferm_coded %>%
  group_by(path) %>%
  nest() %>%
  mutate(w_test = map_dbl(data, function(x) wilcox.test(value ~ site, data = x)$p.value))


fox_coded <- rbind(ox_summary, ferm_summary) %>%
  modify_at('path', factor, levels = c('per_EMP', 'per_ED', 'per_etOH', 'per_lac', 'per_prop', 
                                       'per_but', 'per_AC')) %>%
  mutate(m2 = -m)

### Panel A ----

g_plot <- ggplot(fox_coded, aes(x = path, y = m, fill = site)) +
  geom_errorbar(aes(ymin = m-sd, ymax = m+sd), position = position_dodge()) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  #facet_wrap(~path, nrow = 1, scales = 'free_y') +
  ylab('% MAGs Encoding Pathway') +
  scale_fill_manual(values = pal_site) +
  cowplot::theme_cowplot() 

pal_sxd <- c('#2E603D', '#06B25C', '#402BCA', '#20A9FE')
names(pal_sxd) <- c('bog_ti', 'bog_tf', 'fen_ti', 'fen_tf')

t_plot3 <- ggplot(box_data, aes(y = path, x = value, fill = sxd)) +
  geom_boxplot() +
  scale_fill_manual(values = pal_sxd,
                    name = NULL,
                    labels = c('Bog Initial',
                               'Bog Final',
                               'Fen Initial',
                               'Fen Final')) +
  cowplot::theme_minimal_grid() +
  scale_x_continuous(expand = c(0,0),
                     limits = c(-0.03, 6.3)) +
  theme(legend.position.inside = c(0.7,0.9),
        legend.position = 'inside',
        legend.background = element_rect(fill = 'white'),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(5.5,20.5,5.5,0), 'pt'))  +
  ylab(NULL) +
  guides(y = 'none') +
  xlab('Metatranscriptome Expression (geTMM)')

#t_plot3

g_plot3 <- ggplot(fox_coded, aes(y = path, x = m2, fill = site)) +
  geom_errorbar(aes(xmin = m2-sd, xmax = m2+sd), position = position_dodge()) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  #facet_wrap(~path, nrow = 1, scales = 'free_y') +
  xlab('% MAGs Encoding Pathway') +
  scale_fill_manual(values = pal_site,
                    name = NULL,
                    labels = c('Bog', 'Fen')) +
  scale_x_continuous(expand = c(0,0), 
                     limits = c(-1,0),
                     labels = c('100%', '75%', '50%', '25%', '0%')) +
  cowplot::theme_minimal_grid()  +
  theme(plot.margin = unit(c(5.5,0,5.5,20.5), 'pt'),
        legend.position = 'inside',
        legend.position.inside = c(0.2, 0.9),
        legend.background = element_rect(fill = 'white')) +
  guides(y = 'none') +
  ylab(NULL) 

#g_plot3

lab_df <- fox_coded %>%
  rowwise() %>%
  mutate(label = switch(as.character(path),
                        'per_EMP' = 'EMP \n Glycolysis',
                        'per_ED' = 'ED \n Glycolysis',
                        'per_etOH' = 'Ethanol \n Fermentation',
                        'per_lac' = 'Lactate \n Fermentation',
                        'per_prop' = 'Propionate \n Fermentation',
                        'per_but' = 'Butyrate \n Fermentation',
                        'per_AC' = 'Acetate \n Fermentation'))


m_plot <- ggplot(lab_df, aes(y = path, x = 10)) +
  geom_text(aes(label = label)) +
  theme_void() +
  theme(plot.margin = unit(c(5.5,0,5.5,0), 'pt'),
        legend.position = 'none')

fig2_panelA <- cowplot::plot_grid(g_plot3, m_plot, t_plot3, nrow = 1, axis = 'tblr', align = 'h',
                                  rel_widths = c(10,1.5,10))

### Panel B ----

library(dplyr)
library(ggplot2)
library(ggpubr)

pal_sxd <- c('#2E603D', '#06B25C', '#402BCA', '#20A9FE')
names(pal_sxd) <- c('bog_day0', 'bog_day28', 'fen_day0', 'fen_day35')

# Modify the factor levels and labels for the facets
facet_labels <- c(
  bog_day0 = "Bog Initial",
  fen_day0 = "Fen Initial",
  bog_day28 = "Bog Final",
  fen_day35 = "Fen Final"
)

fig2_panelB <- ggplot(mag_balance, aes(x = ox_full, y = ferm_full, color = sxd)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  stat_cor(method = "pearson", label.x = Inf, label.y = -Inf, 
           hjust = 1.5, vjust = -1, size = 3.5) +
  facet_wrap(~sxd, scales = 'free', labeller = labeller(sxd = facet_labels), nrow = 1) +
  scale_color_manual(values = pal_sxd) +
  labs(x = "Oxidative Pressure", y = "Fermentative Processes", color = "Sample Group") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "gray90", color = NA),
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "none")

### Full figure 2 -----
cowplot::plot_grid(
  fig2_panelA, fig2_panelB,
  ncol = 1,
  labels = c('A', 'B'),
  label_size = 16,
  label_fontface = "bold",
  rel_heights = c(1, 1)  # Adjust this to emphasize one plot over another if needed
)


# LC-MS Analysis - WGNCA: -------------------------------------------------
## WGCNA Analysis ----------------------------------------------------------
library(WGCNA)

gas <- read_csv('Data/Gas/gas_data.csv')

#Plot the Gas data:
gas_plot <- gas %>%
  mutate(tsp = str_split(tube, '_')) %>%
  mutate(rep = map_chr(tsp, function(x) x[3])) %>%
  mutate(Day = paste0('D', day)) %>%
  filter(treatment == 'live_unamended') %>%
  rename('CH4' = `headspace CH4 (ppm)`, "CO2" = `headspace CO2 (ppm)`) %>%
  dplyr::select(site, day, CO2, CH4, rep, Day) %>%
  pivot_longer(cols = c(CO2, CH4), names_to = 'gas', values_to = 'ppm')  %>%
  filter(day != 42)


gas_rat <- gas_plot %>%
  pivot_wider(names_from = 'gas', values_from = 'ppm') %>%
  drop_na() %>%
  mutate(rat = CO2/CH4) 

lab_site <- c('bog' = 'Bog', 'fen' = 'Fen')

ratplot <- ggplot(gas_rat, aes(x = day, y = rat, color = site)) +
  geom_point(size = 2) +
  stat_summary(geom = 'line', fun = mean, linewidth = 1.5) +
  facet_wrap(~site, scales = 'free_y', labeller = as_labeller(lab_site)) +
  scale_color_manual(values = pal_site, labels = c('Bog', 'Fen'), name = NULL) +
  cowplot::theme_cowplot() +
  theme(legend.position = c(0.8, 0.8)) +
  ylab(expression('Ratio ' * CO[2]/CH[4])) +
  geom_hline(yintercept = 1)  +
  xlab('Day')

gas_clean <- gas %>%
  mutate(tsp = str_split(tube, '_')) %>%
  mutate(rep = map_chr(tsp, function(x) x[3])) %>%
  mutate(Day = paste0('D', day)) %>%
  filter(site == 'fen') %>%
  filter(treatment == 'live_unamended') %>%
  filter(rep %in% c('A', 'B', 'C')) %>%
  mutate(dxr = paste0(site, "_", Day, "_", rep)) %>%
  dplyr::select(dxr, `headspace CH4 (ppm)`, `headspace CO2 (ppm)`, site) %>%
  rename('CH4' = `headspace CH4 (ppm)`, "CO2" = `headspace CO2 (ppm)`) %>%
  dplyr::select(dxr, CH4, CO2) 


#Prep the normalized data for analysis
rp_small <- rp_final %>%
  filter(Treatment == 'unamended') %>%
  filter(Status == 'Alive') %>%
  dplyr::select(Sample_Name_Kelly, Site, Day, Rep, starts_with('FT_'))

hilic_small <- hilic_final %>%
  filter(Treatment == 'unamended') %>%
  filter(Status == 'Alive') %>%
  dplyr::select(Sample_Name_Kelly, Site, Day, Rep, starts_with('FT_'))

#Combine positive and negative mode
all_metabs <- left_join(rp_small, hilic_small)
#Check to make sure everything merged properly
sum(complete.cases(all_metabs))/nrow(all_metabs)

#Now combine with the gas data
fen_gas <- all_metabs %>%
  filter(Site == 'fen') %>%
  mutate(dxr = paste0(Site, "_", Day, "_", Rep)) %>%
  left_join(gas_clean) %>%
  column_to_rownames('dxr')

#Isolate just metabolite data
mat_fen <- fen_gas %>%
  dplyr::select(starts_with('FT_'))

#Isolate just gas data
meta_gas <- fen_gas %>%
  dplyr::select(CH4, CO2)

## Model Tuning ------------------------------------------------------------

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=40, by=2))
# Call the network topology analysis function
sft <- pickSoftThreshold(mat_fen, powerVector = powers, verbose = 5, corFnc = 'bicor', networkType = 'signed')

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower <- 14 
adjacency <- adjacency(mat_fen, power = softPower, type = 'signed', corFnc = 'bicor')

# Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency, TOMType = 'signed')
dissTOM <- 1-TOM

# Call the hierarchical clustering function
geneTree <- fastcluster::hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize <- 50
# Module identification using dynamic tree cut:
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList <- moduleEigengenes(mat_fen, colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs)
# Cluster module eigengenes
METree <- fastcluster::hclust(as.dist(MEDiss), method = "average")
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
abline(h=0.25, col = 'blue')
# Call an automatic merging function
merge <- mergeCloseModules(mat_fen, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors <- merge$colors
# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs

sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs

# Define numbers of genes and samples
nGenes <- ncol(mat_fen)
nSamples <- nrow(mat_fen)
# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(mat_fen, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)


moduleTraitCor <- WGCNA::cor(MEs, meta_gas, use = "p")
moduleTraitCor_sp <- WGCNA::cor(MEs, meta_gas, use = "p", method = 'sp')
moduleTraitCor_bi <- WGCNA::bicor(MEs, meta_gas, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitPvalue_sp <- corPvalueStudent(moduleTraitCor_sp, nSamples)
moduleTraitPvalue_bi <- corPvalueStudent(moduleTraitCor_bi, nSamples)
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(meta_gas),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Define variable weight containing the weight column of datTrait
meth <- as.data.frame(meta_gas$CH4);
names(meth) <- "CH4"
# names (colors) of the modules
modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(bicor(mat_fen, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) <- paste("MM", modNames, sep="");
names(MMPvalue) <- paste("p.MM", modNames, sep="");

geneTraitSignificance <- as.data.frame(bicor(mat_fen, meth, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) <- paste("GS.", names(meth), sep="");
names(GSPvalue) <- paste("p.GS.", names(meth), sep="")

module = c("darkolivegreen")
column <- match(module, modNames);
moduleGenes <- moduleColors %in% module;

#sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Metabolite significance for Methane",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

## Compound Metadata -------------------------------------------------------
link_hilic_small <- link_hilic %>%
  left_join(canopus_hilic, by = c('Compound_ID' = 'progenesis_id')) %>%
  dplyr::select(-contains('probability'), -contains("all classifications")) %>%
  dplyr::select(Compound_ID, Name, annotation_type, contains("ClassyFire"), contains('NPC'), MS2)

link_rp_small <- link_rp %>%
  left_join(canopus_rp, by = c('Compound_ID' = 'progenesis_id')) %>%
  dplyr::select(-contains('probability'), -contains("all classifications")) %>%
  dplyr::select(Compound_ID, Name, annotation_type, contains("ClassyFire"), contains('NPC'), MS2)

all_link <- rbind(link_hilic_small, link_rp_small)


#First pull the Dark Magenta Genes which have a positive correlation:
#Compoun names part of dm module
dm_genes <- names(mat_fen)[moduleColors=="darkolivegreen"]

#Now pull the metabolites using the bicor
dm_metabs_bi <- data.frame(bicor(mat_fen, meth, use = "p")) %>%
  rownames_to_column('Compound_ID') %>%
  filter(Compound_ID %in% dm_genes) %>%
  #Pull only the top most correlated genes
  arrange(desc(abs(CH4)))

#Let's plot the top 10:
top_bi <- fen_gas %>%
  #Take only the top 10
  dplyr::select(CH4, Day, dm_metabs_bi$Compound_ID[1:10]) %>%
  pivot_longer(starts_with('FT_'), names_to = 'ft')  %>%
  #Factor them so they're in order
  modify_at('ft', factor, levels = dm_metabs_bi$Compound_ID[1:10])

day_pal2 <- rev(c('#9d0208', "#d00000", "#e85d04", "#f48c06", "#ffba08"))
names(day_pal2) <- c('D0', 'D7', 'D14', 'D21', 'D35')

#Plot the top 10
ggplot(top_bi, aes(x = value, y = CH4, color = Day)) +
  geom_smooth(method = 'lm', se = F, aes(group = 1), color = 'grey') +
  geom_point(size = 2) +
  scale_color_manual(values = day_pal2) +
  facet_wrap(~ft, scales = 'free_x', ncol = 5) +
  cowplot::theme_cowplot() +
  scale_y_continuous(expand = c(0,0))

#Do the smae but based on Spearman's correlation
dm_metabs_sp <- data.frame(cor(mat_fen, meth, use = "p", method = 'sp')) %>%
  rownames_to_column('Compound_ID') %>%
  filter(Compound_ID %in% dm_genes) %>%
  arrange(desc(CH4))

#Let's plot the top 10:
top_sp <- fen_gas %>%
  dplyr::select(CH4, Day, dm_metabs_sp$Compound_ID[1:10]) %>%
  pivot_longer(starts_with('FT_'), names_to = 'ft')  %>%
  modify_at('ft', factor, levels = dm_metabs_sp$Compound_ID[1:10])

ggplot(top_sp, aes(x = value, y = CH4, color = Day)) +
  geom_smooth(method = 'lm', se = F, aes(group = 1), color = 'grey') +
  geom_point(size = 2) +
  scale_color_manual(values = day_pal2) +
  facet_wrap(~ft, scales = 'free_x', ncol = 5) +
  cowplot::theme_cowplot() +
  scale_y_continuous(expand = c(0,0))

## Raw Data ----------------------------------------------------------------

### Data Prep ---------------------------------------------------------------

rp_plot <- rp_raw %>%
  rownames_to_column('Filename') %>%
  left_join(meta_rp) %>%
  filter(str_detect(Filename, 'D'))  %>%
  filter(Day != 'D42') %>%
  modify_at('Day', factor, levels = c('D0', 'D7', 'D14', 'D21', 'D28', 'D35')) %>%
  modify_at('Site', factor, levels = c('bog', 'fen'))

hilic_plot <- hilic_raw %>%
  rownames_to_column('Filename') %>%
  left_join(meta_hilic) %>%
  filter(Day != 'D42') %>%
  filter(str_detect(Filename, 'D'))  %>%
  modify_at('Day', factor, levels = c('D0', 'D7', 'D14', 'D21', 'D28', 'D35')) %>%
  modify_at('Site', factor, levels = c('bog', 'fen'))

rp_plot_small <- rp_plot %>%
  filter(Treatment == 'unamended') %>%
  filter(Status == 'Alive') %>%
  dplyr::select(Sample_Name_Kelly, Site, Day, Rep, starts_with('FT_'))
sum(complete.cases(rp_plot_small))/nrow(rp_plot_small)

hilic_plot_small <- hilic_plot %>%
  filter(Treatment == 'unamended') %>%
  filter(Status == 'Alive') %>%
  dplyr::select(Sample_Name_Kelly, Site, Day, Rep, starts_with('FT_'))

sum(complete.cases(hilic_plot_small))/nrow(hilic_plot_small)

all_plot <- left_join(rp_plot_small, hilic_plot_small)
#Check to make sure everything merged properly
sum(complete.cases(all_plot))/nrow(all_plot)

plot_fen_gas <- all_plot %>%
  filter(Site == 'fen') %>%
  mutate(dxr = paste0(Site, "_", Day, "_", Rep)) %>%
  left_join(gas_clean) %>%
  column_to_rownames('dxr')
sum(complete.cases(plot_fen_gas))/nrow(plot_fen_gas)

mat_fen_raw <- plot_fen_gas %>%
  dplyr::select(starts_with('FT_'))

meta_gas_raw <- plot_fen_gas %>%
  dplyr::select(CH4, CO2)

meth_raw <- as.data.frame(meta_gas_raw$CH4);
names(meth_raw) <- "CH4"

### Evaluation --------------------------------------------------------------

dm_metabs_sp_raw <- data.frame(cor(mat_fen_raw, meth_raw, use = "p", method = 'sp')) %>%
  rownames_to_column('Compound_ID') %>%
  filter(Compound_ID %in% dm_genes) %>%
  arrange(desc(CH4))

top_sp_raw <- plot_fen_gas %>%
  dplyr::select(CH4, Day, dm_metabs_sp_raw$Compound_ID[1:10]) %>%
  pivot_longer(starts_with('FT_'), names_to = 'ft')  %>%
  modify_at('ft', factor, levels = dm_metabs_sp_raw$Compound_ID[1:10])

ggplot(top_sp_raw, aes(x = value, y = CH4, color = Day)) +
  geom_smooth(method = 'lm', se = F, aes(group = 1), color = 'grey') +
  geom_point(size = 2) +
  scale_color_manual(values = day_pal2) +
  facet_wrap(~ft, scales = 'free_x', ncol = 5) +
  cowplot::theme_cowplot() +
  scale_y_continuous(expand = c(0,0))

## Metadata of Correlated Features ----------------------------------------

high_cor <- dm_metabs_sp_raw %>%
  left_join(all_link) %>%
  filter(MS2 != 'No MS2') %>%
  slice_max(CH4, n = 100, with_ties = F)
#79 feautres:
old_mb <- read_csv('highcor_old_pos.cv')
#write_excel_csv(high_cor, file = 'Pos_cor_new.csv')

dm_nosc <- dm_metabs_sp_raw %>%
  left_join(link_hilic)  %>%
  left_join(canopus_hilic, by = c('Compound_ID' = 'progenesis_id')) %>%
  dplyr::select(-contains('probability'), -contains("all classifications")) %>%
  dplyr::select(Compound_ID, molecularFormula, Formula, Name, annotation_type, contains("ClassyFire"), contains('NPC')) %>%
  mutate(chemFormula = ifelse(annotation_type %in% c('L1', 'L2') | is.na(molecularFormula), Formula, molecularFormula)) %>%
  dplyr::select(Compound_ID, chemFormula, Name, annotation_type, contains("ClassyFire"), contains('NPC')) %>%
  filter(!is.na(chemFormula)) %>%
  dplyr::select(Compound_ID, chemFormula) %>%
  #filter(Compound_ID %in% c(d0sig_bog_rp, d28sig_bog_rp, d0sig_fen_rp, d28sig_fen_rp)) %>%
  mutate(breakdown = purrr::map(chemFormula, function(x) data.frame(OrgMassSpecR::ListFormula(x)))) %>%
  unnest('breakdown') %>%
  mutate(NOSC = 4 - ((4*C + H - 2*O - 3*N - 2*S + 5*P)/C)) %>%
  magrittr::inset('Correlation', value = 'positive')

high_cor %>%
  group_by(`ClassyFire#most specific class`) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

high_cor %>%
  group_by(`ClassyFire#level 5`) %>%
  summarise(n = n()) %>%
  arrange(desc(n))


sum(high_cor$Compound_ID %in% old_mb$Compound_ID)

ugh <- high_cor %>%
  filter(Compound_ID %in% old_mb$Compound_ID)

## Light Cyan - Negative Correlation ---------------------------------------

### Transformed Data --------------------------------------------------------

lc_genes <- names(mat_fen)[moduleColors=="brown"]

lc_metabs_bi <- data.frame(bicor(mat_fen, meth, use = "p")) %>%
  rownames_to_column('Compound_ID') %>%
  filter(Compound_ID %in% lc_genes) %>%
  arrange(desc(abs(CH4)))

lc_metabs_sp <- data.frame(cor(mat_fen, meth, use = "p", method = 'sp')) %>%
  rownames_to_column('Compound_ID') %>%
  filter(Compound_ID %in% lc_genes) %>%
  arrange(desc(abs(CH4)))

#Let's plot the top 10:
top_bi_lc <- fen_gas %>%
  dplyr::select(CH4, Day, lc_metabs_bi$Compound_ID[1:10]) %>%
  pivot_longer(starts_with('FT_'), names_to = 'ft')  %>%
  modify_at('ft', factor, levels = lc_metabs_bi$Compound_ID[1:10])

ggplot(top_bi_lc, aes(x = value, y = CH4, color = Day)) +
  geom_smooth(method = 'lm', se = F, aes(group = 1), color = 'grey') +
  geom_point(size = 2) +
  scale_color_manual(values = day_pal2) +
  facet_wrap(~ft, scales = 'free_x', ncol = 5) +
  cowplot::theme_cowplot() +
  scale_y_continuous(expand = c(0,0))

#Let's plot the top 10:
top_sp_lc <- fen_gas %>%
  dplyr::select(CH4, Day, lc_metabs_sp$Compound_ID[1:10]) %>%
  pivot_longer(starts_with('FT_'), names_to = 'ft')  %>%
  modify_at('ft', factor, levels = lc_metabs_sp$Compound_ID[1:10])

ggplot(top_sp_lc, aes(x = value, y = CH4, color = Day)) +
  geom_smooth(method = 'lm', se = F, aes(group = 1), color = 'grey') +
  geom_point(size = 2) +
  scale_color_manual(values = day_pal2) +
  facet_wrap(~ft, scales = 'free_x', ncol = 5) +
  cowplot::theme_cowplot() +
  scale_y_continuous(expand = c(0,0))

### Raw Data ----------------------------------------------------------------

lc_metabs_sp_raw <- data.frame(cor(mat_fen_raw, meth_raw, use = "p", method = 'sp')) %>%
  rownames_to_column('Compound_ID') %>%
  filter(Compound_ID %in% lc_genes) %>%
  arrange(desc(abs(CH4)))

top_sp_raw_lc <- plot_fen_gas %>%
  dplyr::select(CH4, Day, lc_metabs_sp_raw$Compound_ID[1:10]) %>%
  pivot_longer(starts_with('FT_'), names_to = 'ft')  %>%
  modify_at('ft', factor, levels = lc_metabs_sp_raw$Compound_ID[1:10])

ggplot(top_sp_raw_lc, aes(x = value, y = CH4, color = Day)) +
  geom_smooth(method = 'lm', se = F, aes(group = 1), color = 'grey') +
  geom_point(size = 2) +
  scale_color_manual(values = day_pal2) +
  facet_wrap(~ft, scales = 'free_x', ncol = 5) +
  cowplot::theme_cowplot() #+

lc_metabs_bi_raw <- data.frame(bicor(mat_fen_raw, meth_raw, use = "p")) %>%
  rownames_to_column('Compound_ID') %>%
  filter(Compound_ID %in% lc_genes) %>%
  arrange(desc(abs(CH4)))

top_bi_raw_lc <- plot_fen_gas %>%
  dplyr::select(CH4, Day, lc_metabs_bi_raw$Compound_ID[1:10]) %>%
  pivot_longer(starts_with('FT_'), names_to = 'ft')  %>%
  modify_at('ft', factor, levels = lc_metabs_bi_raw$Compound_ID[1:10])

ggplot(top_bi_raw_lc, aes(x = value, y = CH4, color = Day)) +
  geom_smooth(method = 'lm', se = F, aes(group = 1), color = 'grey') +
  geom_point(size = 2) +
  scale_color_manual(values = day_pal2) +
  facet_wrap(~ft, scales = 'free_x', ncol = 5) +
  cowplot::theme_cowplot()

### Metabolite Metadata -----------------------------------------------------

high_cor_lc <- lc_metabs_sp_raw %>%
  left_join(all_link) %>%
  filter(MS2 != 'No MS2') %>%
  slice_min(CH4, n = 100, with_ties = F)

#write_excel_csv(high_cor_lc, file = 'Neg_cor_new.csv')

classes <- high_cor_lc %>%
  group_by(`ClassyFire#level 5`) %>%
  summarise(n_class = n()) %>%
  arrange(desc(n_class))

classes_pos <- high_cor %>%
  group_by(`ClassyFire#level 5`) %>%
  summarise(n_class = n()) %>%
  arrange(desc(n_class))

#Compare back against sig testing:
lm_compare <- rbind(lm_rp_contra_adj, lm_hilic_contra_adj) %>%
  left_join(rbind(hilic_fold, rp_fold)) %>%
  dplyr::select(ft, T0_BvF_adj, T0_l2fc)

#Negative Correlation
lm_compare %>%
  filter(ft %in% high_cor_lc$Compound_ID) %>%
  mutate(sig = ifelse(T0_BvF_adj < 0.05, 'sig', 'ns')) %>%
  mutate(higher = ifelse(T0_l2fc > 0, 'bog', 'fen')) %>%
  group_by(sig, higher) %>%
  summarise(n = n())

#Positive Correlation
lm_compare %>%
  filter(ft %in% high_cor$Compound_ID) %>%
  mutate(sig = ifelse(T0_BvF_adj < 0.05, 'sig', 'ns')) %>%
  mutate(higher = ifelse(T0_l2fc > 0, 'bog', 'fen')) %>%
  group_by(sig, higher) %>%
  summarise(n = n())


#old_lc <- read_csv('high_cor_old_neg.csv')

all_formula <- rbind(dplyr::select(link_hilic, Compound_ID, Name, Formula), 
                     dplyr::select(link_rp, Compound_ID, Name, Formula))

lc_nosc <- high_cor_lc %>%
  left_join(all_formula) %>%
  filter(!is.na(Formula)) %>%
  mutate(breakdown = purrr::map(Formula, function(x) data.frame(OrgMassSpecR::ListFormula(x)))) %>%
  unnest('breakdown') %>%
  mutate(NOSC = 4 - ((4*C + H - 2*O - 3*N - 2*S + 5*P)/C)) %>%
  magrittr::inset('Correlation', value = 'negative')

mb_nosc <- high_cor %>%
  left_join(all_formula) %>%
  filter(!is.na(Formula)) %>%
  mutate(breakdown = purrr::map(Formula, function(x) data.frame(OrgMassSpecR::ListFormula(x)))) %>%
  unnest('breakdown') %>%
  mutate(NOSC = 4 - ((4*C + H - 2*O - 3*N - 2*S + 5*P)/C)) %>%
  magrittr::inset('Correlation', value = 'positive')


high_cor_pos_summary <- high_cor %>%
  group_by(`ClassyFire#level 5`) %>%
  summarise(n = n()) %>%
  magrittr::inset('correlation', value = 'positive')

high_cor_neg_summary <- high_cor_lc %>%
  group_by(`ClassyFire#level 5`) %>%
  summarise(n = n()) %>%
  magrittr::inset('correlation', value = 'negative')

all_cor_figure <- rbind(high_cor_pos_summary, high_cor_neg_summary)

all_nosc <- rbind(lc_nosc, mb_nosc) %>%
  filter(abs(NOSC) < 4) 

ggplot(all_nosc, aes(x = Correlation, y = NOSC, fill = Correlation)) +
  geom_boxplot()

all_nosc %>%
  group_by(Correlation) %>%
  summarise(m = mean(NOSC))

wilcox.test(NOSC ~ Correlation, data = all_nosc)

lm_output <- rbind(lm_rp_cr, lm_hilic_cr) %>%
  mutate(higher = ifelse(sign(effect_size) == 1, 'Bog', 'Fen')) %>%
  mutate(sig = ifelse(pvalue < 0.05, 'significant', 'non-signficicant')) %>%
  dplyr::select(compound_name, higher, sig) %>%
  dplyr::rename('Compound_ID' = compound_name, 'Site More Abundant' = higher, 'Significantly Different' = sig)


out_neg <- high_cor_lc %>%
  left_join(lm_output) %>%
  dplyr::select(Compound_ID, CH4, Name, annotation_type, `Site More Abundant`, `Significantly Different`, 
                starts_with('ClassyFire'), starts_with('NPC'), MS2)
out_pos <- high_cor %>%
  left_join(lm_output) %>%
  dplyr::select(Compound_ID, CH4, Name, annotation_type, `Site More Abundant`, `Significantly Different`, 
                starts_with('ClassyFire'), starts_with('NPC'), MS2)



# Figure 4 ----------------------------------------------------------------

## Panel A: ------
gas_plot_mkd <- gas_plot %>%
  mutate(gas_mkd = ifelse(gas == 'CO2', 'CO<sub>2', 'CH<sub>4'))

#Gas Plot to Start:
gplot <- ggplot(gas_plot_mkd, aes(x = day, y = ppm, color = site, group = site, fill = site)) +
  #geom_point(size = 2) +
  geom_smooth(se = T, linewidth = 2) +
  facet_wrap(~gas_mkd, scales = 'free_y' )  +
  scale_color_manual(values = pal_site, labels = c('Bog','Fen'), name = NULL)  +
  scale_fill_manual(values = pal_site, labels = c('Bog','Fen'), name = NULL)  +
  cowplot::theme_minimal_grid() +
  theme(legend.position = 'none')  +
  xlab('Day')  +
  scale_x_continuous(breaks = c(0, 7, 14, 21, 28, 35),
                     labels = c(c(0,7,14,21,28,35))) +
  theme(legend.position = 'inside',
        strip.text = ggtext::element_markdown(face = 'bold', size = 18),
        strip.background = NULL,
        legend.position.inside = c(0.05, 0.85),
        legend.background = element_rect(fill = 'white'))
gplot

## Panels D E F: Individual Metabolites -----
#Next each individual metabolite plot:
#First set up the palette for each individual day in the fen:
day_pal_cor <- c( '#01267A', '#0238B6', '#0255B7', '#018ADF', '#5DC0FE')

#First Threonic Acid
metab_plot <- ggplot(all_plot, aes(x = Day, y = FT_10.169_135.02995_HILIC, color = Site, group = Site)) +
  geom_point(size = 4) +
  stat_summary(fun = 'mean', geom = 'line', linewidth = 2) +
  scale_color_manual(values = pal_site, name = NULL, labels = c('Bog', 'Fen')) +
  scale_x_discrete(labels = seq(0, 35, by = 7)) +
  cowplot::theme_minimal_grid() +
  ggtitle('Threonic Acid') +
  ylab('Intensity') +
  theme(legend.position = 'none')

cor_plot <- ggplot(plot_fen_gas, aes(x = FT_10.169_135.02995_HILIC, y = CH4, color = Day)) +
  geom_smooth(method = 'lm', se = F, linewidth = 2, aes(group = 1),
              formula = y ~ log(x-min(plot_fen_gas$FT_10.169_135.02995_HILIC)-1), color = 'darkgrey') +
  geom_point(size = 4) +
  scale_color_manual(values = day_pal_cor,
                     name = NULL,
                     labels = c(paste('Day', c(0, 7, 14, 21, 35)))) +
  cowplot::theme_cowplot() +
  theme(legend.position = 'none',
        plot.margin = unit(c(5.5, 18.5, 5.5, 5.5), 'pt'))+
  ylab(expression(CH[4]*" Concentration (ppm)")) +
  xlab('Threonic Acid Intensity') +
  annotate(geom = 'text', x = 3e8, y = 10000, label = expression(rho * ' = -0.82'), size = 5) 

thr_panel <- cowplot::plot_grid(metab_plot, cor_plot, nrow = 1, align = 'h', axis = 'tb')

# Second Phenoxy Compound

ph_metab <- ggplot(all_plot, aes(x = Day, y = FT_2.717_174.02417_HILIC, color = Site, group = Site)) +
  geom_point(size = 4) +
  stat_summary(fun = 'mean', geom = 'line', linewidth = 2) +
  scale_color_manual(values = pal_site) +
  scale_x_discrete(labels = seq(0, 35, by = 7)) +
  cowplot::theme_minimal_grid() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold', size = 18),
        strip.background = NULL) +
  ggtitle('Phenoxy Compound') +
  ylab('Intensity')

ph_cor <- ggplot(plot_fen_gas, aes(x = FT_2.717_174.02417_HILIC, y = CH4, color = Day)) +
  geom_smooth(method = 'lm', se = F, linewidth = 2, aes(group = 1),
              formula = y ~ log(x-min(plot_fen_gas$FT_2.717_174.02417_HILIC)-1), color = 'darkgrey') +
  geom_point(size = 4) +
  scale_color_manual(values = day_pal_cor) +
  cowplot::theme_cowplot() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold', size = 18),
        strip.background = NULL) +
  xlab('Pheonxy Compound Intensity') +
  annotate(geom = 'text', x = 9.8e7, y = 10000, label = expression(rho * ' = -0.87'), size = 5) +
  theme(legend.position = 'none',
        #strip.text = element_text(face = 'bold', size = 18),
        plot.margin = unit(c(5.5, 18.5, 5.5, 5.5), 'pt'),
        strip.background = NULL) +
  # axis.text = element_text(size = 18),
  # axis.title = element_text(size = 20)) +
  ylab(expression(CH[4]*" Concentration (ppm)")) 

ph_panel <- cowplot::plot_grid(ph_metab, ph_cor, nrow = 1, align = 'h', axis = 'tb')

# third: Biflavanoid Compound

bf_metab <- ggplot(all_plot, aes(x = Day, y = FT_6.267_593.1301_HILIC, color = Site, group = Site)) +
  geom_point(size = 4) +
  stat_summary(fun = 'mean', geom = 'line', linewidth = 2) +
  scale_color_manual(values = pal_site) +
  scale_x_discrete(labels = seq(0, 35, by = 7)) +
  cowplot::theme_minimal_grid() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold', size = 18),
        strip.background = NULL) +
  ggtitle('Biflavanoid Compound') +
  ylab('Intensity')

bf_cor <- ggplot(plot_fen_gas, aes(x = FT_6.267_593.1301_HILIC, y = CH4, color = Day)) +
  stat_smooth(method = 'lm', formula = y ~ poly(x,2), se= FALSE, aes(group = 1), color = 'darkgrey', linewidth = 2) +
  geom_point(size = 4) +
  scale_color_manual(values = day_pal_cor) +
  cowplot::theme_cowplot() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold', size = 18),
        strip.background = NULL) +
  xlab('Biflavonoid Compound Intensity') +
  scale_x_continuous(labels = scales::label_scientific()) +
  annotate(geom = 'text', x = 3.1e6, y = 10000, label = expression(rho * ' = 0.95'), size = 5) +
  theme(legend.position = 'none',
        strip.background = NULL,
        plot.margin = unit(c(5.5, 18.5, 5.5, 5.5), 'pt')) +
  ylab(expression(CH[4]*" Concentration (ppm)")) 

bf_panel <- cowplot::plot_grid(bf_metab, bf_cor, nrow = 1, align = 'h', axis = 'tb')

#Now plotting Transcripts for EMP and Methanogenesis:

## Panel B ----
m_plot <- arc_methane %>%
  group_by(sample) %>%
  summarise(ch4_exp = sum(exp)) %>%
  left_join(rbind(ox, ferm)) %>%
  pivot_wider(names_from = 'path', values_from = 'exp') %>%
  filter(site == 'fen')

egm_plot <- ggplot(m_plot, aes(x = EMP_exp, y = ch4_exp, color = day)) +
  geom_smooth(method = 'lm', se = F, linewidth = 2, aes(group = 1),
              formula = y ~ log(x-min(m_plot$ch4_exp)), color = 'darkgrey') +
  geom_point(size = 4) +
  scale_color_manual(values = day_pal_cor,
                     name = NULL,
                     labels = c(paste('Day', c(0, 7, 14, 21, 35)))) +
  #scale_color_manual(values = pal_site) +
  cowplot::theme_cowplot() +
  ylab('Methanogenesis (geTMM)') +
  xlab('EMP Glycolysis (geTMM)') +
  annotate(geom = 'text', x = 1.4, y = 8.5, label = expression(rho * ' = 0.82'), size = 5) +
  theme(legend.position = 'inside',
        strip.background = NULL,
        legend.position.inside = c(0.5, 0.3),
        legend.direction = 'horizontal',
        legend.background = element_rect(fill = alpha('lightgrey', alpha = 0.45))) +
  guides(col = guide_legend(ncol = 2)) 

## Plot C -----
s_data <- expression %>%
  pivot_longer(cols = ends_with('_exp'), names_to = 'path', values_to = 'exp') %>%
  left_join(metadata) %>%
  mutate(sxd = paste0(site, '_', tcat)) %>%
  modify_at('sxd', factor, levels = c('bog_ti', 'bog_tf', 'fen_ti', 'fen_tf'))

splot <- ggplot(s_data, aes(x = path, y = exp, fill = sxd)) +
  geom_boxplot() +
  scale_fill_manual(values = pal_sxd,
                    labels = c('Bog Initial', 'Bog Final', 'Fen Initial', 'Fen Final'),
                    name = NULL) +
  cowplot::theme_minimal_grid() +
  #xlab('Gene') +
  xlab(NULL) +
  ylab('Expression (geTMM)') +
  scale_x_discrete(labels = c('Catalase', 'Glutathione \n Peroxidase', 'Thioredoxin')) + 
  theme(legend.position = 'top',
        plot.margin = unit(c(2, 5.5, 2, 5.5), 'pt'))

## Final fig 4: ----
cowplot::plot_grid(gplot, thr_panel, egm_plot, ph_panel, splot,  bf_panel, ncol = 2,
                   rel_widths = c(7,10), labels = c('A', 'D', 'B', 'E', 'C', 'F'))

# Supplemental Figure 4: --------------------------------------------------

aro_plot <- aromatics %>%
  filter(day %in% c('day0', 'day28', 'day35')) %>%
  modify_at('path', factor, levels  = c('bzcoa_exp', 'phloro_exp', 'caf_exp'))

ggplot(aro_plot, aes(x = path, y  = exp, fill = sxd)) +
  geom_boxplot() +
  scale_fill_manual(values = pal_sxd,
                    labels = c('Bog Initial', 'Bog Final', 'Fen Initial', 'Fen Final'),
                    name = NULL) +
  cowplot::theme_minimal_grid() +
  xlab(NULL) +
  ylab('Expression (geTMM)') +
  scale_x_discrete(labels = c('Benzoyl-CoA \n Pathway', 'Phloroglucinol \n Pathway', 'Caffeate \n Respiration')) 


# Supplemental Figure 5: --------------------------------------------------

alt_sinks <- bf_dram_tc %>%
  mutate(WL_exp = map_dbl(data, function(x) recurse_calc('K00198+(K05299,K15022,K22015+K25123+K25124)+(K01938+K01491,K01500)+(K00297+K25007+K25008)+K15023+K14138', data = x))) %>%
  mutate(NiFe_Hyd_exp = map_dbl(data, function(x) recurse_calc('K00437+K18008', data = x))) %>%
  mutate(FeFe_Hyd_exp = map_dbl(data, function(x) recurse_calc('K17997+K17998+K17999', data = x))) %>%
  mutate(Fe_Hyd_exp = map_dbl(data, function(x) recurse_calc('K00532+K00533+K00534', data = x))) %>%
  dplyr::select(-data) %>%
  pivot_longer(cols = ends_with('_exp'), names_to = 'path', values_to = 'exp') %>%
  left_join(metadata) %>%
  mutate(sxd = paste0(site, '_', tcat)) %>%
  filter(day %in% c('day0', 'day28', 'day35')) %>%
  modify_at('sxd', factor, levels = c('bog_ti', 'bog_tf', 'fen_ti', 'fen_tf')) %>%
  modify_at('path', factor, levels = c('WL_exp','Fe_Hyd_exp', 'FeFe_Hyd_exp', 'NiFe_Hyd_exp'))

pal_sxd <- c('#2E603D', '#06B25C', '#402BCA', '#20A9FE')
names(pal_sxd) <- c('bog_ti', 'bog_tf', 'fen_ti', 'fen_tf')

ggplot(alt_sinks, aes(x = path, y = exp, fill = sxd)) +
  geom_boxplot() +
  scale_fill_manual(values = pal_sxd,
                    name = NULL,
                    labels = c('Bog Initial',
                               'Bog Final',
                               'Fen Initial',
                               'Fen Final')) +
  cowplot::theme_minimal_grid() +
  xlab(NULL) +
  scale_x_discrete(labels = c('Wood-Ljungdahl Pathway', '[Fe]-Hydrogenase', '[FeFe]-Hydrogenase', '[NiFe]-Hydrogenase')) +
  theme(legend.position.inside = c(0.8,0.9),
        legend.position = 'inside',
        legend.background = element_rect(fill = 'white')) +
  ylab('Metatranscriptome Expression (geTMM)')

##Code to conduct the two sample MR analysis in "Adiposity, metabolites, and colorectal cancer risk: Mendelian randomization study"

rm(list=ls(all=TRUE)) 

library(devtools)
library(TwoSampleMR)
library(ggplot2)

ao <- available_outcomes()

##analysis of adipose traits on colorectal cancer 
##snp exposure associations available from https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files
##snp outcome associations available following an application to the Genetics and Epidemiology of Colorectal Cancer Consortium (GECCO) https://www.fredhutch.org/en/research/divisions/public-health-sciences-division/research/cancer-prevention/genetics-epidemiology-colorectal-cancer-consortium-gecco.html

#bmi overall
bmi_dat <- read_exposure_data(
  filename = "pulit_bmi_combined.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta.combined",
  se_col = "se.combined",
  effect_allele_col = "A1.combined",
  other_allele_col = "A2.combined",
  eaf_col = "frqA1.combined",
  pval_col = "pval.combined",
  phenotype_col = "exposure"
)
bmi_dat <- clump_data(bmi_dat)
bmi_dat$t_stat <- (bmi_dat$beta.exposure/bmi_dat$se.exposure)
bmi_dat$f_stat <- (bmi_dat$t_stat)^2
bmi_fstat <- mean(bmi_dat$f_stat)
bmi_fstat <- as.data.frame(bmi_fstat)

#bmi female
bmi_fem_dat <- read_exposure_data(
  filename = "pulit_bmi_female.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta.combined",
  se_col = "se.combined",
  effect_allele_col = "A1.combined",
  other_allele_col = "A2.combined",
  eaf_col = "frqA1.combined",
  pval_col = "pval.combined",
  phenotype_col = "exposure"
)
bmi_fem_dat <- clump_data(bmi_fem_dat)
bmi_fem_dat$t_stat <- (bmi_fem_dat$beta.exposure/bmi_fem_dat$se.exposure)
bmi_fem_dat$f_stat <- (bmi_fem_dat$t_stat)^2
bmifem_fstat <- mean(bmi_fem_dat$f_stat)
bmifem_fstat <- as.data.frame(bmifem_fstat)

#bmi male
bmi_mal_dat <- read_exposure_data(
  filename = "pulit_bmi_male.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta.combined",
  se_col = "se.combined",
  effect_allele_col = "A1.combined",
  other_allele_col = "A2.combined",
  eaf_col = "frqA1.combined",
  pval_col = "pval.combined",
  phenotype_col = "exposure"
)
bmi_mal_dat <- clump_data(bmi_mal_dat)
bmi_mal_dat$t_stat <- (bmi_mal_dat$beta.exposure/bmi_mal_dat$se.exposure)
bmi_mal_dat$f_stat <- (bmi_mal_dat$t_stat)^2
bmimal_fstat <- mean(bmi_mal_dat$f_stat)
bmimal_fstat <- as.data.frame(bmimal_fstat)

#whr overall
whr_dat <- read_exposure_data(
  filename = "pulit_whr_combined.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta.combined",
  se_col = "se.combined",
  effect_allele_col = "A1.combined",
  other_allele_col = "A2.combined",
  eaf_col = "frqA1.combined",
  pval_col = "pval.combined",
  phenotype_col = "exposure"
)
whr_dat <- clump_data(whr_dat)
whr_dat$t_stat <- (whr_dat$beta.exposure/whr_dat$se.exposure)
whr_dat$f_stat <- (whr_dat$t_stat)^2
whr_fstat <- mean(whr_dat$f_stat)
whr_fstat <- as.data.frame(whr_fstat)

#whr female
whr_fem_dat <- read_exposure_data(
  filename = "pulit_whr_female.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta.combined",
  se_col = "se.combined",
  effect_allele_col = "A1.combined",
  other_allele_col = "A2.combined",
  eaf_col = "frqA1.combined",
  pval_col = "pval.combined",
  phenotype_col = "exposure"
)
whr_fem_dat <- clump_data(whr_fem_dat)
whr_fem_dat$t_stat <- (whr_fem_dat$beta.exposure/whr_fem_dat$se.exposure)
whr_fem_dat$f_stat <- (whr_fem_dat$t_stat)^2
whrfem_fstat <- mean(whr_fem_dat$f_stat)
whrfem_fstat <- as.data.frame(whrfem_fstat)

#whr male
whr_mal_dat <- read_exposure_data(
  filename = "pulit_whr_male.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta.combined",
  se_col = "se.combined",
  effect_allele_col = "A1.combined",
  other_allele_col = "A2.combined",
  eaf_col = "frqA1.combined",
  pval_col = "pval.combined",
  phenotype_col = "exposure"
)
whr_mal_dat <- clump_data(whr_mal_dat)
whr_mal_dat$t_stat <- (whr_mal_dat$beta.exposure/whr_mal_dat$se.exposure)
whr_mal_dat$f_stat <- (whr_mal_dat$t_stat)^2
whrmal_fstat <- mean(whr_mal_dat$f_stat)
whrmal_fstat <- as.data.frame(whrmal_fstat)

#colorectal cancer associations
#overall
crc_gecco_overall_dat <- read_outcome_data(
  filename = "1241_GECCO_overall_combined.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P.value",
  phenotype_col = "outcome"
)

#female
crc_gecco_fem_dat <- read_outcome_data(
  filename = "1241_GECCO_female_combined_030219.txt",
  sep = "\t",
  snp_col = "rs",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P.value",
  phenotype_col = "outcome"
)

#male
crc_gecco_mal_dat <- read_outcome_data(
  filename = "1241_GECCO_male_combined.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P.value",
  phenotype_col = "outcome"
)

#colon
crc_gecco_colon_dat <- read_outcome_data(
  filename = "1241_GECCO_colon_combined.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P.value",
  phenotype_col = "outcome"
)

#proximal
crc_gecco_proximal_dat <- read_outcome_data(
  filename = "1241_GECCO_proximal_combined.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P.value",
  phenotype_col = "outcome"
)

#distal
crc_gecco_distal_dat <- read_outcome_data(
  filename = "1241_GECCO_distal_combined.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P.value",
  phenotype_col = "outcome"
)

#rectal
crc_gecco_rectal_dat <- read_outcome_data(
  filename = "1241_GECCO_rectal_combined.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P.value",
  phenotype_col = "outcome"
)

#harmonise and run mr 
#bmi (combined) on crc
bmi_crc_dat <- harmonise_data(bmi_dat, crc_gecco_overall_dat, action = 2)
mr_results_bmi_crc <- mr(bmi_crc_dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_wald_ratio"))
pleiotropy <- mr_pleiotropy_test(bmi_crc_dat)
heterogeneity <- mr_heterogeneity(bmi_crc_dat, method_list=c("mr_egger_regression", "mr_ivw"))

#female bmi on female crc
bmi_fem_crc_dat <- harmonise_data(bmi_fem_dat, crc_gecco_fem_dat, action = 1)
mr_results_bmi_fem_crc <- mr(bmi_fem_crc_dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median",  "mr_weighted_mode", "mr_wald_ratio"))
pleiotropy <- mr_pleiotropy_test(bmi_fem_crc_dat)
heterogeneity <- mr_heterogeneity(bmi_fem_crc_dat, method_list=c("mr_egger_regression", "mr_ivw"))

#male bmi on male crc
bmi_mal_crc_dat <- harmonise_data(bmi_mal_dat, crc_gecco_mal_dat, action = 1)
mr_results_bmi_mal_crc <- mr(bmi_mal_crc_dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median",  "mr_weighted_mode", "mr_wald_ratio"))
pleiotropy <- mr_pleiotropy_test(bmi_mal_crc_dat)
heterogeneity <- mr_heterogeneity(bmi_mal_crc_dat, method_list=c("mr_egger_regression", "mr_ivw"))

#bmi on colon cancer
bmi_colon_dat <- harmonise_data(bmi_dat, crc_gecco_colon_dat, action = 1)
mr_results_bmi_colon <- mr(bmi_colon_dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median",  "mr_weighted_mode", "mr_wald_ratio"))
pleiotropy <- mr_pleiotropy_test(bmi_colon_dat)
heterogeneity <- mr_heterogeneity(bmi_colon_dat, method_list=c("mr_egger_regression", "mr_ivw"))

#bmi on proximal colon cancer
bmi_proximal_dat <- harmonise_data(bmi_dat, crc_gecco_proximal_dat, action = 1)
mr_results_bmi_proximal <- mr(bmi_proximal_dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median",  "mr_weighted_mode", "mr_wald_ratio"))
pleiotropy <- mr_pleiotropy_test(bmi_proximal_dat)
heterogeneity <- mr_heterogeneity(bmi_proximal_dat, method_list=c("mr_egger_regression", "mr_ivw"))

#bmi on distal colon cancer
bmi_distal_dat <- harmonise_data(bmi_dat, crc_gecco_distal_dat, action = 1)
mr_results_bmi_distal <- mr(bmi_distal_dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median",  "mr_weighted_mode", "mr_wald_ratio"))
pleiotropy <- mr_pleiotropy_test(bmi_distal_dat)
heterogeneity <- mr_heterogeneity(bmi_distal_dat, method_list=c("mr_egger_regression", "mr_ivw"))

#bmi on rectal cancer
bmi_rectal_dat <- harmonise_data(bmi_dat, crc_gecco_rectal_dat, action = 1)
mr_results_bmi_rectal <- mr(bmi_rectal_dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median",  "mr_weighted_mode", "mr_wald_ratio"))
pleiotropy <- mr_pleiotropy_test(bmi_rectal_dat)
heterogeneity <- mr_heterogeneity(bmi_rectal_dat, method_list=c("mr_egger_regression", "mr_ivw"))

#whr (combined) on crc
whr_crc_dat <- harmonise_data(whr_dat, crc_gecco_overall_dat, action = 1)
mr_results_whr_crc <- mr(whr_crc_dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median",  "mr_weighted_mode",  "mr_wald_ratio"))
pleiotropy <- mr_pleiotropy_test(whr_crc_dat)
heterogeneity <- mr_heterogeneity(whr_crc_dat, method_list=c("mr_egger_regression", "mr_ivw"))

#whr female on female crc
whr_fem_crc_dat <- harmonise_data(whr_fem_dat, crc_gecco_fem_dat, action = 1)
mr_results_whr_fem_crc <- mr(whr_fem_crc_dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode",  "mr_wald_ratio"))
pleiotropy <- mr_pleiotropy_test(whr_fem_crc_dat)
heterogeneity <- mr_heterogeneity(whr_fem_crc_dat, method_list=c("mr_egger_regression", "mr_ivw"))

#whr male on male crc
whr_mal_crc_dat <- harmonise_data(whr_mal_dat, crc_gecco_mal_dat, action = 1)
mr_results_whr_mal_crc <- mr(whr_mal_crc_dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode",  "mr_wald_ratio"))
pleiotropy <- mr_pleiotropy_test(whr_mal_crc_dat)
heterogeneity <- mr_heterogeneity(whr_mal_crc_dat, method_list=c("mr_egger_regression", "mr_ivw"))

#whr on colon cancer
whr_colon_dat <- harmonise_data(whr_dat, crc_gecco_colon_dat, action = 1)
mr_results_whr_colon <- mr(whr_colon_dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median",  "mr_weighted_mode", "mr_wald_ratio"))
pleiotropy <- mr_pleiotropy_test(whr_colon_dat)
heterogeneity <- mr_heterogeneity(whr_colon_dat, method_list=c("mr_egger_regression", "mr_ivw"))

#whr on proximal colon cancer
whr_proximal_dat <- harmonise_data(whr_dat, crc_gecco_proximal_dat, action = 1)
mr_results_whr_proximal <- mr(whr_proximal_dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median",  "mr_weighted_mode", "mr_wald_ratio"))
pleiotropy <- mr_pleiotropy_test(whr_proximal_dat)
heterogeneity <- mr_heterogeneity(whr_proximal_dat, method_list=c("mr_egger_regression", "mr_ivw"))

#whr on distal colon cancer
whr_distal_dat <- harmonise_data(whr_dat, crc_gecco_distal_dat, action = 1)
mr_results_whr_distal <- mr(whr_distal_dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median",  "mr_weighted_mode", "mr_wald_ratio"))
pleiotropy <- mr_pleiotropy_test(whr_distal_dat)
heterogeneity <- mr_heterogeneity(whr_distal_dat, method_list=c("mr_egger_regression", "mr_ivw"))

#whr on rectal cancer
whr_rectal_dat <- harmonise_data(whr_dat, crc_gecco_rectal_dat, action = 1)
mr_results_whr_rectal <- mr(whr_rectal_dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median",  "mr_weighted_mode", "mr_wald_ratio"))
pleiotropy <- mr_pleiotropy_test(whr_rectal_dat)
heterogeneity <- mr_heterogeneity(whr_rectal_dat, method_list=c("mr_egger_regression", "mr_ivw"))




##analysis of adipose traits on metabolites 
##uses data on snp exposure and snp outcome associations downloaded from https://github.com/lindgrengroup/fatdistnGWAS (adipose traits) and MRBase (metabolites)

#metabolites
kettunen_ids <- ao[ao$pmid%in%"27005778", "id"]

#harmonise and run mr 
#bmi on metabolites
nmr_dat <- extract_outcome_data(bmi_dat$SNP, kettunen_ids$id, proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(bmi_dat, nmr_dat, action = 2)
mr_results_bmi_nmr <- mr(dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_wald_ratio"))
pleiotropy_bmi_nmr <- mr_pleiotropy_test(dat)
heterogeneity_bmi_nmr <- mr_heterogeneity(dat, method_list=c("mr_egger_regression", "mr_ivw"))
list <- list(mr_results_bmi_nmr, pleiotropy_bmi_nmr, heterogeneity_bmi_nmr)

#whr on metabolites
nmr_dat <- extract_outcome_data(whr_dat$SNP, kettunen_ids$id, proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(whr_dat, nmr_dat, action = 1)
mr_results_whr_nmr <- mr(dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_wald_ratio"))
pleiotropy_whr_nmr <- mr_pleiotropy_test(dat)
heterogeneity_whr_nmr <- mr_heterogeneity(dat, method_list=c("mr_egger_regression", "mr_ivw"))
list <- list(mr_results_whr_nmr, pleiotropy_whr_nmr, heterogeneity_whr_nmr)




##analysis of metabolites on colorectal cancer
##snp exposure associations from MRBase (metabolites) 
##snp outcome associations available following an application to the Genetics and Epidemiology of Colorectal Cancer Consortium (GECCO) https://www.fredhutch.org/en/research/divisions/public-health-sciences-division/research/cancer-prevention/genetics-epidemiology-colorectal-cancer-consortium-gecco.html

#metabolite instruments 
nmr_dat <- extract_instruments(kettunen_ids$id)
nmr_dat$t_stat <- (nmr_dat$beta.exposure/nmr_dat$se.exposure)
nmr_dat$f_stat <- (nmr_dat$t_stat)^2
nmr_fstats <- nmr_dat %>%
        group_by(exposure) %>%
        dplyr::summarize(Mean = mean(f_stat, na.rm=TRUE))

#harmonise and run mr 
#metabolites on crc
nmr_crc_dat <- harmonise_data(nmr_dat, crc_dat, action = 1)
mr_results_nmr_crc <- mr(nmr_crc_dat, method_list=c("mr_two_sample_ml", "mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_wald_ratio"))
pleiotropy_nmr_crc <- mr_pleiotropy_test(nmr_crc_dat)

#metabolites on colon cancer
nmr_colon_dat <- harmonise_data(nmr_dat, colon_dat, action = 1)
mr_results_nmr_colon <- mr(nmr_colon_dat, method_list=c("mr_two_sample_ml", "mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_wald_ratio"))
pleiotropy_nmr_colon <- mr_pleiotropy_test(nmr_colon_dat)

#metabolites on proximal colon cancer
nmr_proximal_dat <- harmonise_data(nmr_dat, proximal_dat, action = 1)
mr_results_nmr_proximal <- mr(nmr_proximal_dat, method_list=c("mr_two_sample_ml", "mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_wald_ratio"))
pleiotropy_nmr_proximal <- mr_pleiotropy_test(nmr_proximal_dat)

#metabolites on distal colon cancer
nmr_distal_dat <- harmonise_data(nmr_dat, distal_dat, action = 1)
mr_results_nmr_distal <- mr(nmr_distal_dat, method_list=c("mr_two_sample_ml", "mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_wald_ratio"))
pleiotropy_nmr_distal <- mr_pleiotropy_test(nmr_distal_dat)

#metabolites on rectal cancer
nmr_rectal_dat <- harmonise_data(nmr_dat, rectal_dat, action = 1)
mr_results_nmr_rectal <- mr(nmr_rectal_dat, method_list=c("mr_two_sample_ml", "mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_wald_ratio"))
pleiotropy_nmr_rectal <- mr_pleiotropy_test(nmr_rectal_dat)




##multivariable mr analysis of adipose traits and metabolites on colorectal cancer (overall)
##snp exposure associations from MRBase (adipose traits and metabolites) 
##snp outcome associations available following an application to the Genetics and Epidemiology of Colorectal Cancer Consortium (GECCO) https://www.fredhutch.org/en/research/divisions/public-health-sciences-division/research/cancer-prevention/genetics-epidemiology-colorectal-cancer-consortium-gecco.html

#univariable estimates for the mvmr models
bmi_locke <- extract_instruments(outcomes=2)
outcome_dat <- read_outcome_data(
  snps = bmi_locke$SNP,
  filename = "1241_GECCO_overall_combined.txt",
  sep = "\t",
  snp_col = "rs",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P.value",
  phenotype_col = "outcome"
)
bmi_univar <- mr(dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))

whr_shungin <- extract_instruments(outcomes=73)
outcome_dat <- read_outcome_data(
  snps = whr_shungin$SNP,
  filename = "1241_GECCO_overall_combined.txt",
  sep = "\t",
  snp_col = "rs",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P.value",
  phenotype_col = "outcome"
)
dat <- harmonise_data(
  exposure_dat = whr_shungin, 
  outcome_dat = outcome_dat
)

whr_univar <- mr(dat, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))

#calculate f statistics for the univariable bmi and whr to crc associations
#the f statistic can be calculated as the average of (beta/se), (t statistic) squared  for each snp (beta = effect of snp on exposure, se is the se of that beta)
bmi_locke$t_stat <- (bmi_locke$beta.exposure/bmi_locke$se.exposure)
bmi_locke$f_stat <- (bmi_locke$t_stat)^2
bmi_fstat = mean(bmi_locke$f_stat)
whr_shungin$t_stat <- (whr_shungin$beta.exposure/whr_shungin$se.exposure)
whr_shungin$f_stat <- (whr_shungin$t_stat)^2
whr_fstat = mean(whr_shungin$f_stat)

#estimate the effect of bmi or whr on crc outcomes, adjusted for adipose traits and metabolite clusters
#run for pairs of exposures
#estimate the effect of bmi adjusted for whr on crc (adipose positive control)
id_exposure <- c(2, 73) #(2 = BMI, 73 = WHR)
id_exposure <- c(2, 932) #S.VLDL.TG ("VLDL subclass lipids" block)
id_exposure <- c(2, 867) #IDL.C ("IDL and LDL subclass lipids" block)
id_exposure <- c(2, 904) #M.LDL.C ("IDL and LDL subclass lipids" block)
id_exposure <- c(2, 948) #XL.HDL.TG ("HDL subclass lipids" block)
id_exposure <- c(2, 917) #Other polyunsaturated fatty acids than 18:2 ("Omega-3 and PUFA" block)
id_exposure <- c(2, 893) #18:2, linoleic acid (LA) ("Omega-6 FA" block)
id_exposure <- c(2, 857) #Omega-7, omega-9 and saturated fatty acids ("Monounsaturated and other FAs" block)
id_exposure <- c(2, 859) #Glucose ("Glycemia" block)
id_exposure <- c(2, 849) #Citrate ("Substrates" block)
id_exposure <- c(2, 897) #Leucine ("BCAAs" block
id_exposure <- c(2, 860) #Glutamine ("Other amino acids" block)

#exposure ids for whr adjusted estimates
id_exposure <- c(73, 2) #(2 = BMI, 73 = WHR)
id_exposure <- c(73, 932) #S.VLDL.TG ("VLDL subclass lipids" block)
id_exposure <- c(73, 867) #IDL.C ("IDL and LDL subclass lipids" block)
id_exposure <- c(73, 904) #M.LDL.C ("IDL and LDL subclass lipids" block)
id_exposure <- c(73, 948) #XL.HDL.TG ("HDL subclass lipids" block)
id_exposure <- c(73, 917) #Other polyunsaturated fatty acids than 18:2 ("Omega-3 and PUFA" block)
id_exposure <- c(73, 893) #18:2, linoleic acid (LA) ("Omega-6 FA" block)
id_exposure <- c(73, 857) #Omega-7, omega-9 and saturated fatty acids ("Monounsaturated and other FAs" block)
id_exposure <- c(73, 859) #Glucose ("Glycemia" block)
id_exposure <- c(73, 849) #Citrate ("Substrates" block)
id_exposure <- c(73, 897) #Leucine ("BCAAs" block
id_exposure <- c(73, 860) #Glutamine ("Other amino acids" block)

exposure_dat <- mv_extract_exposures(id_exposure)

# overall crc
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "1241_GECCO_overall_combined.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P.value",
  phenotype_col = "outcome"
)

mvmrdata <- mv_harmonise_data(exposure_dat, outcome_dat)
mv_results <- mv_multiple(mvdat)




##estimating covariance from phenotypic correlation for mvmr dataframes
#code adapted from https://doi.org/10.1101/2020.04.02.021980, "Testing and Correcting for Weak and Pleiotropic Instruments in Two-Sample Multivariable Mendelian Randomisation"

#rho is correlation between x1 and x2. If this is unknown set cov.x1x2 to zero which will impose the assumption of no covariance in the analysis.
#if x1 and x2 are from different datasets then this can be set to zero.

# cov.x1x2 <- rho*mvmrdata$x1.se*mvmrdata$x2.se
cov.x1x2 <- 0

#Estimating the MVMR - using simple weights

mvmr <- lm(mvmrdata$out.beta ~ -1 + mvmrdata$x1.beta + mvmrdata$x2.beta, weights = (mvmrdata$out.se^(-2)))

estbeta1 <- summary(mvmr)$coefficients[1,1]
sebeta1 <- summary(mvmr)$coefficients[1,2]
estbeta2 <- summary(mvmr)$coefficients[2,1]
sebeta2 <- summary(mvmr)$coefficients[2,2]


#calculating the weak instrument test
delta1 <- summary(lm(mvmrdata$x1.beta ~ -1 + mvmrdata$x2.beta))$coefficients[1,1]
v1 <- mvmrdata$x1.se^2 + delta1^2*(mvmrdata$x2.se^2)
F.x1 <- sum((1/v1)*(mvmrdata$x1.beta - (delta1*mvmrdata$x2.beta))^2)/(70-1)

delta2 <- summary(lm(mvmrdata$x2.beta ~ -1 + mvmrdata$x1.beta))$coefficients[1,1]
v2 <- mvmrdata$x2.se^2 + delta2^2*(mvmrdata$x1.se^2)
F.x2 <- sum((1/v2)*(mvmrdata$x2.beta - (delta2*mvmrdata$x1.beta))^2)/(70-1)

#calculating heterogeneity Q statistic for MVMR

seQ <- (mvmrdata$out.se)^2 + (estbeta1^2)*(mvmrdata$x1.se^2) + (estbeta2^2)*(mvmrdata$x2.se^2) + 2*estbeta1*estbeta2*cov.x1x2
Q <- sum(((1/seQ)*((mvmrdata$out.beta-(estbeta1*mvmrdata$x1.beta+estbeta2*mvmrdata$x2.beta))^2)))
critical_value <- qchisq(0.05,70-2,lower.tail = FALSE)
p_value <- pchisq(Q,70-2,lower.tail = FALSE)
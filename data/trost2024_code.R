# Purpose: reproducible code for Trost et al., 2024: analysis of POM characterization across NEON sites

### Load packages----

# install.packages("neonUtilities")
# install.packages("tidyverse")
# install.packages("lubridate")
# install.packages("hydrostats")
# install.packages("reshape2")
# install.packages("gbm")
# install.packages("plotrix")

library(neonUtilities)
library(tidyverse)
library(lubridate)
library(hydrostats)
library(reshape2)
library(gbm)
library(plotrix)

### 1. Download data----

# Download isotope data, NEON data product DP1.20206.001

iso_HOPB_LEWI_POSE_FLNT_list <- loadByProduct(dpID="DP1.20206.001", 
                      site = c("HOPB", "LEWI", "POSE", "FLNT"),
                      tabl = "asi_POMExternalLabDataPerSample",
                      enddate = "2021-12",
                      check.size = F)
iso_HOPB_LEWI_POSE_FLNT <- iso_HOPB_LEWI_POSE_FLNT_list$asi_POMExternalLabDataPerSample

iso_CUPE_GUIL_KING_MCDI_list <- loadByProduct(dpID="DP1.20206.001", 
                     site = c("CUPE", "GUIL", "KING", "MCDI"),
                     tabl = "asi_POMExternalLabDataPerSample",
                     enddate = "2021-12",
                     check.size = F)
iso_CUPE_GUIL_KING_MCDI <- iso_CUPE_GUIL_KING_MCDI_list$asi_POMExternalLabDataPerSample

iso_LECO_WALK_BLWA_MAYF_list <- loadByProduct(dpID="DP1.20206.001", 
                     site = c("LECO", "WALK", "BLWA", "MAYF"),
                     tabl = "asi_POMExternalLabDataPerSample",
                     enddate = "2021-12",
                     check.size = F)
iso_LECO_WALK_BLWA_MAYF <- iso_LECO_WALK_BLWA_MAYF_list$asi_POMExternalLabDataPerSample

iso_TOMB_ARIK_BLUE_PRIN_list <- loadByProduct(dpID="DP1.20206.001", 
                     site = c("TOMB", "ARIK", "BLUE", "PRIN"),
                     tabl = "asi_POMExternalLabDataPerSample",
                     enddate = "2021-12",
                     check.size = F)
iso_TOMB_ARIK_BLUE_PRIN <- iso_TOMB_ARIK_BLUE_PRIN_list$asi_POMExternalLabDataPerSample

iso_BLDE_COMO_WLOU_SYCA_list <- loadByProduct(dpID="DP1.20206.001", 
                     site = c("BLDE", "COMO", "WLOU", "SYCA"),
                     tabl = "asi_POMExternalLabDataPerSample",
                     enddate = "2021-12",
                     check.size = F)
iso_BLDE_COMO_WLOU_SYCA <- iso_BLDE_COMO_WLOU_SYCA_list$asi_POMExternalLabDataPerSample

iso_REDB_MART_MCRA_BIGC_list <- loadByProduct(dpID="DP1.20206.001", 
                     site = c("REDB", "MART", "MCRA", "BIGC"),
                     tabl = "asi_POMExternalLabDataPerSample",
                     enddate = "2021-12",
                     check.size = F)
iso_REDB_MART_MCRA_BIGC <- iso_REDB_MART_MCRA_BIGC_list$asi_POMExternalLabDataPerSample

iso_TECR_OKSR_CARI_list <- loadByProduct(dpID="DP1.20206.001", 
                     site = c("TECR", "OKSR", "CARI"),
                     tabl = "asi_POMExternalLabDataPerSample",
                     enddate = "2021-12",
                     check.size = F)
iso_TECR_OKSR_CARI <- iso_TECR_OKSR_CARI_list$asi_POMExternalLabDataPerSample

iso_data_total <- rbind(iso_HOPB_LEWI_POSE_FLNT,
                        iso_CUPE_GUIL_KING_MCDI,
                        iso_LECO_WALK_BLWA_MAYF,
                        iso_TOMB_ARIK_BLUE_PRIN,
                        iso_BLDE_COMO_WLOU_SYCA,
                        iso_REDB_MART_MCRA_BIGC,
                        iso_TECR_OKSR_CARI)
 
iso_error_list <- loadByProduct(dpID="DP1.20206.001", 
                     site = "POSE", #site with the complete date range
                     package = "expanded",
                     enddate = "2021-12",
                     check.size = F)
iso_error <- iso_error_list$asi_externalLabPOMSummaryData
  

# Download grab data, NEON data product DP1.20093.001

grab_HOPB_LEWI_POSE_FLNT_list <- loadByProduct(dpID="DP1.20093.001", 
                                         site = c("HOPB", "LEWI", "POSE", "FLNT"),
                                         tabl = "swc_externalLabDataByAnalyte",
                                         enddate = "2022-06",
                                         check.size = F)
grab_HOPB_LEWI_POSE_FLNT <- grab_HOPB_LEWI_POSE_FLNT_list$swc_externalLabDataByAnalyte

grab_CUPE_GUIL_KING_MCDI_list <- loadByProduct(dpID="DP1.20093.001", 
                                         site = c("CUPE", "GUIL", "KING", "MCDI"),
                                         tabl = "swc_externalLabDataByAnalyte",
                                         enddate = "2022-06",
                                         check.size = F)
grab_CUPE_GUIL_KING_MCDI <- grab_CUPE_GUIL_KING_MCDI_list$swc_externalLabDataByAnalyte

grab_LECO_WALK_BLWA_MAYF_list <- loadByProduct(dpID="DP1.20093.001", 
                                         site = c("LECO", "WALK", "BLWA", "MAYF"),
                                         tabl = "swc_externalLabDataByAnalyte",
                                         enddate = "2022-06",
                                         check.size = F)
grab_LECO_WALK_BLWA_MAYF <- grab_LECO_WALK_BLWA_MAYF_list$swc_externalLabDataByAnalyte

grab_TOMB_ARIK_BLUE_PRIN_list <- loadByProduct(dpID="DP1.20093.001", 
                                         site = c("TOMB", "ARIK", "BLUE", "PRIN"),
                                         tabl = "swc_externalLabDataByAnalyte",
                                         enddate = "2022-06",
                                         check.size = F)
grab_TOMB_ARIK_BLUE_PRIN <- grab_TOMB_ARIK_BLUE_PRIN_list$swc_externalLabDataByAnalyte

grab_BLDE_COMO_WLOU_SYCA_list <- loadByProduct(dpID="DP1.20093.001", 
                                         site = c("BLDE", "COMO", "WLOU", "SYCA"),
                                         tabl = "swc_externalLabDataByAnalyte",
                                         enddate = "2022-06",
                                         check.size = F)
grab_BLDE_COMO_WLOU_SYCA <- grab_BLDE_COMO_WLOU_SYCA_list$swc_externalLabDataByAnalyte


grab_REDB_MART_MCRA_BIGC_list <- loadByProduct(dpID="DP1.20093.001", 
                                         site = c("REDB", "MART", "MCRA", "BIGC"),
                                         tabl = "swc_externalLabDataByAnalyte",
                                         enddate = "2022-06",
                                         check.size = F)
grab_REDB_MART_MCRA_BIGC <- grab_REDB_MART_MCRA_BIGC_list$swc_externalLabDataByAnalyte

grab_TECR_OKSR_CARI_list <- loadByProduct(dpID="DP1.20093.001", 
                                    site = c("TECR", "OKSR", "CARI"),
                                    tabl = "swc_externalLabDataByAnalyte",
                                    enddate = "2022-06",
                                    check.size = F)
grab_TECR_OKSR_CARI <- grab_TECR_OKSR_CARI_list$swc_externalLabDataByAnalyte

grab_data_total <- rbind(grab_HOPB_LEWI_POSE_FLNT,
                        grab_CUPE_GUIL_KING_MCDI,
                        grab_LECO_WALK_BLWA_MAYF,
                        grab_TOMB_ARIK_BLUE_PRIN,
                        grab_BLDE_COMO_WLOU_SYCA,
                        grab_REDB_MART_MCRA_BIGC,
                        grab_TECR_OKSR_CARI)

grab_error_list <- loadByProduct(dpID="DP1.20093.001", 
                                site = "MAYF", #site with the complete date range
                                package = "expanded",
                                enddate = "2022-06",
                                check.size = F)
grab_error <- grab_error_list$swc_externalLabSummaryData

# Download discharge data, NEON data product DP4.00130.001

discharge_HOPB_LEWI_POSE_FLNT_list <- loadByProduct(dpID="DP4.00130.001", 
                                          site = c("HOPB", "LEWI", "POSE", "FLNT"),
                                          enddate = "2022-10",
                                          check.size=F)
discharge_HOPB_LEWI_POSE_FLNT <- discharge_HOPB_LEWI_POSE_FLNT_list$csd_continuousDischarge

discharge_CUPE_GUIL_KING_MCDI_list <- loadByProduct(dpID="DP4.00130.001", 
                                                    site = c("CUPE", "GUIL", "KING", "MCDI"), 
                                                    enddate = "2022-10",
                                                    check.size=F)
discharge_CUPE_GUIL_KING_MCDI <- discharge_CUPE_GUIL_KING_MCDI_list$csd_continuousDischarge

discharge_LECO_WALK_BLWA_MAYF_list <- loadByProduct(dpID="DP4.00130.001", 
                                                    site = c("LECO", "WALK", "BLWA", "MAYF"), 
                                                    enddate = "2022-10",
                                                    check.size=F)
discharge_LECO_WALK_BLWA_MAYF <- discharge_LECO_WALK_BLWA_MAYF_list$csd_continuousDischarge

discharge_ARIK_BLUE_PRIN_list <- loadByProduct(dpID="DP4.00130.001", 
                                                    site = c("ARIK", "BLUE", "PRIN"),
                                                    enddate = "2022-10",
                                                    check.size=F)
discharge_ARIK_BLUE_PRIN <- discharge_ARIK_BLUE_PRIN_list$csd_continuousDischarge

discharge_TOMB_list <- loadByProduct(dpID="DP4.00130.001", 
                                               site = "TOMB", 
                                               enddate = "2022-10",
                                               check.size=F)
discharge_TOMB <- discharge_TOMB_list$csd_continuousDischargeUSGS

discharge_BLDE_COMO_WLOU_SYCA_list <- loadByProduct(dpID="DP4.00130.001", 
                                               site = c("BLDE", "COMO", "WLOU", "SYCA"),
                                               enddate = "2022-10",
                                               check.size=F)
discharge_BLDE_COMO_WLOU_SYCA <- discharge_BLDE_COMO_WLOU_SYCA_list$csd_continuousDischarge

discharge_REDB_MART_MCRA_BIGC_list <- loadByProduct(dpID="DP4.00130.001", 
                                                    site = c("REDB", "MART", "MCRA", "BIGC"),
                                                    enddate = "2022-10",
                                                    check.size=F)
discharge_REDB_MART_MCRA_BIGC <- discharge_REDB_MART_MCRA_BIGC_list$csd_continuousDischarge

discharge_TECR_OKSR_CARI_list <- loadByProduct(dpID="DP4.00130.001", 
                                                    site = c("TECR", "OKSR", "CARI"),
                                                    enddate = "2022-10",
                                                    check.size=F)
discharge_TECR_OKSR_CARI <- discharge_TECR_OKSR_CARI_list$csd_continuousDischarge

df_list <- list(discharge_HOPB_LEWI_POSE_FLNT, 
                discharge_CUPE_GUIL_KING_MCDI, 
                discharge_LECO_WALK_BLWA_MAYF, 
                discharge_ARIK_BLUE_PRIN, 
                discharge_BLDE_COMO_WLOU_SYCA, 
                discharge_REDB_MART_MCRA_BIGC, 
                discharge_TECR_OKSR_CARI)
df_list_short <- lapply(df_list, function(x) select(x, c(siteID, endDate, maxpostDischarge))) 
q_neon <- df_list_short %>% 
  bind_rows() %>%  
  rename(Q_Ls = maxpostDischarge)

q_usgs <- discharge_TOMB %>% 
  select(siteID, endDate, usgsDischarge) %>% 
  rename(Q_Ls = usgsDischarge)

q_long <- rbind(q_neon, q_usgs)
q_long <- rename(q_long, 
                 collectDatetime = endDate)

### 2. Clean data ----

### a. Remove flagged data points

## Grab

grab_noflag_clean <- grab_data_total %>% 
filter(analyte == "DOC" | 
         analyte == "NO3+NO2 - N" | 
         analyte == "TDN" |
         analyte == "TP" | 
         analyte == "TDP") %>% 
  filter(((sampleCondition == "GOOD" |
            sampleCondition == "OK") &
           is.na(remarks)) |
  (sampleCondition == "Other" &
     (remarks == "File replces original, which had wrong year" |
        remarks == "Low NH4 method and Low NO2+NO3 method" |
        remarks == "Low NO2+NO3 method")))

## Isotopes

iso_noflag_clean <- iso_data_total %>% 
  filter(!(sampleCondition == "other")) %>% 
  filter(is.na(externalRemarks) |
           (externalRemarks == "Below calibration range" |
              externalRemarks == "Below detection limit"))

### b. Set concentrations below MDL to  value of MDL / sqrt(2)

## Grab

# Isolating data NEON has flagged as below MDL with a tag of -1 (will be cleaned below)

belowdet_qf <- grab_noflag_clean %>% 
  filter(!is.na(belowDetectionQF))
belowdet_qf$analyteConcentration <- -1
wo_belowdet_qf <- anti_join(grab_noflag_clean, belowdet_qf, by = ("uid"))
grab_noflags_clean <- rbind(belowdet_qf, wo_belowdet_qf) 

# Adding a column of data specific method detection limits (MDLs) 
# found in swc_externalLabSummaryData data frame available in the 
# “extended” package download for DP1.20093.001

grab_DOC_mdl_1 <- grab_noflags_clean %>% 
  filter(analyte == "DOC" &
           (collectDate < "2016-06-15")) %>% 
  mutate(mdl = 0.1) 
grab_DOC_mdl_2 <- grab_noflags_clean %>% 
  filter(analyte == "DOC" &
           (collectDate >= "2016-06-15" &
              collectDate < "2019-05-28")) %>% 
  mutate(mdl = 0.1)
grab_DOC_mdl_3 <- grab_noflags_clean %>% 
  filter(analyte == "DOC" &
           (collectDate >= "2019-05-28" &
              collectDate < "2021-04-15")) %>% 
  mutate(mdl = 0.097)
grab_DOC_mdl_4 <- grab_noflags_clean %>% 
  filter(analyte == "DOC" &
           (collectDate >= "2021-04-15" &
              collectDate < "2022-02-24")) %>% 
  mutate(mdl = 0.097)
grab_DOC_mdl_5 <- grab_noflags_clean %>% 
  filter(analyte == "DOC" &
           (collectDate >= "2022-02-24")) %>% 
  mutate(mdl = 0.099)


grab_nitrate_mdl_1 <- grab_noflags_clean %>% 
  filter(analyte == "NO3+NO2 - N" &
           (collectDate < "2016-05-02")) %>% 
  mutate(mdl = 0.05)
grab_nitrate_mdl_2 <- grab_noflags_clean %>% 
  filter(analyte == "NO3+NO2 - N" &
           (collectDate >= "2016-05-02" &
              collectDate < "2019-07-03")) %>% 
  mutate(mdl = 0.05)
grab_nitrate_mdl_3 <- grab_noflags_clean %>% 
  filter(analyte == "NO3+NO2 - N" &
           (collectDate >= "2019-07-03" &
              collectDate < "2020-02-05")) %>% 
  mutate(mdl = 0)
grab_nitrate_mdl_4 <- grab_noflags_clean %>% 
  filter(analyte == "NO3+NO2 - N" &
           (collectDate >= "2020-02-05" &
              collectDate < "2021-03-31")) %>% 
  mutate(mdl = 0)
grab_nitrate_mdl_5 <- grab_noflags_clean %>% 
  filter(analyte == "NO3+NO2 - N" &
           (collectDate >= "2021-03-31")) %>% 
  mutate(mdl = 0.001)


grab_TDN_mdl_1 <- grab_noflags_clean %>% 
  filter(analyte == "TDN" &
           (collectDate < "2016-06-15")) %>% 
  mutate(mdl = 0.1)
grab_TDN_mdl_2 <- grab_noflags_clean %>% 
  filter(analyte == "TDN" &
           (collectDate >= "2016-06-15" &
              collectDate < "2019-05-28")) %>% 
  mutate(mdl = 0.1)
grab_TDN_mdl_3 <- grab_noflags_clean %>% 
  filter(analyte == "TDN" &
           (collectDate >= "2019-05-28" &
              collectDate < "2021-04-15")) %>% 
  mutate(mdl = 0.035)
grab_TDN_mdl_4 <- grab_noflags_clean %>% 
  filter(analyte == "TDN" &
           (collectDate >= "2021-04-15" &
              collectDate < "2022-02-24")) %>% 
  mutate(mdl = 0.035)
grab_TDN_mdl_5 <- grab_noflags_clean %>% 
  filter(analyte == "TDN" &
           (collectDate >= "2022-02-24")) %>% 
  mutate(mdl = 0.028)

grab_TDP_mdl_1 <- grab_noflags_clean %>% 
  filter(analyte == "TDP" &
           (collectDate < "2020-03-20")) %>% 
  mutate(mdl = 0.001)
grab_TDP_mdl_2 <- grab_noflags_clean %>% 
  filter(analyte == "TDP" &
           (collectDate >= "2020-03-20")) %>% 
  mutate(mdl = 0.001)

grab_TP_mdl_1 <- grab_noflags_clean %>% 
  filter(analyte == "TP" &
           (collectDate < "2020-03-20")) %>% 
  mutate(mdl = 0.001)
grab_TP_mdl_2 <- grab_noflags_clean %>% 
  filter(analyte == "TP" &
           (collectDate >= "2020-03-20")) %>% 
  mutate(mdl = 0.001)

grab_mdl_bound <- rbind(grab_DOC_mdl_1,
                        grab_DOC_mdl_2,
                        grab_DOC_mdl_3,
                        grab_DOC_mdl_4,
                        grab_DOC_mdl_5,
                        grab_nitrate_mdl_1,
                        grab_nitrate_mdl_2,
                        grab_nitrate_mdl_3,
                        grab_nitrate_mdl_4,
                        grab_nitrate_mdl_5,
                        grab_TDN_mdl_1,
                        grab_TDN_mdl_2,
                        grab_TDN_mdl_3,
                        grab_TDN_mdl_4,
                        grab_TDN_mdl_5,
                        grab_TDP_mdl_1,
                        grab_TDP_mdl_2,
                        grab_TP_mdl_1,
                        grab_TP_mdl_2)

# Replace concentrations below MDL with value of MDL / sqrt(2)

grab_belowdet <- grab_mdl_bound %>% 
  filter(analyteConcentration <= mdl) %>% 
  mutate(below_det = "YES")
grab_belowdet$analyteConcentration <- grab_belowdet$mdl/sqrt(2) 

grab_abovedet <- grab_mdl_bound %>% 
  filter(analyteConcentration > mdl) %>% 
  mutate(below_det = "NO")

grab_mdl_clean <- rbind(grab_belowdet, grab_abovedet) 

## Isotopes

# There are no listed MDLs for this data package, therefore we chose to remove 
# both values NEON flagged as below MDL and negative values

iso_mdl <- iso_noflag_clean %>% 
  filter((externalRemarks == "Below calibration range" |
     externalRemarks == "Below detection limit") |
       ((analyte == "nitrogen" | 
          analyte == "carbon") &
          analyteConcentration < 0)) 

iso_mdl_clean <- anti_join(iso_noflag_clean, iso_mdl, by = ("uid")) %>% 
  mutate(below_det = NA)

# Remove NA values (Grab data package has no NA values after cleaning)

iso_clean_nona <- iso_mdl_clean %>% 
  filter(!is.na(analyteConcentration)) 

## Extra step for Isotope data package to add times to sample collection dates 

iso_fieldData <- loadByProduct(dpID="DP1.20206.001", 
                              site = c("HOPB", "LEWI", "POSE", "FLNT",
                                       "CUPE", "GUIL", "KING", "MCDI",
                                       "LECO", "WALK", "BLWA", "MAYF",
                                       "TOMB", "ARIK", "BLUE", "PRIN",
                                       "BLDE", "COMO", "WLOU", "SYCA",
                                       "REDB", "MART", "MCRA", "BIGC",
                                       "TECR", "OKSR", "CARI"),
                              tabl = "asi_fieldData",
                              check.size = F, 
                              enddate = "2022-01") 

iso_field_data <- iso_fieldData$asi_fieldData

iso_field_select <- select(iso_field_data, collectDate, isotopePOMRep2SampleID,
                           isotopePOMSampleID)
iso_field_id <- iso_field_select %>% 
  pivot_longer(!collectDate, values_to = "sampleID") %>% 
  rename(collectDatetime = collectDate) %>% 
  select(-name) %>% 
  filter(!is.na(sampleID)) 

iso_dateTime <- merge(iso_clean_nona, iso_field_id)

### 3. Wrangle data into long and wide formats----

## Long format 

# Select relevant columns and change column names to match between packages 

iso_dateTime_long <- iso_dateTime %>% 
  select(-collectDate) %>% 
  rename(analyteUnits = plantAlgaeLabUnits, 
         collectDate = collectDatetime)
iso_dateTime_long$mdl <- NA
iso_dateTime_long$NEONpackage <- "DP1.20206.001"

iso_dateTime_long <- select(iso_dateTime_long, 
                            sampleID, domainID, siteID,  collectDate, 
                            analyte, analyteConcentration, analyteUnits, 
                            mdl, below_det, NEONpackage)

iso_dateTime_long$analyte[iso_dateTime_long$analyte == "carbon"] <- "TPC"
iso_dateTime_long$analyte[iso_dateTime_long$analyte == "nitrogen"] <- "TPN"

grab_mdl_clean_long <- grab_mdl_clean

grab_mdl_clean_long$NEONpackage <- "DP1.20093.001"
grab_mdl_clean_long <- select(grab_mdl_clean_long, 
                              sampleID, domainID, siteID,  collectDate, 
                              analyte, analyteConcentration, analyteUnits, 
                              mdl, below_det, NEONpackage)

NEON_chem_long <- rbind(iso_dateTime_long, grab_mdl_clean_long)

## Wide format

wide_grab <- pivot_wider(grab_mdl_clean,
                             id_cols = c("siteID", "collectDate", "domainID"),
                             names_from = analyte, 
                             values_from = analyteConcentration,
                             values_fn = mean) 

wide_iso <- pivot_wider(iso_dateTime,
                             id_cols = c("siteID", "collectDatetime", "domainID"),
                        names_from = analyte, 
                        values_from = analyteConcentration,
                        values_fn = mean) 

wide_merge <- merge(wide_grab, wide_iso, by.x = c("collectDate", "siteID", "domainID"), by.y = c("collectDatetime", "siteID", "domainID"), all = TRUE) 

# Tidy names and convert concentrations to mg/L

wide_merge_tidy <- wide_merge %>% 
  mutate(TPN_mgL = nitrogen/10^3, #in isotope data package, TPN and TPC recorded in micrograms/L
         TPC_mgL = carbon/10^3,
         TPP_mgL = TP-TDP, na.rm = T,
         month = format(collectDate, "%m"),
         year = format(collectDate, "%Y")) %>% 
  rename(DOC_mgL = DOC,
         TDN_mgL = TDN, 
         nitrate_mgL = "NO3+NO2 - N",
         TDP_mgL = TDP, 
         TP_mgL = TP,
         d13C_perMill = d13C, 
         d15N_perMill = d15N) %>% 
  select(siteID, domainID, collectDate, year, month, 
         DOC_mgL, TDN_mgL, nitrate_mgL, TDP_mgL, TPC_mgL, TPN_mgL, TPP_mgL, TP_mgL, d13C_perMill, d15N_perMill)

# Remove meaningless TPP values

wide_merge_tidy$TPP_mgL[wide_merge_tidy$TPP_mgL <= 0 & !is.na(wide_merge_tidy$TPP_mgL)] <- NA

# Calculate molar conversions and nutrient ratios. Unit conversion formula: 1 ug/L = 1/MW*umol/L

wide_merge_calcs <- wide_merge_tidy %>% 
  mutate(DOC_to_TDN_mgL = DOC_mgL/TDN_mgL,
         DOC_to_nitrate_mgL = DOC_mgL/nitrate_mgL,
         DOC_to_TDP_mgL = DOC_mgL/TDP_mgL, 
         TDN_to_TDP_mgL = TDN_mgL/TDP_mgL,
         TPC_to_TPN_mgL = TPC_mgL/TPN_mgL,
         TPC_to_TPP_mgL = TPC_mgL/TPP_mgL, 
         TPN_to_TPP_mgL = TPN_mgL/TPP_mgL,
         DOC_mmolL = DOC_mgL*(1/12.011),
         TDN_mmolL = TDN_mgL*(1/14.007),
         TDP_mmolL = TDP_mgL*(1/30.974),
         TPC_mmolL = TPC_mgL*(1/12.011),
         TPN_mmolL = TPN_mgL*(1/14.007),
         TPP_mmolL = TPP_mgL*(1/30.974),
         DOC_to_TDN_mw = DOC_mmolL/TDN_mmolL,
         DOC_to_TDP_mw = DOC_mmolL/TDP_mmolL,
         TDN_to_TDP_mw = TDN_mmolL/TDP_mmolL,
         TPC_to_TPN_mw = TPC_mmolL/TPN_mmolL,
         TPC_to_TPP_mw = TPC_mmolL/TPP_mmolL,
         TPN_to_TPP_mw = TPN_mmolL/TPP_mmolL)

### 4. Summary stats for chemistry data----

analyte_sum <- wide_merge_calcs %>% 
  dplyr::group_by(siteID) %>% 
  dplyr::summarise(mean_DOC_mgL = mean(DOC_mgL, na.rm = TRUE),
                   error_DOC_mgL = std.error(DOC_mgL, na.rm= T), 
                   mean_TDN_mgL = mean(TDN_mgL, na.rm = TRUE),
                   error_TDN_mgL = std.error(TDN_mgL, na.rm= T),
                   mean_TDP_mgL = mean(TDP_mgL, na.rm = TRUE),
                   error_TDP_mgL = std.error(TDP_mgL, na.rm= T),
                   mean_TPC_mgL = mean(TPC_mgL, na.rm = TRUE),
                   error_TPC_mgL = std.error(TPC_mgL, na.rm= T),
                   mean_TPN_mgL = mean(TPN_mgL, na.rm = TRUE),
                   error_TPN_mgL = std.error(TPN_mgL, na.rm= T),
                   mean_TPP_mgL = mean(TPP_mgL, na.rm = TRUE),
                   error_TPP_mgL = std.error(TPP_mgL, na.rm= T),
                   mean_DOC_mmolL = mean(DOC_mmolL, na.rm = TRUE),
                   error_DOC_mmolL = std.error(DOC_mmolL, na.rm= T), 
                   mean_TDN_mmolL = mean(TDN_mmolL, na.rm = TRUE),
                   error_TDN_mmolL = std.error(TDN_mmolL, na.rm= T), 
                   mean_TDP_mmolL = mean(TDP_mmolL, na.rm = TRUE),
                   error_TDP_mmolL = std.error(TDP_mmolL, na.rm= T), 
                   mean_TPC_mmolL = mean(TPC_mmolL, na.rm = TRUE),
                   error_TPC_mmolL = std.error(TPC_mmolL, na.rm= T),
                   mean_TPN_mmolL = mean(TPN_mmolL, na.rm = TRUE),
                   error_TPN_mmolL = std.error(TPN_mmolL, na.rm= T),
                   mean_TPP_mmolL = mean(TPP_mmolL, na.rm = TRUE),
                   error_TPP_mmolL = std.error(TPP_mmolL, na.rm= T),
                   mean_DOC_to_TDN_mw = mean(DOC_to_TDN_mw, na.rm = TRUE),
                   error_DOC_to_TDN_mw = std.error(DOC_to_TDN_mw, na.rm= T),
                   mean_DOC_to_TDP_mw = mean(DOC_to_TDP_mw, na.rm = TRUE),
                   error_DOC_to_TDP_mw = std.error(DOC_to_TDP_mw, na.rm= T),
                   mean_TDN_to_TDP_mw = mean(TDN_to_TDP_mw, na.rm = TRUE),
                   error_TDN_to_TDP_mw = std.error(TDN_to_TDP_mw, na.rm= T),
                   mean_TPC_to_TPN_mw = mean(TPC_to_TPN_mw, na.rm = TRUE),
                   error_TPC_to_TPN_mw = std.error(TPC_to_TPN_mw, na.rm= T),
                   mean_TPC_to_TPP_mw = mean(TPC_to_TPP_mw, na.rm = TRUE), 
                   error_TPC_to_TPP_mw = std.error(TPC_to_TPP_mw, na.rm= T), 
                   mean_TPN_to_TPP_mw = mean(TPN_to_TPP_mw, na.rm = TRUE), 
                   error_TPN_to_TPP_mw = std.error(TPN_to_TPP_mw, na.rm= T),
                   mean_d13C_perMill = mean(d13C_perMill, na.rm = TRUE), 
                   error_d13C_perMill = std.error(d13C_perMill, na.rm= T),
                   mean_d15N_perMill = mean(d15N_perMill, na.rm = TRUE),
                   error_d15N_perMill = std.error(d15N_perMill, na.rm= T))

### 5. Assemble predictor variable data frames----

## Calculate seasonality index with "hydrostats" package

HOPB_seas_index <- discharge_HOPB_LEWI_POSE_FLNT  %>% 
  filter(siteID == "HOPB") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
HOPB_seas_index_num <- seasonality(HOPB_seas_index)
HOPB_seas_index_df <- data.frame(siteID = "HOPB", seasonality = HOPB_seas_index_num)

LEWI_seas_index <- discharge_HOPB_LEWI_POSE_FLNT  %>% 
  filter(siteID == "LEWI") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
LEWI_seas_index_num <- seasonality(LEWI_seas_index)
LEWI_seas_index_df <- data.frame(siteID = "LEWI", seasonality = LEWI_seas_index_num)

POSE_seas_index <- discharge_HOPB_LEWI_POSE_FLNT  %>% 
  filter(siteID == "POSE") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
POSE_seas_index_num <- seasonality(POSE_seas_index)
POSE_seas_index_df <- data.frame(siteID = "POSE", seasonality = POSE_seas_index_num)

FLNT_seas_index <- discharge_HOPB_LEWI_POSE_FLNT  %>% 
  filter(siteID == "FLNT") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
FLNT_seas_index_num <- seasonality(FLNT_seas_index)
FLNT_seas_index_df <- data.frame(siteID = "FLNT", seasonality = FLNT_seas_index_num)

CUPE_seas_index <- discharge_CUPE_GUIL_KING_MCDI  %>% 
  filter(siteID == "CUPE") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
CUPE_seas_index_num <- seasonality(CUPE_seas_index)
CUPE_seas_index_df <- data.frame(siteID = "CUPE", seasonality = CUPE_seas_index_num)

GUIL_seas_index <- discharge_CUPE_GUIL_KING_MCDI  %>% 
  filter(siteID == "GUIL") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
GUIL_seas_index_num <- seasonality(GUIL_seas_index)
GUIL_seas_index_df <- data.frame(siteID = "GUIL", seasonality = GUIL_seas_index_num)

KING_seas_index <- discharge_CUPE_GUIL_KING_MCDI  %>% 
  filter(siteID == "KING") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
KING_seas_index_num <- seasonality(KING_seas_index)
KING_seas_index_df <- data.frame(siteID = "KING", seasonality = KING_seas_index_num)

MCDI_seas_index <- discharge_CUPE_GUIL_KING_MCDI  %>% 
  filter(siteID == "MCDI") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
MCDI_seas_index_num <- seasonality(MCDI_seas_index)
MCDI_seas_index_df <- data.frame(siteID = "MCDI", seasonality = MCDI_seas_index_num)

LECO_seas_index <- discharge_LECO_WALK_BLWA_MAYF  %>% 
  filter(siteID == "LECO") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
LECO_seas_index_num <- seasonality(LECO_seas_index)
LECO_seas_index_df <- data.frame(siteID = "LECO", seasonality = LECO_seas_index_num)

WALK_seas_index <- discharge_LECO_WALK_BLWA_MAYF  %>% 
  filter(siteID == "WALK") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
WALK_seas_index_num <- seasonality(WALK_seas_index)
WALK_seas_index_df <- data.frame(siteID = "WALK", seasonality = WALK_seas_index_num)

BLWA_seas_index <- discharge_LECO_WALK_BLWA_MAYF  %>% 
  filter(siteID == "BLWA") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
BLWA_seas_index_num <- seasonality(BLWA_seas_index)
BLWA_seas_index_df <- data.frame(siteID = "BLWA", seasonality = BLWA_seas_index_num)

MAYF_seas_index <- discharge_LECO_WALK_BLWA_MAYF  %>% 
  filter(siteID == "MAYF") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
MAYF_seas_index_num <- seasonality(MAYF_seas_index)
MAYF_seas_index_df <- data.frame(siteID = "MAYF", seasonality = MAYF_seas_index_num)

TOMB_seas_index <- discharge_TOMB  %>% 
  filter(siteID == "TOMB") %>% 
  rename(Date = endDate, 
         Q = usgsDischarge) %>% 
  select(Date, Q)
TOMB_seas_index_num <- seasonality(TOMB_seas_index)
TOMB_seas_index_df <- data.frame(siteID = "TOMB", seasonality = TOMB_seas_index_num)

ARIK_seas_index <- discharge_ARIK_BLUE_PRIN  %>% 
  filter(siteID == "ARIK") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
ARIK_seas_index_num <- seasonality(ARIK_seas_index)
ARIK_seas_index_df <- data.frame(siteID = "ARIK", seasonality = ARIK_seas_index_num)

BLUE_seas_index <- discharge_ARIK_BLUE_PRIN  %>% 
  filter(siteID == "BLUE") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
BLUE_seas_index_num <- seasonality(BLUE_seas_index)
BLUE_seas_index_df <- data.frame(siteID = "BLUE", seasonality = BLUE_seas_index_num)

PRIN_seas_index <- discharge_ARIK_BLUE_PRIN  %>% 
  filter(siteID == "PRIN") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
PRIN_seas_index_num <- seasonality(PRIN_seas_index)
PRIN_seas_index_df <- data.frame(siteID = "PRIN", seasonality = PRIN_seas_index_num)

BLDE_seas_index <- discharge_BLDE_COMO_WLOU_SYCA  %>% 
  filter(siteID == "BLDE") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
BLDE_seas_index_num <- seasonality(BLDE_seas_index)
BLDE_seas_index_df <- data.frame(siteID = "BLDE", seasonality = BLDE_seas_index_num)

COMO_seas_index <- discharge_BLDE_COMO_WLOU_SYCA  %>% 
  filter(siteID == "COMO") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
COMO_seas_index_num <- seasonality(COMO_seas_index)
COMO_seas_index_df <- data.frame(siteID = "COMO", seasonality = COMO_seas_index_num)

WLOU_seas_index <- discharge_BLDE_COMO_WLOU_SYCA  %>% 
  filter(siteID == "WLOU") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
WLOU_seas_index_num <- seasonality(WLOU_seas_index)
WLOU_seas_index_df <- data.frame(siteID = "WLOU", seasonality = WLOU_seas_index_num)

SYCA_seas_index <- discharge_BLDE_COMO_WLOU_SYCA  %>% 
  filter(siteID == "SYCA") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
SYCA_seas_index_num <- seasonality(SYCA_seas_index)
SYCA_seas_index_df <- data.frame(siteID = "SYCA", seasonality = SYCA_seas_index_num)

REDB_seas_index <- discharge_REDB_MART_MCRA_BIGC  %>% 
  filter(siteID == "REDB") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
REDB_seas_index_num <- seasonality(REDB_seas_index)
REDB_seas_index_df <- data.frame(siteID = "REDB", seasonality = REDB_seas_index_num)

MART_seas_index <- discharge_REDB_MART_MCRA_BIGC  %>% 
  filter(siteID == "MART") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
MART_seas_index_num <- seasonality(MART_seas_index)
MART_seas_index_df <- data.frame(siteID = "MART", seasonality = MART_seas_index_num)

MCRA_seas_index <- discharge_REDB_MART_MCRA_BIGC  %>% 
  filter(siteID == "MCRA") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
MCRA_seas_index_num <- seasonality(MCRA_seas_index)
MCRA_seas_index_df <- data.frame(siteID = "MCRA", seasonality = MCRA_seas_index_num)

BIGC_seas_index <- discharge_REDB_MART_MCRA_BIGC  %>% 
  filter(siteID == "BIGC") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
BIGC_seas_index_num <- seasonality(BIGC_seas_index)
BIGC_seas_index_df <- data.frame(siteID = "BIGC", seasonality = BIGC_seas_index_num)

TECR_seas_index <- discharge_TECR_OKSR_CARI  %>% 
  filter(siteID == "TECR") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
TECR_seas_index_num <- seasonality(TECR_seas_index)
TECR_seas_index_df <- data.frame(siteID = "TECR", seasonality = TECR_seas_index_num)

OKSR_seas_index <- discharge_TECR_OKSR_CARI  %>% 
  filter(siteID == "OKSR") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
OKSR_seas_index_num <- seasonality(OKSR_seas_index)
OKSR_seas_index_df <- data.frame(siteID = "OKSR", seasonality = OKSR_seas_index_num)

CARI_seas_index <- discharge_TECR_OKSR_CARI  %>% 
  filter(siteID == "CARI") %>% 
  rename(Date = endDate, 
         Q = maxpostDischarge) %>% 
  select(Date, Q)
CARI_seas_index_num <- seasonality(CARI_seas_index)
CARI_seas_index_df <- data.frame(siteID = "CARI", seasonality = CARI_seas_index_num)

seas_index_total <- rbind(HOPB_seas_index_df, LEWI_seas_index_df, POSE_seas_index_df, FLNT_seas_index_df,
                          CUPE_seas_index_df, GUIL_seas_index_df, KING_seas_index_df, MCDI_seas_index_df,
                          LECO_seas_index_df, WALK_seas_index_df, BLWA_seas_index_df, MAYF_seas_index_df,
                          TOMB_seas_index_df, ARIK_seas_index_df, BLUE_seas_index_df, PRIN_seas_index_df,
                          BLDE_seas_index_df, COMO_seas_index_df, WLOU_seas_index_df, SYCA_seas_index_df,
                          REDB_seas_index_df, MART_seas_index_df, MCRA_seas_index_df, BIGC_seas_index_df,
                          TECR_seas_index_df, OKSR_seas_index_df, CARI_seas_index_df)

## Wrangle data for watershed area and NLCD land use percentages
# (NEON stores watershed land cover data here: https://www.neonscience.org/data-samples/data/spatial-data-maps 
# in .csv files within the "Aquatic Sites Watersheds" section download folder)

# Watershed area. NEON file name is "NEON_Aquatic_Watershed.csv"

NEON_site_info <- read_csv("NEON_Aquatic_Watershed.csv") # must create a directory a directory

NEON_stream_sites <- c("HOPB", "LEWI", "POSE", "BLWA", "FLNT", "MAYF", "PRIN", "TOMB",
                       "CUPE", "GUIL", "LECO", "WALK", "ARIK", "BLDE", "BLUE", "COMO", "KING",
                       "MCDI", "SYCA", "BIGC", "MART", "MCRA", "REDB", "TECR", "WLOU", "CARI", "OKSR") # List of NEON aquatic streams

NEON_site_info <- filter(NEON_site_info, SiteID %in% NEON_stream_sites) %>% 
                  select(SiteID, WSAreaKm2, Area_NED)

NEON_site_ned_subset <- filter(NEON_site_info, WSAreaKm2 == 0)
NEON_site_ned_subset$WSAreaKm2 <- NEON_site_ned_subset$Area_NED
NEON_site_ned_subset <- select(NEON_site_ned_subset, SiteID, WSAreaKm2)

NEON_site_info <- select(NEON_site_info, SiteID, WSAreaKm2) %>% 
  filter(!(SiteID %in% NEON_site_ned_subset$SiteID))

NEON_ws_area <- rbind(NEON_site_ned_subset, NEON_site_info) %>% 
             rename(siteID = SiteID,
                    area_km2 = WSAreaKm2)

# NLCD land use percentages. NEON file name is "NEON_Aquatic_Watersheds_nlcdArea.csv"

NEON_NLCD_area <- read_csv("NEON_Aquatic_Watersheds_nlcdArea.csv") # must create a directory

NEON_NLCD <- NEON_NLCD_area %>% 
  select(site, nlcdClass, percent) %>% 
  filter(site %in% NEON_stream_sites) %>% 
  pivot_wider(names_from = "nlcdClass", values_from = "percent") %>% 
  rename(siteID = site)
NEON_NLCD[is.na(NEON_NLCD)] <- 0

### 6. BRT modeling process----

### Combine chemistry and predictor data frames to create model data set

predictors <- merge(NEON_ws_area, NEON_NLCD)
predictors <- merge(predictors, seas_index_total)
brt_data <- merge(wide_merge_calcs, predictors)

### BRT for DOC

## Re-sample data so equal number of samples per site 

count_list <- brt_data[!is.na(brt_data$DOC_mgL),] %>% 
  count(siteID) # Identify which site has least number of samples
min <- min(count_list$n)

# Recorded these values for SI table
# min
# count_list[count_list$n == min, 1]
# 27*min/sum(count_list$n)*100 # Proportion of re-sampled of total data as % 


brt_nona <- brt_data[!is.na(brt_data$DOC_mgL) &
                       !is.infinite(brt_data$DOC_mgL) ,] 

brt_data_sample <- brt_nona %>%
  group_by(siteID) %>% 
  slice_sample(n = min) %>% 
  ungroup() # re-sample data based on "min" value attained in previous code block 

# Select only numeric data from resampled frame

analyte_mat_num <- brt_data_sample  %>% 
  select(DOC_mgL, barrenLand, cultivatedCrops, deciduousForest, developedHighIntensity, developedLowIntensity,      
         developedMediumIntensity, developedOpenSpace, dwarfScrub, emergentHerbaceousWetlands,evergreenForest,           
         grasslandHerbaceous, mixedForest, openWater, pastureHay, perennialIceSnow, sedgeHerbaceous, shrubScrub,                 
         woodyWetlands, area_km2, seasonality) 
analyte_mat_num <- na.omit(analyte_mat_num)
analyte_mat_num_shuf <-  analyte_mat_num[sample(1:nrow(analyte_mat_num)), ] # shuffle data to randomize splitting in next step

## Run BRT with the GBM package

fit <- gbm::gbm(DOC_mgL ~., 
                data = analyte_mat_num_shuf,
                verbose = F, 
                shrinkage = 0.005, # learning rate, user must change depending on number of trees model produces (want highest model performance when n trees > 1000)
                interaction.depth = 3, 
                #n.minobsinnode = 5,
                train.fraction = 0.7, 
                # NOTE: model does not split data randomly. It will do train.fraction*nrow(data) to model. 
                # Therefore need to shuffle when creating data set. 
                n.trees = 5000, # numebr of iterations
                cv.folds = 10, # Number of cross-validation folds 
                bag.fraction = 0.5
)

## Model Performance

perf = gbm::gbm.perf(fit, method = "cv") # number of trees where model performs the best
vari <- var(analyte_mat_num_shuf$DOC_mgL) 
mse_cv <- fit$cv.error[perf] #this is error as MSE
mse_test <- fit$valid.error[perf]

trainR2 <- 1-(mse_cv/vari) # calculate pseudo-R2 values
testR2 <- 1-(mse_test/vari)

DOC <- data.frame(perf, vari, mse_cv, mse_test, trainR2, testR2)

## Relative influence

rel_influence_doc <- tibble::as_tibble(gbm::summary.gbm(fit, 
                                                        plotit = F))


### BRT for TDN

## Re-sample data so equal number of samples per site 

count_list <- brt_data[!is.na(brt_data$TDN_mgL),] %>% 
  count(siteID) # Identify which site has least number of samples
min <- min(count_list$n)

# Recorded these values for SI table
# min
# count_list[count_list$n == min, 1]
# 27*min/sum(count_list$n)*100 # Proportion of re-sampled of total data as % 


brt_nona <- brt_data[!is.na(brt_data$TDN_mgL) &
                     !is.infinite(brt_data$TDN_mgL) ,] 

brt_data_sample <- brt_nona %>%
  group_by(siteID) %>% 
  slice_sample(n = min) %>% 
  ungroup() # re-sample data based on "min" value attained in previous code block 

# Select only numeric data from resampled frame

analyte_mat_num <- brt_data_sample  %>% 
  select(TDN_mgL, barrenLand, cultivatedCrops, deciduousForest, developedHighIntensity, developedLowIntensity,      
                developedMediumIntensity, developedOpenSpace, dwarfScrub, emergentHerbaceousWetlands,evergreenForest,           
                grasslandHerbaceous, mixedForest, openWater, pastureHay, perennialIceSnow, sedgeHerbaceous, shrubScrub,                 
                woodyWetlands, area_km2, seasonality) 
analyte_mat_num <- na.omit(analyte_mat_num)
analyte_mat_num_shuf <-  analyte_mat_num[sample(1:nrow(analyte_mat_num)), ] # shuffle data to randomize splitting in next step

## Run BRT with the GBM package

fit <- gbm::gbm(TDN_mgL ~., 
         data = analyte_mat_num_shuf,
         verbose = F, 
         shrinkage = 0.005, # learning rate, user must change depending on number of trees model produces (want highest model performance when n trees > 1000)
         interaction.depth = 3, 
         #n.minobsinnode = 5,
         train.fraction = 0.7, 
         #NOTE: this is not random. It will do train.fraction×nrow(data) to model. 
         # Therefore need to shuffle when creating data set. 
         n.trees = 5000, # numebr of iterations
         cv.folds = 10, # Number of cross-validation folds 
         bag.fraction = 0.5
)

## Model Performance

perf = gbm::gbm.perf(fit, method = "cv") # number of trees where model performs the best
vari <- var(analyte_mat_num_shuf$TDN_mgL) 
mse_cv <- fit$cv.error[perf] #this is error as MSE
mse_test <- fit$valid.error[perf]

trainR2 <- 1-(mse_cv/vari)
testR2 <- 1-(mse_test/vari)

TDN <- data.frame(perf, vari, mse_cv, mse_test, trainR2, testR2)

## Relative influence

rel_influence_tdn <- tibble::as_tibble(gbm::summary.gbm(fit, 
                                                    plotit = F))

### BRT for TDP

## Re-sample data so equal number of samples per site 

count_list <- brt_data[!is.na(brt_data$TDP_mgL),] %>% 
  count(siteID) # Identify which site has least number of samples
min <- min(count_list$n)

# Recorded these values for SI table
# min
# count_list[count_list$n == min, 1]
# 27*min/sum(count_list$n)*100 # Proportion of re-sampled of total data as % 


brt_nona <- brt_data[!is.na(brt_data$TDP_mgL) &
                       !is.infinite(brt_data$TDP_mgL) ,] 

brt_data_sample <- brt_nona %>%
  group_by(siteID) %>% 
  slice_sample(n = min) %>% 
  ungroup() # re-sample data based on "min" value attained in previous code block 

# Select only numeric data from resampled frame

analyte_mat_num <- brt_data_sample  %>% 
  select(TDP_mgL, barrenLand, cultivatedCrops, deciduousForest, developedHighIntensity, developedLowIntensity,      
         developedMediumIntensity, developedOpenSpace, dwarfScrub, emergentHerbaceousWetlands,evergreenForest,           
         grasslandHerbaceous, mixedForest, openWater, pastureHay, perennialIceSnow, sedgeHerbaceous, shrubScrub,                 
         woodyWetlands, area_km2, seasonality) 
analyte_mat_num <- na.omit(analyte_mat_num)
analyte_mat_num_shuf <-  analyte_mat_num[sample(1:nrow(analyte_mat_num)), ] # shuffle data to randomize splitting in next step

## Run BRT with the GBM package

fit <- gbm::gbm(TDP_mgL ~., 
                data = analyte_mat_num_shuf,
                verbose = F, 
                shrinkage = 0.001, # learning rate, user must change depending on number of trees model produces (want highest model performance when n trees > 1000)
                interaction.depth = 3, 
                #n.minobsinnode = 5,
                train.fraction = 0.7, 
                #NOTE: this is not random. It will do train.fraction×nrow(data) to model. 
                # Therefore need to shuffle when creating data set. 
                n.trees = 5000, # numebr of iterations
                cv.folds = 10, # Number of cross-validation folds 
                bag.fraction = 0.5
)

## Model Performance

perf = gbm::gbm.perf(fit, method = "cv") # number of trees where model performs the best
vari <- var(analyte_mat_num_shuf$TDP_mgL) 
mse_cv <- fit$cv.error[perf] #this is error as MSE
mse_test <- fit$valid.error[perf]

trainR2 <- 1-(mse_cv/vari)
testR2 <- 1-(mse_test/vari)

TDP <- data.frame(perf, vari, mse_cv, mse_test, trainR2, testR2)

## Relative influence

rel_influence_tdp <- tibble::as_tibble(gbm::summary.gbm(fit, 
                                                        plotit = F))

### BRT for DOC_to_TDN_mgL

## Re-sample data so equal number of samples per site 

count_list <- brt_data[!is.na(brt_data$DOC_to_TDN_mgL),] %>% 
  count(siteID) # Identify which site has least number of samples
min <- min(count_list$n)

# Recorded these values for SI table
# min
# count_list[count_list$n == min, 1]
# 27*min/sum(count_list$n)*100 # Proportion of re-sampled of total data as % 


brt_nona <- brt_data[!is.na(brt_data$DOC_to_TDN_mgL) &
                       !is.infinite(brt_data$DOC_to_TDN_mgL) ,] 

brt_data_sample <- brt_nona %>%
  group_by(siteID) %>% 
  slice_sample(n = min) %>% 
  ungroup() # re-sample data based on "min" value attained in previous code block 

# Select only numeric data from resampled frame

analyte_mat_num <- brt_data_sample  %>% 
  select(DOC_to_TDN_mgL, barrenLand, cultivatedCrops, deciduousForest, developedHighIntensity, developedLowIntensity,      
         developedMediumIntensity, developedOpenSpace, dwarfScrub, emergentHerbaceousWetlands,evergreenForest,           
         grasslandHerbaceous, mixedForest, openWater, pastureHay, perennialIceSnow, sedgeHerbaceous, shrubScrub,                 
         woodyWetlands, area_km2, seasonality) 
analyte_mat_num <- na.omit(analyte_mat_num)
analyte_mat_num_shuf <-  analyte_mat_num[sample(1:nrow(analyte_mat_num)), ] # shuffle data to randomize splitting in next step

## Run BRT with the GBM package

fit <- gbm::gbm(DOC_to_TDN_mgL ~., 
                data = analyte_mat_num_shuf,
                verbose = F, 
                shrinkage = 0.01, # learning rate, user must change depending on number of trees model produces (want highest model performance when n trees > 1000)
                interaction.depth = 3, 
                #n.minobsinnode = 5,
                train.fraction = 0.7, 
                #NOTE: this is not random. It will do train.fraction×nrow(data) to model. 
                # Therefore need to shuffle when creating data set. 
                n.trees = 5000, # numebr of iterations
                cv.folds = 10, # Number of cross-validation folds 
                bag.fraction = 0.5
)

## Model Performance

perf = gbm::gbm.perf(fit, method = "cv") # number of trees where model performs the best
vari <- var(analyte_mat_num_shuf$DOC_to_TDN_mgL) 
mse_cv <- fit$cv.error[perf] #this is error as MSE
mse_test <- fit$valid.error[perf]

trainR2 <- 1-(mse_cv/vari)
testR2 <- 1-(mse_test/vari)

DOC_to_TDN_mgL <- data.frame(perf, vari, mse_cv, mse_test, trainR2, testR2)

## Relative influence

rel_influence_doc_to_tdn_mgL <- tibble::as_tibble(gbm::summary.gbm(fit, 
                                                        plotit = F))

### BRT for DOC_to_TDP_mgL

## Re-sample data so equal number of samples per site 

count_list <- brt_data[!is.na(brt_data$DOC_to_TDP_mgL),] %>% 
  count(siteID) # Identify which site has least number of samples
min <- min(count_list$n)

# Recorded these values for SI table
# min
# count_list[count_list$n == min, 1]
# 27*min/sum(count_list$n)*100 # Proportion of re-sampled of total data as % 


brt_nona <- brt_data[!is.na(brt_data$DOC_to_TDP_mgL) &
                       !is.infinite(brt_data$DOC_to_TDP_mgL) ,] 

brt_data_sample <- brt_nona %>%
  group_by(siteID) %>% 
  slice_sample(n = min) %>% 
  ungroup() # re-sample data based on "min" value attained in previous code block 

# Select only numeric data from resampled frame

analyte_mat_num <- brt_data_sample  %>% 
  select(DOC_to_TDP_mgL, barrenLand, cultivatedCrops, deciduousForest, developedHighIntensity, developedLowIntensity,      
         developedMediumIntensity, developedOpenSpace, dwarfScrub, emergentHerbaceousWetlands,evergreenForest,           
         grasslandHerbaceous, mixedForest, openWater, pastureHay, perennialIceSnow, sedgeHerbaceous, shrubScrub,                 
         woodyWetlands, area_km2, seasonality) 
analyte_mat_num <- na.omit(analyte_mat_num)
analyte_mat_num_shuf <-  analyte_mat_num[sample(1:nrow(analyte_mat_num)), ] # shuffle data to randomize splitting in next step

## Run BRT with the GBM package

fit <- gbm::gbm(DOC_to_TDP_mgL ~., 
                data = analyte_mat_num_shuf,
                verbose = F, 
                shrinkage = 0.001, # learning rate, user must change depending on number of trees model produces (want highest model performance when n trees > 1000)
                interaction.depth = 3, 
                #n.minobsinnode = 5,
                train.fraction = 0.7, 
                #NOTE: this is not random. It will do train.fraction×nrow(data) to model. 
                # Therefore need to shuffle when creating data set. 
                n.trees = 5000, # numebr of iterations
                cv.folds = 10, # Number of cross-validation folds 
                bag.fraction = 0.5
)

## Model Performance

perf = gbm::gbm.perf(fit, method = "cv") # number of trees where model performs the best
vari <- var(analyte_mat_num_shuf$DOC_to_TDP_mgL) 
mse_cv <- fit$cv.error[perf] #this is error as MSE
mse_test <- fit$valid.error[perf]

trainR2 <- 1-(mse_cv/vari)
testR2 <- 1-(mse_test/vari)

DOC_to_TDP_mgL <- data.frame(perf, vari, mse_cv, mse_test, trainR2, testR2)

## Relative influence

rel_influence_doc_to_tdp_mgL <- tibble::as_tibble(gbm::summary.gbm(fit, 
                                                               plotit = F))

### BRT for TDN_to_TDP_mgL

## Re-sample data so equal number of samples per site 

count_list <- brt_data[!is.na(brt_data$TDN_to_TDP_mgL),] %>% 
  count(siteID) # Identify which site has least number of samples
min <- min(count_list$n)

# Recorded these values for SI table
# min
# count_list[count_list$n == min, 1]
# 27*min/sum(count_list$n)*100 # Proportion of re-sampled of total data as % 


brt_nona <- brt_data[!is.na(brt_data$TDN_to_TDP_mgL) &
                       !is.infinite(brt_data$TDN_to_TDP_mgL) ,] 

brt_data_sample <- brt_nona %>%
  group_by(siteID) %>% 
  slice_sample(n = min) %>% 
  ungroup() # re-sample data based on "min" value attained in previous code block 

# Select only numeric data from resampled frame

analyte_mat_num <- brt_data_sample  %>% 
  select(TDN_to_TDP_mgL, barrenLand, cultivatedCrops, deciduousForest, developedHighIntensity, developedLowIntensity,      
         developedMediumIntensity, developedOpenSpace, dwarfScrub, emergentHerbaceousWetlands,evergreenForest,           
         grasslandHerbaceous, mixedForest, openWater, pastureHay, perennialIceSnow, sedgeHerbaceous, shrubScrub,                 
         woodyWetlands, area_km2, seasonality) 
analyte_mat_num <- na.omit(analyte_mat_num)
analyte_mat_num_shuf <-  analyte_mat_num[sample(1:nrow(analyte_mat_num)), ] # shuffle data to randomize splitting in next step

## Run BRT with the GBM package

fit <- gbm::gbm(TDN_to_TDP_mgL ~., 
                data = analyte_mat_num_shuf,
                verbose = F, 
                shrinkage = 0.001, # learning rate, user must change depending on number of trees model produces (want highest model performance when n trees > 1000)
                interaction.depth = 3, 
                #n.minobsinnode = 5,
                train.fraction = 0.7, 
                #NOTE: this is not random. It will do train.fraction×nrow(data) to model. 
                # Therefore need to shuffle when creating data set. 
                n.trees = 5000, # numebr of iterations
                cv.folds = 10, # Number of cross-validation folds 
                bag.fraction = 0.5
)

## Model Performance

perf = gbm::gbm.perf(fit, method = "cv") # number of trees where model performs the best
vari <- var(analyte_mat_num_shuf$TDN_to_TDP_mgL) 
mse_cv <- fit$cv.error[perf] #this is error as MSE
mse_test <- fit$valid.error[perf]

trainR2 <- 1-(mse_cv/vari)
testR2 <- 1-(mse_test/vari)

TDN_to_TDP_mgL <- data.frame(perf, vari, mse_cv, mse_test, trainR2, testR2)

## Relative influence

rel_influence_tdn_to_tdp_mgL <- tibble::as_tibble(gbm::summary.gbm(fit, 
                                                               plotit = F))

### BRT for TPC_mgL

## Re-sample data so equal number of samples per site 

count_list <- brt_data[!is.na(brt_data$TPC_mgL),] %>% 
  count(siteID) # Identify which site has least number of samples
min <- min(count_list$n)

# Recorded these values for SI table
# min
# count_list[count_list$n == min, 1]
# 27*min/sum(count_list$n)*100 # Proportion of re-sampled of total data as % 


brt_nona <- brt_data[!is.na(brt_data$TPC_mgL) &
                       !is.infinite(brt_data$TPC_mgL) ,] 

brt_data_sample <- brt_nona %>%
  group_by(siteID) %>% 
  slice_sample(n = min) %>% 
  ungroup() # re-sample data based on "min" value attained in previous code block 

# Select only numeric data from resampled frame

analyte_mat_num <- brt_data_sample  %>% 
  select(TPC_mgL, barrenLand, cultivatedCrops, deciduousForest, developedHighIntensity, developedLowIntensity,      
         developedMediumIntensity, developedOpenSpace, dwarfScrub, emergentHerbaceousWetlands,evergreenForest,           
         grasslandHerbaceous, mixedForest, openWater, pastureHay, perennialIceSnow, sedgeHerbaceous, shrubScrub,                 
         woodyWetlands, area_km2, seasonality) 
analyte_mat_num <- na.omit(analyte_mat_num)
analyte_mat_num_shuf <-  analyte_mat_num[sample(1:nrow(analyte_mat_num)), ] # shuffle data to randomize splitting in next step

## Run BRT with the GBM package

fit <- gbm::gbm(TPC_mgL ~., 
                data = analyte_mat_num_shuf,
                verbose = F, 
                shrinkage = 0.01, # learning rate, user must change depending on number of trees model produces (want highest model performance when n trees > 1000)
                interaction.depth = 3, 
                #n.minobsinnode = 5,
                train.fraction = 0.7, 
                #NOTE: this is not random. It will do train.fraction×nrow(data) to model. 
                # Therefore need to shuffle when creating data set. 
                n.trees = 5000, # numebr of iterations
                cv.folds = 10, # Number of cross-validation folds 
                bag.fraction = 0.5
)

## Model Performance

perf = gbm::gbm.perf(fit, method = "cv") # number of trees where model performs the best
vari <- var(analyte_mat_num_shuf$TPC_mgL) 
mse_cv <- fit$cv.error[perf] #this is error as MSE
mse_test <- fit$valid.error[perf]

trainR2 <- 1-(mse_cv/vari)
testR2 <- 1-(mse_test/vari)

TPC <- data.frame(perf, vari, mse_cv, mse_test, trainR2, testR2)

## Relative influence

rel_influence_tpc <- tibble::as_tibble(gbm::summary.gbm(fit, 
                                                        plotit = F))

### BRT for TPN_mgL

## Re-sample data so equal number of samples per site 

count_list <- brt_data[!is.na(brt_data$TPN_mgL),] %>% 
  count(siteID) # Identify which site has least number of samples
min <- min(count_list$n)

# Recorded these values for SI table
# min
# count_list[count_list$n == min, 1]
# 27*min/sum(count_list$n)*100 # Proportion of re-sampled of total data as % 


brt_nona <- brt_data[!is.na(brt_data$TPN_mgL) &
                       !is.infinite(brt_data$TPN_mgL) ,] 

brt_data_sample <- brt_nona %>%
  group_by(siteID) %>% 
  slice_sample(n = min) %>% 
  ungroup() # re-sample data based on "min" value attained in previous code block 

# Select only numeric data from resampled frame

analyte_mat_num <- brt_data_sample  %>% 
  select(TPN_mgL, barrenLand, cultivatedCrops, deciduousForest, developedHighIntensity, developedLowIntensity,      
         developedMediumIntensity, developedOpenSpace, dwarfScrub, emergentHerbaceousWetlands,evergreenForest,           
         grasslandHerbaceous, mixedForest, openWater, pastureHay, perennialIceSnow, sedgeHerbaceous, shrubScrub,                 
         woodyWetlands, area_km2, seasonality) 
analyte_mat_num <- na.omit(analyte_mat_num)
analyte_mat_num_shuf <-  analyte_mat_num[sample(1:nrow(analyte_mat_num)), ] # shuffle data to randomize splitting in next step

## Run BRT with the GBM package

fit <- gbm::gbm(TPN_mgL ~., 
                data = analyte_mat_num_shuf,
                verbose = F, 
                shrinkage = 0.0005, # learning rate, user must change depending on number of trees model produces (want highest model performance when n trees > 1000)
                interaction.depth = 3, 
                #n.minobsinnode = 5,
                train.fraction = 0.7, 
                #NOTE: this is not random. It will do train.fraction×nrow(data) to model. 
                # Therefore need to shuffle when creating data set. 
                n.trees = 5000, # numebr of iterations
                cv.folds = 10, # Number of cross-validation folds 
                bag.fraction = 0.5
)

## Model Performance

perf = gbm::gbm.perf(fit, method = "cv") # number of trees where model performs the best
vari <- var(analyte_mat_num_shuf$TPN_mgL) 
mse_cv <- fit$cv.error[perf] #this is error as MSE
mse_test <- fit$valid.error[perf]

trainR2 <- 1-(mse_cv/vari)
testR2 <- 1-(mse_test/vari)

TPN <- data.frame(perf, vari, mse_cv, mse_test, trainR2, testR2)

## Relative influence

rel_influence_tpn <- tibble::as_tibble(gbm::summary.gbm(fit, 
                                                            plotit = F))

### BRT for TPP_mgL

## Re-sample data so equal number of samples per site 

count_list <- brt_data[!is.na(brt_data$TPP_mgL),] %>% 
  count(siteID) # Identify which site has least number of samples
min <- min(count_list$n)

# Recorded these values for SI table
# min
# count_list[count_list$n == min, 1]
# 27*min/sum(count_list$n)*100 # Proportion of re-sampled of total data as % 


brt_nona <- brt_data[!is.na(brt_data$TPP_mgL) &
                       !is.infinite(brt_data$TPP_mgL) ,] 

brt_data_sample <- brt_nona %>%
  group_by(siteID) %>% 
  slice_sample(n = min) %>% 
  ungroup() # re-sample data based on "min" value attained in previous code block 

# Select only numeric data from resampled frame

analyte_mat_num <- brt_data_sample  %>% 
  select(TPP_mgL, barrenLand, cultivatedCrops, deciduousForest, developedHighIntensity, developedLowIntensity,      
         developedMediumIntensity, developedOpenSpace, dwarfScrub, emergentHerbaceousWetlands,evergreenForest,           
         grasslandHerbaceous, mixedForest, openWater, pastureHay, perennialIceSnow, sedgeHerbaceous, shrubScrub,                 
         woodyWetlands, area_km2, seasonality) 
analyte_mat_num <- na.omit(analyte_mat_num)
analyte_mat_num_shuf <-  analyte_mat_num[sample(1:nrow(analyte_mat_num)), ] # shuffle data to randomize splitting in next step

## Run BRT with the GBM package

fit <- gbm::gbm(TPP_mgL ~., 
                data = analyte_mat_num_shuf,
                verbose = F, 
                shrinkage = 0.0001, # learning rate, user must change depending on number of trees model produces (want highest model performance when n trees > 1000)
                interaction.depth = 3, 
                #n.minobsinnode = 5,
                train.fraction = 0.7, 
                #NOTE: this is not random. It will do train.fraction×nrow(data) to model. 
                # Therefore need to shuffle when creating data set. 
                n.trees = 5000, # numebr of iterations
                cv.folds = 10, # Number of cross-validation folds 
                bag.fraction = 0.5
)

## Model Performance

perf = gbm::gbm.perf(fit, method = "cv") # number of trees where model performs the best
vari <- var(analyte_mat_num_shuf$TPP_mgL) 
mse_cv <- fit$cv.error[perf] #this is error as MSE
mse_test <- fit$valid.error[perf]

trainR2 <- 1-(mse_cv/vari)
testR2 <- 1-(mse_test/vari)

TPP <- data.frame(perf, vari, mse_cv, mse_test, trainR2, testR2)

## Relative influence

rel_influence_tpp <- tibble::as_tibble(gbm::summary.gbm(fit, 
                                                            plotit = F))

### BRT for TPC_to_TPN_mgL

## Re-sample data so equal number of samples per site 

count_list <- brt_data[!is.na(brt_data$TPC_to_TPN_mgL),] %>% 
  count(siteID) # Identify which site has least number of samples
min <- min(count_list$n)

# Recorded these values for SI table
# min
# count_list[count_list$n == min, 1]
# 27*min/sum(count_list$n)*100 # Proportion of re-sampled of total data as % 


brt_nona <- brt_data[!is.na(brt_data$TPC_to_TPN_mgL) &
                       !is.infinite(brt_data$TPC_to_TPN_mgL) ,] 

brt_data_sample <- brt_nona %>%
  group_by(siteID) %>% 
  slice_sample(n = min) %>% 
  ungroup() # re-sample data based on "min" value attained in previous code block 

# Select only numeric data from resampled frame

analyte_mat_num <- brt_data_sample  %>% 
  select(TPC_to_TPN_mgL, barrenLand, cultivatedCrops, deciduousForest, developedHighIntensity, developedLowIntensity,      
         developedMediumIntensity, developedOpenSpace, dwarfScrub, emergentHerbaceousWetlands,evergreenForest,           
         grasslandHerbaceous, mixedForest, openWater, pastureHay, perennialIceSnow, sedgeHerbaceous, shrubScrub,                 
         woodyWetlands, area_km2, seasonality) 
analyte_mat_num <- na.omit(analyte_mat_num)
analyte_mat_num_shuf <-  analyte_mat_num[sample(1:nrow(analyte_mat_num)), ] # shuffle data to randomize splitting in next step

## Run BRT with the GBM package

fit <- gbm::gbm(TPC_to_TPN_mgL ~., 
                data = analyte_mat_num_shuf,
                verbose = F, 
                shrinkage = 0.001, # learning rate, user must change depending on number of trees model produces (want highest model performance when n trees > 1000)
                interaction.depth = 3, 
                #n.minobsinnode = 5,
                train.fraction = 0.7, 
                #NOTE: this is not random. It will do train.fraction×nrow(data) to model. 
                # Therefore need to shuffle when creating data set. 
                n.trees = 5000, # numebr of iterations
                cv.folds = 10, # Number of cross-validation folds 
                bag.fraction = 0.5
)

## Model Performance

perf = gbm::gbm.perf(fit, method = "cv") # number of trees where model performs the best
vari <- var(analyte_mat_num_shuf$TPC_to_TPN_mgL) 
mse_cv <- fit$cv.error[perf] #this is error as MSE
mse_test <- fit$valid.error[perf]

trainR2 <- 1-(mse_cv/vari)
testR2 <- 1-(mse_test/vari)

TPC_to_TPN_mgL <- data.frame(perf, vari, mse_cv, mse_test, trainR2, testR2)

## Relative influence

rel_influence_tpc_to_tpn_mgL <- tibble::as_tibble(gbm::summary.gbm(fit, 
                                                            plotit = F))


### BRT for TPC_to_TPP_mgL

## Re-sample data so equal number of samples per site 

count_list <- brt_data[!is.na(brt_data$TPC_to_TPP_mgL),] %>% 
  count(siteID) # Identify which site has least number of samples
min <- min(count_list$n)

# Recorded these values for SI table
# min
# count_list[count_list$n == min, 1]
# 27*min/sum(count_list$n)*100 # Proportion of re-sampled of total data as % 


brt_nona <- brt_data[!is.na(brt_data$TPC_to_TPP_mgL) &
                       !is.infinite(brt_data$TPC_to_TPP_mgL) ,] 

brt_data_sample <- brt_nona %>%
  group_by(siteID) %>% 
  slice_sample(n = min) %>% 
  ungroup() # re-sample data based on "min" value attained in previous code block 

# Select only numeric data from resampled frame

analyte_mat_num <- brt_data_sample  %>% 
  select(TPC_to_TPP_mgL, barrenLand, cultivatedCrops, deciduousForest, developedHighIntensity, developedLowIntensity,      
         developedMediumIntensity, developedOpenSpace, dwarfScrub, emergentHerbaceousWetlands,evergreenForest,           
         grasslandHerbaceous, mixedForest, openWater, pastureHay, perennialIceSnow, sedgeHerbaceous, shrubScrub,                 
         woodyWetlands, area_km2, seasonality) 
analyte_mat_num <- na.omit(analyte_mat_num)
analyte_mat_num_shuf <-  analyte_mat_num[sample(1:nrow(analyte_mat_num)), ] # shuffle data to randomize splitting in next step

## Run BRT with the GBM package

fit <- gbm::gbm(TPC_to_TPP_mgL ~., 
                data = analyte_mat_num_shuf,
                verbose = F, 
                shrinkage = 0.0005, # learning rate, user must change depending on number of trees model produces (want highest model performance when n trees > 1000)
                interaction.depth = 3, 
                #n.minobsinnode = 5,
                train.fraction = 0.7, 
                #NOTE: this is not random. It will do train.fraction×nrow(data) to model. 
                # Therefore need to shuffle when creating data set. 
                n.trees = 5000, # numebr of iterations
                cv.folds = 10, # Number of cross-validation folds 
                bag.fraction = 0.5
)

## Model Performance

perf = gbm::gbm.perf(fit, method = "cv") # number of trees where model performs the best
vari <- var(analyte_mat_num_shuf$TPC_to_TPP_mgL) 
mse_cv <- fit$cv.error[perf] #this is error as MSE
mse_test <- fit$valid.error[perf]

trainR2 <- 1-(mse_cv/vari)
testR2 <- 1-(mse_test/vari)

TPC_to_TPP_mgL <- data.frame(perf, vari, mse_cv, mse_test, trainR2, testR2)

## Relative influence

rel_influence_tpc_to_tpp_mgL <- tibble::as_tibble(gbm::summary.gbm(fit, 
                                                                   plotit = F))

### BRT for TPN_to_TPP_mgL

## Re-sample data so equal number of samples per site 

count_list <- brt_data[!is.na(brt_data$TPN_to_TPP_mgL),] %>% 
  count(siteID) # Identify which site has least number of samples
min <- min(count_list$n)

# Recorded these values for SI table
# min
# count_list[count_list$n == min, 1]
# 27*min/sum(count_list$n)*100 # Proportion of re-sampled of total data as % 


brt_nona <- brt_data[!is.na(brt_data$TPN_to_TPP_mgL) &
                       !is.infinite(brt_data$TPN_to_TPP_mgL) ,] 

brt_data_sample <- brt_nona %>%
  group_by(siteID) %>% 
  slice_sample(n = min) %>% 
  ungroup() # re-sample data based on "min" value attained in previous code block 

# Select only numeric data from resampled frame

analyte_mat_num <- brt_data_sample  %>% 
  select(TPN_to_TPP_mgL, barrenLand, cultivatedCrops, deciduousForest, developedHighIntensity, developedLowIntensity,      
         developedMediumIntensity, developedOpenSpace, dwarfScrub, emergentHerbaceousWetlands,evergreenForest,           
         grasslandHerbaceous, mixedForest, openWater, pastureHay, perennialIceSnow, sedgeHerbaceous, shrubScrub,                 
         woodyWetlands, area_km2, seasonality) 
analyte_mat_num <- na.omit(analyte_mat_num)
analyte_mat_num_shuf <-  analyte_mat_num[sample(1:nrow(analyte_mat_num)), ] # shuffle data to randomize splitting in next step

## Run BRT with the GBM package

fit <- gbm::gbm(TPN_to_TPP_mgL ~., 
                data = analyte_mat_num_shuf,
                verbose = F, 
                shrinkage = 0.0005, # learning rate, user must change depending on number of trees model produces (want highest model performance when n trees > 1000)
                interaction.depth = 3, 
                #n.minobsinnode = 5,
                train.fraction = 0.7, 
                #NOTE: this is not random. It will do train.fraction×nrow(data) to model. 
                # Therefore need to shuffle when creating data set. 
                n.trees = 5000, # numebr of iterations
                cv.folds = 10, # Number of cross-validation folds 
                bag.fraction = 0.5
)

## Model Performance

perf = gbm::gbm.perf(fit, method = "cv") # number of trees where model performs the best
vari <- var(analyte_mat_num_shuf$TPN_to_TPP_mgL) 
mse_cv <- fit$cv.error[perf] #this is error as MSE
mse_test <- fit$valid.error[perf]

trainR2 <- 1-(mse_cv/vari)
testR2 <- 1-(mse_test/vari)

TPN_to_TPP_mgL <- data.frame(perf, vari, mse_cv, mse_test, trainR2, testR2)

## Relative influence

rel_influence_tpn_to_tpp_mgL <- tibble::as_tibble(gbm::summary.gbm(fit, 
                                                                   plotit = F))

### BRT for d13C_perMill

## Re-sample data so equal number of samples per site 

count_list <- brt_data[!is.na(brt_data$d13C_perMill),] %>% 
  count(siteID) # Identify which site has least number of samples
min <- min(count_list$n)

# Recorded these values for SI table
# min
# count_list[count_list$n == min, 1]
# 27*min/sum(count_list$n)*100 # Proportion of re-sampled of total data as % 


brt_nona <- brt_data[!is.na(brt_data$d13C_perMill) &
                       !is.infinite(brt_data$d13C_perMill) ,] 

brt_data_sample <- brt_nona %>%
  group_by(siteID) %>% 
  slice_sample(n = min) %>% 
  ungroup() # re-sample data based on "min" value attained in previous code block 

# Select only numeric data from resampled frame

analyte_mat_num <- brt_data_sample  %>% 
  select(d13C_perMill, barrenLand, cultivatedCrops, deciduousForest, developedHighIntensity, developedLowIntensity,      
         developedMediumIntensity, developedOpenSpace, dwarfScrub, emergentHerbaceousWetlands,evergreenForest,           
         grasslandHerbaceous, mixedForest, openWater, pastureHay, perennialIceSnow, sedgeHerbaceous, shrubScrub,                 
         woodyWetlands, area_km2, seasonality) 
analyte_mat_num <- na.omit(analyte_mat_num)
analyte_mat_num_shuf <-  analyte_mat_num[sample(1:nrow(analyte_mat_num)), ] # shuffle data to randomize splitting in next step

min_d13C <- min(analyte_mat_num_shuf$d13C_perMill) # have to add this constant to make all values positive, BRT cannot use negative values
analyte_mat_num_shuf <- analyte_mat_num_shuf %>% 
  mutate(d13C_perMill_pos = d13C_perMill + min) %>% 
  select(-d13C_perMill)

## Run BRT with the GBM package

fit <- gbm::gbm(d13C_perMill_pos ~., 
                data = analyte_mat_num_shuf,
                verbose = F, 
                shrinkage = 0.001, # learning rate, user must change depending on number of trees model produces (want highest model performance when n trees > 1000)
                interaction.depth = 3, 
                #n.minobsinnode = 5,
                train.fraction = 0.7, 
                #NOTE: this is not random. It will do train.fraction×nrow(data) to model. 
                # Therefore need to shuffle when creating data set. 
                n.trees = 5000, # numebr of iterations
                cv.folds = 10, # Number of cross-validation folds 
                bag.fraction = 0.5
)

## Model Performance

perf = gbm::gbm.perf(fit, method = "cv") # number of trees where model performs the best
vari <- var(analyte_mat_num_shuf$d13C_perMill_pos) 
mse_cv <- fit$cv.error[perf] #this is error as MSE
mse_test <- fit$valid.error[perf]

trainR2 <- 1-(mse_cv/vari)
testR2 <- 1-(mse_test/vari)

d13C_perMill <- data.frame(perf, vari, mse_cv, mse_test, trainR2, testR2)

## Relative influence

rel_influence_d13C_perMill <- tibble::as_tibble(gbm::summary.gbm(fit, 
                                                                   plotit = F))

### BRT for d15N_perMill

## Re-sample data so equal number of samples per site 

count_list <- brt_data[!is.na(brt_data$d15N_perMill),] %>% 
  count(siteID) # Identify which site has least number of samples
min <- min(count_list$n)

# Recorded these values for SI table
# min
# count_list[count_list$n == min, 1]
# 27*min/sum(count_list$n)*100 # Proportion of re-sampled of total data as % 


brt_nona <- brt_data[!is.na(brt_data$d15N_perMill) &
                       !is.infinite(brt_data$d15N_perMill) ,] 

brt_data_sample <- brt_nona %>%
  group_by(siteID) %>% 
  slice_sample(n = min) %>% 
  ungroup() # re-sample data based on "min" value attained in previous code block 

# Select only numeric data from resampled frame

analyte_mat_num <- brt_data_sample  %>% 
  select(d15N_perMill, barrenLand, cultivatedCrops, deciduousForest, developedHighIntensity, developedLowIntensity,      
         developedMediumIntensity, developedOpenSpace, dwarfScrub, emergentHerbaceousWetlands,evergreenForest,           
         grasslandHerbaceous, mixedForest, openWater, pastureHay, perennialIceSnow, sedgeHerbaceous, shrubScrub,                 
         woodyWetlands, area_km2, seasonality) 
analyte_mat_num <- na.omit(analyte_mat_num)
analyte_mat_num_shuf <-  analyte_mat_num[sample(1:nrow(analyte_mat_num)), ] # shuffle data to randomize splitting in next step

min_d15N <- min(analyte_mat_num_shuf$d15N_perMill) # have to add this constant to make all values positive, BRT cannot use negative values
analyte_mat_num_shuf <- analyte_mat_num_shuf %>% 
  mutate(d15N_perMill_pos = d15N_perMill + min) %>% 
  select(-d15N_perMill)

## Run BRT with the GBM package

fit <- gbm::gbm(d15N_perMill_pos ~., 
                data = analyte_mat_num_shuf,
                verbose = F, 
                shrinkage = 0.001, # learning rate, user must change depending on number of trees model produces (want highest model performance when n trees > 1000)
                interaction.depth = 3, 
                #n.minobsinnode = 5,
                train.fraction = 0.7, 
                #NOTE: this is not random. It will do train.fraction×nrow(data) to model. 
                # Therefore need to shuffle when creating data set. 
                n.trees = 5000, # numebr of iterations
                cv.folds = 10, # Number of cross-validation folds 
                bag.fraction = 0.5
)

## Model Performance

perf = gbm::gbm.perf(fit, method = "cv") # number of trees where model performs the best
vari <- var(analyte_mat_num_shuf$d15N_perMill_pos) 
mse_cv <- fit$cv.error[perf] #this is error as MSE
mse_test <- fit$valid.error[perf]

trainR2 <- 1-(mse_cv/vari)
testR2 <- 1-(mse_test/vari)

d15N_perMill <- data.frame(perf, vari, mse_cv, mse_test, trainR2, testR2)

## Relative influence

rel_influence_d15N_perMill <- tibble::as_tibble(gbm::summary.gbm(fit, 
                                                                 plotit = F))

### 7. BRT performance metrics, relative influence values, and partial dependency plots ----

## Performance metrics table

analyte_names <- c("DOC", "TDN", "TDP", "DOC_to_TDN_mgL", "DOC_to_TDP_mgL", "TDN_to_TDP_mgL", 
                   "TPC", "TPN", "TPP", "TPC_to_TPN_mgL", "TPC_to_TPP_mgL", "TPN_to_TPP_mgL", 
                   "d13C_perMill", "d15N_perMill")
analyte_names <-  as.data.frame(analyte_names)
             
perf_table <- list(DOC, TDN, TDP, DOC_to_TDN_mgL, DOC_to_TDP_mgL, TDN_to_TDP_mgL, 
             TPC, TPN, TPP, TPC_to_TPN_mgL, TPC_to_TPP_mgL, TPN_to_TPP_mgL, 
             d13C_perMill, d15N_perMill) # how to have these names so they carry over
perf_table <- bind_rows(perf_table)
perf_table <- cbind(perf_table, analyte_names)
perf_table <- column_to_rownames(perf_table, var = "analyte_names") %>% 
  rename(n_trees = perf)

## Relative influence values table

rel_inf_table <- list(rel_influence_doc, rel_influence_tdn, rel_influence_tdp, rel_influence_doc_to_tdn_mgL, 
                   rel_influence_doc_to_tdn_mgL, rel_influence_tdn_to_tdp_mgL, rel_influence_tpc, rel_influence_tpn, 
                   rel_influence_tpp, rel_influence_tpc_to_tpn_mgL, rel_influence_tpc_to_tpp_mgL, rel_influence_tpn_to_tpp_mgL, 
                  rel_influence_d13C_perMill, rel_influence_d15N_perMill)

rel_inf_table <- reduce(rel_inf_table, left_join, by = "var")

colnames(rel_inf_table) <- c("predictor", "DOC", "TDN", "TDP", "DOC_to_TDN_mgL", "DOC_to_TDP_mgL", "TDN_to_TDP_mgL", 
                   "TPC", "TPN", "TPP", "TPC_to_TPN_mgL", "TPC_to_TPP_mgL", "TPN_to_TPP_mgL", 
                   "d13C_perMill", "d15N_perMill")

## Example code to visualize partial dependency plots

which(colnames(analyte_mat_num)=="pastureHay") # input n-1 into i.var below to choose desired predictor variable
partial_dep <- gbm::plot.gbm(fit, i.var = c(14)) 

partial_dep$xlab = "Pasture/Hay" # Change labels manually depending on what is graphed
partial_dep$ylab = "TDN (mg/L)" 
partial_dep

### 8. Flow regime hydrology summary stats----

## Summary statistics of the flow regime for each site with hydrostats::baseflows

q_split <- list()
q_split <- split(q_long, q_long$siteID)
q_baseflows <- lapply(q_split, function(x) select(x, collectDatetime, Q_Ls))
q_baseflows <- lapply(q_baseflows, function(x) rename(x, Date = collectDatetime,  
                                                        Q = Q_Ls))
q_baseflows <- lapply(q_baseflows, function(x) baseflows(x))
q_baseflows <- bind_rows(q_baseflows, .id = "siteID")

## Calculate Q10 and Q90 for each site 

q_quantile <- q_long %>% 
  group_by(siteID) %>% 
  summarise(q_q90 = quantile(Q_Ls, probs = 0.90, na.rm = T),
            q_q10 = quantile(Q_Ls, probs = 0.10, na.rm = T)) %>% 
  merge(q_baseflows)

## Count number of particulate chemistry samples taken at very high/low flow events

# Make a data frame in long format of TPC, TPN, and TPP

tpp_subset <- wide_merge_tidy %>% 
  select(siteID, collectDate, TPP_mgL) %>% 
  rename(analyteConcentration = TPP_mgL) %>% 
  na.omit()
tpp_subset$analyte <- "TPP"

tpc_tpn_subset <- NEON_chem_long %>% 
  filter(analyte == "TPC" |
         analyte == "TPN") %>% 
  select(siteID, collectDate, analyte, analyteConcentration)
tpc_tpn_subset$analyteConcentration <- tpc_tpn_subset$analyteConcentration / 10^3 # Convert to mg/L

part_subset <- rbind(tpp_subset, tpc_tpn_subset)

# Merge frames to determine the sampling distribution above/below Q90 and Q10

cq_merge <- merge(part_subset, q_long, by.x = c("collectDate", "siteID"), by.y = c("collectDatetime", "siteID"), all.x = T)

q_quantile_merge <- merge(cq_merge, q_quantile)
q_quantile_sample_count <-  q_quantile_merge %>% 
  group_by(siteID) %>% 
  summarise(n_abv90 = sum(Q_Ls > q_q90, na.rm = T),
            n_blw10 = sum(Q_Ls < q_q10, na.rm = T),
            n_total = sum(!is.na(Q_Ls)))

## Calculate Q10, Q50, and Q90 for specific discharge

q_baseflows_area <- merge(q_quantile, NEON_ws_area)
q_specific <- q_baseflows_area %>% 
  mutate(Q50_spec = Q50 / area_km2,
         q10_spec = q_q10 / area_km2,
         q90_spec = q_q90 / area_km2)

## Summary table of hydrology by site 

hydro_sum_stats <- merge(q_quantile, q_specific)
hydro_sum_stats <- merge(hydro_sum_stats, q_quantile_sample_count)
hydro_sum_stats <- merge(hydro_sum_stats, seas_index_total) %>% 
  select(-area_km2)

### 9. Process NEON benthic chemistry data to compare with NEON suspended seston data----

## Download algal data, NEON data product DP1.20166.001

algae_neon_list <- loadByProduct(dpID="DP1.20163.001", 
                            site = c("HOPB", "LEWI", "POSE", "FLNT",
                                     "CUPE", "GUIL", "KING", "MCDI",
                                     "LECO", "WALK", "BLWA", "MAYF",
                                     "TOMB", "ARIK", "BLUE", "PRIN",
                                     "BLDE", "COMO", "WLOU", "SYCA",
                                     "REDB", "MART", "MCRA", "BIGC",
                                     "TECR", "OKSR", "CARI"),
                                              check.size = F, 
                            startdate = '2014-12-08', #min(iso_data_total$collectDate)
                            enddate =  "2021-12-01", #max(iso_data_total$collectDate)
                            include.provisional=TRUE)

algae_neon <- algae_neon_list$alg_algaeExternalLabDataPerSample

## Filter data to select only benthic samples

# table(algae_neon_list$alg_fieldData$algalSampleType) ## this returns epilithon, epilithon_largeSubstrate, 
# epipelon, epiphyton, epipsammon, epixylon, epixylon_largeSubstrate, phytoplankton, and seston. 
# Want to filter for benthic samples so remove data labeled phytoplankton or seston.

algae_seston <- grepl("seston", algae_neon$sampleID, ignore.case = T)
algae_phyto <- grepl("phytoplankton", algae_neon$sampleID, ignore.case = T)
algae_rm <- Reduce("|", list(algae_seston, algae_phyto))
algae_benthic <- algae_neon[!algae_rm, ]

## Filter for C, N and P 

algae_benthic_clean <- algae_benthic %>% 
  filter(analyte == "carbon" |
           analyte == "nitrogen" |
           analyte == "phosphorus") %>% 
  mutate(conc_mgL = analyteConcentration/10^3) # convert from ug/L to mg/L

## Calculate molar stoichiometric ratios

algae_benthic_calcs_bysite <- algae_benthic_clean %>% 
select(siteID, analyte, conc_mgL) %>% 
  group_by(siteID, analyte) %>% 
  summarise(mean_conc_mgL = mean(conc_mgL, na.rm = T)) %>%
  pivot_wider(names_from = "analyte", values_from = "mean_conc_mgL") %>% 
  rename(mean_c_mgL = carbon,
         mean_n_mgL = nitrogen, 
         mean_p_mgL = phosphorus) %>% 
  mutate(mean_c_mmolL = mean_c_mgL*(1/12.011), # Unit conversion formula: 1 ug/L = 1/MW*umol/L
         mean_n_mmolL = mean_n_mgL*(1/14.007),
         mean_p_mmolL = mean_p_mgL*(1/30.974),
         cn_ratio = mean_c_mmolL/mean_n_mmolL,
         cp_ratio = mean_c_mmolL/mean_p_mmolL,
         np_ratio = mean_n_mmolL/mean_p_mmolL)

################### Magnesium Analysis ####################

########### Libraries, working drive, functions ###########
setwd("C:/Users/chris/OneDrive/Documents/Research/Electrolyte Replacement CaMPPhIRE/Instrumental Variable SHN LHO/Magnesium/Analysis")

library(tidyverse)
library(brms)
library(mice)
library(xgboost)
library(rsample)
library(recipes)
library(purrr)
library(grf)
library(fastshap)
library(ggplot2)
library(yardstick)
library(pROC)
library(patchwork)

# 3-class YlGnBu
green <- "#7fcdbb"
blue  <- "#2c7fb8"

# 3-class YlOrRd
lightorange <- "#feb24c"
burntorange <- "#f03b20"


########## Functions ##########################

# Define the standardized difference function, handling NA values
standardized_diff <- function(mean1, mean2, sd1, sd2) {
  mean_diff <- abs(mean1 - mean2)
  pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
  mean_diff / pooled_sd
}

# Mean and 95% credible interval from logodds
mean_95CrI_OR <- function(x){
  c(
    "mean" = exp(mean(x)),
    "lb95" = exp(quantile(x, 0.025)),
    "ub95" = exp(quantile(x, 0.975)),
    "P_leq_1" = mean(x < 0)
  )
}

# Make ARR from control group event rate
meanARR_95CrI <- function(x, rate){
  post_arr = rate - inv_logit_scaled(logit_scaled(rate)+x)
  c(mean = mean(post_arr),
    lb95 = quantile(post_arr, 0.025),
    ub95 = quantile(post_arr, 0.975))
}

# median odds ratio
mor <- function(v) exp(sqrt(2)*v*qnorm(0.75))


########## Load baseline data ############################


df_baseline <- 
  read_csv("Deidentified Data/Mag_Baseline_2025-06-14.csv") %>%
   # dichotomize preferred languages
  mutate(Value = ifelse(
    Characteristic %in% "Preferred Language",
    ifelse(Value == "English", Value, "Not English"),
    Value)) %>%
  # remove hospital names
  mutate(Value = ifelse(
    Characteristic %in% "Hospital",
    case_when(Value == "Scarborough Health Network - G" ~ "A",
              Value == "Scarborough Health Network - B" ~ "B",
              Value == "Scarborough Health Network - C" ~ "C",
              Value == "Lakeridge Health - Oshawa" ~ "D",
              Value == "Lakeridge Health - Ajax/Picker" ~ "E"),
    Value
  )) %>% # categorize principal problem diagnoses
  left_join(
    select(read_csv("PrincipalProblems2.csv"), Value, Category),
    by = "Value"
  ) %>%
  mutate(Value = ifelse(
    Characteristic == "Principal Problem",
    Category, Value
  )) %>%
  select(-Category)

# 6 patients less than 18 years old
teens <- filter(df_baseline, Characteristic == "Age", as.numeric(Value) < 18)

df_baseline <- filter(df_baseline, 
                      !(Patient %in% teens$Patient))

###### Load outcomes data ############

df_outcome <- 
  read_csv("Deidentified Data/Mag_Outcomes_2025-06-14.csv") %>%
  filter(!(Patient %in% teens$Patient))

patients_with_rhythms <- filter(df_outcome, Characteristic == "Rhythm")$Patient

patients_without_rhythms <- filter(df_outcome,
                                   !(Patient %in% 
                                       patients_with_rhythms)) %>%
  left_join(df_baseline %>%
              filter(Characteristic == "Magnesium") %>%
              mutate(Replacement = as.numeric(Value) < 0.96) %>%
              select(Patient, Replacement),
            by = "Patient")

df_full_baseline <- df_baseline
df_full_outcome <- df_outcome

df_baseline <- filter(df_baseline,
                      Patient %in% patients_with_rhythms)

df_outcome <- filter(df_outcome,
                     Patient %in% df_baseline$Patient)

N_total <- nrow(count(df_baseline, Patient))

N_group <- 
  df_baseline %>%
  filter(Characteristic == "Magnesium") %>%
  summarise(replacement = sum(Value <= 0.95),
            noreplacement = sum(Value > 0.95))


# Fix blood pressures
SBP_DBP<-
  df_baseline %>%
  filter(Characteristic == "BLOOD PRESSURE") %>%
  separate(Value, into = c("SBP", "DBP"), sep = "/", convert = TRUE) %>%
  select(Patient, SBP, DBP) %>%
  pivot_longer(SBP:DBP,
               names_to = "Characteristic",
               values_to = "Value") %>%
  mutate(Value = as.character(Value))

df_baseline2 <- 
  df_baseline %>% 
    filter(Characteristic != "BLOOD PRESSURE") %>%
  filter(Characteristic == "Magnesium") %>%
  mutate(Exposure = ifelse(Value <= 0.95,
                           "Replacement",
                           "Noreplacement")) %>% 
  select(Patient, Exposure) %>%
  right_join(bind_rows(
    df_baseline,
    SBP_DBP),
    by = "Patient") %>%
  filter(Characteristic != "BLOOD PRESSURE")

################### Table 1 with Standard Diffs #####################

characteristics_category <- c("Comorbidity",
                              "Hospital",
                              "Preferred Language",
                              "Principal Problem",
                              "Sex", 
                              "CRRT",
                              "R CARDIAC RHYTHM",
                              "Intubated")

characteristics_ignore <- c("BLOOD PRESSURE",
                            "HCO3 VEN",
                            "HCO3, MIXED VENOUS",
                            "TOTAL HEMOGLOBIN")

# --- Continuous variables with Total column ---
StandardDiffTableA <-
  df_baseline2 %>%
  filter(!(Characteristic %in% characteristics_category
           | Characteristic %in% characteristics_ignore)) %>%
  mutate(Value = as.numeric(Value)) %>%
  group_by(Characteristic, Exposure) %>%
  summarise(Mean = mean(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE),
            N = n(),
            .groups = "drop") %>%
  pivot_longer(cols = c(Mean, SD, N), names_to = "Measure", values_to = "Value") %>%
  pivot_wider(names_from = Exposure:Measure, values_from = Value) %>%
  mutate(
    Total_Mean = (Replacement_Mean * Replacement_N + Noreplacement_Mean * Noreplacement_N) /
      (Replacement_N + Noreplacement_N),
    Total_SD = sqrt((
      ((Replacement_N - 1) * Replacement_SD^2 +
         (Noreplacement_N - 1) * Noreplacement_SD^2) /
        (Replacement_N + Noreplacement_N - 2)
    )),
    `Standardized Difference` = standardized_diff(
      Replacement_Mean,
      Noreplacement_Mean,
      Replacement_SD,
      Noreplacement_SD
    ),
    Total = paste0(round(Total_Mean, 2), " (", round(Total_SD, 2), ")"),
    Replacement = paste0(round(Replacement_Mean, 2), " (", round(Replacement_SD, 2), ")"),
    `No replacement` = paste0(round(Noreplacement_Mean, 2), " (", round(Noreplacement_SD, 2), ")")
  ) %>%
  select(
    Characteristic,
    Total,
    `No replacement`,
    Replacement,
    `Standardized Difference`
  ) %>%
  arrange(Characteristic)

saveRDS(StandardDiffTableA, "Table1A.RDS")

# --- Categorical variables with Total column ---
Table1B_with_diff <- 
  df_baseline2 %>%
  mutate(Value = ifelse(
    Characteristic == "R CARDIAC RHYTHM",
    ifelse(Value %in% c("A-fib", "A-flutter"),
           "Afib/flutter", "Not AFib/flutter"),
    Value
  )) %>%
  filter(Characteristic %in% characteristics_category) %>%
  group_by(Characteristic, Value, Exposure) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(N = ifelse(Exposure == "Replacement", N_group[[1]], N_group[[2]])) %>%
  mutate(prop = count / N) %>%
  select(Characteristic, Value, Exposure, prop) %>%
  pivot_wider(names_from = Exposure, values_from = prop, values_fill = 0) %>%
  mutate(
    pooled = (Replacement * N_group[[1]] + `Noreplacement` * N_group[[2]]) / (N_group[[1]] + N_group[[2]]),
    `Standardized Difference` = abs(`Replacement` - `Noreplacement`) / sqrt(pooled * (1 - pooled))
  ) %>%
  rename(Noreplacement_Proportion = Noreplacement,
         Replacement_Proportion = Replacement,
         Pooled_Proportion = pooled) %>%
  left_join(
    df_baseline2 %>%
      mutate(Value = ifelse(
        Characteristic == "R CARDIAC RHYTHM",
        ifelse(Value %in% c("A-fib", "A-flutter"),
               "Afib/flutter", "Not AFib/flutter"),
        Value
      )) %>%
      filter(Characteristic %in% characteristics_category) %>%
      group_by(Characteristic, Value, Exposure) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(N = ifelse(Exposure == "Replacement", N_group[[1]], N_group[[2]])) %>%
      mutate(perc = round(100 * count / N),
             label = paste0(count, " (", perc, "%)")) %>%
      select(Characteristic, Value, Exposure, label) %>%
      pivot_wider(names_from = Exposure, values_from = label),
    by = c("Characteristic", "Value")
  ) %>%
  mutate(
    total_count = (N_group[[1]] + N_group[[2]]) * Pooled_Proportion,
    total_perc = round(100 * Pooled_Proportion),
    Total = paste0(round(total_count), " (", total_perc, "%)")
  ) %>%
  select(Characteristic, Value, Total, `No replacement` = Noreplacement, Replacement, `Standardized Difference`) %>%
  mutate(Characteristic = ifelse(Characteristic == "R CARDIAC RHYTHM",
                                 "Baseline cardiac rhythm",
                                 Characteristic)) %>%
  arrange(Characteristic, Value)

saveRDS(Table1B_with_diff, "Table1B.RDS")


# Add 'Section' column to each table
StandardDiffTableA <- StandardDiffTableA %>%
  mutate(Section = "Continuous Variables",
         Value = "") %>%
  relocate(Section, Characteristic, Value)

Table1B_with_diff <- Table1B_with_diff %>%
  mutate(Section = "Categorical Variables") %>%
  relocate(Section)

# Ensure both are numeric for combining
StandardDiffTableA <- StandardDiffTableA %>%
  mutate(`Standardized Difference` = as.numeric(`Standardized Difference`),
         Section = "Continuous Variables")

Table1B_with_diff <- Table1B_with_diff %>%
  mutate(`Standardized Difference` = as.numeric(`Standardized Difference`),
         Section = "Categorical Variables")

# Combine
CombinedTable1 <- bind_rows(StandardDiffTableA, Table1B_with_diff) %>%
  arrange(Section, Characteristic)

# Optional: format SD column with rounding if you prefer
CombinedTable1 <- CombinedTable1 %>%
  mutate(`Standardized Difference` = round(`Standardized Difference`, 2))

# Export
write.csv(CombinedTable1, "Tables/Table1_Combined.csv", row.names = FALSE)

# median IQRs for some variables
df_baseline2 %>% filter(Characteristic %in% c("Age",
                                              "Hospital Length of Stay (Days)",
                                              "ICU Length of Stay (Hours)",
                                              "Creatinine",
                                              "Norepi Equivalents Per Hour")) %>%
  bind_rows(df_baseline2 %>% 
              filter(Characteristic %in% c("Age",
                                                          "Hospital Length of Stay (Days)",
                                                          "ICU Length of Stay (Hours)",
                                           "Creatinine",
                                           "Norepi Equivalents Per Hour")) %>%
              mutate(Exposure = "Total")
              ) %>%
  group_by(Exposure, Characteristic) %>%
  mutate(Value = as.numeric(Value)) %>%
  mutate(Value = ifelse(Characteristic == "Norepi Equivalents Per Hour",
                        Value/60, Value)) %>%
  summarise(median = median(Value, na.rm = T),
            IQR25 = quantile(Value, 0.25, na.rm = T),
            IQR75 = quantile(Value, 0.75, na.rm = T)) %>%
  arrange(Characteristic)


############ Table 2 ####################

df_outcomes_wide <-
  df_outcome %>%
  mutate(AfibFlutter = ifelse(Characteristic == "Rhythm",
                              ifelse(Value %in% c("A-fib",
                                                  "A-flutter"), 
                                     1, 0),
                              0),
         MagReplaced = ifelse(Characteristic == "NamedMedication", 
                              1,0),
         Amiodarone = ifelse(Characteristic == "GroupedMedication" &
                               Value == "Amiodarone",
                             1,0),
         BloodDraws = ifelse(Characteristic == "BloodDraws_24h",
                             as.numeric(Value), 0),
         MagLevel = ifelse(Characteristic == "Magnesium",
                           as.numeric(Value), NA),
         Tachyarrhythmia = ifelse(Characteristic == "Rhythm",
                                  ifelse(Value %in% c("VF","VT", "SVT"), 1, 0),
                                  0),
         NEqDose = ifelse(Characteristic == "NorepiEquivalents_24h",
                          as.numeric(Value),
                          0),
         Death24h = ifelse(Characteristic == "TimeToDeath_h" &
                             Value != "NULL",
                           ifelse(as.numeric(Value) <= 24,1,0),
                           0),
         RRT72h = ifelse(Characteristic == "CRRT_72h",
                         1,
                         0),
         Ett72h = ifelse(Characteristic == "Intubation_72h",
                         1,
                         0),
         Ext72h = ifelse(Characteristic == "Extubation_72h",
                         1,
                         0)#,
#         Vent72h = ifelse(Characteristic == "VentDuration_72h",
 #                         as.numeric(Value),0)
                        ) %>% 
  group_by(Patient) %>%
  summarise(AfibFlutter = max(AfibFlutter, na.rm = T),
            MagReplaced = max(MagReplaced, na.rm = T),
            BloodDraws = max(BloodDraws, na.rm = T),
            MagLevel = max(MagLevel, na.rm = T),
            Amiodarone = max(Amiodarone, na.rm = T),
            Tachyarrhythmia = max(Tachyarrhythmia, na.rm = T),
            NEqDose = ifelse(is.na(max(NEqDose, na.rm = T)), 0,
                             max(NEqDose)),
            Death24h= max(Death24h),
            RRT72h = max(RRT72h),
            Ett72h = max(Ett72h),
            Ext72h = max(Ext72h)
#            Vent72h = max(Vent72h)
            ) %>%
  mutate(MagLevel = ifelse(MagLevel > 0, MagLevel, NA)) %>%
  left_join(distinct(df_baseline2, Patient, Exposure),
            by = "Patient") %>%
  relocate(Patient, Exposure)

table2 <- 
df_outcomes_wide %>% 
  mutate(NorEpi = NEqDose > 0,
 #        Vent72h = Vent72h>0,
         BloodDraws24h_None = BloodDraws == 0,
         BloodDraws24h_One = BloodDraws == 1,
         BloodDraws24h_TwoOrMore = BloodDraws > 1,
         MagLevel_leq_0.7 = MagLevel <= 0.7,
         MagLevel_btw_0.7to0.95 = MagLevel > 0.7 & MagLevel <= 0.95,
         MagLevel_gthan_0.95 = MagLevel > 0.95) %>%
  select(-NEqDose, -BloodDraws, -MagLevel) %>%
  pivot_longer(AfibFlutter:MagLevel_gthan_0.95) %>%
  group_by(Exposure, name) %>%
  summarise(Count = sum(value, na.rm = T),
            Proportion = mean(value, na.rm = T)) %>%
  mutate(combined = paste0(Count, " (",
                           round(Proportion*100,1),
                           "%)")) %>%
  arrange(name) %>%
  select(Exposure, name, combined) %>%
  pivot_wider(names_from = Exposure,
              values_from = combined) %>%
  write_csv("Tables/table2.csv")

# magnesium levels
df_outcomes_wide %>%
  select(Exposure, MagLevel) %>%
  group_by(Exposure) %>%
  summarise(Count = sum(!is.na(MagLevel)),
            Mean = mean(MagLevel, na.rm=T),
            Median = median(MagLevel, na.rm=T),
            q25 = quantile(MagLevel, 0.25, na.rm=T),
            q75 = quantile(MagLevel, 0.75, na.rm=T),
            leq0.7 = mean(MagLevel <= 0.7, na.rm = T),
            geq1.0 = mean(MagLevel >= 1.0, na.rm = T)) %>%
  write_csv("Tables/magnesium_levels.csv")


# mag sulfate vs other

negate = function(x) !x

df_baseline2 %>%
  distinct(Patient, Exposure) %>%
  left_join(
    df_outcome %>%
      filter(Characteristic == "NamedMedication"),
    by = "Patient"
  ) %>% 
  mutate(Value = ifelse(is.na(Value),"None", Value)) %>%
  mutate(Given = 1) %>%
  pivot_wider(names_from = Value, values_from = Given) %>%
  group_by(Exposure) %>% 
  filter(Characteristic == "NamedMedication") %>%
  mutate(across(None:`magnesium citrate`, is.na)) %>%
  mutate(across(None:`magnesium citrate`, negate)) %>%
  summarise(Count = across(None:`magnesium citrate`, sum)) %>%
  t() %>%
  write_csv("Tables/magnesium_agents_within8h.csv")
  

############ Figure 1 ####################

#df_outcomes_fig1 <-
  df_outcome %>%
  mutate(AfibFlutter = ifelse(Characteristic == "Rhythm",
                              ifelse(Value %in% c("A-fib",
                                                  "A-flutter"), 
                                     1, 0),
                              0),
         MagReplaced = ifelse(Characteristic == "NamedMedication", 
                              1,0),
         Tachyarrhythmia = ifelse(Characteristic == "Rhythm",
                                  ifelse(Value %in% c("VF","VT", "SVT"), 1, 0),
                                  0),
         Death24h = ifelse(Characteristic == "TimeToDeath_h" &
                             Value != "NULL",
                           ifelse(as.numeric(Value) <= 24,1,0),
                           0)) %>% 
  group_by(Patient) %>%
  summarise(AfibFlutter = max(AfibFlutter, na.rm = T),
            MagReplaced = max(MagReplaced, na.rm = T),
            Tachyarrhythmia = max(Tachyarrhythmia, na.rm = T),
            Death24h= max(Death24h)) %>%
  left_join(distinct(df_baseline2, Patient, Exposure),
            by = "Patient") %>%
  pivot_longer(AfibFlutter:Death24h,
               names_to = "Outcome",
               values_to = "Value") %>%
  mutate(Outcome = factor(Outcome,
                          levels = c("MagReplaced",
                                     "AfibFlutter",
                                     "Tachyarrhythmia",
                                     "Death24h"),
                          labels = c("Magnesium replaced within 8h",
                                     "Atrial fibrillation within 24h",
                                     "VT, VF, or SVT within 24h",
                                     "Death within 24h"),
                          ordered = T),
         Exposure = factor(Exposure,
                           levels = c("Noreplacement",
                                      "Replacement"),
                           labels = c("No replacement",
                                      "Replacement"))) %>%
  left_join(filter(
    df_baseline, Characteristic == "Magnesium") %>%
      select(-Characteristic) %>% rename(Magnesium = Value) %>%
      mutate(Magnesium = as.numeric(Magnesium)), by = "Patient") %>%
  group_by(Exposure, Magnesium, Outcome) %>%
  summarise(
    Count = sum(Value),
    Proportion = round(100 * mean(Value), 1),
    Label = paste0(Count, " (", round(Proportion), "%)"),
    .groups = "drop"
  ) %>%
  filter(grepl("Atrial", Outcome) | grepl("Magnesium",Outcome)) %>%
  ggplot(aes(fill = Exposure, x = Magnesium, y = Proportion)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(
    aes(label = Label),
#    position = position_dodge(width = 0.9),
    vjust = -1,
    color = "black",
    size = 2
  )  +
  facet_wrap(Outcome ~ .,
             ncol = 2,
             scales = "free_x",
             strip.position = "top") +
  scale_fill_manual(values = c(green, blue),
                    guide = guide_legend(reverse = TRUE)) +
  theme_minimal() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  labs(x = "Magnesium (mmol/L)",
       y = "Proportion (%)") +
  lims(y = c(0,100))

ggsave("Figures/Figure1.svg",
       width = 7, height = 4)

ggsave("Figures/Figure1.jpg",
       width = 7, height = 4,
       dpi = 300)

############ Primary Analysis ####################


primary_outcome <-
  df_outcome %>%
  filter(Characteristic == "Rhythm") %>%
  mutate(Value = as.numeric(ifelse(Value %in% c("A-fib",
                                     "A-flutter"),1,0))) %>%
  group_by(Patient) %>%
  summarise(Outcome = max(Value))

# primary outcome with hospital only
df_brms_primary <-
  df_baseline2 %>% 
    left_join(primary_outcome,
              by = "Patient") %>%
    filter(Characteristic == "Hospital") %>%
  rename(Hospital = Value) %>%
  select(-Characteristic)

# simple count by exposure status and hospital
df_brms_primary %>%
  count(Exposure, Hospital, Outcome)

# primary analysis

fit_primary <- brm(
  Outcome ~ Exposure + (1 | Hospital),
  data = df_brms_primary,
  prior =  c(set_prior("normal(0, 1)", class = "b"),
             set_prior("cauchy(0, 25)", class = "sd")),
  family = "bernoulli",
  cores = 4,
  control = list(adapt_delta = 0.9)
  )

saveRDS(fit_primary, file = "Fits/fit_primary.RDS")

fit_primary <- readRDS("Fits/fit_primary.RDS")

post_primary <- as_draws_df(fit_primary)

primary_OR <- mean_95CrI_OR(post_primary$b_ExposureReplacement)

# Control group event rate
control_rate_primary <- df_brms_primary %>% 
  filter(Exposure == "Noreplacement") %>% 
  summarise(control = mean(Outcome)) %>%
  unlist()

# ARR using control group rate
primary_ARR <- 
  meanARR_95CrI(x = post_primary$b_ExposureReplacement,
                rate = control_rate_primary)

############# Secondary Outcomes ######################
# Outcomes: magnesium replacement, VT, VF, SVT, antiarrhythmic use, 
#           vasopressor dose, mortality, RRT, intubation, 
#           extubation, ventilation duration

# Process secondary outcomes from the df_outcomes database

df_outcomes_secondary <-
  df_outcome %>%
  mutate(AfibFlutter = ifelse(Characteristic == "Rhythm",
                              ifelse(Value %in% c("A-fib",
                                                  "A-flutter"), 
                                     1, 0),
                              0),
         MagReplaced = ifelse(Characteristic == "NamedMedication", 
                              1,0),
         Amiodarone = ifelse(Characteristic == "GroupedMedication" &
                               Value == "Amiodarone",
                             1,0),
         BloodDraws = ifelse(Characteristic == "BloodDraws_24h",
                             as.numeric(Value), 0),
         MagLevel = ifelse(Characteristic == "Magnesium",
                           as.numeric(Value), NA),
         Tachyarrhythmia = ifelse(Characteristic == "Rhythm",
                                  ifelse(Value %in% c("VF","VT", "SVT"), 1, 0),
                                  0),
         NEqDose = ifelse(Characteristic == "NorepiEquivalents_24h",
                          as.numeric(Value),
                          0),
         Death24h = ifelse(Characteristic == "TimeToDeath_h" &
                             Value != "NULL",
                           ifelse(as.numeric(Value) <= 24,1,0),
                           0),
         RRT72h = ifelse(Characteristic == "CRRT_72h",
                         1,
                         0),
         Ett72h = ifelse(Characteristic == "Intubation_72h",
                         1,
                         0),
         Ext72h = ifelse(Characteristic == "Extubation_72h",
                         1,
                         0),
         Vent72h = ifelse(Characteristic == "VentDuration_72h",
                          as.numeric(Value),
                          0)) %>% 
  group_by(Patient) %>%
  summarise(AfibFlutter = max(AfibFlutter, na.rm = T),
            MagReplaced = max(MagReplaced, na.rm = T),
            BloodDraws = max(BloodDraws, na.rm = T),
            MagLevel = max(MagLevel, na.rm = T),
            Amiodarone = max(Amiodarone, na.rm = T),
            Tachyarrhythmia = max(Tachyarrhythmia, na.rm = T),
            NEqDose = ifelse(is.na(max(NEqDose, na.rm = T)), 0,
                             max(NEqDose)),
            Death24h= max(Death24h),
            RRT72h = max(RRT72h),
            Ett72h = max(Ett72h),
            Ext72h = max(Ext72h),
            Vent72h = max(Vent72h)) %>%
  mutate(MagLevel = ifelse(MagLevel > 0, MagLevel, NA)) %>%
  left_join(select(
    filter(
      df_baseline2,
      Characteristic == "Hospital"),
    Patient, Exposure, Value),
            by = "Patient") %>%
  relocate(Patient, Exposure, Value) %>%
  rename(Hospital = Value)


names(df_outcomes_secondary)

binary_outcomes <- c("MagReplaced",
                     "Amiodarone",
                     "Tachyarrhythmia",
                     "Death24h",
                     "RRT72h",
                     "Ett72h",
                     "Ext72h")

secondary_results <- list()

df_outcomes_secondary_long <-
df_outcomes_secondary %>%
  pivot_longer(AfibFlutter:Vent72h) %>%
  rename(Outcome = value,
         Variable = name)


for (outcome in binary_outcomes) {
  message("Running: ", outcome)
  df_tmp <- filter(df_outcomes_secondary_long, Variable == outcome)
  
  tryCatch({
    fit <- brm(data = df_tmp,
               Outcome ~ Exposure + (1 | Hospital),
               prior = c(set_prior("normal(0, 1)", class = "b"),
                         set_prior("cauchy(0, 25)", class = "sd")),
               family = "bernoulli", 
               cores = 4, 
               control = list(adapt_delta = 0.9))
    saveRDS(fit, paste0("Fits/fit_secondary_", gsub(" ", "_", outcome), ".RDS"))
    post <- as_draws_df(fit)
    control_rate <- df_tmp %>% 
      filter(Exposure == "Noreplacement") %>% 
      summarise(control = mean(Outcome)) %>% 
      unlist()
    arr <- meanARR_95CrI(post$b_ExposureReplacement, control_rate)
    secondary_results[[outcome]] <- list(
      OR = mean_95CrI_OR(post$b_ExposureReplacement), 
      ARR = arr,
      control_rate = control_rate
    )
  }, error = function(e) {
    message("Error in outcome: ", outcome)
    message("Error: ", e$message)
  })
}

saveRDS(secondary_results, "Fits/secondary_results.RDS")

secondary_results <- readRDS("Fits/secondary_results.RDS")

# model for death including those with no cardiac rhythms

df_death_composite <-
df_full_outcome %>%
  left_join(
    df_full_baseline %>%
              filter(Characteristic == "Magnesium") %>%
      mutate(Exposure = ifelse(as.numeric(Value > 0.95),
                               "Noreplacement",
                               "Replacement")) %>%
      select(Patient, Exposure),
    by = "Patient"
    ) %>%
  left_join(
    df_full_baseline %>%
      filter(Characteristic == "Hospital") %>%
      mutate(Hospital = Value) %>%
      select(Patient, Hospital),
    by = "Patient"
  ) %>%
  filter(Characteristic == "TimeToDeath_h" |
           Characteristic == "Rhythm") %>%
  mutate(Value = case_when(Characteristic == "TimeToDeath_h" & 
                             Value != "NULL" & as.numeric(Value) <= 24 ~ 1,
                          Characteristic == "Rhythm" &
                            Value %in% c("A-fib",
                                         "A-flutter",
                                         "VF","VT", "SVT") ~ 1,
                          TRUE ~ 0
                            )) %>%
  group_by(Patient, Characteristic, Exposure, Hospital) %>%
  summarise(Value = max(Value)
            ) %>%
  pivot_wider(names_from = "Characteristic",
              values_from = "Value") %>%
  mutate(Composite = TimeToDeath_h == 1 | (!is.na(Rhythm) & Rhythm == 1)) %>%
  rename(Death24h = TimeToDeath_h)


table(df_death_composite$Exposure,
      df_death_composite$Death24h)

fit_death <- brm(Outcome ~ Exposure + (1 | Hospital), data = df_death_composite %>% rename(Outcome = Death24h),
                     prior = c(set_prior("normal(0, 1)", class = "b"), set_prior("cauchy(0, 25)", class = "sd")),
                     family = "bernoulli", cores = 4, control = list(adapt_delta = 0.9))
saveRDS(fit_death, "Fits/fit_death.RDS")

fit_death <- readRDS("Fits/fit_death.RDS")

post_death <- as_draws_df(fit_death)

death_OR <- mean_95CrI_OR(post_death$b_ExposureReplacement)

# Control group event rate
control_rate_death <- df_death_composite %>%
  ungroup() %>%
  filter(Exposure == "Noreplacement") %>% 
  summarise(control = mean(Death24h)) %>%
  unlist()

# ARR using control group rate
death_ARR <- 
  meanARR_95CrI(x = post_death$b_ExposureReplacement,
                rate = control_rate_death)



# Estimate the increment with 2g mag sulfate

mag_increment_df <-
  df_outcomes_secondary %>%
  select(Patient, Exposure, MagLevel) %>%
  left_join(
    filter(
      df_outcome,
      Characteristic == "NamedMedication",
      Value == "magnesium sulfate"),
              by = "Patient") %>%
  filter(!is.na(MagLevel)) %>%
  mutate(MagReplaced = ifelse(is.na(Value), 0, 1)) %>%
  select(-Characteristic, -Value) %>%
  filter(Exposure == "Replacement" & MagReplaced == 1
         | Exposure == "Noreplacement" & MagReplaced == 0) %>%
  left_join(filter(df_baseline,
                   Characteristic == "Magnesium")) %>%
  select(-Characteristic) %>%
  mutate(MagDiff = MagLevel-as.numeric(Value)) %>%
  rename(MagBaseline = Value) %>%
  left_join(filter(df_baseline,
                   Characteristic == "Hospital")) %>%
  select(-Characteristic) %>%
  rename(Hospital = Value)
  

mag_increment_df %>%
  group_by(Exposure) %>%
  summarise(MeanDiff = mean(MagDiff))

fit_magdiff <- brm(
  MagDiff ~ Exposure + (1|Hospital),
  data = mag_increment_df,
  cores = 4,
  prior = c(set_prior("normal(0, 1)", class = "b"),
            set_prior("cauchy(0, 25)", class = "sd")),
  control = list(adapt_delta = 0.9)
)

saveRDS(fit_magdiff, "Fits/fit_magdiff.RDS")

fit_magdiff <- readRDS("Fits/fit_magdiff.RDS")

############# Sensitivity Analyses ####################

# Sensitivity: Composite Outcome #
# use df_death_composite from above to include even those without cardiac rhythms charted

fit_composite <- brm(Outcome ~ Exposure + (1 | Hospital), data = df_death_composite %>% rename(Outcome = Composite),
                     prior = c(set_prior("normal(0, 1)", class = "b"), set_prior("cauchy(0, 25)", class = "sd")),
                     family = "bernoulli", cores = 4, control = list(adapt_delta = 0.9))
saveRDS(fit_composite, "Fits/fit_composite.RDS")

fit_composite <- readRDS("Fits/fit_composite.RDS")

post_composite <- as_draws_df(fit_composite)

composite_OR <- mean_95CrI_OR(post_composite$b_ExposureReplacement)

# Control group event rate
control_rate_composite <- df_death_composite %>% 
  filter(Exposure == "Noreplacement") %>% 
  ungroup() %>%
  summarise(control = mean(Composite)) %>%
  unlist()

# ARR using control group rate
composite_ARR <- 
  meanARR_95CrI(x = post_composite$b_ExposureReplacement,
                rate = control_rate_composite)


# Sensitivity: Random Slope 
fit_randomslope <- brm(Outcome ~ Exposure + (Exposure | Hospital), data = df_brms_primary,
                       prior = c(set_prior("normal(0, 1)", class = "b"), set_prior("cauchy(0, 25)", class = "sd")),
                       family = "bernoulli", cores = 4, control = list(adapt_delta = 0.9))
saveRDS(fit_randomslope, "Fits/fit_randomslope.RDS")

fit_randomslope <- readRDS("Fits/fit_randomslope.RDS")

post_randomslope <- as_draws_df(fit_randomslope)

rs_OR <- mean_95CrI_OR(post_randomslope$b_ExposureReplacement)

# Control group event rate = primary outcome control rate

# ARR using control group rate
rs_ARR <- 
  meanARR_95CrI(x = post_randomslope$b_ExposureReplacement,
                rate = control_rate_primary)

############# Figure: Secondary Outcomes ################################

secondary_results <- readRDS("Fits/secondary_results.RDS")

colnames <- names(purrr::map(secondary_results, 1))

values <- 
  bind_rows(map(secondary_results, 2)) %>%
  bind_cols(bind_rows(map(secondary_results,1)) %>%
                select(P_leq_1)) %>%
  bind_rows(c(primary_ARR,primary_OR[4])) %>%
  bind_rows(c(death_ARR, death_OR[4])) %>%
  bind_rows(c(composite_ARR, composite_OR[4])) %>%
  rename(p_leq_1 = P_leq_1,
         lb95 = `lb95.2.5%`,
         ub95 = `ub95.97.5%`) %>%
  bind_rows(
    tibble(
      mean = NA,
      lb95 = NA,
      ub95 = NA,
    ))


# Prepare data
plot_df_all_secondary_outcomes <- 
  values %>%
  mutate(Outcome = c(colnames, "Afib","Death", "Composite", "Outcome")) %>%
  filter(Outcome != "Death24h") %>%
  mutate(mean = -mean*100,
         tmp = ub95,
         ub95 = -lb95*100,
         lb95 = -tmp*100) %>%
  select(-tmp) %>% 
  filter(Outcome != "MagReplaced") %>%
  mutate(
    Outcome = factor(
      Outcome,
      levels = rev(c("Outcome","Afib", "Death", "RRT72h", "Ett72h", "Tachyarrhythmia", "Composite",  "Ext72h", "Amiodarone")),
      labels = rev(c("Outcome",
                     "A-Fib or A-flutter within 24h", "Death within 24h", "RRT within 72h",
                     "Intubation within 72h", "VF, VT, or SVT within 24h",
                     "Arrhythmia* or death within 24h",
                     "Extubation within 72h", "Amiodarone used within 24h"
      )),
      ordered = TRUE
    ),
    Label = sprintf("%.2f (%.2f to %.2f)", mean, lb95, ub95)) %>%
  mutate(p_leq_1 = as.character(round(p_leq_1,2))) %>%
  mutate(Label = ifelse(Outcome == "Outcome",
                        "ARD (95% CrI)",
                        Label),
         p_leq_1 = ifelse(Outcome == "Outcome",
                          "P(ARD < 0)",
                          p_leq_1)) 

plot_df <- plot_df_all_secondary_outcomes %>%
  filter(!grepl("72h", Outcome)) %>%
  mutate(Outcome = as.character(Outcome)) %>%
  mutate(Outcome = factor(Outcome,
                          levels = 
                            c("Outcome",
                              "A-Fib or A-flutter within 24h", "VF, VT, or SVT within 24h",
                              "Amiodarone used within 24h",
                              "Death within 24h", "Arrhythmia* or death within 24h"
                            ),
                          labels = c(
                            "Outcome",
                              "Atrial fibrillation within 24h", "VF, VT, or SVT within 24h",
                              "Antiarrhythmics used within 24h",
                              "Death within 24h", "Arrhythmia* or death within 24h"
                            ),
                          ordered = T
    )
  )

# Combine header and data, assign manual row numbers
plot_df_with_header <- bind_rows(
  tibble::tibble(
    Outcome = "Outcome",
    Label = "ARD (95% CrI)",
    p_leq_1 = "P(ARD < 0)",
    mean = NA, lb95 = NA, ub95 = NA,
    is_header = TRUE,
    row = 1
  ),
  plot_df %>%
    filter(!(Outcome == "Outcome")) %>%
    mutate(is_header = FALSE) %>%
    arrange(Outcome) %>%
    mutate(row = row_number() + 1)  # push all down by 1
)

ggplot(plot_df_with_header, aes(y = row, x = mean)) +
  geom_vline(xintercept = 0, color = "gray80") +
  
  geom_pointrange(
    data = subset(plot_df_with_header, !is_header),
    aes(xmin = lb95, xmax = ub95),
    size = 0.6, na.rm = TRUE
  ) +
  geom_point(
    data = subset(plot_df_with_header, !is_header),
    size = 2, na.rm = TRUE, color = "grey85"
  ) +
  
  # Outcome column (left-aligned)
  geom_text(aes(x = -25, label = Outcome,
                fontface = ifelse(is_header, "bold", "plain")),
            hjust = 0, size = 3, na.rm = TRUE) +
  
  # OR (95% CrI) column
  geom_text(aes(x = -8, label = Label,
                fontface = ifelse(is_header, "bold", "plain")),
            hjust = 1, size = 3, na.rm = TRUE) +
  
  # P(OR < 1) column
  geom_text(aes(x = 6, label = p_leq_1,
                fontface = ifelse(is_header, "bold", "plain")),
            hjust = 1, size = 3, na.rm = TRUE) +
  
  scale_x_continuous(
    name = "",
    breaks = c(-8, 0, 8),
    limits = c(-25, 10),
    labels = c("<<< Replacement better",
               "ARD = 0",
               "No replacement better >>>")
  ) +
  scale_y_reverse(breaks = NULL) +  # reverse to have row 1 on top
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )


ggsave("Figures/Figure_ARR_PrimarySecondaryOutcomes.svg",
       width = 7, height = 3)

ggsave("Figures/Figure2.jpg",
       width = 7, height = 3, dpi = 300)

############# Subgroup Analyses #######################

# outcomes database with subgroup info 
  # AFib/Flutter
  # Sex, 
  # Age quintile 
  # History of AFib, 
  # Principal problem of acute coronary syndrome
  # Vasopressor use at baseline
  # Hypokalemia  
  # Sepsis
  # CRRT dropped due to small numbers
  
# Make wide baseline dataframe

df_brms_subgroup <-
  # has Patient, Exposure, Hospital, primary outcome
  df_brms_primary %>% 
  left_join(filter(df_baseline,
                   Characteristic %in% c(
                     "R CARDIAC RHYTHM",
                     "Sex",
                     "Age",
                     "Comorbidity",
                     "Principal Problem",
                     "Norepi Equivalents Per Hour",
                     "Potassium")),
                   by = "Patient") %>%
  mutate(
    Afib_Baseline = ifelse(
      Characteristic == "R CARDIAC RHYTHM" &
      (Value == "A-fib" | Value == "A-flutter"),
      1,0),
    Sex_Female = ifelse(Characteristic == "Sex" & Value == "Female",1,0),
    # rough age quintiles: 18-50, 50-65, 65-73, 73-80, >80 
    Age_51to65 = ifelse(Characteristic == "Age" & 
                          as.numeric(Value) > 50 &
                          as.numeric(Value) <= 65,1,0),
    Age_66to73 = ifelse(Characteristic == "Age" & 
                          as.numeric(Value) > 65 &
                          as.numeric(Value) <= 73,1,0),
    Age_74to80 = ifelse(Characteristic == "Age" & 
                          as.numeric(Value) > 73 &
                          as.numeric(Value) <= 80,1,0),
    Age_81orMore = ifelse(Characteristic == "Age" & 
                          as.numeric(Value) > 80,1,0),
    AFib_Hx = ifelse(Characteristic == "Comorbidity" &
                     Value == "AFib/flutter",1,0),
    ACS = ifelse(Characteristic == "Principal Problem" &
                   Value == "Acute coronary syndromes",1,0),
    Sepsis = ifelse(Characteristic == "Principal Problem" &
                      Value == "Sepsis and septic shock",1,0),
    HypoK = ifelse(Characteristic == "Potassium" &
                     as.numeric(Value) < 4.0, 1, 0),
    Norepi = ifelse(Characteristic == "Norepi Equivalents Per Hour" &
                      as.numeric(Value) > 0, 1, 0)
  ) %>%
  # dichotomous structure treats missing as "0"
  mutate(Value = ifelse(is.na(Value),0,Value)) %>%
  select(-Characteristic, -Value) %>%
  group_by(Patient, Exposure, Hospital, Outcome) %>%
  summarise(across(Afib_Baseline:Norepi,~max(.x, na.rm = T))) %>%
  pivot_longer(Afib_Baseline:Norepi) %>%
  mutate(value = ifelse(value < 0, 0, value)) %>%
  rename(Subgroup = name, Covariate = value)
  
# analysis

# age first, separately

df_age <-
  df_brms_subgroup %>%
    filter(grepl("Age",Subgroup)) %>%
    pivot_wider(names_from = Subgroup,
                values_from = Covariate)

fit_age <- brm(
  Outcome ~ Exposure*(Age_51to65 + Age_66to73 + Age_74to80 + Age_81orMore) + (1 | Hospital),
  data = df_age,
  prior =  c(set_prior("normal(0, 1)", class = "b"),
             set_prior("cauchy(0, 25)", class = "sd")),
  family = "bernoulli",
  cores = 4,
  control = list(adapt_delta = 0.9)
)

saveRDS(fit_age, "Fits/fit_age.RDS")

fit_age <- readRDS("Fits/fit_age.RDS")

post_age <- as_draws_df(fit_age)

control_rate_age <- 
  df_age %>%
  ungroup() %>%
  group_by(Age_51to65, Age_66to73, Age_74to80, Age_81orMore) %>%
    filter(Exposure == "Noreplacement") %>%
    summarise(control = mean(Outcome))

age_ARRs <- 
  bind_rows(meanARR_95CrI(x = post_age$b_ExposureReplacement,
                         rate = unlist(control_rate_age[1,"control"])),
            meanARR_95CrI(x = post_age$b_ExposureReplacement + post_age$`b_ExposureReplacement:Age_51to65`,
                          rate = unlist(control_rate_age[5,"control"]))) %>%
    bind_rows(meanARR_95CrI(x = post_age$b_ExposureReplacement + post_age$`b_ExposureReplacement:Age_66to73`,
                            rate = unlist(control_rate_age[4,"control"]))) %>%
    bind_rows(meanARR_95CrI(x = post_age$b_ExposureReplacement + post_age$`b_ExposureReplacement:Age_74to80`,
                            rate = unlist(control_rate_age[3,"control"]))) %>%
    bind_rows(meanARR_95CrI(x = post_age$b_ExposureReplacement + post_age$`b_ExposureReplacement:Age_81orMore`,
                            rate = unlist(control_rate_age[2,"control"]))) %>%
    mutate(covariate = c("<51",
                    "51-65",
                    "66-73",
                    "74-80",
                    ">80")) %>%
    rename(lb95 = `lb95.2.5%`,
           ub95 = `ub95.97.5%`) %>%
    mutate(subgroup = "Age") %>%
    relocate(mean, lb95, ub95, subgroup, covariate)

saveRDS(age_ARRs, "Fits/age_ARRs.rds")


# all the rest, in a loop

subgroups <- 
df_brms_subgroup %>%
  ungroup() %>%
  distinct(Subgroup) %>%
  filter(!grepl("Age",Subgroup)) %>%
  unlist()

subgroup_results <- list()

for (group in subgroups) {
  print(paste0("Running ", group))
  df_tmp <- 
    filter(df_brms_subgroup,
           Subgroup == group)
  fit <- brm(
    Outcome ~ Exposure*Covariate + (1 | Hospital),
    data = df_tmp,
    prior =  c(set_prior("normal(0, 1)", class = "b"),
               set_prior("cauchy(0, 25)", class = "sd")),
    family = "bernoulli",
    cores = 4,
    control = list(adapt_delta = 0.9)
  )
  saveRDS(fit, paste0("Fits/fit_subgroup_", gsub(" ", "_", group), ".RDS"))
  
  print(paste0("Saved ", group))
  
  post <- as_draws_df(fit)
  
  # subgroup cov = 0 control rate
  control_rate_subgroup0 <- df_tmp %>%
    ungroup() %>%
    filter(Exposure == "Noreplacement" &
             Covariate == 0) %>% 
    summarise(control = mean(Outcome)) %>% 
    unlist()

  # subgroup cov = 1 control rate
  control_rate_subgroup1 <- df_tmp %>%
    ungroup() %>%
    filter(Exposure == "Noreplacement" &
             Covariate == 1) %>% 
    summarise(control = mean(Outcome)) %>% 
    unlist()
  
    arr_subgroup0 <- meanARR_95CrI(post$b_ExposureReplacement, control_rate_subgroup0)
    arr_subgroup1 <- meanARR_95CrI(post$b_ExposureReplacement+post$`b_ExposureReplacement:Covariate`,
                                   control_rate_subgroup1)

    subgroup_results[[group]] <- list(
    OR_subgroup0 = mean_95CrI_OR(post$b_ExposureReplacement), 
    OR_subgroup1 = mean_95CrI_OR(post$b_ExposureReplacement + post$`b_ExposureReplacement:Covariate`),
    ARR_subgroup0 = arr_subgroup0,
    ARR_subgroup1 = arr_subgroup1,
    control_rate_subgroup0 = control_rate_subgroup0,
    control_rate_subgroup1 = control_rate_subgroup1
  )
}

saveRDS(subgroup_results, "Fits/subgroup_results.RDS")

subgroup_results <- readRDS("Fits/subgroup_results.RDS")

post <- as_draws_df(readRDS("Fits/fit_subgroup_Afib_Baseline.rds"))
mean(post$`b_ExposureReplacement:Covariate` < 0)

post <- as_draws_df(readRDS("Fits/fit_subgroup_Sex_Female.rds"))
mean(post$`b_ExposureReplacement:Covariate` < 0)

post <- as_draws_df(readRDS("Fits/fit_subgroup_AFib_Hx.rds"))
mean(post$`b_ExposureReplacement:Covariate` < 0)

post <- as_draws_df(readRDS("Fits/fit_subgroup_ACS.rds"))
mean(post$`b_ExposureReplacement:Covariate` < 0)

post <- as_draws_df(readRDS("Fits/fit_subgroup_Sepsis.rds"))
mean(post$`b_ExposureReplacement:Covariate` < 0)

post <- as_draws_df(readRDS("Fits/fit_subgroup_HypoK.rds"))
mean(post$`b_ExposureReplacement:Covariate` < 0)

post <- as_draws_df(readRDS("Fits/fit_subgroup_Norepi.rds"))
mean(post$`b_ExposureReplacement:Covariate` < 0)


# average heart rate in patients with AFib at time 0

df_hr <-
df_outcome %>%
  filter(Characteristic == "HeartRate_24h") %>%
  select(-Characteristic) %>%
  rename(Outcome = Value) %>%
  mutate(Outcome = as.numeric(Outcome)) %>%
  left_join(filter(df_baseline2,
                   Characteristic == "R CARDIAC RHYTHM"),
            by = "Patient") %>%
  rename(BaselineRhythm = Value) %>%
  select(-Characteristic) %>%
  filter(BaselineRhythm %in% c("A-fib",
                               "A-flutter")) %>%
  left_join(filter(df_baseline,
                   Characteristic == "Hospital"), 
            by = "Patient") %>%
  select(-Characteristic, -BaselineRhythm) %>%
  rename(Hospital = Value) %>%
  relocate(Patient, Exposure, Hospital, Outcome)
  
# center and scale HR
hr_mean <- mean(df_hr$Outcome)
hr_sd <- sd(df_hr$Outcome)

df_hr$Outcome <- (df_hr$Outcome-hr_mean)/hr_sd
  
fit_HRifAfib <- brm(
  Outcome ~ Exposure + (1 | Hospital),
  data = df_hr,
  prior =  c(set_prior("normal(0, 1)", class = "b"),
             set_prior("cauchy(0, 25)", class = "sd")),
  family = "gaussian",
  cores = 4,
  control = list(adapt_delta = 0.9)
)

saveRDS(fit_HRifAfib, "Fits/fit_HRifAfib.RDS")

post_hr <- as_draws_df(readRDS("Fits/fit_HRifAfib.RDS"))

control_hr <- df_hr %>%
  filter(Exposure == "Noreplacement") %>%
  summarise(mean(Outcome)) %>%
  unlist()

mean(post_hr$b_ExposureReplacement)*hr_sd
quantile(post_hr$b_ExposureReplacement, 0.025)*hr_sd
quantile(post_hr$b_ExposureReplacement, 0.975)*hr_sd

############ Figure of Subgroup Analyses ########

# load primary analysis and get ARR
fit_primary <- readRDS("Fits/fit_primary.RDS")
post_primary <- as_draws_df(fit_primary)
primary_OR <- mean_95CrI_OR(post_primary$b_ExposureReplacement)

# Control group event rate
control_rate_primary <- df_brms_primary %>% 
  filter(Exposure == "Noreplacement") %>% 
  summarise(control = mean(Outcome)) %>%
  unlist()

# ARR using control group rate
primary_ARR <- 
  meanARR_95CrI(x = post_primary$b_ExposureReplacement,
                rate = control_rate_primary)

# load age subgroup ARRs

age_ARRs <- readRDS("Fits/age_ARRs.rds")

subgroup_results <- readRDS("Fits/subgroup_results.RDS")

subgroups <- names(subgroup_results)


names(subgroup_results[[1]])

data_for_plot <- 
  bind_rows(map(subgroup_results, 3)) %>%
  mutate(subgroup = subgroups,
         covariate = 0) %>%
  bind_rows(
    bind_rows(map(subgroup_results, 4)) %>%
      mutate(subgroup = subgroups,
             covariate = 1)
  ) %>%
  rename(lb95 = `lb95.2.5%`,
         ub95 = `ub95.97.5%`) %>%
  mutate(covariate = ifelse(covariate == 1,
                            "Yes","No")) %>%
  bind_rows(
    age_ARRs 
  ) %>%
  bind_rows(as.data.frame(t(primary_ARR)) %>%
              mutate(subgroup = "Overall",
                     covariate = "") %>%
              rename(lb95 = `lb95.2.5%`,
                     ub95 = `ub95.97.5%`)) %>%  
mutate(across(mean:ub95, ~.x * -100),
         label = paste0(covariate,
                        ifelse(subgroup == "Overall","", ": "),
                        format(round(mean, 1), nsmall = 1),
                        " (",
                        format(round(ub95, 1), nsmall = 1),
                        " to ",
                        format(round(lb95, 1), nsmall = 1),
                        ")"))  %>%
  mutate(covariate = factor(covariate,
                            levels = c("No",
                                       "Yes",
                                       "",
                                       ">80",
                                       "74-80",
                                       "66-73",
                                       "51-65",
                                       "<51"),
                            ordered = T)) %>%
  arrange(desc(mean)) %>%
  mutate(subgroup = factor(
    subgroup,
    levels = rev(c(
      "Afib_Baseline",
      "ACS",
      "HypoK",
      "Sex_Female",
      "Norepi",
      "AFib_Hx",
      "Sepsis",
      "Age",
      "Overall"
    )),
    ordered = T,
    labels = rev(c("Baseline atrial fibrillation",
               "Acute coronary syndrome",
               "Hypokalemia (K < 4)",
               "Female sex",
               "Receiving vasopressor",
               "History of atrial fibrillation",
               "Sepsis or septic shock",
               "Age group (years)",
               "Overall"
               ))
  )) %>%
  mutate(group_colour = case_when(subgroup == "Overall" ~ "black",
                                    subgroup == "Age group (years)" ~ "grey75",
                                  covariate == "Yes" ~ "grey55",
                                  covariate == "No" ~ "grey95"))



data_with_headers <- 
  data_for_plot %>%
  group_by(subgroup) %>%
  arrange(desc(covariate), .by_group = TRUE) %>%
  mutate(row_type = "data") %>%
  ungroup() %>%
  group_split(subgroup) %>%
  purrr::imap_dfr(function(df, i) {
    subgroup_name <- as.character(df$subgroup[1])
    
    # Header row
    header_row <- df[1, ]
    header_row[] <- NA
    header_row$label <- subgroup_name
    header_row$row_type <- "header"
    
    # Spacer row
    spacer_row <- df[1, ]
    spacer_row[] <- NA
    spacer_row$label <- ""
    spacer_row$row_type <- "spacer"
    
    bind_rows(header_row, df, spacer_row)
  }) %>%
  mutate(row = dplyr::row_number()) %>%
  mutate(row = max(row) - row + 1)  # reverse row order for coord_flip()



ggplot(data_with_headers, aes(y = row, x = mean)) +
  # Reference line at ARR = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  
  # CI bars and points only for data rows
  geom_pointrange(
    data = subset(data_with_headers, row_type == "data"),
    aes(xmin = ub95, xmax = lb95),
    size = 0.3, na.rm = TRUE
  ) +
  geom_point(
    data = subset(data_with_headers, row_type == "data"),
    aes(color = group_colour),
    size = 1, na.rm = TRUE
  ) +
  
  # Text labels
  geom_text(
    data = subset(data_with_headers, row_type == "header"),
    aes(x = -36, label = label),  # left-align header
    hjust = 0, fontface = "bold", size = 3.5, na.rm = TRUE
  ) +
  geom_text(
    data = subset(data_with_headers, row_type == "data"),
    aes(label = label),
    x = -18, hjust = 1, size = 3, na.rm = TRUE
  ) +
  
  scale_color_manual(
    values = c("grey55" = "grey55",
               "grey75" = "grey75", 
               "black" = "black", 
               "grey95" = "grey95"),
    guide = "none"
  ) +
  scale_x_continuous(
    name = "Absolute Risk Difference (%)",
    limits = c(-40, 15),
    breaks = c(-10,0,10),
    labels = c("<<< Replacement better",
               "0",
               "No replacement better >>>")
  ) +
  scale_y_continuous(
    name = "",
    breaks = NULL  # hide tick labels
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5.5, 40, 5.5, 5.5)
  ) +
  labs(
    title = "Absolute difference in 24-hr risk of A-Fib"  ) 


ggsave("Figures/Figure_ARR.svg",
       width = 6, height = 6)

ggsave("Figures/Figure_ARR.jpeg",
       width = 6, height = 6,
       dpi = 300)

############ Multiple Imputation for HTE, CATE, adjusted analyses ###############

# Step 1: wide-format baseline dataframe

characteristics_ordered <- c("Age",
                             "TotalCalcium",
                             "CorrectedCalcium",
                             "IonizedCalcium",
                             "Hospital Length of Stay (Days)",
                             "ICU Length of Stay (Hours)",
                             "Phosphate",
                             "PULSE",
                             "RESPIRATIONS",
                             "TEMPERATURE",
                             "SBP",
                             "Albumin",
                             "Bicarbonate",
                             "Hemoglobin",
                             "Norepi Equivalents Per Hour",
                             "Lactate",
                             "Potassium",
                             "WEIGHT/SCALE",
                             "DBP",
                             "Chloride",
                             "Creatinine",
                             "HEIGHT",
                             "Sodium",
                             "UrineOutput_24h",
                             "Glucose")

characteristics_categorical <- c("Preferred Language",
                                 "Comorbidity",
                                 "ICU Admission Month",
                                 "Principal Problem",
                                 #"BaselineAFib/flutter", # treat this separate because has missing data
                                 "Sex",
                                 "Intubated",
                                 "CRRT")

df_categorical <- 
  df_baseline2 %>%
  filter(Characteristic %in% characteristics_categorical) %>%
  mutate(OneHot = gsub(" ", "_", paste(Characteristic, Value, sep = "_"))) %>%
  distinct(Patient, OneHot) %>% # all duplicates are duplicate intubated variables
  mutate(Value = 1.0) %>%
  pivot_wider(names_from = OneHot,
              values_from = Value,
              values_fill = 0) %>% 
#values_fill = 0 because all variables above have
# a default value if missing except age, 
# which has no missing values
# however, baseline rhythm has about 5% missing, so:
left_join(df_baseline %>%
            # convert baseline cardiac rhythm into Afib/flutter
            filter(Characteristic == "R CARDIAC RHYTHM") %>%
            mutate(Characteristic = "BaselineAfib/flutter",
                   Value = ifelse(Value %in% c("A-fib","A-flutter"),1,0)) %>%
            pivot_wider(names_from = Characteristic, 
                        values_from = Value),
          by = "Patient"
)

df_ordinal <- 
  df_baseline2 %>%
  filter(Characteristic %in% characteristics_ordered) %>%
  mutate(Value = as.numeric(Value),
         Characteristic = gsub(" ", "_", Characteristic)) %>%
  select(-Exposure) %>%
  pivot_wider(names_from = "Characteristic",
              values_from = "Value")

df_risk <- 
  df_categorical %>%
  left_join(df_ordinal,
            by = "Patient") %>%
  left_join(primary_outcome,
            by = "Patient") %>% #primary outcome comes from primary analysis section
  relocate(Patient, Outcome)


## missingness table
missingness_table <- 
  df_risk %>%
  left_join(
    distinct(df_baseline2, Patient, Exposure),
    by = "Patient"
  ) %>%
  pivot_longer(cols = -c(Patient, Exposure), names_to = "Variable", values_to = "Value") %>%
  mutate(Missing = is.na(Value)) %>%
  group_by(Variable, Exposure) %>%
  summarise(count = sum(Missing), 
            prop = mean(Missing), 
            .groups = "drop") %>%
  mutate(label = paste0(count, " (", round(100 * prop), "%)")) %>%
  select(Variable, Exposure, label) %>%
  pivot_wider(names_from = Exposure, values_from = label) %>%
  left_join(
    df_risk %>%
      pivot_longer(cols = -Patient, names_to = "Variable", values_to = "Value") %>%
      mutate(Missing = is.na(Value)) %>%
      group_by(Variable) %>%
      summarise(count = sum(Missing),
                prop = mean(Missing),
                .groups = "drop") %>%
      mutate(Total = paste0(count, " (", round(100 * prop), "%)")) %>%
      select(Variable, Total),
    by = "Variable"
  ) %>%
  select(Variable, Total, `Noreplacement`, Replacement)

write_csv(missingness_table, "Tables/missingness.csv")


# HEIGHT and WEIGHT have too much missingness so we remove them
# the Calcium values also have high missingness for corrected
# evidence for formulas connecting these is weak, so we will drop them

df_risk_for_mice <- 
  df_risk %>%
    select(-HEIGHT, -`WEIGHT/SCALE`,
           -CorrectedCalcium) %>%
# Step 2: Remove identifier columns before imputation
    select(-Patient)  

# Step 3: Run MICE
# ignore outcome so that cross-validation is not contaminated

imp <- mice(select(df_risk_for_mice, -Outcome), 
            m = 10, method = "pmm", maxit = 5, seed = 123)
            

############# Risk-based HTE with 5-fold CV ######################

all_preds <- vector("list", length = 10)
metrics_list <- vector("list", length = 10)

for (i in 1:10) {
  message("Imputation dataset ", i)
  imp_data <- complete(imp, i)
  imp_data$Patient <- df_risk$Patient
  imp_data$Outcome <- df_risk$Outcome
  
  # 5-fold CV
  folds <- vfold_cv(imp_data, v = 5)
  
  fold_preds <- map_dfr(folds$splits, function(split) {
    train <- analysis(split)
    test <- assessment(split)
    
    test_patients <- test$Patient
    test_truth <- test$Outcome  # for metric calc
    
    train <- train %>% select(-Patient)
    test <- test %>% select(-Patient)
    
    rec <- recipe(Outcome ~ ., data = train) %>%
      step_dummy(all_nominal_predictors(), -all_outcomes()) %>%
      prep()
    
    X_train <- bake(rec, new_data = train) %>%
      select(-Outcome) %>%
      mutate(across(everything(), as.numeric))
    y_train <- train$Outcome
    
    model <- xgboost(data = as.matrix(X_train), label = y_train,
                     objective = "binary:logistic", nrounds = 100, verbose = 0)
    
    X_test <- bake(rec, new_data = test) %>%
      select(-Outcome) %>%
      mutate(across(everything(), as.numeric))
    
    preds <- predict(model, as.matrix(X_test))
    
    tibble(Patient = test_patients, pred = preds, truth = test_truth)
  })
  
  # Save predictions
  all_preds[[i]] <- fold_preds %>% select(Patient, pred)
  
  # === METRICS ===
  
  # AUC
  auc_val <- roc(fold_preds$truth, fold_preds$pred)$auc
  
  # Calibration: decile binning
  cal_df <- fold_preds %>%
    mutate(decile = ntile(pred, 10)) %>%
    group_by(decile) %>%
    summarise(
      mean_pred = mean(pred, na.rm = TRUE),
      obs_rate = mean(truth, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  metrics_list[[i]] <- list(
    auc = auc_val,
    calibration = cal_df
  )
}

saveRDS(metrics_list, "Fits/metrics_list_riskhte.RDS")
saveRDS(all_preds, "Fits/allpreds_riskhte.RDS")

metrics_list <- readRDS("Fits/metrics_list_riskhte.RDS")

# Average AUC
avg_auc <- mean(sapply(metrics_list, function(x) x$auc))
print(avg_auc)

# Combine calibration data from all 10 models
calibration_df <- bind_rows(
  lapply(seq_along(metrics_list), function(i) {
    metrics_list[[i]]$calibration %>%
      mutate(model = i)
  })
)

# Plot all 100 points
ggplot(calibration_df, aes(x = mean_pred, y = obs_rate)) +
  geom_point(alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted Risk", y = "Observed Event Rate",
       title = "Calibration of A-Fib Risk Model") +
  theme_minimal()

ggsave("Figures/riskhte_calibration.svg",
       width = 5, height = 5)

all_preds <- readRDS("Fits/allpreds_riskhte.RDS")

# Step 4: Average predictions across imputations
risk_preds_combined <- 
  bind_rows(all_preds, .id = "imputation") %>%
  group_by(Patient) %>%
  summarise(mean_pred = mean(pred), .groups = "drop") %>% 
  mutate(RiskGroup = ntile(mean_pred, 10)) %>%
  left_join(primary_outcome,
            by = "Patient") 

df_riskhte <- df_brms_primary %>% left_join(risk_preds_combined %>% select(Patient, RiskGroup), by = "Patient")

fit_riskhte <- brm(Outcome ~ Exposure * as.factor(RiskGroup) + (1 | Hospital), data = df_riskhte,
  prior = c(set_prior("normal(0, 1)", class = "b"), set_prior("cauchy(0, 25)", class = "sd")),
  family = "bernoulli", cores = 4, control = list(adapt_delta = 0.9))

saveRDS(fit_riskhte, "Fits/fit_riskhte.RDS")

fit_riskhte <- readRDS("Fits/fit_riskhte.rds")

post_riskhte <- as_draws_df(fit_riskhte)

bind_rows(
  list(Group1 = mean_95CrI_OR(post_riskhte$b_ExposureReplacement)
    ,Group2 = mean_95CrI_OR(post_riskhte$b_ExposureReplacement + post_riskhte$`b_ExposureReplacement:as.factorRiskGroup2`)
    ,Group3 = mean_95CrI_OR(post_riskhte$b_ExposureReplacement + post_riskhte$`b_ExposureReplacement:as.factorRiskGroup3`)
    ,Group4 = mean_95CrI_OR(post_riskhte$b_ExposureReplacement + post_riskhte$`b_ExposureReplacement:as.factorRiskGroup4`)
    ,Group5 = mean_95CrI_OR(post_riskhte$b_ExposureReplacement + post_riskhte$`b_ExposureReplacement:as.factorRiskGroup5`)
    ,Group6 = mean_95CrI_OR(post_riskhte$b_ExposureReplacement + post_riskhte$`b_ExposureReplacement:as.factorRiskGroup6`)
    ,Group7 = mean_95CrI_OR(post_riskhte$b_ExposureReplacement + post_riskhte$`b_ExposureReplacement:as.factorRiskGroup7`)
    ,Group8 = mean_95CrI_OR(post_riskhte$b_ExposureReplacement + post_riskhte$`b_ExposureReplacement:as.factorRiskGroup8`)
    ,Group9 = mean_95CrI_OR(post_riskhte$b_ExposureReplacement + post_riskhte$`b_ExposureReplacement:as.factorRiskGroup9`)
    ,Group10 = mean_95CrI_OR(post_riskhte$b_ExposureReplacement + post_riskhte$`b_ExposureReplacement:as.factorRiskGroup10`)
), .id = "Group") %>%
  rename(lb95 = `lb95.2.5%`,
         ub95 = `ub95.97.5%`) %>%
  mutate(Group = factor(
    Group,
    levels = c(
      "Group10"
    , "Group9"
    , "Group8"
    , "Group7"
    , "Group6"
    , "Group5"
    , "Group4"
    , "Group3"
    , "Group2"
    , "Group1"
    ),
    labels = c(
      "Highest risk"
      , "Group 9"
      , "Group 8"
      , "Group 7"
      , "Group 6"
      , "Group 5"
      , "Group 4"
      , "Group 3"
      , "Group 2"
      , "Lowest risk"
    ),
    ordered = T),
    label = paste0(
      format(round(mean,2),2)
    , " ("
    , format(round(lb95,2),2)
    , " to "
    , format(round(ub95,2),nsmall = 2)
    , ")")) %>%
  ggplot(aes(x = Group,
             y = mean,
             ymin = lb95,
             ymax = ub95)) +
  geom_hline(yintercept = 1) +
  geom_pointrange() +
  theme_minimal() +
  scale_y_continuous(trans= "log",
                     breaks = c(0.35, 1, 3),
                     limits = c(0.01,4),
                     labels = c(
                       "<<< Replacement better",
                       "OR = 1.0",
                       "No replacement better >>>"
                     )) +
  theme(panel.grid = element_blank()) + 
  coord_flip() +
  geom_text(aes(label = label),
            y = log(0.05),
            size = 3) +
  labs(x = "",
       y = "",
       title = "Risk of atrial fibrillation over next 24 hours",
       subtitle = "By predicted risk decile of subsequent atrial fibrillation")

ggsave("Figures/benefit_riskhte.svg",
       width = 8, height = 5)

# just look at interaction between predicted risk and outcome


df_riskhte_continuous <- 
df_brms_primary %>%
  left_join(
    risk_preds_combined %>%
      select(Patient, mean_pred) %>%
      rename(pred_risk = mean_pred),
            by = "Patient")
  
fit_riskhte_continuous <- 
  brm(data = df_riskhte_continuous,
      sample_prior = "yes",
      Outcome ~ Exposure * pred_risk + (1 | Hospital),
      prior = c(set_prior("normal(0, 1)", class = "b"), set_prior("cauchy(0, 25)", class = "sd")),
      family = "bernoulli", cores = 4, control = list(adapt_delta = 0.9))

saveRDS(fit_riskhte_continuous, "Fits/fit_riskhte_continuous.RDS")

fit_riskhte_continuous <- readRDS("Fits/fit_riskhte_continuous.RDS")

hypothesis(fit_riskhte_continuous,
           "ExposureReplacement:pred_risk = 0")

# Extract posterior draws
draws <- as_draws_df(fit_riskhte_continuous)

mean(draws$`b_ExposureReplacement:pred_risk`>0)

# Adjust variable names as needed
draws <- draws %>%
  rename(
    beta_rep = `b_ExposureReplacement`,
    beta_interaction = `b_ExposureReplacement:pred_risk`
  )

# Create a grid of pred_risk values
pred_grid <- tibble(pred_risk = seq(0, 1, length.out = 100))

# Compute the exponentiated linear combination for each draw and each pred_risk
plot_data <- pred_grid %>%
  crossing(draw = 1:nrow(draws)) %>%
  mutate(
    beta_rep = draws$beta_rep[draw],
    beta_interaction = draws$beta_interaction[draw],
    log_or = beta_rep + beta_interaction * pred_risk,
    or = exp(log_or)
  )

# Summarize across posterior draws
plot_summary <- plot_data %>%
  group_by(pred_risk) %>%
  summarise(
    median = median(or),
    lower = quantile(or, 0.025),
    upper = quantile(or, 0.975)
  )

# Plot
ggplot(plot_summary, aes(x = pred_risk, y = median)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_y_continuous(trans = "log",
                     breaks = c(0.5, 1, 2)) +
  labs(
    x = "Baseline predicted risk of atrial fibrillation over next 24h",
    y = "Odds Ratio for Replacement vs No Replacement",
    title = "Effect Modification by Baseline Risk"
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

ggsave("Figures/RiskHTE_Continuous.svg",
       width = 6, height = 5)

############# Effect-based HTE with 100x CV and SHAP ######################

# run mice again, but include assignment
df_risk_for_mice2 <- 
  df_risk %>%
  left_join(
    distinct(df_baseline2,
           Patient, Exposure),
    by = "Patient"
  ) %>%
  # 1 = replacement, 0 = no replacement
  mutate(treat = ifelse(Exposure == "Replacement", 1, 0))%>%
  select(-HEIGHT, -`WEIGHT/SCALE`,
         -CorrectedCalcium,
         -Exposure) %>%
  # Step 2: Remove identifier and outcome columns before imputation
  select(-Patient, -Outcome)  

# Step 3: Run MICE
# ignore outcome so that cross-validation is not contaminated
# note that the original plan (imputation within the training splits)
# was not feasible because of computation and low-prevalence categorical
# variables leading to singular training data matrices

imp2 <- mice(df_risk_for_mice2, 
            m = 10, method = "pmm", maxit = 5, seed = 123)

df_forest_raw <- 
  primary_outcome %>%
  left_join(distinct(df_baseline2,
                   Patient, Exposure),
            by = "Patient") %>%
  # 1 = replacement, 0 = no replacement
    mutate(treat = ifelse(Exposure == "Replacement", 1, 0))%>%
  select(treat, Outcome)

# - imp2: imputed object of 10 datasets (covariates + Exposure)
# - df_forest_raw: has Outcome and treat (Exposure), but not Patient

results <- vector("list", length = 100)
patient_cate_list <- vector("list", length = 100)

set.seed(123)

for (i in 1:100) {
  message("Iteration ", i)
  
  # 1. Create train/test split
  idx <- sample(1:nrow(df_forest_raw), floor(2/3 * nrow(df_forest_raw)))
  train_raw <- df_forest_raw[idx, ]
  test_raw  <- df_forest_raw[-idx, ]
  
  # 2. Subset imputed covariates to match train/test rows
  train_list <- map(1:10, ~ complete(imp2, .x)[idx, ])
  test_list  <- map(1:10, ~ complete(imp2, .x)[-idx, ])
  
  # 3. Add Outcome and treat from raw data
  train_list <- map(train_list, ~ mutate(.x,
                                         Outcome = train_raw$Outcome,
                                         treat = train_raw$treat))
  test_list  <- map(test_list,  ~ mutate(.x,
                                         Outcome = test_raw$Outcome,
                                         treat = test_raw$treat))
  
  # 4. Fit causal forest on training sets
  forest_list <- map(train_list, ~ causal_forest(
    X = as.matrix(select(.x, -Outcome, -treat)),
    Y = .x$Outcome,
    W = .x$treat
  ))
  
  # 5. Predict CATEs on test sets
  preds <- map2(test_list, forest_list, function(test, model) {
    predict(model, as.matrix(select(test, -Outcome, -treat)))$predictions
  })
  
  # 6. Average predicted treatment effects
  avg_preds <- reduce(preds, `+`) / length(preds)
  
  # Track which patients are in this test set
  test_patients <- df_risk$Patient[-idx]  # match order of test_raw
  
  # Store CATEs for this fold with patient IDs
  patient_cate_list[[i]] <- tibble(
    iteration = i,
    Patient = test_patients,
    predicted_cate = avg_preds
  )
  
  # 7. Create deciles and summarize observed effect
  decile <- ntile(avg_preds, 10)
  test_outcomes <- test_raw$Outcome
  test_treats   <- test_raw$treat
  
  observed <- map_dfr(1:10, function(d) {
    idxs <- which(decile == d)
    tibble(
      decile = d,
      treat_mean = mean(test_outcomes[idxs][test_treats[idxs] == 1], na.rm = TRUE),
      control_mean = mean(test_outcomes[idxs][test_treats[idxs] == 0], na.rm = TRUE)
    )
  }) %>% mutate(diff = treat_mean - control_mean)
  
  # 8. Compute expected benefit if treatment were assigned by predicted CATE
  gain <- ifelse(
    (test_treats == 0 & avg_preds < 0) | (test_treats == 1 & avg_preds > 0),
    abs(avg_preds),
    0
  )
  
  expected_gain <- mean(gain, na.rm = TRUE)
  # calculate absolute risk difference quantiles,
  # and retain the median and 95% CI ARD within each quantile
  
  # Store both decile-wise results and policy gain
  results[[i]] <- list(
    diff_by_decile = observed$diff,
    expected_gain = expected_gain
  )
  
}


saveRDS(results, "Fits/effecthte_results.RDS")

results <- readRDS("Fits/effecthte_results.RDS")


all_patient_cates <- bind_rows(patient_cate_list)

saveRDS(all_patient_cates, "Fits/all_patient_cates.RDS")

all_patient_cates <- readRDS("Fits/all_patient_cates.RDS")

# Best Linear Predictor of heterogeneity
# this is the coefficient gamma in regression
# Y = alpha + gamma * cate_i * treat_i
# where Y is outcome, alpha is ate, 
# cate_i is CENTERED cate for patient i,
# (centered to avoid collinearity)
# (scaled because log-odds space, and to help fitting)
# and treat_i is treatment assignment for patient i
# calculate by iteration
# pool posterior distributions

# find ATE to center CATEs
df_ate <- 
  all_patient_cates %>%
  group_by(iteration) %>%
  summarise(ate = mean(predicted_cate))

# blp dataframe
df_blp <- 
  all_patient_cates %>% 
    left_join(df_baseline2 %>% 
                distinct(Patient, Exposure), 
              by= "Patient") %>%
    left_join(df_risk %>%
                select(Patient, Outcome),
              by = "Patient") %>%
    left_join(df_ate, 
              by = "iteration") %>%
  mutate(treat = ifelse(Exposure == "Replacement",1,0)) %>%
  mutate(cate_centered = predicted_cate - ate) %>%
  mutate(cate_z = cate_centered/(sd(cate_centered)))
  

# BAYESIAN version to calculate probability of BLP beta2 > 0

# Use first iteration's data to compile model
init_data <- df_blp %>% filter(iteration == min(iteration))

brm_model <- brm(
  formula = Outcome ~ treat + treat:cate_z,
  data = init_data,
  family = bernoulli(),
  chains = 2, iter = 1000, warmup = 500,
  refresh = 0,
  prior = prior(normal(0, 5), class = "b"),  # weak prior
  seed = 123
)


posterior_list <- df_blp %>%
  group_split(iteration) %>%
  map2_dfr(
    .x = .,
    .y = map_int(., ~ unique(.x$iteration)),
    .f = function(df_iter, iter_id) {
      fit_iter <- update(
        brm_model,
        newdata = df_iter,
        recompile = FALSE,
        refresh = 0,
        chains = 2, iter = 1000, warmup = 500,
        seed = 123 + iter_id
      )
      
      as_draws_df(fit_iter, variable = c("b_treat","b_treat:cate_z")) %>%
        transmute(
          iteration = iter_id,
          beta1 = b_treat,
          beta2 = `b_treat:cate_z`
        )
    }
  )

saveRDS(posterior_list, "fits/blp_bayesian.rds")

posterior_list <- readRDS("fits/blp_bayesian.rds")

# plot

# 1. Create grid of beta2 values (z-scores)
z_seq <- seq(-2, 2, length.out = 200)

# 2. For each z, calculate OR = exp(beta1 + beta2 * z)
# Across all posterior samples

# Reframe version — much faster and avoids dplyr warnings
or_posterior <- map_dfr(z_seq, function(z) {
  logit <- posterior_list$beta1 + posterior_list$beta2 * z
  tibble(
    z = z,
    OR_median = median(exp(logit)),
    OR_low = quantile(exp(logit), 0.025),
    OR_high = quantile(exp(logit), 0.975)
  )
})

ate <- data.frame(mean = exp(mean(posterior_list$beta1)),
                  lb95 = exp(quantile(posterior_list$beta1, 0.025)),
                  ub95 = exp(quantile(posterior_list$beta1, 0.975)))


p1 <- 
  ggplot(or_posterior, aes(x = z, y = OR_median)) +
  # Main line and ribbon
  geom_line(color = "steelblue") +
  geom_ribbon(aes(ymin = OR_low, ymax = OR_high), fill = "steelblue", alpha = 0.2) +
  
  # ATE ribbon and line
  annotate("rect", xmin = -2, xmax = 2,
           ymin = ate$lb95, 
           ymax = ate$ub95,
           fill = "red", alpha = 0.1) +  # faint red ribbon
  
  geom_hline(yintercept = ate$mean, 
             color = "red", linewidth = 0.8) +
  annotate("text", x = 1, y = exp(mean(posterior_list$beta1)), 
           label = paste0("ATE = ",
                          format(round(ate$mean,2), nsmall = 2),
                          " (",
                          format(round(ate$lb95,2), nsmall = 2),
                          " to ",
                          format(round(ate$ub95,2), nsmall = 2),
                          ")"
                          ), color = "red", vjust = 1.5) +
  
  # Axis and theme
  scale_y_continuous(
    trans = "log",
    breaks = c(0.25, 0.5, 1.0, 2, 4),
    labels = c("0.25", "0.5", "1", "2", "4")
  ) +
  labs(
    x = "",
    y = "Odds Ratio"
  ) +
  geom_text(data = data.frame(x = 0, 
                              y = 2,
                              label = paste0(
                                "P(coefficient for interaction > 0) = ",
                              round(mean(posterior_list$beta2>0),2))),
            aes(x = x, y = y, label = label),
            color = "steelblue") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

# 4. Histogram of beta2 samples
p2 <- 
  ggplot(df_blp, aes(x = cate_z)) +
  geom_histogram(fill = "steelblue", adjust = 1.2, color = NA,
               alpha = 0.5) +
  coord_cartesian(xlim = c(-2, 2)) +
  labs(
    x = "Standardized CATE",
    y = "Observed frequency"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.grid.minor = element_blank())

# 5. Combine with shared x-axis using patchwork
p1 / p2 + plot_layout(heights = c(2, 1))

ggsave("Figures/blp_bayesian.svg",
       height = 6, width = 6)

  # medians of quantiles with 95% intervals 
  # again by iteration
  
# quantile_stats <-

quantile_df <- 
  all_patient_cates %>% 
  left_join(df_baseline2 %>% 
              distinct(Patient, Exposure), 
            by= "Patient") %>%
  left_join(df_risk %>%
              select(Patient, Outcome),
            by = "Patient") %>%
  group_by(iteration) %>%
  mutate(quintile = ntile(predicted_cate, 3)) %>%
  mutate(treat = ifelse(Exposure == "Replacement",1,0)) %>%
  select(-Exposure) %>%
  group_by(iteration, treat, quintile) %>%
  summarise(mean = mean(Outcome),
            n = n()) %>%
  pivot_wider(names_from = treat,
              values_from = c(mean, n)) %>%
  mutate(treatment_effect = mean_1 - mean_0,
         se = sqrt(
           mean_1*(1-mean_1)/n_1 + 
             mean_0*(1-mean_0)/n_0),
         conf.low = treatment_effect - 1.96 * se,
         conf.high = treatment_effect + 1.96 * se,
         mean_control = mean_0
  )

quantile_df %>%
  group_by(quintile) %>%
  summarise(across(.cols = c(treatment_effect, conf.low, conf.high),
                   median)) %>%
  mutate(label = paste0(
    format(round(treatment_effect*100,1), nsmall = 1),
    " (",
    format(round(conf.low*100,1), nsmall = 1),
    " to ",
    format(round(conf.high*100,1), nsmall = 1),
    ")"
  ),
  quintile_label = paste0("Q", quintile)
) %>%
  ggplot(aes(x = factor(quintile_label), y = treatment_effect)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(-0.25,0.25)) +
  geom_text(aes(label = label),
            y = 0.25,
            hjust = 1, size = 3.5) +
  labs(
    title = "Treatment Effect by Predicted CATE Tertile",
    x = "CATE Tertile",
    y = "Treatment Effect (Risk Difference)"
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave("Figures/CATEterttiles.svg",
       width = 7, height = 5)

# Qini score function

compute_qini_score <- function(df_iter) {
  df <- df_iter[order(df_iter$predicted_cate), ] # benefit is negative
  n <- nrow(df)
  
  df$treated_y <- df$Outcome * df$treat
  df$control_y <- df$Outcome * (1 - df$treat)
  
  df$cum_treated <- cumsum(df$treat)
  df$cum_control <- cumsum(1 - df$treat)
  
  df$cum_y_treated <- cumsum(df$treated_y)
  df$cum_y_control <- cumsum(df$control_y)
  
  df$rate_treated <- ifelse(df$cum_treated > 0, df$cum_y_treated / df$cum_treated, 0)
  df$rate_control <- ifelse(df$cum_control > 0, df$cum_y_control / df$cum_control, 0)
  
  df$uplift <- df$rate_treated - df$rate_control
  df$cumulative_uplift <- cumsum(df$uplift)
  df$percent <- seq_len(n) / n
  
  # Area under Qini curve (trapezoidal rule)
  area_model <- sum(diff(df$percent) * head(df$cumulative_uplift, -1))
  
  # Random line baseline
  total_uplift <- tail(df$cumulative_uplift, 1)
  area_random <- 0.5 * total_uplift
  
  qini_score <- area_model - area_random
  return(qini_score)
}

# C-for-benefit function

compute_c_for_benefit <- function(df_iter) {
  y <- df_iter$Outcome
  w <- df_iter$treat
  tau <- -df_iter$predicted_cate # more negative cate suggests more benefit
  
  n <- length(y)
  concordant <- 0
  tied <- 0
  discordant <- 0
  total <- 0
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (w[i] != w[j]) {
        obs_diff <- (y[i] - y[j]) * (w[i] - w[j])
        pred_diff <- tau[i] - tau[j]
        
        score <- obs_diff * pred_diff
        
        if (score > 0) {
          concordant <- concordant + 1
        } else if (score == 0) {
          tied <- tied + 1
        } else {
          discordant <- discordant + 1
        }
        
        total <- total + 1
      }
    }
  }
  
  cfb <- if (total == 0) NA else (concordant + 0.5 * tied) / total
  
  return(list(
    c_for_benefit = cfb,
    concordant = concordant,
    tied = tied,
    discordant = discordant,
    total_pairs = total
  ))
}

# Qini and c for benefit by iteration


cfb_results <- df_blp %>%
  group_by(iteration) %>%
  group_map(~ {
    res <- compute_c_for_benefit(.x)
    tibble(
      iteration = .y$iteration,
      c_for_benefit = res$c_for_benefit,
      concordant = res$concordant,
      tied = res$tied,
      discordant = res$discordant,
      total_pairs = res$total_pairs
    )
  }) %>%
  bind_rows()




qini_results <- 
  df_blp %>%
  group_by(iteration) %>%
  group_map(~ tibble(
    iteration = .y$iteration,
    qini_score = compute_qini_score(.x) 
  )) %>%
  bind_rows()

# View summary
summary(qini_results)
summary(cfb_results)


# Vector of expected gains
expected_gain_vec <- sapply(results, `[[`, "expected_gain")

# Summary of expected gain
median(expected_gain_vec)
quantile(expected_gain_vec, probs = c(0.025, 0.975))


# Wide to long
results_long <- 
  df %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(-iteration, names_to = "decile", values_to = "observed_diff") %>%
  mutate(decile = as.integer(gsub("V", "", decile)))

# "The observed ARD in each quantile is estimated as the median of 
# estimated ARDs across the 100 iterations, and 90% confidence bounds 
# are calculated as the medians of the one hundred 95% confidence bounds."



# Plot
ggplot(results_long, aes(x = decile, y = observed_diff, group = iteration)) +
  geom_line(alpha = 0.2, color = "gray") +
#  stat_summary(fun = mean, geom = "line", color = "grey90", size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Observed Treatment Effect by Decile of Predicted CATE",
       x = "Decile of Predicted CATE",
       y = "Observed Treatment Effect (Treated - Control)") +
  theme_minimal()


# Shap values
cf_final <- causal_forest(as.matrix(complete(imp, 1)),
                          df_forest_raw$Outcome,
                          df_forest_raw$treat)


saveRDS(cf_final, "Fits/fit_causalforest_final.RDS")

cf_final <- readRDS("Fits/fit_causalforest_final.RDS")



# Step 1: Use the original covariates (must be a data frame)
X <- as.data.frame(cf_final$X.orig)

# Step 2: Define the prediction function (returns CATEs)
predict_cate <- function(object, newdata) {
  predict(object, newdata = as.matrix(newdata))$predictions
}

# Step 3: Estimate SHAP values using fastshap::explain()
# Use a small sample first (e.g. 100 obs, 10 features) to test

N_explain <- 2053

set.seed(123)
X_sample <- X[1:N_explain, ]  # or however many you want to explain

shap_values <- explain(
  object = cf_final,
  X = X_sample,
  pred_wrapper = predict_cate,
  nsim = 100,        # number of permutations for approximation
  adjust = TRUE,     # use marginal adjustment
  parallel = FALSE  # set to TRUE if you want to parallelize
)

saveRDS(shap_values,
        "Fits/shap_values.RDS")

shap_values <- readRDS("Fits/shap_values.RDS")

# Step 4: Visualize SHAP values

# Summarize mean absolute SHAP value per feature
shap_summary <- as.data.frame(shap_values) %>%
  summarise(across(everything(), ~ mean(abs(.), na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "mean_abs_shap") %>%
  arrange(desc(mean_abs_shap))

# Plot
shap_summary %>%
  arrange(desc(mean_abs_shap)) %>%
  mutate(row = row_number()) %>%
  filter(row <= 25) %>%
ggplot(aes(x = reorder(variable, mean_abs_shap), y = mean_abs_shap)) +
  geom_col(fill = green) +
  coord_flip() +
  labs(
    title = "Mean Absolute SHAP Value per Variable",
    x = "Variable (top 25)",
    y = "Mean |SHAP|"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

ggsave("Figures/shap_global.svg",
       width = 8, height = 7)

# 1. Get the top 10 features by mean absolute SHAP value
top_features <- shap_summary %>%
  arrange(desc(mean_abs_shap)) %>%
  slice_head(n = 10) %>%
  pull(variable)

custom_limits <- list(
  Creatinine = list(x = c(0, 750)),
  Age = list(x = c(18,90)),
  Phosphate = list(x = c(0, 4)),
  IonizedCalcium = list(x = c(0.7, 1.7)),
  Norepi_Equivalents_Per_Hour = list(x = c(0,3000)),
  Glucose = list(x = c(0, 35))
  # Add more features as needed
)

# 2. Plotting function
plot_shap_dependence <- function(feature_name, shap_values, X_sample, shap_summary, custom_limits = NULL) {
  mean_shap <- shap_summary %>%
    filter(variable == feature_name) %>%
    pull(mean_abs_shap)
  
  df <- tibble(
    shap = shap_values[, feature_name],
    value = X_sample[[feature_name]]
  )
  
  p <- ggplot(df, aes(x = value, y = shap)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    labs(
      title = paste("SHAP Dependence:", feature_name),
      subtitle = paste("Mean |SHAP| =", formatC(mean_shap, format = "f", digits = 4)),
      x = feature_name,
      y = "SHAP Value"
    ) +
    theme_minimal()
  
  # Apply custom limits if defined
  if (!is.null(custom_limits[[feature_name]])) {
    if (!is.null(custom_limits[[feature_name]]$x)) {
      p <- p + xlim(custom_limits[[feature_name]]$x)
    }
    if (!is.null(custom_limits[[feature_name]]$y)) {
      p <- p + ylim(custom_limits[[feature_name]]$y)
    }
  }
  
  return(p)
}

shap_plots <- lapply(top_features, plot_shap_dependence, shap_values, X_sample, shap_summary,
                     custom_limits = custom_limits)

combined_plot <- wrap_plots(shap_plots, ncol = 2)

ggsave("Figures/SHAP_dependence_grid_top10.svg", combined_plot, width = 12, height = 15)



# Assume you already have:
# - shap_values: SHAP matrix (rows = patients, cols = features)
# - avg_preds: vector of predicted CATEs for each patient
# - X_sample: original feature data used in SHAP (same rows as shap_values)

all_patient_cates %>% group_by(Patient) %>% summarise(avg_pred = mean(predicted_cate))

avg_preds <- all_patient_cates$predicted_cate[1:N_explain]

# Step 1: Identify patients
cate_df <- tibble(
  row_id = 1:N_explain,
  cate = avg_preds
)

top_benefit <- cate_df %>% arrange(cate) %>% slice_head(n = 3)
top_harm    <- cate_df %>% arrange(desc(cate)) %>% slice_head(n = 2)

selected_rows <- bind_rows(top_benefit, top_harm)$row_id

# Step 2: Waterfall plot function
plot_waterfall_shap <- function(row_index, shap_matrix, feature_matrix, top_n = 15) {
  shap_row <- shap_matrix[row_index, ]
  baseline <- mean(avg_preds)  # base value (optional)
  patient_features <- feature_matrix[row_index, ]
  
  df <- tibble(
    feature = names(shap_row),
    shap_value = as.numeric(shap_row),
    value = unlist(patient_features)
  ) %>%
    arrange(desc(abs(shap_value))) %>%
    slice_head(n = top_n) %>%
    mutate(
      direction = ifelse(shap_value > 0, "Increase", "Decrease"),
      feature_label = paste0(feature, " = ", round(value, 2))
    )
  
  ggplot(df, aes(x = reorder(feature_label, shap_value), y = shap_value, fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Increase" = "red", "Decrease" = "blue")) +
    labs(
      title = paste("SHAP Waterfall Plot for Patient", row_index),
      subtitle = paste("Predicted CATE:", round(avg_preds[row_index], 3)),
      x = "Feature (with value)",
      y = "SHAP Value (Impact on Prediction)"
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top")
}

# Step 3: Plot for selected patients
plots <- lapply(selected_rows, function(i) {
  plot_waterfall_shap(i, shap_values, X_sample)
})

# Optional: Print plots
for (p in plots) {
  print(p)
}

############## Post Hoc Investigations ######################

# distribution of creatinines

filter(df_baseline2, Characteristic == "Creatinine") %>%
  ggplot(aes(y = as.numeric(Value),
             fill = Exposure,
             x = Exposure)) +
  geom_violin() +
  labs(y = "Creatinine",
       x = "") +
  theme_minimal()

# mean creatinines after excluding all those with Cr > 400
filter(df_baseline2, Characteristic == "Creatinine") %>%
  mutate(Value = as.numeric(Value)) %>%
  filter(Value <= 400) %>%
  group_by(Exposure) %>%
  summarise(mean = mean(Value, na.rm = T),
            sd = sd(Value, na.rm = T))


# Oncology diagnosis subgroup

df_oncodx <- 
  df_baseline2 %>%
    mutate(OncologicDx = ifelse(
      Characteristic == "Principal Problem" &
      Value == "Cancer", 1, 0
    )) %>%
    group_by(Patient, Exposure) %>%
    summarise(OncologicDx = max(OncologicDx)) %>%
    left_join(df_brms_primary,
              by = c("Patient", "Exposure")) 

fit_onco <- brm(
  Outcome ~ Exposure*OncologicDx + (1 | Hospital),
  data = df_oncodx,
  prior =  c(set_prior("normal(0, 1)", class = "b"),
             set_prior("cauchy(0, 25)", class = "sd")),
  family = "bernoulli",
  cores = 4,
  control = list(adapt_delta = 0.9)
)

saveRDS(fit_onco, paste0("Fits/fit_onco.RDS"))

post_onco <- as_draws_df(fit_onco)

# subgroup cov = 0 control rate
control_rate_OncoDx0 <- df_oncodx %>%
  ungroup() %>%
  filter(Exposure == "Noreplacement" &
           OncologicDx == 0) %>% 
  summarise(control = mean(Outcome)) %>% 
  unlist()

# subgroup cov = 1 control rate
control_rate_OncoDx1 <- df_oncodx %>%
  ungroup() %>%
  filter(Exposure == "Noreplacement" &
           OncologicDx == 1) %>% 
  summarise(control = mean(Outcome)) %>% 
  unlist()

arr_OncoDx0 <- meanARR_95CrI(post_onco$b_ExposureReplacement, control_rate_OncoDx0)
arr_OncoDx1 <- meanARR_95CrI(post_onco$b_ExposureReplacement+post_onco$`b_ExposureReplacement:OncologicDx`,
                               control_rate_OncoDx1)

############# Post Hoc Adjusted Analysis ##################

# Factors associated with afib and magnesium level

# Strong confounders (strong association with magnesium level and AFib)
  # Potassium
  # Creatinine
  # Urine output
  # septic shock
# Weak confounders (weak association with one of magnesium level or AFib, strong association with the other)
  # Baseline AFib (ie AFib at time 0)
  # History of AFib
  # History of CHF or ACS
  # Age
  # lactate
  # Diagnosis of acute coronary syndrome or CHF
  # vasopressor use

# using df_risk from HTE analysis

# center and scale continuous variables
mean_sd <- 
  df_risk %>%
  mutate(Creatinine = log(Creatinine)) %>%
  select(Potassium,
         Creatinine,
         UrineOutput_24h,
         Age,
         Lactate,
         Norepi_Equivalents_Per_Hour) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(
    mean = mean(value, na.rm = T),
    sd = sd(value, na.rm = T)
  ) %>%
  ungroup() %>%
  column_to_rownames(var = "name")


df_risk_for_mice3 <- 
  df_risk %>%
  mutate(Creatinine = log(Creatinine)) %>%
  mutate(Potassium = (Potassium - mean_sd["Potassium", "mean"])/(mean_sd["Potassium","sd"]),
         Creatinine = (Creatinine - mean_sd["Creatinine", "mean"])/(mean_sd["Creatinine","sd"]),
         UrineOutput_24h = (UrineOutput_24h - mean_sd["UrineOutput_24h", "mean"])/(mean_sd["UrineOutput_24h","sd"]),
         Age = (Age - mean_sd["Age", "mean"])/(mean_sd["Age","sd"]),
         Lactate = (Lactate - mean_sd["Lactate", "mean"])/(mean_sd["Lactate","sd"]),
         # don't center norepi_equivalents at the mean
         Norepi_Equivalents_Per_Hour = (Norepi_Equivalents_Per_Hour)/(mean_sd["Norepi_Equivalents_Per_Hour","sd"])) %>%
  left_join(
    distinct(df_baseline2,
             Patient, Exposure),
    by = "Patient"
  ) %>%
  select(-HEIGHT, -`WEIGHT/SCALE`,
         -CorrectedCalcium)

imp3 <- mice(df_risk_for_mice3 %>% select(-Patient), 
             m = 10, method = "pmm", maxit = 5, seed = 123)

# isolate relevant confounders

library(Hmisc)

# Helper to make ordered Q1–Q5 factors safely
to_quintile <- function(x) {
  factor(
    as.integer(cut2(x, g = 5)),
    levels = 1:5,
    labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
    ordered = TRUE
  )
}

data_list <- list()

for (i in 1:imp3$m) {
  data_list[[i]] <-
    complete(imp3, i) |>
    select(
      Exposure, Outcome,
      Sex_Male,
      Potassium,
      Creatinine,
      UrineOutput_24h,
      Principal_Problem_Sepsis_and_septic_shock,
      `Comorbidity_AFib/flutter`,
      `BaselineAfib/flutter`,
      Principal_Problem_Congestive_heart_failure_or_other_cardiac_problem_not_due_to_acute_ischemia,
      Principal_Problem_Acute_coronary_syndromes,
      Age,
      Lactate,
      Norepi_Equivalents_Per_Hour
    ) |>
    mutate(
      CHF_or_ACS = Principal_Problem_Congestive_heart_failure_or_other_cardiac_problem_not_due_to_acute_ischemia > 0 |
        Principal_Problem_Acute_coronary_syndromes > 0
    ) |>
    select(
      -Principal_Problem_Congestive_heart_failure_or_other_cardiac_problem_not_due_to_acute_ischemia,
      -Principal_Problem_Acute_coronary_syndromes
    ) |>
    rename(
      Afib_Hx = `Comorbidity_AFib/flutter`,
      Afib_Baseline = `BaselineAfib/flutter`
    ) |>
    mutate(
      Age = to_quintile(Age),
      Creatinine = to_quintile(Creatinine),
      UrineOutput_24h = to_quintile(UrineOutput_24h),
      Norepi_Equivalents_Per_Hour = to_quintile(Norepi_Equivalents_Per_Hour)
    ) |>
    bind_cols(
      df_baseline |>
        filter(Characteristic == "Hospital") |>
        select(Value) |>
        rename(Hospital = Value)
    )
}



# adjusted_fit <- 
#   brm_multiple(data = data_list,
#       formula = Outcome ~ Exposure + (1|Hospital) +
#         Age + Sex_Male +
#         Potassium + Creatinine + UrineOutput_24h +
#         Principal_Problem_Sepsis_and_septic_shock +
#         Afib_Hx + Afib_Baseline +
#         CHF_or_ACS + Lactate + Norepi_Equivalents_Per_Hour,
#       cores = 12,
#       chains = 4,
#       prior =  c(set_prior("normal(0, 1)", class = "b"),
#                  set_prior("cauchy(0, 25)", class = "sd")),
#       family = "bernoulli",
#       control = list(adapt_delta = 0.9)
#     )

library(brms)

# Compile model only once on the first dataset
base_fit <- brm(
  formula = Outcome ~ Exposure + (1 | Hospital) +
    Age + Sex_Male +
    Potassium + Creatinine + UrineOutput_24h +
    Principal_Problem_Sepsis_and_septic_shock +
    Afib_Hx + Afib_Baseline +
    CHF_or_ACS + Lactate + Norepi_Equivalents_Per_Hour,
  data = data_list[[1]],
  family = bernoulli(),
  prior = c(
    set_prior("normal(0, 1)", class = "b"),
    set_prior("cauchy(0, 25)", class = "sd")
  ),
  cores = 12,
  chains = 4,
  control = list(adapt_delta = 0.9),
  save_model = TRUE  # speeds up reuse
)

# Fit to remaining datasets without recompiling
fit_list <- vector("list", length(data_list))
fit_list[[1]] <- base_fit

for (i in 2:length(data_list)) {
  fit_list[[i]] <- update(
    base_fit,
    newdata = data_list[[i]],
    recompile = FALSE  # <- this is key
  )
}


saveRDS(fit_list, "Fits/adjusted_fit_list.RDS")
fit_list <- readRDS("Fits/adjusted_fit_list.RDS")

exposure_draws <- map_dfr(
  fit_list,
  ~ as_draws_df(.x, variable = "b_ExposureReplacement") |>
    select(b_ExposureReplacement),
  .id = "imputation"
)

exp(mean(exposure_draws$b_ExposureReplacement))
exp(quantile(exposure_draws$b_ExposureReplacement, c(0.025, 0.975)))


post_adjusted <- as_draws_df(adjusted_fit)

adjusted_OR <- mean_95CrI_OR(exposure_draws$b_ExposureReplacement)

# Control group event rate
control_rate_adjusted <- df_brms_primary %>% 
  filter(Exposure == "Noreplacement") %>% 
  summarise(control = mean(Outcome)) %>%
  unlist()

# ARR using control group rate
adjusted_ARR <- 
  meanARR_95CrI(x = post_adjusted$b_ExposureReplacement,
                rate = control_rate_adjusted)


# exclude those with baseline Afib

data_list <- list()


for (i in 1:imp3$m) {
  data_list[[i]] <- 
    complete(imp3, i) %>%
    select(Exposure, Outcome
           , Sex_Male
           , Potassium
           , Creatinine
           , UrineOutput_24h
           , Principal_Problem_Sepsis_and_septic_shock
           , 'Comorbidity_AFib/flutter'
           , 'BaselineAfib/flutter'
           , Principal_Problem_Congestive_heart_failure_or_other_cardiac_problem_not_due_to_acute_ischemia
           , Principal_Problem_Acute_coronary_syndromes
           , Age
           , Lactate
           , Norepi_Equivalents_Per_Hour) %>%
    mutate(CHF_or_ACS = Principal_Problem_Congestive_heart_failure_or_other_cardiac_problem_not_due_to_acute_ischemia > 0 |
             Principal_Problem_Acute_coronary_syndromes > 0) %>%
    select(-Principal_Problem_Congestive_heart_failure_or_other_cardiac_problem_not_due_to_acute_ischemia,
           -Principal_Problem_Acute_coronary_syndromes) %>%
    rename(Afib_Hx = `Comorbidity_AFib/flutter`,
           Afib_Baseline = `BaselineAfib/flutter`) %>%
    bind_cols(df_baseline %>% 
                filter(Characteristic == "Hospital") %>% 
                select(Value) %>% 
                rename(Hospital = Value)) %>%
    filter(Afib_Baseline == 0)
}

adjusted_fit_noAfib <- 
  brm_multiple(data = data_list,
               formula = Outcome ~ Exposure + (1|Hospital) +
                 Age + Sex_Male +
                 Potassium + Creatinine + UrineOutput_24h +
                 Principal_Problem_Sepsis_and_septic_shock +
                 Afib_Hx + 
                 CHF_or_ACS + Lactate + Norepi_Equivalents_Per_Hour,
               cores = 12,
               chains = 4,
               prior =  c(set_prior("normal(0, 1)", class = "b"),
                          set_prior("cauchy(0, 25)", class = "sd")),
               family = "bernoulli",
               control = list(adapt_delta = 0.9)
  )

saveRDS(adjusted_fit_noAfib, "Fits/adjusted_fit_noAfib.RDS")

########### Baseline characteristics by magnesium level ###############

df_baseline3 <- 
  df_baseline2 %>%
  left_join(
    df_baseline %>%
    filter(Characteristic == "Magnesium") %>%
    mutate(
      Magnesium = case_when(
        Value >= 0.92 & Value <= 0.93 ~ "0.92–0.93",
        Value >= 0.94 & Value <= 0.95 ~ "0.94–0.95",
        Value >= 0.96 & Value <= 0.97 ~ "0.96–0.97",
        Value >= 0.98 & Value <= 0.99 ~ "0.98–0.99",
        TRUE ~ NA_character_
      )
    ) %>%
    select(-Characteristic, -Value),
    by = "Patient") %>%
  select(-Exposure) %>%
  rename(Exposure = Magnesium)

TableA <-
  df_baseline3 %>%
  filter(!(Characteristic %in% characteristics_category
           | Characteristic %in% characteristics_ignore)) %>%
  mutate(Value = as.numeric(Value)) %>%
  group_by(Characteristic, Exposure) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Label = paste0(round(Mean, 2), " (", round(SD, 2), ")")) %>%
  select(Characteristic, Exposure, Label) %>%
  pivot_wider(names_from = Exposure, values_from = Label) %>%
  arrange(Characteristic)

TableB <- 
  df_baseline3 %>%
  mutate(Value = ifelse(
    Characteristic == "R CARDIAC RHYTHM",
    ifelse(Value %in% c("A-fib", "A-flutter"),
           "Afib/flutter", "Not AFib/flutter"),
    Value
  )) %>%
  filter(Characteristic %in% characteristics_category) %>%
  group_by(Characteristic, Value, Exposure) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Characteristic, Exposure) %>%
  mutate(N = sum(count)) %>%
  ungroup() %>%
  mutate(perc = round(100 * count / N),
         label = paste0(count, " (", perc, "%)")) %>%
  select(Characteristic, Value, Exposure, label) %>%
  pivot_wider(names_from = Exposure, values_from = label) %>%
  arrange(Characteristic, Value)

TableA <- TableA %>%
  mutate(Section = "Continuous Variables", Value = "") %>%
  relocate(Section, Characteristic, Value)

TableB <- TableB %>%
  mutate(Section = "Categorical Variables") %>%
  relocate(Section)

CombinedTable <- bind_rows(TableA, TableB) %>%
  arrange(Section, Characteristic)

write.csv(CombinedTable, "Tables/Baseline_by_Mg.csv", row.names = FALSE)

############## Per protocol analysis ##################

df_per_prot <- 
  df_brms_primary %>%
  left_join(
    df_outcome %>%
      mutate(MagReplaced = ifelse(Characteristic == "NamedMedication", 1,0)) %>%
      group_by(Patient) %>%
      summarise(MagReplaced = max(MagReplaced)),
    by = "Patient"
  ) %>%
  filter(Exposure == "Noreplacement" & MagReplaced == 0 | 
           Exposure == "Replacement" & MagReplaced == 1)

df_per_prot %>%
  ungroup() %>%
  group_by(Exposure) %>%
  summarise(n(),
            sum(Outcome))

fit_per_prot <- brm(
  Outcome ~ Exposure + (1 | Hospital),
  data = df_per_prot,
  prior =  c(set_prior("normal(0, 1)", class = "b"),
             set_prior("cauchy(0, 25)", class = "sd")),
  family = "bernoulli",
  cores = 4,
  control = list(adapt_delta = 0.9)
)

saveRDS(fit_per_prot, file = "Fits/fit_per_prot.RDS")

fit_per_prot <- readRDS("Fits/fit_per_prot.RDS")

post_per_prot <- as_draws_df(fit_per_prot)

per_prot_OR <- mean_95CrI_OR(post_per_prot$b_ExposureReplacement)

# Control group event rate
control_rate_per_prot <- df_per_prot %>% 
  filter(Exposure == "Noreplacement") %>% 
  summarise(control = mean(Outcome)) %>%
  unlist()

# ARR using control group rate
per_prot_ARR <- 
  meanARR_95CrI(x = post_per_prot$b_ExposureReplacement,
                rate = control_rate_per_prot)

###################### e-value ############

# use 1/RR + sqrt(1/RR*(1/RR-1)) from Chung et al

# 1/RR = 20.1/(20.1-3.1) = 1.18

eval <- function(x){x + sqrt(x*(x-1))}

eval(1.18)

eval(20.1/(20.1-1.7))


##################### SMDs of the log ####################

Cr <- 
  df_baseline2 %>%
  filter(Characteristic == "Creatinine") %>%
  rename(Creatinine = Value) %>%
  select(-Characteristic, -Patient) %>%
  mutate(Creatinine = as.numeric(Creatinine))

Cr %>%
  group_by(Exposure) %>%
  summarise(quantile(Creatinine, 0.25, na.rm = T),
            median(Creatinine, na.rm = T),
            quantile(Creatinine, 0.75, na.rm = T))


# Example calculation of SMD
calculate_smd <- function(df, value_col = "value", group_col = "group", group1 = "R", group2 = "N") {
  x1 <- df[[value_col]][df[[group_col]] == group1]
  x2 <- df[[value_col]][df[[group_col]] == group2]
  
  m1 <- mean(x1, na.rm = TRUE)
  m2 <- mean(x2, na.rm = TRUE)
  sd1 <- sd(x1, na.rm = TRUE)
  sd2 <- sd(x2, na.rm = TRUE)
  n1 <- sum(!is.na(x1))
  n2 <- sum(!is.na(x2))
  
  # Pooled standard deviation
  s_pooled <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  
  smd <- (m1 - m2) / s_pooled
  return(smd)
}

smd_result <- calculate_smd(Cr %>% mutate(Creatinine = log(Creatinine)), value_col = "Creatinine", group_col = "Exposure",
                            group1 = "Replacement", group2 = "Noreplacement")
print(smd_result)


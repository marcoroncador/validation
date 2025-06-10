setwd("~/Desktop/PROJECTS/MayoValidation/")

# Load required libraries
library(survival)
library(survminer)
library(pROC)
library(AICcmodavg)

# Load the dataset
patients <- read.csv("training_6k_postreview.csv", header = TRUE)
external <- read.table("NCRI_paper_aml_prognosis_updated.tsv", header = T)

# Select patients with ID starting with "PD" and Treatment == 1
patients <- patients[grep("^PD", patients$ID), ]
patients <- subset(patients, Treatment == 1)

external <- external[external$ID %in% patients$ID,]

# Define your columns
num_cols <- c("bm_blasts", "hb", "plt", "wbc", "age")
fac_cols <- c("Treatment", "SEX",
              "WHO_2016", "WHO_2022", "ICC",
              "IPSSR_ELN", "ELN2022_IPSSM")
library(dplyr)

external <- external %>%
  mutate(hb = case_when(
    hb > 20  ~ hb / 10,   # anything over 20 → assume it was 10× too big
    hb < 2   ~ hb * 10,   # anything under 2  → assume it was 10× too small
    TRUE     ~ hb         # otherwise leave it alone
  ))


# 1a) Quick overall summary (numeric → min/1Q/median/mean/3Q/max; factor → counts)
summary(external[, c(num_cols)])

# 1b) If you want detail on each factor separately:
lapply(patients[fac_cols], table, useNA = "ifany")


# Load additional data
NCRI <- read.table("NCRI_paper_aml_prognosis_updated.tsv", header = TRUE)
AMLSG <- read.csv("full_data_validation.csv", header = TRUE)

# Summarize NRAS in NCRI dataset
NCRI$NRAS <- ifelse(NCRI$NRAS_other == 1 | NCRI$NRAS_p.G12_13 == 1 | NCRI$NRAS_p.Q61_62 == 1, 1, 0)



# Restore original cytogenetic information 
cytofix <- AMLSG[,c("ID","complex","del_5", "del_7","del_17")]
colnames(cytofix) <- c("ID","complex","X.5","X.7", "X.17")
cytofix <- rbind(cytofix,NCRI[,c("ID","complex","X.5","X.7", "X.17")])

patients[,c("complex","X.5","X.7", "X.17")] <- NULL
patients <- merge(patients, cytofix, by = "ID", all.x = TRUE)

ITDvec <- rbind(NCRI[, c("ID", "ITD")], AMLSG[, c("ID", "ITD")])
patients <- merge(patients, ITDvec, by = "ID", all.x = TRUE)

# Calculate Mayo Risk Score
patients$Mayo <- (patients$complex == 1 | patients$X.5 ==1 | patients$X.7 ==1 | patients$X.17 ==1)*1 +
  (patients$IDH2 == 0)*1 +  
  (patients$TP53 == 1)*1 +  
  (patients$KRAS == 1)*1 +  
  ((patients$t.9.11. == 1 | patients$t.v.11. == 1) * 2)

patients$Mayo_risk <- cut(patients$Mayo, 
                          breaks = c(-Inf, 0, 1, Inf),
                          labels = c("favourable", "intermediate", "adverse"))

patients$OS <- patients$OS*12

# Calculate ELN2024 Risk Score
patients$ELN2024 <- with(patients, ifelse(TP53 == 1, "adverse", 
                                          ifelse(ITD == 1 | KRAS == 1 | NRAS == 1, "intermediate", "favourable")))


patients$ELN2024 <- as.factor(patients$ELN2024)
patients$ELN2024 <- factor(
  patients$ELN2024,
  levels = c("favourable", "intermediate", "adverse")
)
# Calculate VEN-PRS Risk Score (OS and PFS)
patients$VEN_PRS_OS <- with(patients, 
                            PTPN11 * 0.8 + FLT3 * 0.49 + TP53 * 0.39 + NF1 * 0.94 + SF3B1 * 0.98)
patients$VEN_PRS_PFS <- with(patients, 
                             PTPN11 * 0.73 + FLT3 * 0.6 + TP53 * 0 + NF1 * 0.81 + SF3B1 * 0.87)

patients$VEN_PRS_OS_risk <- cut(patients$VEN_PRS_OS, 
                                breaks = c(-Inf, 0.25, 0.75, Inf),
                                labels = c("favourable", "intermediate", "adverse"))

patients$VEN_PRS_PFS_risk <- cut(patients$VEN_PRS_PFS, 
                                 breaks = c(-Inf, 0.25, 0.75, Inf),
                                 labels = c("favourable", "intermediate", "adverse"))

patients$ELN2022_IPSSM[patients$ELN2022_IPSSM =="ELN2022_favorable"] <- "favourable"
patients$ELN2022_IPSSM[patients$ELN2022_IPSSM =="ELN2022_adverse"] <- "adverse"
patients$ELN2022_IPSSM[patients$ELN2022_IPSSM =="ELN2022_intermediate"] <- "intermediate"
patients$ELN2022_IPSSM <- as.factor(patients$ELN2022_IPSSM)

patients$ELN2022_IPSSM <- factor(
  patients$ELN2022_IPSSM,
  levels = c("favourable", "intermediate", "adverse")
)

patients$Mayo_risk <- factor(
  patients$Mayo_risk,
  levels = c("favourable", "intermediate", "adverse")
)

patients$ELN2024 <- factor(
  patients$ELN2024,
  levels = c("favourable", "intermediate", "adverse")
)





### COUNT MUTATIONS

# pick out all mutation columns from ASXL1 through FLT3
mut_cols <- names(patients)[2:which(names(patients) == "FLT3")]

# count how many 1’s (mutations) each gene has
mutation_counts <- colSums(patients[ , mut_cols])

# turn into a nice two‐column table
counts_table <- data.frame(
  Mutation = names(mutation_counts),
  Patients = as.integer(mutation_counts),
  row.names = NULL
)

print(counts_table)



# your raw lines
lines <- c(
  "TP53 (394)\t101 (26)",
  "RUNX1 (392)\t76 (19)",
  "TET2 (387)\t74 (19)",
  "SRSF2 (387)\t70 (18)",
  "ASXL1 (387)\t70 (18)",
  "DNMT3A (392)\t57 (15)",
  "NPM1 (394)\t49 (12)",
  "IDH2 (393)\t48 (12)",
  "IDH1 (393)\t26 (7)",
  "FLT3-ITD (397)\t39 (10)",
  "NRAS (392)\t34 (9)",
  "BCOR (387)\t32 (8)",
  "U2AF1 (372)\t22 (6)",
  "CEBPA (394)\t22 (6)",
  "CEBPA bZIP (393)\t12 (3)",
  "STAG2 (376)\t22 (8)",
  "SF3B1 (376)\t19 (5)",
  "KRAS (392)\t15 (4)",
  "DDX41 (376)\t14 (4)",
  "PTPN11 (376)\t12 (3)",
  "EZH2 (387)\t11 (3)",
  "WT1 (387)\t11 (3)",
  "SETBP1 (387)\t11 (3)",
  "PHF6 (376)\t10 (3)",
  "JAK2 (376)\t9 (2)",
  "CBL (375)\t7 (2)"
)

# remove any space+parenthesized-number
clean <- gsub("\\s*\\(\\d+\\)", "", lines)

# now split on the tab (or whitespace) into two columns
mayo_mut <- read.table(text = clean,
                 sep   = "\t",
                 header= FALSE,
                 stringsAsFactors = FALSE,
                 col.names = c("Mutation", "Patients"))


# assume mayo_mut and counts_table have Mutation & Patients columns
N1 <- 400
N2 <- 352

# join the two tables
df <- merge(mayo_mut, counts_table, by="Mutation", suffixes=c("_mayo","_ct"))

# for each gene, build 2×2 and choose test
results <- lapply(seq_len(nrow(df)), function(i) {
  x1 <- df$Patients_mayo[i]
  x2 <- df$Patients_ct[i]
  tbl <- matrix(c(x1, N1 - x1, x2, N2 - x2),
                nrow=2,
                dimnames=list(Cohort=c("Mayo","Counts"),
                              Status=c("Mut","WT")))
  # decide test
  if(any(chisq.test(tbl, simulate.p.value=FALSE)$expected < 5)) {
    res <- fisher.test(tbl)
    test <- "Fisher"
  } else {
    res <- chisq.test(tbl, correct=FALSE)
    test <- "Chi‐square"
  }
  c(Test      = test,
    p.value   = res$p.value,
    OddsRatio = if(!is.null(res$estimate)) unname(res$estimate) else NA)
})

res_df <- do.call(rbind, results)
res_df <- data.frame(Mutation=df$Mutation, res_df,
                     p.adj = p.adjust(res_df[,"p.value"], method="BH"),
                     row.names=NULL)
print(res_df)




#### ist ELN assigning more patients to high risk?

# Ensure necessary package is installed
# Set levels explicitly to ensure order matches for Bowker's test
patients$ELN2024 <- factor(patients$ELN2024, levels = c("favourable", "intermediate", "adverse"))
patients$Mayo_risk <- factor(patients$Mayo_risk, levels = c("favourable", "intermediate", "adverse"))


# Create 3x3 contingency table from your dataframe
# Replace df with your actual dataframe name
table3x3 <- table(patients$ELN2024, patients$Mayo_risk)

# View the table to confirm structure
print(table3x3)

# Perform symmetry test (Bowker-type) via coin
mcnemar.test(table3x3)




# Kaplan-Meier analysis
fit_mayo <- survfit(Surv(OS, OS_STATUS) ~ Mayo_risk, data = patients)
fit_eln <- survfit(Surv(OS, OS_STATUS) ~ ELN2024, data = patients)
fit_ven_prs <- survfit(Surv(OS, OS_STATUS) ~ VEN_PRS_OS_risk, data = patients)
fit_eln2022 <- survfit(Surv(OS, OS_STATUS) ~ ELN2022_IPSSM, data = patients)


# Plot Kaplan-Meier curves
png(filename = "validation_mayo.png", width = 2400, height = 2000, res = 300)
ggsurvplot(fit_mayo, data = patients, pval = TRUE, conf.int = TRUE, 
           title = paste0("Kaplan-Meier Curve for Mayo Risk Score - n=", dim(patients)[1]),
           risk.table = TRUE, ggtheme = theme_minimal(), palette =c( "#9ACD32",  "#FFA500","#8B0000" ) )
dev.off()

png(filename = "validation_ELN.png", width = 2400, height = 2000, res = 300)
ggsurvplot(fit_eln, data = patients, pval = TRUE, conf.int = TRUE, 
           title = paste0("Kaplan-Meier Curve for ELN2024 Risk Score - n=", dim(patients)[1]),
           risk.table = TRUE, ggtheme = theme_minimal(), palette =c( "#9ACD32",  "#FFA500","#8B0000" ) )
dev.off()


ggsurvplot(fit_eln2022, data = patients, pval = TRUE, conf.int = TRUE, xlim =c(0,36),
           title = paste0("Kaplan-Meier Curve for ELN2022 Risk Score - n=", dim(patients)[1]),
           risk.table = TRUE, ggtheme = theme_minimal(), palette =c( "#e2B48C", "#7B0000", "black") )
dev.off()


plot_mayo <- ggsurvplot(fit_mayo, data = patients, pval = TRUE, conf.int = TRUE, 
                        title = paste0("Mayo Risk Score - n=", dim(patients)[1]), 
                        xlab = "Time (months)", risk.table.fontsize = 4, 
                        risk.table = TRUE, ggtheme = theme_minimal(), break.time.by = 12, pval.coord = c(30, 0.5), xlim =c(0,36),
                        palette = c("#9ACD32", "#FFA500", "#8B0000"), legend.labs  = c("favourable","intermediate","adverse") )

plot_eln <- ggsurvplot(fit_eln, data = patients, pval = TRUE, conf.int = TRUE,  xlab = "Time (months)", risk.table.fontsize = 4, 
                       title = paste0("ELN2024 Risk Score - n=", dim(patients)[1]),
                       risk.table = TRUE, ggtheme = theme_minimal(), break.time.by = 12, pval.coord = c(30, 0.5), xlim =c(0,36),
                       palette = c("#9ACD32", "#FFA500", "#8B0000"), legend.labs  = c("favourable","intermediate","adverse")  )

plot_eln2022 <- ggsurvplot(fit_eln2022, data = patients, pval = TRUE, conf.int = TRUE,  xlab = "Time (months)", risk.table.fontsize = 4, 
           title = paste0("ELN2022 Risk Score - n=", dim(patients)[1]),
           risk.table = TRUE, ggtheme = theme_minimal(), break.time.by = 12, pval.coord = c(30, 0.5), xlim =c(0,36),
           palette = c("#e2B48C","#7B0000", "black" ) )

plot_eln2022$plot <- plot_eln2022$plot +
  theme(
    plot.title   = element_text(face = "bold", size = 16),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12, color = "black"),
    legend.text  = element_text(size = 10),    # <<<<< Increase legend font
    legend.title = element_text(size = 10)     # (optional: also legend title bigger)
  )

# Now fix the p-value (assuming it's a "subtitle" or "caption" or annotation)
plot_eln2022$plot <- plot_eln2022$plot +
  theme(
    plot.subtitle = element_text(size = 10)    # <<<<< Smaller p-value if in subtitle
  )

png(filename = "validation_ELN2022.png", width = 2400, height = 2000, res = 300)
print(plot_eln2022)
dev.off()


 plot_mayo$plot <- plot_mayo$plot +
  theme(
    plot.title   = element_text(face = "bold", size = 16),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12, color = "black")
  )

plot_eln$plot <- plot_eln$plot +
  theme(
    plot.title   = element_text(face = "bold", size = 16),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12, color = "black")
  )


plot_eln$plot <- plot_eln$plot +
  theme(
    plot.title   = element_text(face = "bold", size = 16),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12, color = "black")
  )


big_legend <- theme(
  legend.title      = element_text(size = 14),
  legend.text       = element_text(size = 14),
  legend.key.size   = unit(0.6, "cm")      # make the little boxes bigger
)

# apply it to both plots
plot_mayo$plot <- plot_mayo$plot + big_legend
plot_eln$plot  <- plot_eln$plot  + big_legend

# Arrange the plots using ggarrange
combined_plot <- ggarrange(plot_mayo$plot, plot_eln$plot, 
                           plot_mayo$table, plot_eln$table,
                           ncol = 2, nrow = 2, heights = c(3,1,3,1),common.legend = T,          # pull the legends into one
                           legend        = "bottom",
                           labels = c("", "", "", "")); combined_plot

# Save the arranged plot
ggsave("validation_combined.png", combined_plot, width = 16, height = 8, dpi = 450)
ggsave("Fig_2_validation_combined.png", combined_plot, width = 12, height = 7, dpi = 300)



png(filename = "validation_VEN_PRS.png", width = 2400, height = 2000, res = 300)
ggsurvplot(fit_ven_prs, data = patients, pval = TRUE, conf.int = TRUE, 
           title = paste0("Kaplan-Meier Curve for VEN-PRS OS - n=", dim(patients)[1]),
           risk.table = TRUE, ggtheme = theme_minimal(), palette =c( "#9ACD32",  "#FFA500","#8B0000" ) )
dev.off()

write.csv(file = "patients_3scores.csv",patients,row.names = F,quote = T)

# Compare AUC
roc_mayo <- roc(patients$OS_STATUS, as.numeric(patients$Mayo))
roc_eln <- roc(patients$OS_STATUS, as.numeric(factor(patients$ELN2024)))
roc_ven_prs <- roc(patients$OS_STATUS, as.numeric(factor(patients$VEN_PRS_OS_risk)))

auc_mayo <- auc(roc_mayo)
auc_eln <- auc(roc_eln)
auc_ven_prs <- auc(roc_ven_prs)

# Compare AICc values
fit_mayo_aic <- coxph(Surv(OS, OS_STATUS) ~ Mayo_risk, data = patients)
fit_eln_aic <- coxph(Surv(OS, OS_STATUS) ~ ELN2024, data = patients)
fit_ven_prs_aic <- coxph(Surv(OS, OS_STATUS) ~ VEN_PRS_OS_risk, data = patients)

aicc_mayo <- AICc(fit_mayo_aic)
aicc_eln <- AICc(fit_eln_aic)
aicc_ven_prs <- AICc(fit_ven_prs_aic)

# Print results
cat("AUC Mayo:", auc_mayo, "\n")
cat("AUC ELN2024:", auc_eln, "\n")
cat("AUC VEN-PRS:", auc_ven_prs, "\n")
cat("AICc Mayo:", aicc_mayo, "\n")
cat("AICc ELN2024:", aicc_eln, "\n")
cat("AICc VEN-PRS:", aicc_ven_prs, "\n")


##### MORE STATISTICAL INDICATORS:

library(survival)
library(pROC)
library(MuMIn)  # for AICc
library(dplyr)

# Fit models
fit_mayo <- coxph(Surv(OS, OS_STATUS) ~ Mayo_risk, data = patients)
fit_eln <- coxph(Surv(OS, OS_STATUS) ~ ELN2024, data = patients)
fit_ven_prs <- coxph(Surv(OS, OS_STATUS) ~ VEN_PRS_OS_risk, data = patients)

# Compute ROC and AUC
roc_mayo <- roc(patients$OS_STATUS, as.numeric(patients$Mayo))
roc_eln <- roc(patients$OS_STATUS, as.numeric(factor(patients$ELN2024)))
roc_ven_prs <- roc(patients$OS_STATUS, as.numeric(factor(patients$VEN_PRS_OS_risk)))

# Store results
results <- data.frame(
  Model = c("Mayo_risk", "ELN2024", "VEN_PRS_OS_risk"),
  AUC = c(auc(roc_mayo), auc(roc_eln), auc(roc_ven_prs)),
  AICc = c(AICc(fit_mayo), AICc(fit_eln), AICc(fit_ven_prs)),
  BIC = c(BIC(fit_mayo), BIC(fit_eln), BIC(fit_ven_prs)),
  C_index = c(summary(fit_mayo)$concordance[1],
              summary(fit_eln)$concordance[1],
              summary(fit_ven_prs)$concordance[1]),
  LR_pvalue = c(summary(fit_mayo)$logtest["pvalue"],
                summary(fit_eln)$logtest["pvalue"],
                summary(fit_ven_prs)$logtest["pvalue"]),
  Wald_pvalue = c(summary(fit_mayo)$waldtest["pvalue"],
                  summary(fit_eln)$waldtest["pvalue"],
                  summary(fit_ven_prs)$waldtest["pvalue"]),
  Score_pvalue = c(summary(fit_mayo)$sctest["pvalue"],
                   summary(fit_eln)$sctest["pvalue"],
                   summary(fit_ven_prs)$sctest["pvalue"])
)

# Round for nicer output
results <- results %>%
  mutate(across(where(is.numeric), round, digits = 6))

# Write to CSV
write.csv(results, "risk_model_comparison.csv", row.names = FALSE)

# Print to console
print(results)







### Sankey 
# install.packages("remotes")
# remotes::install_github("fbreitwieser/sankeyD3")

# install.packages("remotes")
# remotes::install_github("fbreitwieser/sankeyD3")

library(dplyr)
library(sankeyD3)

patients$ELN2022_IPSSM[patients$ELN2022_IPSSM =="ELN2022_favorable"] <- "favourable"
patients$ELN2022_IPSSM[patients$ELN2022_IPSSM =="ELN2022_adverse"] <- "adverse"
patients$ELN2022_IPSSM[patients$ELN2022_IPSSM =="ELN2022_intermediate"] <- "intermediate"
patients$ELN2022_IPSSM <- as.factor(patients$ELN2022_IPSSM)

patients$ELN2022_IPSSM <- factor(
  patients$ELN2022_IPSSM,
  levels = c("favourable", "intermediate", "adverse")
)

patients$Mayo_risk <- factor(
  patients$Mayo_risk,
  levels = c("favourable", "intermediate", "adverse")
)

patients$ELN2024 <- factor(
  patients$ELN2024,
  levels = c("favourable", "intermediate", "adverse")
)

# Step 1: Standardize text values
patients <- patients %>%
  mutate(
    Mayo_risk         = factor(Mayo_risk,         levels = c("favourable", "intermediate", "adverse")),
    ELN2022_IPSSM     = factor(ELN2022_IPSSM,      levels = c("favourable", "intermediate", "adverse")),
    ELN2024           = factor(ELN2024,           levels = c("favourable", "intermediate", "adverse"))
  )

# Step 2: Create long format edge list
flow1 <- patients %>%
  count(
    source = paste0("Mayo: ", Mayo_risk),
    target = paste0("ELN2022: ", ELN2022_IPSSM),
    name   = "value"
  )

flow2 <- patients %>%
  count(
    source = paste0("ELN2022: ", ELN2022_IPSSM),
    target = paste0("ELN2024: ", ELN2024),
    name   = "value"
  )

links <- bind_rows(flow1, flow2)

# Step 3: Create node list from unique sources and targets
nodes <- data.frame(name = unique(c(links$source, links$target)), stringsAsFactors = FALSE)

# Step 4: Map node names to indices
links$source <- match(links$source, nodes$name) - 1
links$target <- match(links$target, nodes$name) - 1

# Step 5: Assign colors per node
color_map <- c(
  "Mayo: intermediate"     = adjustcolor("darkorange",  alpha.f = 0.75),
  "Mayo: favourable"       = adjustcolor("forestgreen", alpha.f = 0.75),
  "Mayo: adverse"          = adjustcolor("firebrick3",  alpha.f = 0.75),
  
  
  "ELN2022: intermediate" = "#7B0000",   # mid brown
  "ELN2022: favourable"   = "#e2B48C",   # light brown
  "ELN2022: adverse"      = "black",
  
  "ELN2024: intermediate"   = adjustcolor("darkorange",  alpha.f = 0.75),
  "ELN2024: favourable"     = adjustcolor("forestgreen", alpha.f = 0.75),
  "ELN2024: adverse"        = adjustcolor("firebrick3",  alpha.f = 0.75)
)

node_groups <- color_map[nodes$name]


nodes <- data.frame(
  name = c(
    
    "Mayo: adverse",
    "Mayo: favourable",
    "Mayo: intermediate",
    
    
    
    "ELN2022: adverse",
    "ELN2022: favourable",
    "ELN2022: intermediate",
    
    "ELN2024: adverse",
    "ELN2024: favourable",
    "ELN2024: intermediate"
  ),
  stringsAsFactors = FALSE
)

links <- data.frame(
  source = c(0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5),
  target = c(3, 4, 5, 3, 4, 5, 3, 4, 5, 6, 7, 8, 7, 8, 7, 8),
  value  = c(95, 1, 2, 26, 4, 7, 133, 29, 55, 52, 121, 81, 28, 6, 32, 32)
)

font_family <- "Arial"
font_size   <- 19

# Step 6: Plot the Sankey
sankey1 <- sankeyNetwork(
  Links        = links,
  Nodes        = nodes,
  Source       = "source",
  Target       = "target",
  Value        = "value",
  NodeID       = "name",
  NodeGroup    = "name",  # use name to assign color
  colourScale  = htmlwidgets::JS(sprintf(
    "d3.scaleOrdinal().domain([\"%s\"]).range([\"%s\"])",
    paste(names(color_map), collapse = "\",\""),
    paste(color_map, collapse = "\",\"")
  )),
  numberFormat = ".0f",  # zero decimals
  nodeWidth    = 20,
  nodePadding  = 22,
  width        = 1300,
  height       = 800, 
  fontSize     = font_size,
  fontFamily   = font_family, iterations = 1
)

sankey1


nodes <- nodes[c(1, 3, 2, 4, 6, 5, 7, 9, 8), , drop = FALSE]


# 1. Build a map from old to new
old_to_new <- setNames(0:(nrow(nodes)-1), nodes$name)

# 2. But wait: your links$source and links$target are numbers, not names
# => you need to know which name was associated to which original index

# Here was the original node list BEFORE reordering
original_node_order <- c(
  "Mayo: adverse",
  "Mayo: favourable",
  "Mayo: intermediate",
  "ELN2022: adverse",
  "ELN2022: favourable",
  "ELN2022: intermediate",
  "ELN2024: adverse",
  "ELN2024: favourable",
  "ELN2024: intermediate"
)

# Now create a dataframe linking old indices to names
lookup_table <- data.frame(
  old_index = 0:(length(original_node_order)-1),
  name = original_node_order,
  stringsAsFactors = FALSE
)

# Now: map old index -> name -> new index
links$source <- old_to_new[lookup_table$name[links$source + 1]]
links$target <- old_to_new[lookup_table$name[links$target + 1]]

# Final correction: links$source and links$target must be numeric
links$source <- as.numeric(links$source)
links$target <- as.numeric(links$target)



# Step 6: Plot the Sankey
sankey2 <- sankeyNetwork(
  Links        = links,
  Nodes        = nodes,
  Source       = "source",
  Target       = "target",
  Value        = "value",
  NodeID       = "name",
  NodeGroup    = "name",  # use name to assign color
  colourScale  = htmlwidgets::JS(sprintf(
    "d3.scaleOrdinal().domain([\"%s\"]).range([\"%s\"])",
    paste(names(color_map), collapse = "\",\""),
    paste(color_map, collapse = "\",\"")
  )),
  numberFormat = ".0f",  # zero decimals
  nodeWidth    = 20,
  nodePadding  = 22,
  width        = 640,
  height       = 1300, 
  fontSize     = font_size,
  fontFamily   = font_family, iterations = 0, nodeAlign = "center"
)

sankey2






















# COMBINE THE PLOTS

library(webshot2)
library(magick)
library(cowplot)
library(ggplot2)
library(grid)

# 1. Save your sankey widget to HTML
htmlwidgets::saveWidget(
  sankey1,
  file = "sankey_temp.html",
  selfcontained = TRUE, 
  background = "transparent" 
)

htmlwidgets::saveWidget(
  sankey2,
  file = "sankey_temp2.html",
  selfcontained = TRUE, 
  background = "transparent" 
)

# 2. Render it to PNG with webshot2
webshot("sankey_temp.html", file = "sankey_temp.png",
        vwidth = 1330, vheight = 800, zoom = 4, cliprect = "viewport",   delay  = 0.5) # only capture the plot area)

sankey_img <- image_read("sankey_temp.png") %>%
  image_trim()  # auto-crops any solid border

webshot("sankey_temp2.html", file = "sankey_temp2.png",
        vwidth = 640, vheight = 1330, zoom = 4, cliprect = "viewport",   delay  = 0.5) # only capture the plot area)

sankey_img2 <- image_read("sankey_temp2.png") %>%
  image_trim()  # auto-crops any solid border

# convert to grob
sankey_grob <- rasterGrob(sankey_img, interpolate = TRUE)
sankey_grob2 <- rasterGrob(sankey_img2, interpolate = TRUE)

# 5. Combine ggarrange output with sankey image

final_plot <- plot_grid(
  combined_plot,
  ggplot() +
    annotation_custom(sankey_grob) +
    theme_void() +
    theme(
      plot.margin      = margin(0, 0, 0, 0),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background  = element_rect(fill = "transparent", color = NA)
    ),
  ncol             = 1,
  rel_heights      = c(1, 1.1),
  labels           = c("A", "B"),
  label_size       = 16,
  label_fontfamily = "Arial",
  label_y          = c(1, 1.025)  
); final_plot


plot_mayo$plot <- plot_mayo$plot +
  theme(
    plot.title   = element_text(face = "bold", size = 18),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 14, color = "black")
  )

plot_eln$plot <- plot_eln$plot +
  theme(
    plot.title   = element_text(face = "bold", size = 18),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 14, color = "black")
  )


plot_mayo$table <- plot_mayo$table +
  theme(text = element_text(size = 16))

plot_eln$table <- plot_eln$table +
  theme(text = element_text(size = 16))


big_legend <- theme(
  legend.title      = element_text(size = 16),
  legend.text       = element_text(size = 16),
  legend.key.size   = unit(0.6, "cm")      # make the little boxes bigger
)

# apply it to both plots
plot_mayo$plot <- plot_mayo$plot + big_legend
plot_eln$plot  <- plot_eln$plot  + big_legend


combined_plot2 <- ggarrange(plot_mayo$plot, 
                           plot_mayo$table, 
                           plot_eln$plot, 
                           plot_eln$table,
                           ncol = 1, nrow = 4, heights = c(3,1,3,1),common.legend = T,          # pull the legends into one
                           legend        = "bottom",
                           labels = c("", "", "", "")); combined_plot2


final_plot2 <- plot_grid(
  ggplot() +
    annotation_custom(sankey_grob2) +
    theme_void() +
    theme(
      plot.margin      = margin(0, 0, 0, 0),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background  = element_rect(fill = "transparent", color = NA)
    ),
  combined_plot2,
  ncol             = 2,
  rel_widths =    c(1, 0.8),
  labels           = c("A", "B "),
  label_size       = 16,
  label_fontfamily = "Arial",
  label_x          = c(0, -0.02)
)




# 6. Save the final image
ggsave(
  filename = "full_combined_plot.png",
  plot     = final_plot,
  width    = 12,
  height   = 14,
  dpi      = 300
)

ggsave(
  filename = "full_combined_plot2.png",
  plot     = final_plot2,
  width    = 16,
  height   = 14,
  dpi      = 300
)


library(dplyr)
# For ELN2022_IPSSM
median_os_eln2022 <- patients %>%
  group_by(ELN2022_IPSSM) %>%
  summarise(
    median_OS = median(OS, na.rm = TRUE),
    min_OS = min(OS, na.rm = TRUE),
    max_OS = max(OS, na.rm = TRUE)
  ) %>%
  arrange(ELN2022_IPSSM)

# For ELN2024
median_os_eln2024 <- patients %>%
  group_by(ELN2024) %>%
  summarise(
    median_OS = median(OS, na.rm = TRUE),
    min_OS = min(OS, na.rm = TRUE),
    max_OS = max(OS, na.rm = TRUE)
  ) %>%
  arrange(ELN2024)

# For Mayo_risk
median_os_mayo <- patients %>%
  group_by(Mayo_risk) %>%
  summarise(
    median_OS = median(OS, na.rm = TRUE),
    min_OS = min(OS, na.rm = TRUE),
    max_OS = max(OS, na.rm = TRUE)
  ) %>%
  arrange(Mayo_risk)

# View results
median_os_eln2022
median_os_eln2024
median_os_mayo




table(patients[,c("Mayo_risk","ELN2024")])
 
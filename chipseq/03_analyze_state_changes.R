## Analyze fuso-induced state changes.
## Aug 28, 2020
## sbrown

setwd("/projects/sbrown_prj/190505_Fuso-Inv")


if(!"alluvial" %in% installed.packages()){
  install.packages("alluvial")
}
library(alluvial)


all_states <- c("1_TssA", "2_TssFlnkU", "3_TssFlnkD", "4_EnhWk", "5_EnhA1", "6_EnhA2", "7_Tx", "9_ReprPC", "10_ZNF/Rpts/Het", "8_Quies", "X_NoConsistentState")


## Load in regions that have consistent state calls across replicates

#strategy <- "CONCORDANT"
strategy <- "MAJORITY"

if(strategy == "CONCORDANT"){
  message("Reading CONCORDANT regions...")
  hce_reg <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/consistent/HCE_fuso_correlative_concordant_10_dense_ucsc.bed", header=F,
                        stringsAsFactors=F, sep="\t", skip=1)
  hcae_reg <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/consistent/HCAE_fuso_correlative_concordant_10_dense_ucsc.bed", header=F,
                         stringsAsFactors=F, sep="\t", skip=1)
}else if(strategy == "MAJORITY"){
  message("Reading MAJORITY regions...")
  hce_reg <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/consistent/HCE_fuso_correlative_majority_10_dense_ucsc.bed", header=F,
                        stringsAsFactors=F, sep="\t", skip=1)
  hcae_reg <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/consistent/HCAE_fuso_correlative_majority_10_dense_ucsc.bed", header=F,
                         stringsAsFactors=F, sep="\t", skip=1)
}


names(hce_reg) <- c("chr", "start", "end", "state_change", "V5", "V6", "V7", "V8", "V9")
names(hcae_reg) <- c("chr", "start", "end", "state_change", "V5", "V6", "V7", "V8", "V9")

## get region sizes
hce_reg$region_length <- hce_reg$end - hce_reg$start
hcae_reg$region_length <- hcae_reg$end - hcae_reg$start

## count number of 200bp chunks in these regions (for calculating fraction of genome.)
hce_reg$numregions <- hce_reg$region_length/200
hcae_reg$numregions <- hcae_reg$region_length/200

## extract the start state (in unexposed sample) and end state (in Fn exposed sample)
## split the state change by "->". This gives [1] start [2] end
## but, it does that for each row. Unlisting gives [1] row1 start [2] row1 end [3] row2 start [4] row2 end
## so, then take the odds (starts) or evens (ends)
hce_reg$start_state <- unlist(strsplit(hce_reg$state_change, "->"))[c(TRUE,FALSE)]
hce_reg$end_state <- unlist(strsplit(hce_reg$state_change, "->"))[c(FALSE,TRUE)]
hcae_reg$start_state <- unlist(strsplit(hcae_reg$state_change, "->"))[c(TRUE,FALSE)]
hcae_reg$end_state <- unlist(strsplit(hcae_reg$state_change, "->"))[c(FALSE,TRUE)]


num_regions <- sum(hce_reg$numregions)

## summary of counts
hce_reg$condition <- "HCE"
hcae_reg$condition <- "HCAE"


## HCE
sb_hce <- do.call(rbind, by(hce_reg, interaction(hce_reg$start_state, hce_reg$end_state), function(x){
  return(data.frame(control_state = x[["start_state"]][1],
                    fuso_state = x[["end_state"]][1],
                    numregions = sum(as.numeric(x[["numregions"]]))))
}))
sb_hce$numControlRegions <- unlist(lapply(sb_hce$control_state, function(x){
  return(sum(subset(sb_hce, control_state==x)$numregions))
}))
sb_hce$percentFuso <- sb_hce$numregions / sb_hce$numControlRegions


sb_hce$control_state <- factor(sb_hce$control_state, levels = all_states)
sb_hce$fuso_state <- factor(sb_hce$fuso_state, levels = all_states)

## HCAE
sb_hcae <- do.call(rbind, by(hcae_reg, interaction(hcae_reg$start_state, hcae_reg$end_state), function(x){
  return(data.frame(control_state = x[["start_state"]][1],
                    fuso_state = x[["end_state"]][1],
                    numregions = sum(as.numeric(x[["numregions"]]))))
}))
sb_hcae$numControlRegions <- unlist(lapply(sb_hcae$control_state, function(x){
  return(sum(subset(sb_hcae, control_state==x)$numregions))
}))
sb_hcae$percentFuso <- sb_hcae$numregions / sb_hcae$numControlRegions

sb_hcae$control_state <- factor(sb_hcae$control_state, levels = all_states)
sb_hcae$fuso_state <- factor(sb_hcae$fuso_state, levels = all_states)

#### Remove state X (no consistent state) ####
### HCE
sb_hce_consistent <- subset(sb_hce, control_state != "X_NoConsistentState" & fuso_state != "X_NoConsistentState")
sb_hce_consistent$numControlRegions <- unlist(lapply(sb_hce_consistent$control_state, function(x){
  return(sum(subset(sb_hce_consistent, control_state==x)$numregions))
}))
sb_hce_consistent$percentFuso <- sb_hce_consistent$numregions / sb_hce_consistent$numControlRegions

### HCAE
sb_hcae_consistent <- subset(sb_hcae, control_state != "X_NoConsistentState" & fuso_state != "X_NoConsistentState")
sb_hcae_consistent$numControlRegions <- unlist(lapply(sb_hcae_consistent$control_state, function(x){
  return(sum(subset(sb_hcae_consistent, control_state==x)$numregions))
}))
sb_hcae_consistent$percentFuso <- sb_hcae_consistent$numregions / sb_hcae_consistent$numControlRegions








### Sankey diagrams

state_colours = list("1_TssA" = "#006d2c",
                     "2_TssFlnkU" = "#41ae76",
                     "3_TssFlnkD" = "#99d8c9",
                     "4_EnhWk" = "#8c6bb1",
                     "5_EnhA1" = "#c51b7d",
                     "6_EnhA2" = "#f1b6da",
                     "7_Tx" = "#b30000",
                     "8_Quies" = "#bdbdbd",
                     "9_ReprPC" = "#2171b5",
                     "10_ZNF/Rpts/Het" = "#6baed6",
                     "X_NoConsistentState" = "#737373")

all_states_rev <- rev(c("1_TssA", "2_TssFlnkU", "3_TssFlnkD", "4_EnhWk", "5_EnhA1", "6_EnhA2", "7_Tx", "9_ReprPC", "10_ZNF/Rpts/Het", "8_Quies", "X_NoConsistentState"))

sb_hce_consistent$color <- unlist(apply(sb_hce_consistent, 1, function(x){
  return(state_colours[x[["control_state"]][1]])
}))

sb_hce_consistent$control_state <- factor(sb_hce_consistent$control_state, levels=all_states_rev)
sb_hce_consistent$fuso_state <- factor(sb_hce_consistent$fuso_state, levels=all_states_rev)


## Weighted block sizes.
hcae_control_reg <- sum(unique(sb_hcae_consistent$numControlRegions))
hce_control_reg <- sum(unique(sb_hce_consistent$numControlRegions))

sb_hcae_consistent$percentAll <- sb_hcae_consistent$numregions / hcae_control_reg
sb_hce_consistent$percentAll <- sb_hce_consistent$numregions / hce_control_reg

## calculate fraction of genome the plots represent
hcae_control_reg / num_regions
hce_control_reg / num_regions

pdf(paste("figures/manuscript/HCAE_stateDynamics_noQuies_",strategy,"_alluvial_v3_scaled.pdf",sep=""), width=2, height = 8)
alluvial(sb_hcae_consistent[,c("control_state", "fuso_state")],
         freq=sb_hcae_consistent$percentAll,
         #hide=sb_hcae_consistent$percentFuso < 0.0001,
         col = sb_hcae_consistent$color,
         border = sb_hcae_consistent$color,
         gap.width=0.1,
         xw=0.2,
         cw=0.02,
         alpha=0.6,
         blocks=T)
dev.off()

pdf(paste("figures/manuscript/HCE_stateDynamics_noQuies_",strategy,"_alluvial_v3_scaled.pdf",sep=""), width=2, height = 8)
alluvial(sb_hce_consistent[,c("control_state", "fuso_state")],
         freq=sb_hce_consistent$percentAll,
         #hide=sb_hcae_consistent$percentFuso < 0.0001,
         col = sb_hce_consistent$color,
         border = sb_hce_consistent$color,
         gap.width=0.1,
         xw=0.2,
         cw=0.02,
         alpha=0.6,
         blocks=T)
dev.off()



#### Stats for Manuscript ####

## percent of genome with a consistent state called
## hce unexposed
sum(subset(hce_reg, start_state != "X_NoConsistentState")$numregions)/num_regions
## hce exposed
sum(subset(hce_reg, end_state != "X_NoConsistentState")$numregions)/num_regions
## hcae unexposed
sum(subset(hcae_reg, start_state != "X_NoConsistentState")$numregions)/num_regions
## hcae exposed
sum(subset(hcae_reg, end_state != "X_NoConsistentState")$numregions)/num_regions

## percent of genome with a consistent state called in both unexposed and exposed
## hce
sum(subset(hce_reg, start_state != "X_NoConsistentState" & end_state != "X_NoConsistentState")$numregions)/num_regions
## hcae
sum(subset(hcae_reg, start_state != "X_NoConsistentState" & end_state != "X_NoConsistentState")$numregions)/num_regions

## percent of genome with state change in HCAE
sum(subset(sb_hcae_consistent, control_state != fuso_state)$numregions)/num_regions

## numbers of windows
sum(subset(sb_hce_consistent, control_state != fuso_state)$numregions)
sum(subset(sb_hcae_consistent, control_state != fuso_state)$numregions)

## percent of genome without state change in HCE
sum(subset(sb_hce_consistent, control_state == fuso_state)$numregions)/num_regions










## select subset where there are changes, for annotating to genes.

hce_reg_sub <- subset(hce_reg, start_state != end_state & !grepl("X_NoConsistentState",state_change, fixed=T))
hcae_reg_sub <- subset(hcae_reg, start_state != end_state & !grepl("X_NoConsistentState",state_change, fixed=T))


hce_reg_sub$condition <- "HCE"
hcae_reg_sub$condition <- "HCAE"

reg_sub <- rbind(hce_reg_sub, hcae_reg_sub)

ggplot(reg_sub, aes(x=region_length, color=condition)) +
  geom_density() +
  facet_wrap(~state_change, scales="free") +
  theme_bw()

## number of regions with state change, period.
by(reg_sub, reg_sub$condition, nrow)

## write these

write.table(subset(reg_sub, condition=="HCE")[,1:9], paste("data/processed/chip-seq/ChromHMM/v4_hce-hcae/statechanges_v2/HCE_fusoConcordant_",strategy,"_withQuies_210301.bed",sep=""), col.names = F, sep="\t", row.names=F, quote = F)
write.table(subset(reg_sub, condition=="HCAE")[,1:9], paste("data/processed/chip-seq/ChromHMM/v4_hce-hcae/statechanges_v2/HCAE_fusoConcordant_",strategy,"_withQuies_210301.bed",sep=""), col.names = F, sep="\t", row.names=F, quote = F)


## POST PROCESS THE ABOVE USING 
# python 03a_breakup_chromhmm_regions_into_200bp_beds.py statechanges_v2/HCE_fusoConcordant_MAJORITY_withQuies_210301.bed statechanges_v3/HCE_fusoConcordant_MAJORITY_withQuies_210301_windows_set



#### Run GREAT externally, and then filter further ####


hce_mapping <- data.frame()
for(i in c(1:5)){
  hce_mapping <- rbind(hce_mapping, read.table(paste("data/processed/chip-seq/ChromHMM/v4_hce-hcae/statechanges_v3/HCE_fusoConcordant_MAJORITY_withQuies_210301_windows_set_", i, ".genes.txt", sep=""), sep="\t", stringsAsFactors = F))
}
names(hce_mapping) <- c("gene", "regions")

hcae_mapping <- data.frame()
for(i in c(1:12)){
  hcae_mapping <- rbind(hcae_mapping, read.table(paste("data/processed/chip-seq/ChromHMM/v4_hce-hcae/statechanges_v3/HCAE_fusoConcordant_MAJORITY_withQuies_210301_windows_set_", i, ".genes.txt", sep=""), sep="\t", stringsAsFactors = F))
}
names(hcae_mapping) <- c("gene", "regions")



state_scores = list("1_TssA" = 2,
                     "2_TssFlnkU" = 1,
                     "3_TssFlnkD" = 1,
                     "4_EnhWk" = 1,
                     "5_EnhA1" = 1,
                     "6_EnhA2" = 1,
                     "7_Tx" = 2,
                     "8_Quies" = 0,
                     "9_ReprPC" = -1,
                     "10_ZNF/Rpts/Het" = -1)

## calculate summed score for a gene.
summedStateChangeScore <- function(rs){
  score <- 0
  for(r in strsplit(rs, ", ")[[1]]){
    rsplit <- strsplit(r, " ")[[1]]
    sc <- strsplit(rsplit[1], "-chr")[[1]][1]   ## state change
    sc <- strsplit(sc, "->")[[1]]
    
    score <- score + (as.numeric(state_scores[sc[2]]) - as.numeric(state_scores[sc[1]]))
  }
  return(score)
}



### Set 1: All state changes
hce_mapping$set1_nes <- unlist(lapply(hce_mapping$regions, summedStateChangeScore))
hcae_mapping$set1_nes <- unlist(lapply(hcae_mapping$regions, summedStateChangeScore))

write.table(subset(hce_mapping, regions != "")[,c("gene","regions","set1_nes")], paste("data/processed/chip-seq/ChromHMM/v4_hce-hcae/statechanges_v3/HCE_fusoConcordant_",strategy,"_withQuies_set1.genes.txt",sep=""), col.names = F, sep="\t", row.names=F, quote = F)
write.table(subset(hcae_mapping, regions != "")[,c("gene","regions","set1_nes")], paste("data/processed/chip-seq/ChromHMM/v4_hce-hcae/statechanges_v3/HCAE_fusoConcordant_",strategy,"_withQuies_set1.genes.txt",sep=""), col.names = F, sep="\t", row.names=F, quote = F)


### Set 1proximal: Subset of proximal state changes
filterProximal <- function(rs){
  rsn <- unlist(lapply(strsplit(rs, ", ")[[1]], function(r){
    rsplit <- strsplit(r, " ")[[1]]
    sc <- strsplit(rsplit[1], "-chr")[[1]][1]
    dst <- rsplit[2]
    dst <- as.numeric(substr(dst, 2, nchar(dst)-1))
    if(dst > -5000 & dst < 1000){
      return(r)
    }else{
      return("")
    }
  }))
  rsn <- rsn[rsn!=""]
  return(paste(rsn, collapse=", "))
}

hce_mapping$set1prox <- unlist(lapply(hce_mapping$regions, filterProximal))
hcae_mapping$set1prox <- unlist(lapply(hcae_mapping$regions, filterProximal))

hce_mapping$set1prox_nes <- unlist(lapply(hce_mapping$set1prox, summedStateChangeScore))
hcae_mapping$set1prox_nes <- unlist(lapply(hcae_mapping$set1prox, summedStateChangeScore))

write.table(subset(hce_mapping, set1prox != "")[,c("gene","set1prox","set1prox_nes")], paste("data/processed/chip-seq/ChromHMM/v4_hce-hcae/statechanges_v3/HCE_fusoConcordant_",strategy,"_withQuies_set1prox.genes.txt",sep=""), col.names = F, sep="\t", row.names=F, quote = F)
write.table(subset(hcae_mapping, set1prox != "")[,c("gene","set1prox","set1prox_nes")], paste("data/processed/chip-seq/ChromHMM/v4_hce-hcae/statechanges_v3/HCAE_fusoConcordant_",strategy,"_withQuies_set1prox.genes.txt",sep=""), col.names = F, sep="\t", row.names=F, quote = F)



## check correlations with log2FC

hce_nes <- hce_mapping[,c("gene", "set1_nes", "set1prox_nes")]
hcae_nes <- hcae_mapping[,c("gene", "set1_nes", "set1prox_nes")]

## load expression data
hce_ge <- read.table("data/processed/DEseq/v4/HCE_7.1.withNames.deseq.tsv", sep="\t", stringsAsFactors = F, header=T)
hcae_ge <- read.table("data/processed/DEseq/v4/HCAE_7.1.withNames.deseq.tsv", sep="\t", stringsAsFactors = F, header=T)

## merge with epigenetic scores
hce_nes_ge <- merge(hce_nes, hce_ge, by.x="gene", by.y="human_gene_name", all=T)
hcae_nes_ge <- merge(hcae_nes, hcae_ge, by.x="gene", by.y="human_gene_name", all=T)

## set missing values (because no states associated) to 0
hce_nes_ge$set1_nes[is.na(hce_nes_ge$set1_nes)] <- 0
hcae_nes_ge$set1_nes[is.na(hcae_nes_ge$set1_nes)] <- 0

hce_nes_ge$set1prox_nes[is.na(hce_nes_ge$set1prox_nes)] <- 0
hcae_nes_ge$set1prox_nes[is.na(hcae_nes_ge$set1prox_nes)] <- 0


## select subset of genes, and select random subset of genes (same number)

## given a number of genes, do repeated random resampling to report p-value
random_test <- function(y_col, x_col, num_runs, obs_val, num_genes){
  test_dat <- data.frame(y=y_col, x=x_col)
  metric_res <- c()
  for(i in c(1:num_runs)){
    # do random subset
    sub_dat <- test_dat[sample(c(1:nrow(test_dat)), num_genes),]
    # Calculate metric
    m <- cor(sub_dat$x, sub_dat$y, use="complete", method="pearson")
    # add metric to metric_res
    metric_res <- c(metric_res, m)
  }
  # calculate proportion of runs with a metric more extreme than obs_val
  if(obs_val > 0){
    #print(length(metric_res[metric_res > obs_val]))
    #print(length(metric_res))
    #print(metric_res)
    p <- length(metric_res[metric_res > obs_val])/length(metric_res)
  }else{
    p <- length(metric_res[metric_res < obs_val])/length(metric_res)
  }
  # return that proportion
  return(p)
}

## load genes from HCE×HCAE scatterplot.
red_genes <- read.table("data/interim/genelists/CONSERVED-UP-red_genelist.csv", header=T, sep=",", stringsAsFactors = F)
names(red_genes) <- c("X", "gene")
orange_genes <- read.table("data/interim/genelists/HCAE-UP-orange_genelist.csv", header=T, sep=",", stringsAsFactors = F)
names(orange_genes) <- c("X", "gene")
purple_genes <- read.table("data/interim/genelists/HCE-DOWN-purple_genelist.csv", header=T, sep=",", stringsAsFactors = F)
names(purple_genes) <- c("X", "gene")
blue_genes <- read.table("data/interim/genelists/CONSERVED-DOWN-blue_genelist.csv", header=T, sep=",", stringsAsFactors = F)
names(blue_genes) <- c("X", "gene")
green_genes <- read.table("data/interim/genelists/HCAE-DOWN-green_genelist.csv", header=T, sep=",", stringsAsFactors = F)
names(green_genes) <- c("X", "gene")
pink_genes <- read.table("data/interim/genelists/HCE-UP-pink_genelist.csv", header=T, sep=",", stringsAsFactors = F)
names(pink_genes) <- c("X", "gene")



## Red genes
red_hce_prox_cor <- with(subset(hce_nes_ge, gene %in% red_genes$gene), cor(set1prox_nes, log2FoldChange, use = "complete", method="pearson"))
red_hce_prox_bootp <- random_test(hce_nes_ge$log2FoldChange[hce_nes_ge$gene_biotype=="protein_coding"],
                                  hce_nes_ge$set1prox_nes[hce_nes_ge$gene_biotype=="protein_coding"],
                                  1000,
                                  red_hce_prox_cor,
                                  nrow(subset(hce_nes_ge, gene %in% red_genes$gene)))
red_hce_all_cor <- with(subset(hce_nes_ge, gene %in% red_genes$gene), cor(set1_nes, log2FoldChange, use = "complete", method="pearson"))
red_hce_all_bootp <- random_test(hce_nes_ge$log2FoldChange[hce_nes_ge$gene_biotype=="protein_coding"],
                                  hce_nes_ge$set1_nes[hce_nes_ge$gene_biotype=="protein_coding"],
                                  1000,
                                  red_hce_all_cor,
                                  nrow(subset(hce_nes_ge, gene %in% red_genes$gene)))
red_hcae_prox_cor <- with(subset(hcae_nes_ge, gene %in% red_genes$gene), cor(set1prox_nes, log2FoldChange, use = "complete", method="pearson"))
red_hcae_prox_bootp <- random_test(hcae_nes_ge$log2FoldChange[hcae_nes_ge$gene_biotype=="protein_coding"],
                                  hcae_nes_ge$set1prox_nes[hcae_nes_ge$gene_biotype=="protein_coding"],
                                  1000,
                                  red_hcae_prox_cor,
                                  nrow(subset(hcae_nes_ge, gene %in% red_genes$gene)))
red_hcae_all_cor <- with(subset(hcae_nes_ge, gene %in% red_genes$gene), cor(set1_nes, log2FoldChange, use = "complete", method="pearson"))
red_hcae_all_bootp <- random_test(hcae_nes_ge$log2FoldChange[hcae_nes_ge$gene_biotype=="protein_coding"],
                                 hcae_nes_ge$set1_nes[hcae_nes_ge$gene_biotype=="protein_coding"],
                                 1000,
                                 red_hcae_all_cor,
                                 nrow(subset(hcae_nes_ge, gene %in% red_genes$gene)))

## Orange genes
orange_hce_prox_cor <- with(subset(hce_nes_ge, gene %in% orange_genes$gene), cor(set1prox_nes, log2FoldChange, use = "complete", method="pearson"))
orange_hce_prox_bootp <- random_test(hce_nes_ge$log2FoldChange[hce_nes_ge$gene_biotype=="protein_coding"],
                                  hce_nes_ge$set1prox_nes[hce_nes_ge$gene_biotype=="protein_coding"],
                                  1000,
                                  orange_hce_prox_cor,
                                  nrow(subset(hce_nes_ge, gene %in% orange_genes$gene)))
orange_hce_all_cor <- with(subset(hce_nes_ge, gene %in% orange_genes$gene), cor(set1_nes, log2FoldChange, use = "complete", method="pearson"))
orange_hce_all_bootp <- random_test(hce_nes_ge$log2FoldChange[hce_nes_ge$gene_biotype=="protein_coding"],
                                 hce_nes_ge$set1_nes[hce_nes_ge$gene_biotype=="protein_coding"],
                                 1000,
                                 orange_hce_all_cor,
                                 nrow(subset(hce_nes_ge, gene %in% orange_genes$gene)))
orange_hcae_prox_cor <- with(subset(hcae_nes_ge, gene %in% orange_genes$gene), cor(set1prox_nes, log2FoldChange, use = "complete", method="pearson"))
orange_hcae_prox_bootp <- random_test(hcae_nes_ge$log2FoldChange[hcae_nes_ge$gene_biotype=="protein_coding"],
                                   hcae_nes_ge$set1prox_nes[hcae_nes_ge$gene_biotype=="protein_coding"],
                                   1000,
                                   orange_hcae_prox_cor,
                                   nrow(subset(hcae_nes_ge, gene %in% orange_genes$gene)))
orange_hcae_all_cor <- with(subset(hcae_nes_ge, gene %in% orange_genes$gene), cor(set1_nes, log2FoldChange, use = "complete", method="pearson"))
orange_hcae_all_bootp <- random_test(hcae_nes_ge$log2FoldChange[hcae_nes_ge$gene_biotype=="protein_coding"],
                                  hcae_nes_ge$set1_nes[hcae_nes_ge$gene_biotype=="protein_coding"],
                                  1000,
                                  orange_hcae_all_cor,
                                  nrow(subset(hcae_nes_ge, gene %in% orange_genes$gene)))

## Purple genes
purple_hce_prox_cor <- with(subset(hce_nes_ge, gene %in% purple_genes$gene), cor(set1prox_nes, log2FoldChange, use = "complete", method="pearson"))
purple_hce_prox_bootp <- random_test(hce_nes_ge$log2FoldChange[hce_nes_ge$gene_biotype=="protein_coding"],
                                     hce_nes_ge$set1prox_nes[hce_nes_ge$gene_biotype=="protein_coding"],
                                     1000,
                                     purple_hce_prox_cor,
                                     nrow(subset(hce_nes_ge, gene %in% purple_genes$gene)))
purple_hce_all_cor <- with(subset(hce_nes_ge, gene %in% purple_genes$gene), cor(set1_nes, log2FoldChange, use = "complete", method="pearson"))
purple_hce_all_bootp <- random_test(hce_nes_ge$log2FoldChange[hce_nes_ge$gene_biotype=="protein_coding"],
                                    hce_nes_ge$set1_nes[hce_nes_ge$gene_biotype=="protein_coding"],
                                    1000,
                                    purple_hce_all_cor,
                                    nrow(subset(hce_nes_ge, gene %in% purple_genes$gene)))
purple_hcae_prox_cor <- with(subset(hcae_nes_ge, gene %in% purple_genes$gene), cor(set1prox_nes, log2FoldChange, use = "complete", method="pearson"))
purple_hcae_prox_bootp <- random_test(hcae_nes_ge$log2FoldChange[hcae_nes_ge$gene_biotype=="protein_coding"],
                                      hcae_nes_ge$set1prox_nes[hcae_nes_ge$gene_biotype=="protein_coding"],
                                      1000,
                                      purple_hcae_prox_cor,
                                      nrow(subset(hcae_nes_ge, gene %in% purple_genes$gene)))
purple_hcae_all_cor <- with(subset(hcae_nes_ge, gene %in% purple_genes$gene), cor(set1_nes, log2FoldChange, use = "complete", method="pearson"))
purple_hcae_all_bootp <- random_test(hcae_nes_ge$log2FoldChange[hcae_nes_ge$gene_biotype=="protein_coding"],
                                     hcae_nes_ge$set1_nes[hcae_nes_ge$gene_biotype=="protein_coding"],
                                     1000,
                                     purple_hcae_all_cor,
                                     nrow(subset(hcae_nes_ge, gene %in% purple_genes$gene)))

## Blue genes
blue_hce_prox_cor <- with(subset(hce_nes_ge, gene %in% blue_genes$gene), cor(set1prox_nes, log2FoldChange, use = "complete", method="pearson"))
blue_hce_prox_bootp <- random_test(hce_nes_ge$log2FoldChange[hce_nes_ge$gene_biotype=="protein_coding"],
                                     hce_nes_ge$set1prox_nes[hce_nes_ge$gene_biotype=="protein_coding"],
                                     1000,
                                     blue_hce_prox_cor,
                                     nrow(subset(hce_nes_ge, gene %in% blue_genes$gene)))
blue_hce_all_cor <- with(subset(hce_nes_ge, gene %in% blue_genes$gene), cor(set1_nes, log2FoldChange, use = "complete", method="pearson"))
blue_hce_all_bootp <- random_test(hce_nes_ge$log2FoldChange[hce_nes_ge$gene_biotype=="protein_coding"],
                                    hce_nes_ge$set1_nes[hce_nes_ge$gene_biotype=="protein_coding"],
                                    1000,
                                    blue_hce_all_cor,
                                    nrow(subset(hce_nes_ge, gene %in% blue_genes$gene)))
blue_hcae_prox_cor <- with(subset(hcae_nes_ge, gene %in% blue_genes$gene), cor(set1prox_nes, log2FoldChange, use = "complete", method="pearson"))
blue_hcae_prox_bootp <- random_test(hcae_nes_ge$log2FoldChange[hcae_nes_ge$gene_biotype=="protein_coding"],
                                      hcae_nes_ge$set1prox_nes[hcae_nes_ge$gene_biotype=="protein_coding"],
                                      1000,
                                      blue_hcae_prox_cor,
                                      nrow(subset(hcae_nes_ge, gene %in% blue_genes$gene)))
blue_hcae_all_cor <- with(subset(hcae_nes_ge, gene %in% blue_genes$gene), cor(set1_nes, log2FoldChange, use = "complete", method="pearson"))
blue_hcae_all_bootp <- random_test(hcae_nes_ge$log2FoldChange[hcae_nes_ge$gene_biotype=="protein_coding"],
                                     hcae_nes_ge$set1_nes[hcae_nes_ge$gene_biotype=="protein_coding"],
                                     1000,
                                     blue_hcae_all_cor,
                                     nrow(subset(hcae_nes_ge, gene %in% blue_genes$gene)))

## Green genes
green_hce_prox_cor <- with(subset(hce_nes_ge, gene %in% green_genes$gene), cor(set1prox_nes, log2FoldChange, use = "complete", method="pearson"))
green_hce_prox_bootp <- random_test(hce_nes_ge$log2FoldChange[hce_nes_ge$gene_biotype=="protein_coding"],
                                     hce_nes_ge$set1prox_nes[hce_nes_ge$gene_biotype=="protein_coding"],
                                     1000,
                                     green_hce_prox_cor,
                                     nrow(subset(hce_nes_ge, gene %in% green_genes$gene)))
green_hce_all_cor <- with(subset(hce_nes_ge, gene %in% green_genes$gene), cor(set1_nes, log2FoldChange, use = "complete", method="pearson"))
green_hce_all_bootp <- random_test(hce_nes_ge$log2FoldChange[hce_nes_ge$gene_biotype=="protein_coding"],
                                    hce_nes_ge$set1_nes[hce_nes_ge$gene_biotype=="protein_coding"],
                                    1000,
                                    green_hce_all_cor,
                                    nrow(subset(hce_nes_ge, gene %in% green_genes$gene)))
green_hcae_prox_cor <- with(subset(hcae_nes_ge, gene %in% green_genes$gene), cor(set1prox_nes, log2FoldChange, use = "complete", method="pearson"))
green_hcae_prox_bootp <- random_test(hcae_nes_ge$log2FoldChange[hcae_nes_ge$gene_biotype=="protein_coding"],
                                      hcae_nes_ge$set1prox_nes[hcae_nes_ge$gene_biotype=="protein_coding"],
                                      1000,
                                      green_hcae_prox_cor,
                                      nrow(subset(hcae_nes_ge, gene %in% green_genes$gene)))
green_hcae_all_cor <- with(subset(hcae_nes_ge, gene %in% green_genes$gene), cor(set1_nes, log2FoldChange, use = "complete", method="pearson"))
green_hcae_all_bootp <- random_test(hcae_nes_ge$log2FoldChange[hcae_nes_ge$gene_biotype=="protein_coding"],
                                     hcae_nes_ge$set1_nes[hcae_nes_ge$gene_biotype=="protein_coding"],
                                     1000,
                                     green_hcae_all_cor,
                                     nrow(subset(hcae_nes_ge, gene %in% green_genes$gene)))

## Pink genes
pink_hce_prox_cor <- with(subset(hce_nes_ge, gene %in% pink_genes$gene), cor(set1prox_nes, log2FoldChange, use = "complete", method="pearson"))
pink_hce_prox_bootp <- random_test(hce_nes_ge$log2FoldChange[hce_nes_ge$gene_biotype=="protein_coding"],
                                     hce_nes_ge$set1prox_nes[hce_nes_ge$gene_biotype=="protein_coding"],
                                     1000,
                                     pink_hce_prox_cor,
                                     nrow(subset(hce_nes_ge, gene %in% pink_genes$gene)))
pink_hce_all_cor <- with(subset(hce_nes_ge, gene %in% pink_genes$gene), cor(set1_nes, log2FoldChange, use = "complete", method="pearson"))
pink_hce_all_bootp <- random_test(hce_nes_ge$log2FoldChange[hce_nes_ge$gene_biotype=="protein_coding"],
                                    hce_nes_ge$set1_nes[hce_nes_ge$gene_biotype=="protein_coding"],
                                    1000,
                                    pink_hce_all_cor,
                                    nrow(subset(hce_nes_ge, gene %in% pink_genes$gene)))
pink_hcae_prox_cor <- with(subset(hcae_nes_ge, gene %in% pink_genes$gene), cor(set1prox_nes, log2FoldChange, use = "complete", method="pearson"))
pink_hcae_prox_bootp <- random_test(hcae_nes_ge$log2FoldChange[hcae_nes_ge$gene_biotype=="protein_coding"],
                                      hcae_nes_ge$set1prox_nes[hcae_nes_ge$gene_biotype=="protein_coding"],
                                      1000,
                                      pink_hcae_prox_cor,
                                      nrow(subset(hcae_nes_ge, gene %in% pink_genes$gene)))
pink_hcae_all_cor <- with(subset(hcae_nes_ge, gene %in% pink_genes$gene), cor(set1_nes, log2FoldChange, use = "complete", method="pearson"))
pink_hcae_all_bootp <- random_test(hcae_nes_ge$log2FoldChange[hcae_nes_ge$gene_biotype=="protein_coding"],
                                     hcae_nes_ge$set1_nes[hcae_nes_ge$gene_biotype=="protein_coding"],
                                     1000,
                                     pink_hcae_all_cor,
                                     nrow(subset(hcae_nes_ge, gene %in% pink_genes$gene)))



## red
red_hce_prox_cor
red_hce_prox_bootp
red_hcae_prox_cor
red_hcae_prox_bootp

red_hce_all_cor
red_hce_all_bootp
red_hcae_all_cor
red_hcae_all_bootp

## orange
orange_hce_prox_cor
orange_hce_prox_bootp
orange_hcae_prox_cor
orange_hcae_prox_bootp

orange_hce_all_cor
orange_hce_all_bootp
orange_hcae_all_cor
orange_hcae_all_bootp

## purple
purple_hce_prox_cor
purple_hce_prox_bootp
purple_hcae_prox_cor
purple_hcae_prox_bootp

purple_hce_all_cor
purple_hce_all_bootp
purple_hcae_all_cor
purple_hcae_all_bootp

## blue
blue_hce_prox_cor
blue_hce_prox_bootp
blue_hcae_prox_cor
blue_hcae_prox_bootp

blue_hce_all_cor
blue_hce_all_bootp
blue_hcae_all_cor
blue_hcae_all_bootp

## green
green_hce_prox_cor
green_hce_prox_bootp
green_hcae_prox_cor
green_hcae_prox_bootp

green_hce_all_cor
green_hce_all_bootp
green_hcae_all_cor
green_hcae_all_bootp

## pink
pink_hce_prox_cor
pink_hce_prox_bootp
pink_hcae_prox_cor
pink_hcae_prox_bootp

pink_hce_all_cor
pink_hce_all_bootp
pink_hcae_all_cor
pink_hcae_all_bootp



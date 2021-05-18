## Genome chromatin characteristics
## Dec 18 2020
## Scott Brown

setwd("/projects/sbrown_prj/190505_Fuso-Inv")


## Load data for Raw chipseq marks across genome

hce_f71_1_bin <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/binarized_merged/HCE_fn71_1_binary.txt", header=T, sep="\t", stringsAsFactors = F)
hce_f71_2_bin <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/binarized_merged/HCE_fn71_2_binary.txt", header=T, sep="\t", stringsAsFactors = F)
hce_f71_3_bin <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/binarized_merged/HCE_fn71_3_binary.txt", header=T, sep="\t", stringsAsFactors = F)
hce_neg_1_bin <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/binarized_merged/HCE_control_1_binary.txt", header=T, sep="\t", stringsAsFactors = F)
hce_neg_2_bin <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/binarized_merged/HCE_control_2_binary.txt", header=T, sep="\t", stringsAsFactors = F)
hce_neg_3_bin <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/binarized_merged/HCE_control_3_binary.txt", header=T, sep="\t", stringsAsFactors = F)
hcae_f71_1_bin <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/binarized_merged/HCAE_fn71_1_binary.txt", header=T, sep="\t", stringsAsFactors = F)
hcae_f71_2_bin <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/binarized_merged/HCAE_fn71_2_binary.txt", header=T, sep="\t", stringsAsFactors = F)
hcae_f71_3_bin <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/binarized_merged/HCAE_fn71_3_binary.txt", header=T, sep="\t", stringsAsFactors = F)
hcae_neg_1_bin <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/binarized_merged/HCAE_control_1_binary.txt", header=T, sep="\t", stringsAsFactors = F)
hcae_neg_2_bin <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/binarized_merged/HCAE_control_2_binary.txt", header=T, sep="\t", stringsAsFactors = F)
hcae_neg_3_bin <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/binarized_merged/HCAE_control_3_binary.txt", header=T, sep="\t", stringsAsFactors = F)

bins <- list(hce_f71_1_bin, hce_f71_2_bin, hce_f71_3_bin,
             hce_neg_1_bin, hce_neg_2_bin, hce_neg_3_bin,
             hcae_f71_1_bin, hcae_f71_2_bin, hcae_f71_3_bin,
             hcae_neg_1_bin, hcae_neg_2_bin, hcae_neg_3_bin)
names(bins) <- c("hce_f71_1", "hce_f71_2", "hce_f71_3",
                 "hce_neg_1", "hce_neg_2", "hce_neg_3",
                 "hcae_f71_1", "hcae_f71_2", "hcae_f71_3",
                 "hcae_neg_1", "hcae_neg_2", "hcae_neg_3")


num_bins <- nrow(hce_f71_1_bin) ## this is the same for all samples, as it is a product of reference genome size

## calculate total number of bins for each mark for each condition
bins_sums <- lapply(names(bins), function(n){
  x <- bins[[n]]
  d <- as.data.frame(t(apply(x, 2, sum)))
  d$condition <- n
  return(d)
})
names(bins_sums) <- names(bins)

## merge into one dataframe
bin_dat <- do.call(rbind, bins_sums)

## reshape for plotting
bin_dat <- melt(bin_dat, id.vars="condition")

## calculate fraction of genome for each
names(bin_dat) <- c("condition", "mark", "fracOfGenome")
bin_dat$fracOfGenome <- bin_dat$fracOfGenome / num_bins

bin_dat$condition <- apply(bin_dat, 1, function(x){
  return(substr(x[["condition"]], 1, nchar(x[["condition"]])-2))
})

bin_dat$cell_line <- apply(bin_dat, 1, function(x){
  return(strsplit(x[["condition"]][1], "_")[[1]][1])
})
bin_dat$fuso <- apply(bin_dat, 1, function(x){
  return(strsplit(x[["condition"]][1], "_")[[1]][2])
})
bin_dat$fuso <- factor(bin_dat$fuso, levels = c("neg", "f71"))
bin_dat$cell_line <- factor(bin_dat$cell_line, levels = c("hce", "hcae"))

## get mean and stdev of reps
bin_dat_agg <- do.call(rbind, by(bin_dat, interaction(bin_dat$condition, bin_dat$mark), function(x){
  return(data.frame(condition = x[["condition"]][1],
                    mark = x[["mark"]][1],
                    mn = mean(x[["fracOfGenome"]]),
                    sdv = sd(x[["fracOfGenome"]]),
                    cell_line = x[["cell_line"]][1],
                    fuso = x[["fuso"]][1]))
}))



(p <- ggplot(bin_dat_agg, aes(x=mark, y=mn, fill=fuso, group=fuso)) +
  geom_errorbar(aes(ymin = mn-sdv, ymax=mn+sdv), width = 0.5, position = position_dodge(width = 0.9)) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~cell_line, ncol=1) +
  theme_bw() +
  scale_fill_manual(values = c("grey", "#404040")) +
  xlab("Histone Mark") +
  ylab("Fraction of genome with each mark") +
  theme(panel.grid.minor = element_blank()))

ggsave("figures/manuscript/genomeWideChipMarks_v4.pdf", p, width=8, height = 6)

## t-test for fuso negative vs fuso positive for each mark
by(bin_dat, bin_dat$cell_line, function(cl){
  cl <- data.frame(cl)
  by(cl, cl$mark, function(x){
    x <- data.frame(x)
    t.test(x$fracOfGenome~x$fuso)$p.val
  })
})


## Load in data for chromatin states across genome

## chromhmm states across genome
hce_f71_1_hmm <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/renamed/HCE_fn71_1_10_dense_ucsc_200bpregions.tsv", sep="\t", header=T, stringsAsFactors = F)
hce_f71_2_hmm <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/renamed/HCE_fn71_2_10_dense_ucsc_200bpregions.tsv", sep="\t", header=T, stringsAsFactors = F)
hce_f71_3_hmm <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/renamed/HCE_fn71_3_10_dense_ucsc_200bpregions.tsv", sep="\t", header=T, stringsAsFactors = F)
hce_neg_1_hmm <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/renamed/HCE_control_1_10_dense_ucsc_200bpregions.tsv", sep="\t", header=T, stringsAsFactors = F)
hce_neg_2_hmm <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/renamed/HCE_control_2_10_dense_ucsc_200bpregions.tsv", sep="\t", header=T, stringsAsFactors = F)
hce_neg_3_hmm <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/renamed/HCE_control_3_10_dense_ucsc_200bpregions.tsv", sep="\t", header=T, stringsAsFactors = F)
hcae_f71_1_hmm <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/renamed/HCAE_fn71_1_10_dense_ucsc_200bpregions.tsv", sep="\t", header=T, stringsAsFactors = F)
hcae_f71_2_hmm <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/renamed/HCAE_fn71_2_10_dense_ucsc_200bpregions.tsv", sep="\t", header=T, stringsAsFactors = F)
hcae_f71_3_hmm <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/renamed/HCAE_fn71_3_10_dense_ucsc_200bpregions.tsv", sep="\t", header=T, stringsAsFactors = F)
hcae_neg_1_hmm <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/renamed/HCAE_control_1_10_dense_ucsc_200bpregions.tsv", sep="\t", header=T, stringsAsFactors = F)
hcae_neg_2_hmm <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/renamed/HCAE_control_2_10_dense_ucsc_200bpregions.tsv", sep="\t", header=T, stringsAsFactors = F)
hcae_neg_3_hmm <- read.table("data/processed/chip-seq/ChromHMM/v4_hce-hcae/renamed/HCAE_control_3_10_dense_ucsc_200bpregions.tsv", sep="\t", header=T, stringsAsFactors = F)


hmms <- list(hce_f71_1_hmm, hce_f71_2_hmm, hce_f71_3_hmm,
             hce_neg_1_hmm, hce_neg_2_hmm, hce_neg_3_hmm,
             hcae_f71_1_hmm, hcae_f71_2_hmm, hcae_f71_3_hmm,
             hcae_neg_1_hmm, hcae_neg_2_hmm, hcae_neg_3_hmm)
names(hmms) <- c("hce_f71_1", "hce_f71_2", "hce_f71_3",
                 "hce_neg_1", "hce_neg_2", "hce_neg_3",
                 "hcae_f71_1", "hcae_f71_2", "hcae_f71_3",
                 "hcae_neg_1", "hcae_neg_2", "hcae_neg_3")

num_bins <- nrow(hce_f71_1_hmm)

hmms_sums <- lapply(names(hmms), function(n){
  x <- hmms[[n]]
  d <- as.data.frame(t(apply(x, 2, sum)))
  d$condition <- n
  return(d)
})
names(hmms_sums) <- names(hmms)


hmm_dat <- do.call(rbind, hmms_sums)

hmm_dat <- melt(hmm_dat, id.vars="condition")
names(hmm_dat) <- c("condition", "state", "fracOfGenome")
hmm_dat$fracOfGenome <- hmm_dat$fracOfGenome / num_bins

hmm_dat$condition <- apply(hmm_dat, 1, function(x){
  return(substr(x[["condition"]], 1, nchar(x[["condition"]])-2))
})

hmm_dat$cell_line <- apply(hmm_dat, 1, function(x){
  return(strsplit(x[["condition"]][1], "_")[[1]][1])
})
hmm_dat$fuso <- apply(hmm_dat, 1, function(x){
  return(strsplit(x[["condition"]][1], "_")[[1]][2])
})
hmm_dat$fuso <- factor(hmm_dat$fuso, levels = c("neg", "f71"))
hmm_dat$cell_line <- factor(hmm_dat$cell_line, levels = c("hce", "hcae"))

## get mean and stdev of reps
hmm_dat_agg <- do.call(rbind, by(hmm_dat, interaction(hmm_dat$condition, hmm_dat$state), function(x){
  return(data.frame(condition = x[["condition"]][1],
                    state = x[["state"]][1],
                    mn = mean(x[["fracOfGenome"]]),
                    sdv = sd(x[["fracOfGenome"]]),
                    cell_line = x[["cell_line"]][1],
                    fuso = x[["fuso"]][1]))
}))

## y-axis scaled for first set of marks
(p1 <- ggplot(hmm_dat_agg, aes(x=state, y=mn, fill=fuso, group=fuso)) +
    geom_errorbar(aes(ymin = mn-sdv, ymax=mn+sdv), width = 0.5, position = position_dodge(width = 0.9)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(cell_line~state=="X8_Quies",scales="free", space="free_x") +
    theme_bw() +
    scale_fill_manual(values = c("grey", "#404040")) +
    xlab("Chromatin State") +
    ylab("Fraction of genome with each State") +
    theme(panel.grid.minor = element_blank()))

## y-axis scaled for Quies
(p2 <- ggplot(hmm_dat_agg, aes(x=state, y=mn, fill=fuso, group=fuso)) +
    geom_errorbar(aes(ymin = mn-sdv, ymax=mn+sdv), width = 0.5, position = position_dodge(width = 0.9)) +
    geom_bar(stat="identity", position="dodge") +
    facet_grid(cell_line~state=="X8_Quies",scales="free", space="free_x") +
    theme_bw() +
    scale_fill_manual(values = c("grey", "#404040")) +
    xlab("Chromatin State") +
    ylab("Fraction of genome with each State") +
    coord_cartesian(ylim=c(0,0.2)) +
    theme(panel.grid.minor = element_blank()))

ggsave("figures/manuscript/genomeWideChromatinStates_v4_full.pdf", p1, width=15, height = 6)
ggsave("figures/manuscript/genomeWideChromatinStates_v4_truncatedaxis.pdf", p2, width=15, height = 6)
## these are manually merged later

## t-tests
by(hmm_dat, hmm_dat$cell_line, function(cl){
  cl <- data.frame(cl)
  by(cl, cl$state, function(x){
    x <- data.frame(x)
    t.test(x$fracOfGenome~x$fuso)$p.val
  })
})

## deltas
by(hmm_dat_agg, hmm_dat_agg$cell_line, function(cl){
  cl <- data.frame(cl)
  by(cl, cl$state, function(x){
    x <- data.frame(x)
    (x$mn[x$fuso=="f71"] - x$mn[x$fuso=="neg"])*100
  })
})
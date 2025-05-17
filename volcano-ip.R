#ferran volcano

library(tidyverse)
library(dplyr)
library(limma)
library(gdata)
library(ggplot2)
library(gplots)
library("Rsamtools")
library("DESeq2")
library(Rsubread)
library(RColorBrewer)
library(edgeR)

setwd("/Users/mperchard/Dropbox/+LMS-db/papers-db/+2025-Garcia-Llagostera/")

## read data
volcano <- read.table("ip_volcano5.txt", header=T, sep="\t", na.strings="", fill=T)
volcano2 <-na.omit(volcano)

vol<-data.frame(logFC=volcano2$welch.diff.log, FDR=volcano2$welch.pval.log, label=(toupper(volcano2$Symbol)))
common <- read.table("common.txt", header=T)
sel <- c("RPA1", "XRCC5", "XRCC6", "PCNA")

com <-subset(vol, vol$label %in% common$common)
sel2 <-subset(vol, vol$label %in% sel)

vol$Color = "grey40"
vol$Color[vol$label == "ORF1"] = "blue4"
vol$Color[vol$label %in% common$common] = "red3"
vol$Color[vol$label %in% sel] = "purple"

# , color=vol$Color
ggplot(data=vol, aes(logFC,FDR)) +       
  geom_point(alpha=0.4, size=1, color=vol$Color) +
  geom_point(data=com, aes(logFC,FDR), color="green4") +
  geom_point(data=sel2, aes(logFC,FDR), color="purple") +
  geom_text(data=sel2, aes(logFC,FDR, label=label), vjust=-1.5, hjust=-1, size=2) +
  #         geom_text(data=nanog, aes(logFC,-log10(P.Value), label="Nanog"), vjust=-0.8, hjust=-0.8,color="red4") +
  #         geom_text(data=otx2, aes(logFC,-log10(P.Value), label="Otx2"), vjust=0, hjust=2.6, color="green4") + 
  #         geom_text(data=zic2, aes(logFC,-log10(P.Value), label="Zic2"), vjust=0, hjust=3.9, color="green4") + 
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(color="black", vjust=-0.35), 
        axis.title.y = element_text(color="black" , vjust=0.35))   +
  geom_hline(yintercept = 0, colour = "grey65") +
  geom_vline(xintercept = -10, colour = "grey65") +
  xlim(c(-10, 12)) + ylim(c(0, 10)) +
  xlab("-log10 Welch's t-test difference") + ylab("-log10 Welch's t-test p value") +
  theme(panel.background=element_rect(fill="white"), panel.grid.major=element_line(color="white")) +
  theme(axis.text.x=element_text(size=10, color="black"), axis.text.y=element_text(size=10, color="black")) +
  ggtitle("Volcano Plot - ORF1 IP-MS") +
  theme(plot.title = element_text(size=18,lineheight=.8, vjust=2))

###nicer one?
BiocManager::install("ggrepel")
library(ggrepel)


common <- read.table("common.txt", header=T)
sel <- c("RPA1", "XRCC5", "XRCC6", "PCNA")

volSP <- vol %>%
mutate(gene_type = case_when(label %in% common$common ~ "common",
                             label %in% sel ~ "ssdna",
                             TRUE ~ "ns"))

cols <- c("common" = "green4", "ssdna" = "purple", "ns" = "grey")
sizes <- c("common" = 2, "ssdna" = 2, "ns" = 1)
alphas <- c("common" = 1, "ssdna" = 1, "ns" = 0.1)

# pdf(paste0(plot_dir, today(), "_",  "Volcano_SENvsProlif.pdf"))


volSP <- volSP %>%
  arrange(factor(gene_type, levels = c("ns", "ssdna", "common")))

ggplot(data = volSP,
       aes(x = logFC,
           y = FDR)) +
  geom_point(aes(colour = gene_type),
             alpha = 0.8,
             shape = 16,
             size = 2) +
  geom_text_repel(aes(label = ifelse(gene_type %in% c("common", "ssdna"), as.character(label), "")),
                  nudge_y = 0.5, color = "black", max.overlaps = 200, size = 2) +
  geom_point(data = volSP[volSP$label=="ORF1",],
             shape = 21,
             size = 2,
             fill = "firebrick",
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(-2, 2),
             linetype = "dashed") +
  geom_label_repel(data = volSP[volSP$label=="ORF1",],
                   aes(label = "ORF1p"),
                   force = 2,
                   nudge_x = -3) +
  scale_colour_manual(values = cols) +
  # scale_x_continuous(breaks = c(seq(-8, 8, 2)),
  #                    limits = c(-8, 8)) +
  labs(title = "ORF1-IP",
       x = "Welch t-test difference",
       y = "-log10(Welch t-test p value)",
       colour = "Enrichment")
dev.off()
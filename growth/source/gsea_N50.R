library(BBmisc)
library(piano)
library(snow)
library(XLConnect)


# From differential gene expression analysis
doublingTime_N <- readRDS("~/Desktop/growth/data/results/diff_gene_exp_doublingTime_N.Rda")


# Top and bottom 50 samples
doublingTime_N50 <- doublingTime_N$`n=50`


# Gene set enrichment analysis
geneStat <- doublingTime_N50$logFC
names(geneStat) <- rownames(doublingTime_N50)

## We selected the C5 Biological Processes, KEGG, and Reactome gene sets for the pathwayDBs
c5bp <- "~/Desktop/growth/data/gene_sets/c5.bp.v6.2.symbols.gmt"
kegg <- "~/Desktop/growth/data/gene_sets/c2.cp.kegg.v6.2.symbols.gmt"
reactome <- "~/Desktop/growth/data/gene_sets/c2.cp.reactome.v6.2.symbols.gmt"

doGSEA <- function(geneStat, pathwayDB, nPerm, ncpus) {
  geneStatDir=NULL
  geneSetStat="gsea"
  
  gsc <- piano::loadGSC(pathwayDB)
  gsaRes<- piano::runGSA(geneLevelStats=geneStat,
                         geneSetStat=geneSetStat,
                         directions=geneStatDir,
                         gsc=gsc, nPerm=nPerm, ncpus=ncpus)
  
  gst <- GSAsummaryTable(gsaRes);rownames(gst) <- NULL
  if("p (non-dir.)" %in% colnames(gst)) {
    gst <- BBmisc::sortByCol(gst, c("p (non-dir.)"), asc = c(TRUE))
  } else {
    gst[, c("p", "fdr", "dir")] <- NA
    for(i in 1:nrow(gst)) {
      if(!is.na(gst[i,"p (dist.dir.up)"])) {
        gst[i, "p"] <- gst[i,"p (dist.dir.up)"]
        gst[i, "dir"] <- "up"
      } else {
        gst[i, "p"] <- gst[i,"p (dist.dir.dn)"]
        gst[i, "dir"] <- "dn"
      }
    }
    gst$fdr <- p.adjust(gst$p, "fdr")
    gst <- gst[, c("Name", "Genes (tot)", "Stat (dist.dir)",
                   "Genes (up)", "Genes (down)", "p", "fdr", "dir")]
    gst <- sortByCol(gst, c("fdr", "p"))
    rownames(gst) <- NULL
  }
  return(gst)
}

## Permute 10,000 times for all gene sets
gsea_N50_c5bp <- doGSEA(geneStat=geneStat, pathwayDB=c5bp, nPerm=10000, ncpus=4)
gsea_N50_kegg <- doGSEA(geneStat=geneStat, pathwayDB=kegg, nPerm=10000, ncpus=4)
gsea_N50_reactome <- doGSEA(geneStat=geneStat, pathwayDB=reactome, nPerm=10000, ncpus=4)

writeWorksheetToFile(file="~/Desktop/growth/data/results/GSEA_doublingTime_N50.xlsx",
                     data=list(gsea_N50_c5bp, gsea_N50_kegg, gsea_N50_reactome), sheet=c("C5 Biological Processes", "KEGG", "Reactome"))

## Get an idea for yourself: 10 most enriched pathways (both directions)
c5bp_up <- gsea_N50_c5bp[gsea_N50_c5bp$dir=="up", "Name"][1:10]
c5bp_dn <- gsea_N50_c5bp[gsea_N50_c5bp$dir=="dn", "Name"][1:10]

kegg_up <- gsea_N50_kegg[gsea_N50_kegg$dir=="up", "Name"][1:10]
kegg_dn <- gsea_N50_kegg[gsea_N50_kegg$dir=="dn", "Name"][1:10]

reactome_up <- gsea_N50_reactome[gsea_N50_reactome$dir=="up", "Name"][1:10]
reactome_dn <- gsea_N50_reactome[gsea_N50_reactome$dir=="dn", "Name"][1:10]

most_enriched_N50 <- list(C5BP=list(Up=c5bp_up, Down=c5bp_dn),
                          KEGG=list(Up=kegg_up, Down=kegg_dn),
                          Reactome=list(Up=reactome_up, Down=reactome_dn))

saveRDS(most_enriched_N50, file="~/Desktop/growth/data/results/most_enriched_N50.Rda")
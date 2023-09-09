#' multiNetEnrich
#'
#' @param multiGene a data.frame of multi-omics gene difference analysis results (pvalue).
#' Each row is a gene, and each column represents an omics dataset. 
#' @param network network
#' @param p restart probability
#' @param threshold threshold
#' @param n number of bins
#' @param combineMethod The method of combining pvalues.
#' @param nperm Number of permutations to do.
#' @param cutoff Pvalue threshold of differentially expressed genes.
#' @param pvalueCutoff Cutoff value of pvalue.
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", 
#' "bonferroni", "BH", "BY", "fdr", "none"
#' @param qvalueCutoff Cutoff of qvalue.
#' @param TERM2GENE user input annotation of TERM TO GENE mapping, 
#' a data.frame of 2 column with term and gene
#' @param TERM2NAME user input of TERM TO NAME mapping, 
#' a data.frame of 2 column with term and name
#' @param combineMethod The method of combining pvalues, one of 
#' "fisher", "edgington", "stouffer" and "Brown"(only used in ActivePathways method).
#' @param stoufferWeights weights of stouffer combine method.
#' @importFrom clusterProfiler compareCluster
#' @export
multiNetEnrich <- function(multiGene, network, p = 0, TERM2GENE = NULL,
                      TERM2NAME = NULL, threshold = 1e-9, n = 10, 
                      nperm = 100,  pvalueCutoff = 0.05, cutoff = 0.05,
                      pAdjustMethod = "BH", qvalueCutoff = 0.2,
                      combineMethod = "fisher",
                      stoufferWeights = NULL) {
    gene_list <- vector("list", ncol(multiGene))
    names(gene_list) <- colnames(multiGene)
    for (i in names(gene_list)) {
        gene_list[[i]] <- rownames(multiGene)[multiGene[, i] < cutoff]
    }
    compareClusterResult <- compareCluster(gene_list, fun = NetEnrich, network = network, 
                      p = p, TERM2GENE = TERM2GENE,
                      TERM2NAME = TERM2NAME, threshold = threshold, n = n, 
                      nperm = nperm,  pvalueCutoff = pvalueCutoff, 
                      pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff)
    enrichResultList <- vector("list", nrow(multiGene))
    names(enrichResultList) <- nrow(multiGene)
    enrichResultList <- lapply(gene_list, function(x) {
        NetEnrich(x, network = network, 
            p = p, TERM2GENE = TERM2GENE,
            TERM2NAME = TERM2NAME, threshold = threshold, n = n, 
            nperm = nperm,  pvalueCutoff = 1, 
            pAdjustMethod = pAdjustMethod, qvalueCutoff = 1)
    })

    # library(metap)   
    # resultdf <- matrix(0, nrow = enrichResultList[[1]], ncol = length(gene_list))
    # rownames(resultdf) <- enrichResultList[[1]]@result$ID
    # colnames(resultdf) <- names(gene_list) 
    # for (i in names(gene_list)) {
    #     ii <- match(rownames(resultdf), enrichResultList[[i]]@result$ID)
    #     resultdf[, i] <- enrichResultList[[i]]@result[ii, "pvalue"]
    # }     
   
    # pathway_meta <- rep(1, nrow(resultdf))
    # for (k in 1:length(pathway_meta)) {
    #     pvalues <- as.numeric(resultdf[k, ])
    #     pvalues2 <- pvalues[!is.na(pvalues)]
    #     if (length(pvalues2) < 2) {
    #         pathway1_meta[k] <- pvalues[1]
    #     }  else {
    #         pathway1_meta[k] <- sumlog(pvalues2)$p
    #     }       
    # }
    em <- combine_enricher(multiEm = enrichResultList, method = combineMethod, 
        stoufferWeights = stoufferWeights, pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)
    em@pAdjustMethod <- pAdjustMethod
    em@pvalueCutoff <- pvalueCutoff
    em@qvalueCutoff <- qvalueCutoff
    enrichResult <- get_enriched2(em)
    return(list(compareClusterResult = compareClusterResult, enrichResult = enrichResult))
}
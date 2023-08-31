#' multiNetEnrich
#'
#' @param genelist genelist
#' @param network network
#' @param p restart probability
#' @param threshold threshold
#' @param TERM2GENE TERM2GENE
#' @param TERM2NAME TERM2NAME
#' @param n n
#' @param combineMethod The method of combining pvalues.
#' @param nperm Number of permutations to do.
#' @export
multiNetEnrich <- function(multiGene, network, p = 0, TERM2GENE = NULL,
                      TERM2NAME = NULL, threshold = 1e-9, n = 10, 
                      nperm = 100,  pvalueCutoff = 0.05, 
                      pAdjustMethod = "BH", qvalueCutoff = 0.2,
                      combineMethod = "fisher",
                      stoufferWeights = NULL) {
    gene_list <- vector("list", nrow(multiGene))
    names(gene_list) <- colnames(multiGene)
    for (i in names(gene_list)) {
        gene_list[[i]] <- rownames(multiGene)[multiGene[, i] < 0.05]
    }
    compareClusterResult <- compareCluster(gene_list, fun = NetEnrich, network = network, 
                      p = p, TERM2GENE = NULL,
                      TERM2NAME = TERM2NAME, threshold = threshold, n = n, 
                      nperm = nperm,  pvalueCutoff = pvalueCutoff, 
                      pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff)
    enrichResultList <- vector("list", nrow(multiGene))
    names(enrichResultList) <- nrow(multiGene)
    enrichResultList <- lapply(gene_list, function(x) {
        NetEnrich(x, network = network, 
            p = p, TERM2GENE = NULL,
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
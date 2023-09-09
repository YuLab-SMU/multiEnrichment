#' Multi-omics ORA enrichment analysis
#'
#' @param multiGene a data.frame of multi-omics gene difference analysis results (pvalue).
#' Each row is a gene, and each column represents an omics dataset. 
#' @param cutoff Pvalue threshold of differentially expressed genes.
#' @param pvalueCutoff Cutoff value of pvalue.
#' @param pAdjustMethod one of "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none"
#' @param universe background genes
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of each geneSet for analyzing
#' @param qvalueCutoff Cutoff of qvalue.
#' @param TERM2GENE user input annotation of TERM TO GENE mapping,
#' a data.frame of 2 column with term and gene
#' @param TERM2NAME user input of TERM TO NAME mapping,
#' a data.frame of 2 column with term and name
#' @param combineMethod The method of combining pvalues.
#' @param ... Other parameters.
#' @noRd
multi_enricher <- function(multiGene,
                           cutoff = 0.05,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           universe = NULL,
                           minGSSize=10,
                           maxGSSize=500,
                           qvalueCutoff = 0.2,
                           TERM2GENE,
                           TERM2NAME = NA,
                           combineMethod = "fisher",
                           stoufferWeights = NULL,
                           ...) {

    if (class(multiGene) == "data.frame") {
        multiGene <- split(multiGene, multiGene$omic)
    }                           
    run_enricher <- function(df, ...) {
        genes <- df[df$pvalue < cutoff, "gene"]
        clusterProfiler::enricher(genes, TERM2GENE = TERM2GENE,
            pvalueCutoff = 1, qvalueCutoff = 1, ...)
    }
    multiEm <- lapply(multiGene, run_enricher, ...)
    em <- combine_enricher(multiEm = multiEm, method = combineMethod, 
        stoufferWeights = stoufferWeights, pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)
    em@pAdjustMethod <- pAdjustMethod
    em@pvalueCutoff <- pvalueCutoff
    em@qvalueCutoff <- qvalueCutoff
    get_enriched2(em)
}

#' Multi-omics GSEA enrichment analysis
#'
#' @param multiGene a list of differential result. 
#' Each is an data.frame of omic data containg pvalue and logFC. 
#' @param pvalueCutoff Cutoff value of pvalue.
#' @param pAdjustMethod one of "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none"
#' @param universe background genes
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of each geneSet for analyzing
#' @param TERM2GENE user input annotation of TERM TO GENE mapping,
#' a data.frame of 2 column with term and gene
#' @param TERM2NAME user input of TERM TO NAME mapping,
#' a data.frame of 2 column with term and name
#' @param combineMethod The method of combining pvalues.
#' @param ... Other parameters.
#' @noRd
multi_GSEA <- function(multiGene,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       universe = NULL,
                       minGSSize=10,
                       maxGSSize=500,
                       qvalueCutoff = 0.2,
                       TERM2GENE,
                       TERM2NAME = NA,
                       combineMethod = "fisher",
                       stoufferWeights = NULL,
                       ...) {
    if (class(multiGene) == "data.frame") {
        multiGene <- split(multiGene, multiGene$omic)
    }   
    run_GSEA <- function(df, ...) {
      genelist <- -sign(df$logFC) * log10(df$pvalue)
      names(genelist) <- df$gene
      genelist <- sort(genelist, decreasing = TRUE)
      clusterProfiler::GSEA(genelist, TERM2GENE = TERM2GENE,
           pvalueCutoff = 1)
    }
    multiEm <- lapply(multiGene, run_GSEA, ...)
    em <- combine_GSEA(multiEm, method = combineMethod, 
        stoufferWeights = stoufferWeights, pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)
    em@params$pAdjustMethod <- pAdjustMethod
    get_enriched2(em)
}



#' Combine multi-omics ORA enrichment analysis results
#'
#' @param multiEm multi-omics ORA enrichment analysis results.
#' @noRd
combine_enricher <- function(multiEm, method, stoufferWeights, 
                             pAdjustMethod, pvalueCutoff,
                             qvalueCutoff) {
    result <- combine_pvalue(multiEm = multiEm, method  = method, 
        stoufferWeights = stoufferWeights, pAdjustMethod = pAdjustMethod)
    # The same part of each omic
    em <- lapply(multiEm, as.data.frame)
    em <- do.call(rbind, em)
    em <- em[, c("ID", "Description", "BgRatio")]
    em <- unique(em)
    rownames(em) <- em$ID
    ID <- result$ID
    result$Description <- em[ID, "Description"]
    result$BgRatio    <- em[ID, "BgRatio"]
    result <- result[!is.na(result$pvalue), ]
    result <- result[result$pvalue < pvalueCutoff, ]
    result <- result[result$qvalue < qvalueCutoff, ]
    result <- result[order(result$pvalue), ]
    # The different part of each omic
    geneID <- vector2list(multiEm, "geneID", termID = result$ID)
    Count <- vector2list(multiEm, "Count", termID = result$ID)
    GeneRatio <- vector2list(multiEm, "GeneRatio", termID = result$ID)
    gene <- slot2list(multiEm, "gene")
    em1 <- multiEm[[1]]
    return(
        new("multiEnrichResult",
            result         = result,
            geneID         = geneID,
            Count          = Count,
            GeneRatio      = GeneRatio,
            pvalueCutoff   = em1@pvalueCutoff,
            pAdjustMethod  = em1@pAdjustMethod,
            qvalueCutoff   = em1@qvalueCutoff,
            organism       = em1@organism,
            ontology       = em1@ontology,
            gene           = gene,
            keytype        = em1@keytype,
            universe       = em1@universe,
            gene2Symbol    = em1@gene2Symbol,
            geneSets       = em1@geneSets,
            readable       = em1@readable,
            termsim        = em1@termsim,
            method         = em1@method
        )
    )
}




#' Combine pvalues
#'
#' @param multiEm multi-omics enrichment analysis results.
#' @param method The method of combining pvalues, one of
#' "fisher", "edgington" and "stouffer".
#' @param stoufferWeights weights of stouffer method.
#' @param pAdjustMethod one of "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none".
#' @noRd
combine_pvalue <- function(multiEm, method = "fisher", 
                           stoufferWeights = NULL, 
                           pAdjustMethod = "BH") {
    extraPvalue <- function(em) {
        em <- as.data.frame(em)
        data.frame(ID = em$ID, pvalue = em$pvalue)
    }
    result <- lapply(multiEm, extraPvalue)
    result2 <- lapply(result, function(x) x$ID)
    categorys <- unique(unlist(result2))
    values <- as.data.frame(matrix(NA, length(categorys), length(result)))
    rownames(values) <- categorys
    for (i in seq_len(length(result))) {
      values[result[[i]]$ID, i] <-  result[[i]]$pvalue
    }
    pvalue <- rep(0, nrow(values))
    for (i in seq_len(nrow(values))) {
        pp <- as.numeric(values[i, ])
        delete <- which(is.na(pp))
        if (length(delete) > 0) pp <- pp[-delete]
        if (method == "fisher") pvalue[i] <- metap::sumlog(pp)$p
        if (method == "edgington") pvalue[i] <- metap::sump(pp)$p
        if (method == "stouffer") pvalue[i] <- metap::sumz(pp, weights = stoufferWeights)
    }
    #names(pvalue) <- categorys
    p.adj <- p.adjust(pvalue, method=pAdjustMethod)
    qobj <- tryCatch(qvalue::qvalue(p=pvalue, lambda=0.05, pi0.method="bootstrap"),
        error=function(e) NULL)
    data.frame(ID = categorys, pvalue = pvalue,
               p.adjust = p.adj, qvalue = qobj$qvalues)
}

#' Combine multi-omics GSEA enrichment analysis results
#'
#' @param multiEm multi-omics GSEA enrichment analysis results.
#' @noRd
combine_GSEA <- function(multiEm, method, stoufferWeights, pAdjustMethod,
                         pvalueCutoff, qvalueCutoff) {
    result <- combine_pvalue(multiEm = multiEm, method  = method,
        stoufferWeights = stoufferWeights, pAdjustMethod = pAdjustMethod)
    # The same part of each omic
    
    em <- lapply(multiEm, as.data.frame)
    em <- do.call(rbind, em)
    em <- em[, c("ID", "Description")]
    em <- unique(em)
    rownames(em) <- em$ID
    
    ID <- result$ID
    result$Description <- em[ID, "Description"]
    result$setSize     <- em[ID, "setSize"]
    result <- result[!is.na(result$pvalue), ]
    result <- result[result$pvalue < pvalueCutoff, ]
    result <- result[result$qvalue < qvalueCutoff, ]
    result <- result[order(result$pvalue), ]
    # The different part of each omic
    enrichmentScore <- vector2list(multiEm, "enrichmentScore", termID = result$ID)
    NES <- vector2list(multiEm, "NES", termID = result$ID)
    rank <- vector2list(multiEm, "rank", termID = result$ID)
    leading_edge <- vector2list(multiEm, "leading_edge", termID = result$ID)
    core_enrichment <- vector2list(multiEm, "core_enrichment", termID = result$ID)
    setSize <- vector2list(multiEm, "setSize", termID = result$ID)
    
    geneList <- slot2list(multiEm, "geneList")
    permScores <- slot2list(multiEm, "permScores")
    em1 <- multiEm[[1]]
    return(
        new("multiGseaResult",
            result          = result,
            enrichmentScore = enrichmentScore,
            NES             = NES,
            rank            = rank,
            leading_edge    = leading_edge,
            core_enrichment = core_enrichment,
            geneList        = geneList,
            setSize         = setSize,
            organism        = em1@organism,
            setType         = em1@setType,
            geneSets        = em1@geneSets,
            keytype         = em1@keytype,
            permScores      = permScores,
            params          = em1@params,
            gene2Symbol     = em1@gene2Symbol,
            readable        = em1@readable,
            termsim         = em1@termsim,
            method          = em1@method
        )
    )
}

#' slot2list
#'
#' @param multiEm multi-omics enrichment analysis results.
#' @param slotName slot name.
#' @importFrom methods slot
#' @noRd
slot2list <- function(multiEm, slotName) {
    slotList <- vector("list", length(multiEm))
    for (i in seq_len(length(slotList))) {
        slotList[[i]] <- slot(multiEm[[i]], slotName)
    }
    slotList <- setNames(slotList, names(multiEm))
    return(slotList)
}

#' vector2list
#'
#' @param multiEm multi-omics enrichment analysis results.
#' @param vectorName vector name.
#' @noRd
vector2list <- function(multiEm, vectorName, termID) {
    slotList <- vector("list", length(multiEm))
    for (i in seq_len(length(slotList))) {
        slotListI <- multiEm[[i]][, vectorName]
        names(slotListI) <- multiEm[[i]][, "ID"]
        termID2 <- intersect(termID, names(slotListI))
        slotList[[i]] <- slotListI[termID2]
    }
    slotList <- setNames(slotList, names(multiEm))
    return(slotList)
}



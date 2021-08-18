#' multiEnrichResult2enrichResult
#'
#' @param em_enricher em_enricher
#' @param method method
#' @export
multiEnrichResult2enrichResult <- function(em_enricher, method = "union") {
    method <- match.arg(method, c("intersection", "union"))
    result <- as.data.frame(em_enricher)

    combine_col <- function(x) {
        terms <- unique(unlist(lapply(x, names)))
        geneDf <- matrix(NA, nrow  = length(terms), ncol = length(x))
        rownames(geneDf) <- terms
        for (i in seq_len(length(x))) {
            geneDf[names(x[[i]]), i] <- x[[i]]
        }
        return(geneDf)
    }
    geneID <- combine_col(em_enricher@geneID)
    geneID <- geneID[result$ID, ]


    combine_ID <- function(ID) {
        if (method == "union") {
            return(unique(unlist(strsplit(ID, "/"))))
        }
    }

    geneID2 <- apply(geneID, 1, combine_ID)
    Count2 <- lapply(geneID2, length)
    gene2 <- unique(unlist(em_enricher@gene))
    GeneRatio2 <- apply(data.frame(a=as.integer(Count2), b=length(gene2)), 1, function(x)
                       paste(x[1], "/", x[2], sep="", collapse="")
                       )
    geneID3 <- sapply(geneID2, function(i) paste(i, collapse="/"))
    result$geneID <- geneID3
    result$Count <- as.integer(Count2)
    result$GeneRatio <- GeneRatio2
    result <- result[, c("ID", "Description", "GeneRatio", "BgRatio",
        "pvalue", "p.adjust", "qvalue", "geneID", "Count")]
    rownames(result) <- result$ID
    enricherResult <- new("enrichResult",
                         result         = result,
                         pvalueCutoff   = em_enricher@pvalueCutoff,
                         pAdjustMethod  = em_enricher@pAdjustMethod,
                         qvalueCutoff   = em_enricher@qvalueCutoff,
                         gene           = gene2,
                         universe       = em_enricher@universe,
                         geneSets       = em_enricher@geneSets,
                         organism       = em_enricher@organism,
                         keytype        = em_enricher@keytype,
                         ontology       = em_enricher@ontology,
                         readable       = em_enricher@readable
             )
}


multiGseaResult2gseaResult <- function(em_gsea, method = "union") {
    
}










#' Multi-omics enrichment analysis
#'
#' @importFrom stats p.adjust
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom utils tail
#' @importFrom qvalue qvalue
#' @importFrom methods new
#' @param multiGene a data.frame of multi-omics gene difference analysis results.
#' @param method enrichment analysis method, one of "enricher"(the default)
#' "GSEA", "mitch", "ActivePathways".
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
#' @param combineMethod The method of combining pvalues, one of 
#' "fisher", "edgington", "stouffer" and "Brown"(only used in ActivePathways method).
#' @param stoufferWeights weights of stouffer combine method.
#' @param ... Other parameters.
#' @export
multiEnrichment <- function(multiGene,
                            method = "enricher",
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

    method <- match.arg(method,
        c("enricher", "GSEA", "mitch", "ActivePathways"))
    if (method == "mitch") {
        em <- mitch_method(multiGene, TERM2GENE, minGSSize = minGSSize, ...)
    }

    if (method == "ActivePathways") {
        em <- ActivePathways_method(multiGene = multiGene,
                                    pvalueCutoff = pvalueCutoff,
                                    pAdjustMethod = pAdjustMethod,
                                    universe = universe,
                                    minGSSize = minGSSize,
                                    maxGSSize = maxGSSize,
                                    qvalueCutoff = qvalueCutoff,
                                    TERM2GENE = TERM2GENE,
                                    TERM2NAME = TERM2NAME,
                                    combineMethod = combineMethod,
                                    ...)
    }

    if (method == "enricher") {
        em <- multi_enricher(multiGene = multiGene,
                             cutoff = cutoff,
                             pvalueCutoff = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             universe = universe,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             qvalueCutoff = qvalueCutoff,
                             TERM2GENE = TERM2GENE,
                             TERM2NAME = TERM2NAME,
                             combineMethod = combineMethod,
                             stoufferWeights = stoufferWeights,
                             ...)
    }

    if (method == "GSEA") {
        em <- multi_GSEA(multiGene = multiGene,
                         pvalueCutoff = pvalueCutoff,
                         pAdjustMethod = pAdjustMethod,
                         universe = universe,
                         minGSSize = minGSSize,
                         maxGSSize = maxGSSize,
                         qvalueCutoff = qvalueCutoff,
                         TERM2GENE = TERM2GENE,
                         TERM2NAME = TERM2NAME,
                         combineMethod = combineMethod,
                         stoufferWeights = stoufferWeights,
                         ...)
    }
    return(em)
}



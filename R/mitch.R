#' mitch method
#'
#' @param multiGene a data.frame of multi-omics gene difference analysis results (-log10 pvalue * sign(logFC)).
#' Each row is a gene, and each column represents an omics dataset. 
#' @param TERM2GENE user input annotation of TERM TO GENE mapping,
#' a data.frame of 2 column with term and gene
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param ... Other parameters.
#' @noRd
mitch_method <- function(multiGene, TERM2GENE, minGSSize, ...) {
    PATHID2EXTID <- split(as.character(TERM2GENE[,2]), as.character(TERM2GENE[,1]))
    mitch::mitch_calc(x = multiGene, genesets = PATHID2EXTID,
        minsetsize = minGSSize, ...)
}

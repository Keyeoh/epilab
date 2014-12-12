#'
#' Batch effects plot
#'
#' A plot showing the association among both qualitative and quantitative variables and a set of
#' precomputed surogate variables. It is useful for evaluating the relationship among variables and 
#' the possible batch effects present in the data. For now, association between a quantitative 
#' variable and a surogate variable is measured by Pearson Correlation. For qualitative variables, 
#' an analysis of variance is performed and a tick is shown for the significant combinations.
#'
#' @param pdata A data.frame containing the phenotype information. Usually constructed by calling
#' the pData method on a MethylSet or GenomicMethylSet, or the colData method on a 
#' SummarizedExperiment.
#' @param nvars A character vector containing the names of the qualitative variables to analyze. If
#' NULL, its value is set to all the variables which happen to be of character or factor type.
#' @param qvars A character vector containing the names of the quantitative variables to analyze. If
#' NULL, its value is set to all the variables of numeric type. Be careful that this may not work as
#' expected, as some variables can be mistaken, such as the Slide number in a 450k experiment.
#' @param sv A matrix containing the precomputed surogate variables.
#' @param alpha Significance level for the qualitative relationship test.
#' @return A data.frame containing the scores used to generate the plot.
#'
#' @importFrom plyr ddply
#' @importFrom scales muted
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 ggplot aes scale_fill_manual geom_tile xlab ylab ggtitle theme 
#' @importFrom ggplot2 scale_fill_gradient2 ggplotGrob
#' @export
#'
batchPlot <- function(pdata, sv, nvars=NULL, qvars=NULL, alpha=0.05) {

  getScore <- function(comb, pdata, sv, alpha) {
    varData <- pdata[, as.character(comb$varname)]
    svData <- sv[, as.character(comb$svname)]
    if (comb$type == 'n') {
      fit <- aov(svData ~ varData)
      afit <- anova(fit)
      pvalue <- afit$'Pr(>F)'[1]
      if (pvalue < alpha) {
        score <- 1 
      } else {
        score <- 0
      }
    } else if (comb$type == 'q') {
      score <- cor(varData, svData, use='complete.obs')
    } else {
      stop('Error')
    }
    comb$score <- score
    return(comb)
  }

  myData <- cbind(DataFrame(pdata), DataFrame(sv))

  if (is.null(colnames(sv))) {
    colnames(sv) <- paste0('SV', 1:ncol(sv))
  }

  if (is.null(nvars)) {
    pdataClasses <- sapply(pdata, class)
    nvars <- names(pdata)[pdataClasses == 'character' || pdataClasses == 'factor']
  }

  if (is.null(qvars)) {
    pdataClasses <- sapply(pdata, class)
    qvars <- names(pdata)[sapply(pdata, class) == 'numeric']
  }

  ncombs <- expand.grid(varname=nvars, svname=colnames(sv))
  ncombs$type <- 'n'
  qcombs <- expand.grid(varname=qvars, svname=colnames(sv))
  qcombs$type <- 'q'
  combs <- rbind(ncombs, qcombs) 

  scores <- ddply(combs, c('varname', 'svname'), getScore, pdata=pdata, sv=sv, alpha=alpha)

  nGraph <- (
             ggplot(scores[scores$type == 'n', ], aes_string(y='varname', x='svname')) 
             + theme_bw() 
             + geom_tile(aes(fill=as.factor(score)), size=3, colour='white') 
             + scale_fill_manual(values=c('white', muted('green')))
             + xlab('Surrogate Variable')
             + ylab('Variable Name')
             + ggtitle(bquote(paste('Qualitative Associations (', alpha, '=', .(alpha), ')')))
             + theme(
                     legend.position='none'
                     )
             )

  qGraph <- (
             ggplot(scores[scores$type == 'q', ], aes_string(y='varname', x='svname')) 
             + theme_bw() 
             + geom_tile(aes(fill=score), size=3, colour='white') 
             + scale_fill_gradient2(low=muted('blue'), high=muted('red'), midpoint=0, 
                                    limits=c(-1, 1))
             + xlab('Surrogate Variable')
             + ylab('Variable Name')
             + ggtitle('Quantitative Associations')
             )

  grid.arrange(nGraph, qGraph, ncol=2)

  return(scores)
}


library(ggplot2)
library(dplyr)
library(data.table)

#import::from(dplyr, "%>%")

#'
#' Function to generate manhattan plot, QQ plot or top SNPs.
#'
#' @param Path Path to summary statistics
#' @param title title on plots
#' @param sugg suggestive threshold, default is 5e-6
#' @param sigg genome-wide significant threshold, default is 5e-8
#' @param Manhattan True/False variable for whether or not to generate the manhattan plot
#' @param pvalColName Name of p value column in summary statistics file
#' @param basePairName Name of basepair column in summary statistics file
#' @param chrName Name of chromosome column in summary statistics file
#' @param chisqCol Name of ChiSQ column in summary statistics file 
#'
#' 
#' @return Generates a manhattan plot or QQ plot. If asked to, it will also store the top snps that are either genomewide significant or suggestive to separate files.
#'
#'
#' @export
#' 
manhattanPlot_noSave = function(Path,
                         title,
                         sugg = 1e-6,
                         sigg = 5e-8,
                         Manhattan  = T,
                         pvalColName = "P_BOLT_LMM_INF",
                         basePairName = "BP",
                         chrName = "CHR",
                         chisqCol = "CHISQ_LINREG",
                         SNP_distance = 100) {
  data = dplyr::as_tibble(data.table::fread(file = Path, header = T))
  data[["BP"]] = as.double(data[[basePairName]]) #needed to avoid integer overflow
  data[["CHR"]] = as.integer(data[[chrName]])
  data[[pvalColName]] = as.numeric(data[[pvalColName]])
  
  #if(sum(data$'A1FREQ' < 0.01) > 0) {
  # ok_maf = data$'A1FREQ' >= .01
  #  data = data[ok_maf,]
  #}
  
  if (sum(data[[pvalColName]] == 1) != 0) {
    cat("Removing", sum(data[[pvalColName]] == 1), "SNPs with a p-value of 1. \n")
    data = data[!(data[[pvalColName]] == 1), ]
  }
  
  if (sum(data[[pvalColName]] == 0) != 0) {
    cat("Removing", sum(data[[pvalColName]] == 0), "SNPs with a p-value of 0. \n")
    data = data[!(data[[pvalColName]] == 0),]
  }
  
  if (any(data$CHR > 22)) {
    cat("More than 22 chromosomes detected, removing chromosome 23 and above. \n")
    data = filter(data, CHR <= 22)
  }
  #this data will form the basis for extrating the top snps and the manhattan plot
  dataplyr = data %>%
    dplyr::group_by(CHR) %>%
    dplyr::summarise(chr_len = max(BP)) %>%
    dplyr::mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    dplyr::left_join(data, ., by = c("CHR" = "CHR")) %>%
    dplyr::arrange(CHR, BP) %>%
    dplyr::mutate(BPcum = BP + tot)
  
  
  if (sigg != 5e-8) {
    cat("Using Bonferoni corrected alpha\n")
    bonferoni_alpha = 0.05 / nrow(data)
  } else {
    cat("Using fixed significance level: 5e-8 \n")
    bonferoni_alpha = sigg
  }
  
  if (bonferoni_alpha > sugg) {
    cat("Bonferoni corrected significance level *higher* than suggested significance level!\n Bonferioni:", bonferoni_alpha, "\n suggeted:", sugg, "\n Terminating.." )
    stop()
  }
  
  if (Manhattan) {
    cat("Working on Manhattan plot. \n")
    
    
    
    
    dataplyr[["sig_P"]]   = ifelse(dataplyr[[pvalColName]] < bonferoni_alpha, T, F)
    dataplyr[["sugg_P"]]  = ifelse((bonferoni_alpha < dataplyr[[pvalColName]]) & (dataplyr[[pvalColName]] < sugg), T, F)
    dataplyr[["non_sig"]] = ifelse((dataplyr[[pvalColName]] >= sugg), T, F)
    
    dataplyr[["aboveunif"]] = runif(n = nrow(dataplyr))
    
    dataplyr = filter(dataplyr, !!as.symbol(pvalColName) < 0.05)  #removing any snps that are not nominally significant
    
    axisdataplyr = dataplyr %>% 
      dplyr::group_by(CHR) %>%
      dplyr::summarize(center = (max(BPcum) + min(BPcum))/ 2)
    
    non_sig_data = subset(dataplyr, non_sig == T)
    sig_P_data = subset(dataplyr, sig_P == T)
    sugg_P_data = subset(dataplyr, sugg_P == T)
    ymax = max(-log10(dataplyr[[pvalColName]]) + 2, -log10(bonferoni_alpha) + 2)
    ylims = c(-log10(0.05), ymax)
    
    column = ensym(pvalColName)
    # find top snps to mark separately
    cat("Distance used:", SNP_distance, "Kb \n")
    distance = SNP_distance * 1000
    
    column = as.name(pvalColName)
    
    sigg_store = data.frame()
    sigg_snps = dataplyr %>% dplyr::filter(bonferoni_alpha > !!column)
    
    while(nrow(sigg_snps) > 0) {
      cur_sigg_snp = sigg_snps[which.max(sigg_snps[[pvalColName]]),]
      sigg_store   = rbind(sigg_store, cur_sigg_snp)
      sigg_snps    = sigg_snps[abs(cur_sigg_snp$BPcum - sigg_snps$BPcum) > distance,]
    }
    
    
    
    manPlot = ggplot2::ggplot(dataplyr, aes(x = BPcum, y = -log10(!!column)) ) +
      ggplot2::geom_point(data = sugg_P_data, color = "orange", size = 2, alpha = .5) +
      ggplot2::geom_point(data = sig_P_data , color = "red", size = 2, alpha = .5) +
      ggplot2::geom_point(data = sigg_store , color = "red", size = 2, alpha = .5) +
      ggplot2::geom_point(data = non_sig_data, size = 2, alpha = .5, aes(color = as.factor(CHR))) +
      
      ggplot2::scale_color_manual(values = rep(c("gray", "darkgray"), 11)) +
      ggplot2::scale_x_continuous(labels = axisdataplyr$CHR, breaks = axisdataplyr$center) +
      ggplot2::scale_y_continuous(expand = c(0,0), limits = ylims) +
      
      ggplot2::ggtitle(title) +
      ggplot2::labs(x = "Chromosome", y = "-log10(P)") +
      
      ggplot2::geom_hline(yintercept = -log10(bonferoni_alpha)) +
      ggplot2::geom_hline(yintercept = -log10(sugg), linetype = "dashed") +
      
      ggplot2::theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(hjust = .5),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
    return(list("manhattanPlot" = manPlot))
  }
}


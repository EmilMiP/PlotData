library(ggplot2)
library(dplyr)
library(data.table)

#import::from(dplyr, "%>%")

#'
#' Function to generate manhattan plot, QQ plot or top SNPs.
#'
#' @param Path Path to summary statistics
#' @param imgName Name of the images to generate.
#' @param title title on plots
#' @param outDist path to the output destination folder for plots and top snps
#' @param sugg suggestive threshold, default is 5e-6
#' @param sigg genome-wide significant threshold, default is 5e-8
#' @param Manhattan True/False variable for whether or not to generate the manhattan plot
#' @param pvalColName Name of p value column in summary statistics file
#' @param basePairName Name of basepair column in summary statistics file
#' @param chrName Name of chromosome column in summary statistics file
#' @param QQ True/False variable for whether or not to generate the QQ plot
#' @param chisqCol Name of ChiSQ column in summary statistics file 
#' @param topSNPs True/False variable for whether or not to generate top snps files
#' @param SNP_distance distance used to determine independent snps
#'
#' 
#' @return Generates a manhattan plot or QQ plot. If asked to, it will also store the top snps that are either genomewide significant or suggestive to separate files.
#'
#'
#' @export
#' 
manhattanPlot = function(Path,
                         imgName,
                         title,
                         outDist,
                         sugg = 1e-6,
                         sigg = 5e-8,
                         Manhattan  = T,
                         pvalColName = "P_BOLT_LMM_INF",
                         basePairName = "BP",
                         chrName = "CHR",
                         QQ = F,
                         chisqCol = "CHISQ_LINREG",
                         topSNPs = F,
                         SNP_distance = 100) {
  data = dplyr::as_tibble(data.table::fread(file = Path, header = T))
  data[["BP"]] = as.double(data[[basePairName]])
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
    filtervec = (dataplyr[[pvalColName]] > .05) & (dataplyr$aboveunif >= .5 ) ### is this valid ? does it do what i think it does?
    
    dataplyr = dataplyr[!filtervec,]
    axisdataplyr = dataplyr %>% 
      dplyr::group_by(CHR) %>%
      dplyr::summarize(center = (max(BPcum) + min(BPcum))/ 2)
    
    non_sig_data = subset(dataplyr, non_sig == T)
    sig_P_data = subset(dataplyr, sig_P == T)
    sugg_P_data = subset(dataplyr, sugg_P == T)
    ymax = max(-log10(dataplyr[[pvalColName]]) + 2, -log10(bonferoni_alpha) + 2)
    ylims = c(0, ymax)
    
    column = ensym(pvalColName)
    
    manPlot = ggplot2::ggplot(dataplyr, aes(x = BPcum, y = -log10(!!column)) ) +
      ggplot2::geom_point(data = sugg_P_data, color = "orange", size = 2, alpha = .8) +
      ggplot2::geom_point(data = sig_P_data , color = "red", size = 2, alpha = .8) +
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
    manOut = paste(outDist, "/", imgName, ".png", sep = "")
    png(manOut, width = 1920, height = 1080)
    print(manPlot)
    dev.off()
  }
  
  if(QQ) {
    cat("Working on QQ plot. \n")
    
    qqDat = dplyr::tibble(obs_chi = sort(dataplyr[[chisqCol]]),
                   teo_chi = sort(rchisq(n = sum(!is.na(obs_chi)), df = 1))
                   )
    QQplot = ggplot2::ggplot(qqDat, aes(x = teo_chi, y = obs_chi)) +
      ggplot2::geom_points() +
      ggplot2::geom_abline(intercept = 0, slope = 1) + 
      
      ggplot2::ggtitle(title) + 
      ggplot2::labs(x = "theoretical Chisquare", y = "Observed Chisquare")
      
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        plot.title = element_text(hjust = .5),
        legend.position = "none",
        panel.border = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
      
    
    qqOut = paste(outDist, "/", imgName, "QQ.png", sep = "")
    
    png(qqOut, width = 1920, height = 1080)
    print(QQplot)
    dev.off()
  }
  
  if (topSNPs) {
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
    
    if (nrow(sigg_store) > 0) {
      sigg_file_dist = paste(outDist, "/", imgName, ".sigSNPs", sep = "")
      cat("Saving list of independent significant SNPs to ", sigg_file_dist, "\n")
      data.table::fwrite(as.data.frame(sigg_store), file = sigg_file_dist,
             sep = " ", 
             quote = F,
             na = "NA")
    }
    
    sugg_store = data.frame()
    sugg_snps = dataplyr %>% dplyr::filter((bonferoni_alpha < !! column) & (!!column < sugg) )
    #remove suggestive snps that are already close to a significant snp.
    
    BPs_to_remove = sigg_store$BPcum
    
    while (length(BPs_to_remove) > 0) {
      sugg_snps = sugg_snps[abs(BPs_to_remove[1] - sugg_snps$BPcum) > distance,]
      BPs_to_remove = BPs_to_remove[-1]
    }
    cat("exited remove significant snps before finding suggestive snps \n")
    while (nrow(sugg_snps) > 0) {
      cur_sugg_snp = sugg_snps[which.max(sugg_snps[[pvalColName]]),]
      sugg_store   = rbind(sugg_store, cur_sugg_snp)
      sugg_snps    = sugg_snps[abs(cur_sugg_snp$BPcum - sugg_snps$BPcum) > distance,]
    }
    cat("exited finding suggestive snps \n")
    
    if (nrow(sugg_store) > 0) {
      sugg_file_dist = paste(outDist, "/", imgName, ".suggSNPs", sep = "")
      cat("Saving list of independent suggestive SNPs to ", sugg_file_dist, "\n")
      data.table::fwrite(as.data.frame(sugg_store), file = sugg_file_dist,
             sep = " ", 
             quote = F,
             na = "NA")
    }
  }
}


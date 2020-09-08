library(ggplot2)
library(dplyr)


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
  data = as_tibble(fread(file = Path, header = T))
  data[["BP"]] = as.double(data[[BasePairName]])
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
    group_by(CHR) %>%
    summarise(chr_len = max(BP)) %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    left_join(data, ., by = c("CHR" = "CHR")) %>%
    arrange(CHR, BP) %>%
    mutate(BPcum = BP + tot)
  
  
  if (sigg != 5e-8) {
    cat("Using Bonferoni corrected alpha\n")
    bonferoni_alpha = 0.05 / nrow(data)
  } else {
    cat("Using fixed significance level: 5e-8 \n")
    bonferoni_alpha = sigg
  }
  
  if (bonferoni_alpha > sugg) {
    cat("Bonferoni corrected significance level *higher* than suggested significance level!\n Bonferioni:", bonferoni_alpha, "\n suggeted:", sugg, "\n Terminating.." )
    break()
  }
  
  if (noManhattan) {
    cat("Working on Manhattan plot. \n")
    
    
    
    
    dataplyr[["sig_P"]]   = ifelse(dataplyr[[pvalColName]] < bonferoni_alpha, T, F)
    dataplyr[["sugg_P"]]  = ifelse((bonferoni_alpha < dataplyr[[pvalColName]]) & (dataplyr[[pvalColName]] < sugg), T, F)
    dataplyr[["non_sig"]] = ifelse((dataplyr[[pvalColName]] >= sugg), T, F)
    
    dataplyr[["aboveunif"]] = runif(n = nrow(dataplyr))
    filtervec = (dataplyr[[pvalColName]] > .05) & (dataplyr$aboveunif >= .5 ) ### is this valid ? does it do what i think it does?
    
    dataplyr = dataplyr[!filtervec,]
    axisdataplyr = dataplyr %>% 
      group_by(CHR) %>%
      summarize(center = (max(BPcum) + min(BPcum))/ 2)
    
    non_sig_data = subset(dataplyr, non_sig == T)
    sig_P_data = subset(dataplyr, sig_P == T)
    sugg_P_data = subset(dataplyr, sugg_P == T)
    ymax = max(-log10(dataplyr[[pvalColName]]) + 2, -log10(bonferoni_alpha) + 2)
    ylims = c(0, ymax)
    
    column = ensym(pvalColName)
    
    manPlot = ggplot(dataplyr, aes(x = BPcum, y = -log10(!!column)) ) +
      geom_point(data = sugg_P_data, color = "orange", size = 2, alpha = .8) +
      geom_point(data = sig_P_data , color = "red", size = 2, alpha = .8) +
      geom_point(data = non_sig_data, size = 2, alpha = .5, aes(color = as.factor(CHR))) +
    
      scale_color_manual(values = rep(c("gray", "darkgray"), 11)) +
      scale_x_continuous(labels = axisdataplyr$CHR, breaks = axisdataplyr$center) +
      scale_y_continuous(expand = c(0,0), limits = ylims) +
      
      ggtitle(title) +
      labs(x = "Chromosome", y = "-log10(P)") +
      
      geom_hline(yintercept = -log10(bonferoni_alpha)) +
      goem_hline(yintercept = -log10(sugg), linetype = "dashed") +
      
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(hjust = .5),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
    manOut = paste(outDist, "/", imgName, ".png", sep = "")
    png(manOut, width = 950, height = 450)
    print(manPlot)
    dev.off()
  }
  
  if(QQ) {
    cat("Working on QQ plot. \n")
    
    qqDat = tibble(obs_chi = sort(dataplyr[[chisqCol]]),
                   teo_chi = sort(rchisq(n = sum(!is.na(obs_chi)), df = 1))
                   )
    QQplot = ggplot(qqDat, aes(x = teo_chi, y = obs_chi)) +
      geom_points() +
      geom_abline(intercept = 0, slope = 1) + 
      
      ggtitle(title) + 
      labs(x = "theoretical Chisquare", y = "Observed Chisquare")
      
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(hjust = .5),
        legend.position = "none",
        panel.border = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
      
    
    qqOut = paste(outDist, "/", imgName, "QQ.png", sep = "")
    
    png(qqOut, width = 950, height = 450)
    print(QQplot)
    dev.off()
  }
  
  if (topSNPs) {
    cat("Distance used:", SNP_distance, "Kb \n")
    distance = SNP_distance * 1000
    
    column = as.name(pvalColName)
    
    sigg_store = data.frame()
    sigg_snps = dataplyr %>% filter(bonferoni_alpha > !!column)
    
    while(nrow(sigg_snps) > 0) {
      cur_sigg_snp = sigg_snps[which.max(sigg_snps[[pvalColName]]),]
      sigg_store    = rbind(sigg_store, cur_sigg_snp)
      sigg_snps    = sigg_snps[abs(cur_sigg_snp$BPcum - sigg_snps$BPcum) > distance,]
    }
    
    if (nrow(sigg_store) > 0) {
      sigg_file_dist = paste(outDist, "/", imgName, ".sigSNPs", sep = "")
      cat("Saving list of independent significant SNPs to ", sigg_file_dist, "\n")
      fwrite(as.data.frame(sigg_store), file = sigg_file_dist,
             sep = " ", 
             quote = F,
             na = "NA")
    }
    
    sugg_store = data.frame()
    sugg_snps = dataplyr %>% filter((bonferoni_alpha < !! column) & (!!column < sugg) )
    #remove suggestive snps that are already close to a significant snp.
    
    BPs_to_remove = sigg_store$BPcum
    
    while( length(BPs_to_remove) > 0) {
      sugg_snps = sugg_snps[abs(BPs_to_remove[1] - sugg_snps$BPcum) > distance,]
      BPs_to_remove = BPs_to_remove[-1]
    }
    
    while (nrow(sugg_snps) > 0){
      cur_sugg_snp = sugg_snps[which.max(sugg_snps[[pvalColName]]),]
      sugg_store   = rbind(sugg_store, cur_sugg_snp)
      sugg_snps    = sugg_snps[abs(BPs_to_remove[1] - sugg_snps$BPcum) > distance,]
    }
    
    if (nrow(sugg_store) > 0) {
      sugg_file_dist = paste(outDist, "/", imgName, ".suggSNPs", sep = "")
      cat("Saving list of independent suggestive SNPs to ", sugg_file_dist, "\n")
      fwrite(as.data.frame(sugg_store), file = sugg_file_dist,
             sep = " ", 
             quote = F,
             na = "NA")
    }
  }
}


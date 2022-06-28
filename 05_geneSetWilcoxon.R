## calculate P-values in gene-set level and plot
library(ggplot2)
library(ggpubr)
library(stringr)


dir<-'GeneSet_wilcoxon'

wilcoxon <- function(geneset) {
  gsfile<-paste0(dir,'/geneset.',geneset,'.wilcoxon')
  
  all<-read.table(gsfile,sep="\t",header = T)
  
  gs<-subset(all,Group=='Gene set')
  bg<-subset(all,Group=='Background')
  
  wilcox.LoF<-wilcox.test(gs$LoF,bg$LoF,alternative = 'less')
  
  all$Group<-as.factor(all$Group)
  all<-all[order(all[,"LoF"],decreasing = T),]
  rank.LoF<-c(1:length(all$LoF))
  all$rank.LoF<-rank.LoF/length(all$LoF)
  
  p<-signif(wilcox.LoF$p.value,2)
  
  if (p<0.001){
    star<-'***'
  }else if (p<0.01){
    star<-'**'
  }else if (p<0.05){
    star<-'*'
  }else{
    star<-'n.s'
  }
  
  if (p<0.05){
    p.value<-parse(text = paste0('bolditalic(P)','~bold(" = ', p,'")','~bold("',star,'")') ) 
  }else{
    p.value<-parse(text = paste0('italic("P")',' == ', p) )
  }
  
  
  if (p<0.05 && length(gs$Gene) <50){ 
    b<-0.5
  }else{
    b<-0
  }
  
  ## plot 
  plot.LoF <- ggplot(all,aes(x = Group, y = rank.LoF, color = Group)) +
    geom_jitter(shape=16,size=3.5,colour="lightgrey",width=0.2,aes(alpha = Group)) + #
    geom_boxplot(linetype = "dashed",size=1,outlier.shape = NA,fill=NA) +
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), size=1,outlier.shape = NA,fill=NA) +
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), size=1, width = 0.4) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), size=1, width = 0.4) +
    scale_alpha_discrete(range = c(0, b))+
    scale_colour_manual(values=c("#3b8686", "#fe4365"))+
    stat_summary(geom="text", fun=quantile,
                 aes(label=sprintf("%1.2f", ..y..)),
                 position=position_nudge(x=-0.4), size=3.5, color='black') +
    scale_y_continuous(name = NULL,breaks = c(0,0.2,0.4,0.6,0.8,1.0)) +
    scale_x_discrete(name = NULL) +
    labs(title = str_wrap(geneset, width=20),subtitle = p.value) + 
    theme_classic() +
    theme(
      legend.position="none",
      plot.title = element_text(size=15,hjust = 0.5,color = 'black', face='bold'),
      plot.subtitle = element_text(size=12,hjust = 0.5,color = 'black'),
      axis.text.x = element_text(size=15, angle = 90, hjust = 0.5,vjust = 0.5,color = 'black'),
      axis.text.y = element_text(size=12, angle = 90, hjust = 0.5,vjust = 0.5,color = 'black'),
      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
      axis.ticks.x = element_line(color="black", size=1, linetype="solid"),
      axis.ticks.y = element_line(color="black", size=1, linetype="solid"),
      axis.ticks.length=unit(0.25, "cm")
    )
  return(plot.LoF)
}

p.Cell_cycle<-wilcoxon('Cell cycle')
p.Meiosis<-wilcoxon('Meiosis')
p.ProphaseI<-wilcoxon('Meiotic prophase I')
p.DNA_replication<-wilcoxon('DNA replication')
p.DNA_repair_total<-wilcoxon('DNA repair')
p.Checkpoint<-wilcoxon('DNA damage checkpoint')
p.BER_S2<-wilcoxon('Base excision repair')
p.NER<-wilcoxon('Nucleotide excision repair')
p.HRR_S2<-wilcoxon('Homologous recombination repair')
p.NHEJ<-wilcoxon('Non-homologous end-joining')
p.MMR<-wilcoxon('Mismatch repair')
p.FA_core_related<-wilcoxon('Fanconi anemia pathway')

plot.1<-ggarrange(p.Cell_cycle,p.Meiosis,p.ProphaseI,p.DNA_replication,p.Checkpoint,p.DNA_repair_total,
                  p.BER_S2,p.NER,p.HRR_S2,p.NHEJ,p.MMR,p.FA_core_related,
                  nrow = 2,ncol=6,widths = 2,heights=1)
file.1<-paste0(dir,"/wilcox.1.pdf")
# ggsave(file.1,plot = plot.1,width = 14, height = 12,limitsize = FALSE)


p.Mitochondrial_disease<-wilcoxon('Mitochondrial function')
p.Autophagy_total<-wilcoxon('Autophagy')
p.Longevity_regulating_pathway<-wilcoxon('Longevity regulating pathway')
p.Cellular_senescence<-wilcoxon('Cellular senescence')
p.GenAge_all<-wilcoxon('GenAge')
p.Oxidoreductase_activity<-wilcoxon('Oxidoreductase activity')
p.Response_to_oxidative_stress<-wilcoxon('Response to oxidative stress')
p.Oxidative_phosphorylation_KEGG<-wilcoxon('Oxidative phosphorylation')
p.Fatty_acid_metabolism<-wilcoxon('Fatty acid metabolism')
p.Ovarian_steroidogenesis<-wilcoxon('Ovarian steroidogenesis')
p.GnRH<-wilcoxon('GnRH signaling pathway')
p.Estrogen<-wilcoxon('Estrogen signaling pathway')

plot.2<-ggarrange(p.Mitochondrial_disease,p.Autophagy_total,p.Longevity_regulating_pathway,p.Cellular_senescence,p.GenAge_all,p.Oxidoreductase_activity,
                  p.Response_to_oxidative_stress,p.Oxidative_phosphorylation_KEGG,p.Fatty_acid_metabolism,p.Ovarian_steroidogenesis,p.GnRH,p.Estrogen,
                  nrow = 2,ncol=6,widths = 2,heights=1)
file.2<-paste0(dir,"/wilcox.2.pdf")
# ggsave(file.2,plot = plot.2,width = 14, height = 12,limitsize = FALSE)



p.PI3K<-wilcoxon('PI3K-Akt signaling pathway')
p.mTOR<-wilcoxon('mTOR signaling pathway')
p.FoxO<-wilcoxon('FoxO signaling pathway')
p.Hippo<-wilcoxon('Hippo signaling pathway')
p.TGF<-wilcoxon('TGF-beta signaling pathway')
p.Hedgehog<-wilcoxon('Hedgehog signaling pathway')
p.Notch<-wilcoxon('Notch signaling pathway')
p.Wnt<-wilcoxon('Wnt signaling pathway')
p.Ras<-wilcoxon('Ras signaling pathway')
p.ErbB<-wilcoxon('ErbB signaling pathway')
p.JAK<-wilcoxon('JAK-STAT signaling pathway')
p.p53<-wilcoxon('p53 signaling pathway')

plot.3<-ggarrange(p.PI3K,p.mTOR,p.FoxO,p.Hippo,p.TGF,p.Hedgehog,
                  p.Notch,p.Wnt,p.Ras,p.ErbB,p.JAK,p.p53,
                  nrow = 2,ncol=6,widths = 2,heights=1)
file.3<-paste0(dir,"/wilcox.3.pdf")
# ggsave(file.3,plot = plot.3,width = 14, height = 12,limitsize = FALSE)

file.combine<-paste0(dir,"/wilcox.combine.pdf")
plot.combine<-ggarrange(plot.1,plot.2,plot.3,
                        nrow = 3,ncol=1)
ggsave(file.combine,plot = plot.combine,width = 18, height = 24,limitsize = FALSE)



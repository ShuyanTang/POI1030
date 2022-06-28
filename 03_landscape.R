# install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(GenVisR)
library(circlize)

dir='.../landscape'

data2<-read.table(paste0(dir,"/matrix.txt"),sep="\t",header=T,row.names=1)
case<-read.table(paste0(dir,"/cases.txt"),sep="\t",header=T,row.names=1)
landscape<-read.table(paste0(dir,"/landscape.txt"),sep = "\t",header = T)

type<-data2[,1:10]
type$gene<-row.names(data2)
data<-merge(landscape,type,by="gene",all.x = T)
data$group<-factor(data$group,levels = c('Gonadogenesis', 'Oogenesis', 'Folliculogenesis and ovulation','Others'))
data<-arrange(data,group,desc(case_num))

out<-
  waterfall(data, 
            out='data',
            fileType = 'Custom',variant_class_order = levels(factor(data$variant_class)),
               geneOrder = factor(unique(data$gene)))

main<-out$main
out<-c()
main$order<-as.numeric(row.names(main))
main$sample<-gsub("-",".",main$sample)
main<-arrange(main,order)
sample_order<-unique(main$sample)


vc<-c("Nonsense","Frameshift indel","Splice","Start-loss","Inframe indel","Missense")
vc2<-c("Nonsense","Frameshift indel","Splice site/region","Start-loss","Inframe indel","Missense")
vc3<-c("Nonsense","Frameshift.indel","Splice","Start.loss","Inframe.indel","Missense")
col_vc<-c('#D25B55','#5DBA77','#A56FAF','#FFC000','#4557A5','#39B9EB')
col_main<-structure(col_vc,names=vc)
col_vc_lgd<-structure(col_vc,names=vc2)
col_vc_mut_load<-structure(col_vc,names=vc3)

lgd_sv = Legend(title = "Mutation Type", legend_gp = gpar(fill= col_vc_lgd), labels = vc2)

rgb2hex=function (a,b,c){
  return(rgb(a/255,b/255,c/255))
}

col_pli = colorRamp2(c(0,1), c("#FFDEDE", "#C80000")) 
col_prec = colorRamp2(c(0,0.9,1), c("#DEFFDE", "#00C8C8","#1fab89"))
col_misz = colorRamp2(c(min(data2$mis_z),0,max(data2$mis_z)), c( rgb2hex(94,124,107),rgb2hex(255,236,199),rgb2hex(208,118,2))) 
col_pvalue = colorRamp2(c(min(-log10(data2$p_value)),1,2,max(-log10(data2$p_value))), c( rgb2hex(230,240,250), rgb2hex(159,203,226), rgb2hex(32,115,183), rgb2hex(6,63,140)))

lgd_pli = Legend(title = "pLI", col_fun = col_pli, at = c(0,0.9, 1))
lgd_prec = Legend(title = "pRec", col_fun = col_prec, at = c(0,0.9, 1))
lgd_misz = Legend(title = "Mis-Z", col_fun = col_misz, at = c(trunc(min(data2$mis_z)),0, round(max(data2$mis_z),0)))
lgd_pvalue = Legend(title = "P-value", col_fun = col_pvalue, at = c(0, 1, 2, 3), labels = c("1", "0.1", "0.01", "0.001"))

lgd_am = Legend(title = "Amenorrhea", legend_gp = gpar(fill= c("#b7955c","#b1d39c")), labels = c("PA", "SA"))
lgd_moi = Legend(title = "MOI", legend_gp = gpar(fill= c("#96483B","#E07462",'#41609C')), labels = c("Biallelic", "Multi-het",'Monoallelic'))

gene_color<-c()
for (i in 1:length(data2$list)){
  if (data2$list[i] == "known"){
    gene_color[i]<-'black'
  }else{
    gene_color[i]<-'red'
  }
}

hm<-
  Heatmap(data2[,12:length(data2)],
        row_names_side = "left",
        show_column_names = FALSE,
        na_col = "#E9E9E9",
        col = col_main,
        
        row_order = rownames(arrange(type,desc(case_num),gene)),
        column_order = sample_order,
        
        row_split = factor( data2$group, levels = c('Gonadogenesis', 'Oogenesis', 'Folliculogenesis and ovulation','Others')),
        row_title_gp = gpar(fontsize = 12,fontface = 'bold'),
        row_names_gp = gpar(
          fontsize = 8,
          col = gene_color,
          fontface = 'italic'
          ),
        row_gap = unit(1, "mm"),

        show_heatmap_legend = FALSE,
        
        left_annotation = rowAnnotation(
          'Cases' = anno_barplot(
            data2$case_num,
            axis_param = list(direction = "reverse"), 
            # border = FALSE,
            gp = gpar(col='black',fill='MidnightBlue'),
            bar_width = 1
            ) 
          ),
        
        right_annotation = rowAnnotation(
          # annotation_name_side  = 'top',
          'pLI' = anno_simple(
            data2$pLI,
            col=col_pli
          ),

          'Mis-Z' = anno_simple(
            data2$mis_z,
            col=col_misz
          ),
          
          'P-value' = anno_simple(
            -log10(data2$p_value),
            col=col_pvalue
          )
          
        )
        ,

        top_annotation = columnAnnotation(
          # annotation_name_side  = 'left',
          'Mutation load' = anno_barplot( 
            case[,3:8],
            gp = gpar(col=col_vc_mut_load,fill=col_vc_mut_load),
            bar_width = 1
            ),
          
          'MOI' = anno_simple(
            case$MOI,
            col =c("Biallelic" = "#96483B", "Multi-het" = "#E07462", 'Monoallelic' ='#41609C')
          ),
          
          'Amenorrhea' = anno_simple(
            case$Amenorrhea,
            col =c("PA" = "#b7955c", "SA" = "#b1d39c")
          )
          
        )
)

pdf(file=paste0(dir,"/landscape.pdf"),width = 14, height=12)

pd = packLegend(lgd_sv,lgd_am,lgd_moi,lgd_pli,lgd_misz,lgd_pvalue,row_gap  = unit(1, "cm")) 
draw(hm,
     annotation_legend_list = pd,
     merge_legend = TRUE)

dev.off()


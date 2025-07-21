###usage: Rscript manqq.R <result.txt> <p-function:p_wald,p_lrt,p_score>
rm(list=ls())
args <- commandArgs(T)

library(qqman)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

manfile<-read.csv(file=args[1],sep = "\t", header=T, stringsAsFactors=F)
print("读取文件成功")
manfile$log10P=-log10(manfile[,args[2]])
pdf(file=paste(args[1],args[2],"_manhattan.pdf",sep=""), width=20, height=6)
manhattan(manfile, chr="chr",bp="ps",p="log10P",snp = "rs",main = "Manhattan Plot", cex = 0.6, cex.axis = 0.9,
          col = getPalette(15), chrlabs = NULL,
          suggestiveline = F,genomewideline = F,
          highlight = NULL, logp = F, annotatePval = NULL,
          annotateTop = TRUE)
dev.off()
pdf(file=paste(args[1],args[2],"_qq.pdf",sep=""), width=6, height=6)
qq(manfile[,args[2]], pch = 18, col = "blue4", cex = 1.5, las = 1)
dev.off()

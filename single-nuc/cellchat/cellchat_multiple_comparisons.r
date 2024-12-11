
suppressMessages({
    library(plyr)
    library(CellChat)
    library(patchwork)
    library(Seurat)
    library(SingleCellExperiment)
    library(ggalluvial)
    library(repr)
    library(ggplot2)
    library(gplots)
    library(RColorBrewer)
    library(gplots)
    library(tidyr)
    library(stringr)
    library(ComplexHeatmap)
 })
options(stringsAsFactors = FALSE)

'%!in%' <- function(x,y)!('%in%'(x,y))

path_data="/cellchat/output_files/"

obese=readRDS(paste0(path_data,"cellchat_obese_no_lymphatic_bcells_kit_scott_only_ct_fine_stringent_subsampled.rds"))
lean=readRDS(paste0(path_data,"cellchat_lean_no_lymphatic_bcells_kit_scott_only_ct_fine_stringent_subsampled.rds"))
wl=readRDS(paste0(path_data,"cellchat_weightloss_no_lymphatic_bcells_kit_scott_only_ct_fine_stringent_subsampled.rds"))

#weightloss does not have Mu5 cluster. Need to add it.
group.new = levels(obese@idents)
wl <- liftCellChat(wl, group.new)
lean <- liftCellChat(lean, group.new)

object.list <- list(Obese = obese, Lean = lean)
obese_lean <- mergeCellChat(object.list, add.names = names(object.list))
object.list <- list(Obese = obese, Weightloss = wl)
obese_wl <- mergeCellChat(object.list, add.names = names(object.list))
object.list <- list(Obese = obese, Lean = lean, Weightloss = wl)
cellchat_all <- mergeCellChat(object.list, add.names = names(object.list))

obese <- netAnalysis_computeCentrality(obese, slot.name = "netP")

lean <- netAnalysis_computeCentrality(lean, slot.name = "netP")

wl <- netAnalysis_computeCentrality(wl, slot.name = "netP")

pathways.show.all <- unique(c(obese@netP$pathways,wl@netP$pathways,lean@netP$pathways))
measure = c("outdeg","indeg","flowbet","info")
measure.name = c("Sender","Receiver","Mediator","Influencer")
mat2=NULL
for(j in pathways.show.all){
    centr=obese@netP$centr[j]
    for(i in 1:length(centr)) {
        centr0 <- centr[[i]]
        if(!is.null(centr0)){
        mat <- matrix(unlist(centr0), ncol = length(centr0), byrow = FALSE)
        mat <- t(mat)
        rownames(mat) <- names(centr0); colnames(mat) <- names(centr0$outdeg)
        if (!is.null(measure)) {
          mat <- mat[measure,]
          if (!is.null(measure.name)) {
            rownames(mat) <- measure.name
          }
        }
    mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
    mat=as.data.frame(mat)
    mat$Pathway=j
    mat$Role=rownames(mat)
        
    if (is.null(mat2)){
        mat2=data.frame(mat)
    }
    else{
        mat2=rbind(mat2,mat)
    }
    }    
    }
    }
mat2$Condition="Obese"
mat_obese=mat2
mat2=NULL
for(j in pathways.show.all){
    centr=lean@netP$centr[j]
    for(i in 1:length(centr)) {
        centr0 <- centr[[i]]
        if(!is.null(centr0)){
        mat <- matrix(unlist(centr0), ncol = length(centr0), byrow = FALSE)
        mat <- t(mat)
        rownames(mat) <- names(centr0); colnames(mat) <- names(centr0$outdeg)
        if (!is.null(measure)) {
          mat <- mat[measure,]
          if (!is.null(measure.name)) {
            rownames(mat) <- measure.name
          }
        }
    mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
    mat=as.data.frame(mat)
    mat$Pathway=j
    mat$Role=rownames(mat)
        
    if (is.null(mat2)){
        mat2=data.frame(mat)
    }
    else{
        mat2=rbind(mat2,mat)
    }
    }    
    }
    }
mat2$Condition="Lean"
mat_lean=mat2
mat2=NULL
for(j in pathways.show.all){
    centr=wl@netP$centr[j]
    for(i in 1:length(centr)) {
        centr0 <- centr[[i]]
        if(!is.null(centr0)){
        mat <- matrix(unlist(centr0), ncol = length(centr0), byrow = FALSE)
        mat <- t(mat)
        rownames(mat) <- names(centr0); colnames(mat) <- names(centr0$outdeg)
        if (!is.null(measure)) {
          mat <- mat[measure,]
          if (!is.null(measure.name)) {
            rownames(mat) <- measure.name
          }
        }
    mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
    mat=as.data.frame(mat)
    mat$Pathway=j
    mat$Role=rownames(mat)
        
    if (is.null(mat2)){
        mat2=data.frame(mat)
    }
    else{
        mat2=rbind(mat2,mat)
    }
    }    
    }
    }
mat2$Condition="Weightloss"
mat_wl=mat2
mat_final=rbind(mat_obese,mat_lean,mat_wl)

head(mat_final)

write.csv(mat_final, paste0(path_data,"cellchat_ct_fine_no_lymphatic_bcells_kit_stringent_subsampled_scott_only_centrality_all.csv"), row.names=FALSE)

df.net <- subsetCommunication(obese)

unique(df.net$pathway_name)

write.csv(df.net, paste0(path_data,"cellchat_ct_fine_no_lymphatic_bcells_kit_stringent_subsampled_scott_only_net_obese.csv"), row.names=FALSE)

df.net <- subsetCommunication(lean)

write.csv(df.net, paste0(path_data,"cellchat_ct_fine_no_lymphatic_bcells_kit_stringent_subsampled_scott_only_net_lean.csv"), row.names=FALSE)

df.net <- subsetCommunication(wl)

write.csv(df.net, paste0(path_data,"cellchat_ct_fine_no_lymphatic_bcells_kit_stringent_subsampled_scott_only_net_wl.csv"), row.names=FALSE)

options(repr.plot.width = 6, repr.plot.height = 6)

heatmap.df <- list(contrib = data.frame(matrix(ncol = 0, nrow = 0)), pval = data.frame(matrix(ncol = 0, nrow = 0)))
cellchat=cellchat_all
genotypes <- names(cellchat@net)
for (i in c(2:3)){
    gg <- rankNet(cellchat, mode = "comparison", comparison=c(1,i)
            , title= genotypes[i]
            , cutoff.pvalue = 0.05
            , stacked = T, show.raw = T, do.stat = TRUE, return.data = FALSE)
    temp <- data.frame(gg$data[gg$data$group == levels(gg$data$group)[1],c("contribution", "pvalues")], 
                       control.contribution=gg$data[gg$data$group == levels(gg$data$group)[2],c("contribution")], 
                       row.names = rownames(gg$data)[gg$data$group == levels(gg$data$group)[2]])
    temp2 <- data.frame((temp[,1]/temp[,3]), row.names = rownames(temp))
    colnames(temp2) <- genotypes[i]
    if(sum(dim(heatmap.df$contrib) == c(0,0))){
        heatmap.df$contrib <- temp2
    } else{
        heatmap.df$contrib <- merge(heatmap.df$contrib, temp2, all = TRUE, by = "row.names")
        rownames(heatmap.df$contrib) <- heatmap.df$contrib$Row.names
        heatmap.df$contrib <- heatmap.df$contrib[,c(-1)]
    }
    
    temp2 <- data.frame(p.adjust(temp[,2], method="bonferroni"), row.names = rownames(temp))
    colnames(temp2) <- genotypes[i]
    if(sum(dim(heatmap.df$pval) == c(0,0))){
        heatmap.df$pval <- temp2
    } else {
        heatmap.df$pval <- merge(heatmap.df$pval, temp2, all = TRUE, by = "row.names")
        rownames(heatmap.df$pval) <- heatmap.df$pval$Row.names
        heatmap.df$pval <- heatmap.df$pval[,c(-1)]
    }
}



heatmap.df$contrib <- log2(heatmap.df$contrib)

heatmap.df$pval <- data.frame(sapply(heatmap.df$pval, as.numeric), row.names=rownames(heatmap.df$pval))
heatmap.df$contrib <- data.frame(sapply(heatmap.df$contrib, as.numeric), row.names=rownames(heatmap.df$contrib))

temp <- heatmap.df$contrib
temp$Pathway <- rownames(temp)
heatmap.df.contrib.wide <- gather(temp, Group, log2FC_contrib, Lean:Weightloss, factor_key=TRUE)
print(dim(heatmap.df.contrib.wide))
head(heatmap.df.contrib.wide)

temp <- heatmap.df$pval
temp$Pathway <- rownames(temp)
heatmap.df.pval.wide <- gather(temp, Group, adjusted_pval,Lean:Weightloss, factor_key=TRUE)
print(dim(heatmap.df.pval.wide))
head(heatmap.df.pval.wide)

heatmap.df.wide <- merge(heatmap.df.contrib.wide, heatmap.df.pval.wide, by=c("Pathway","Group"), all.x = TRUE, all.y = TRUE)

head(heatmap.df.wide)

heatmap.df.wide[which(heatmap.df.wide$Pathway=="NRG"),]

MAX.CONTRIB <- max(abs(heatmap.df$contrib[!(sapply(heatmap.df$contrib, is.infinite) | is.na(heatmap.df$contrib))]))
MAX.CONTRIB <- round(MAX.CONTRIB)
MAX.CONTRIB

INF.VALUE <- MAX.CONTRIB+1
INF.VALUE

heatmap.df.wide$adjusted_pval[heatmap.df.wide$	log2FC_contrib<0.2]=1

heatmap.df.wide$adjusted_pval[heatmap.df.wide$	log2FC_contrib<0.2]=1

options(repr.plot.width = 20, repr.plot.height = 5)

heatmap.df2 <- heatmap.df

heatmap.df2$contrib[sapply(heatmap.df2$contrib, is.infinite) & heatmap.df2$contrib > 0] <- INF.VALUE
heatmap.df2$contrib[sapply(heatmap.df2$contrib, is.infinite) & heatmap.df2$contrib < 0] <- -INF.VALUE
heatmap.df2$pval[(abs(heatmap.df2$contrib) < 0.2)] <- 1
heatmap.df2$pval[is.na(heatmap.df2$pval)] <- 2


heatmap.df2$pval <- heatmap.df2$pval[(apply(is.na(heatmap.df2$contrib), 1, sum) < 2),]
heatmap.df2$contrib <- heatmap.df2$contrib[(apply(is.na(heatmap.df2$contrib), 1, sum) < 2),]

heatmap.df2$contrib <- heatmap.df2$contrib[(apply(heatmap.df2$pval > 0.01, 1, sum) < 2),]
heatmap.df2$pval <- heatmap.df2$pval[(apply(heatmap.df2$pval > 0.01, 1, sum) < 2),]

heatmap.df2$pval[heatmap.df2$pval == 2] <- NA

heatmap.df2$sig <- heatmap.df2$pval
heatmap.df2$sig[] <- ''
heatmap.df2$sig[heatmap.df2$pval <= 0.01] <- '*'

heatmap.df2$sig[is.na(heatmap.df2$pval)] <- '-'

heatmap.df2$contrib[is.na(heatmap.df2$contrib)] <- 0

heatmap.df2$sig <- heatmap.df2$sig[,c("Lean","Weightloss")]
heatmap.df2$pval <- heatmap.df2$pval[,c("Lean","Weightloss")]
heatmap.df2$contrib <- heatmap.df2$contrib[,c("Lean","Weightloss")]

scale.interval.size <- 0.05
num_breaks <- (MAX.CONTRIB * 2 / scale.interval.size)+1
breaks <- c(-(INF.VALUE+1), seq(from=-MAX.CONTRIB, to=MAX.CONTRIB, length.out=num_breaks), INF.VALUE+1)
midpoint <- 0 # the mid of the "real" scale

rampCol2 <- colorRampPalette(c("#6699FF", "white", "#FF6600"))(length(breaks)-1)
mypalette <- c(rampCol2)
mypalette[1] <- "#5689EF" # just a random extreme color for -inf
mypalette[num_breaks+1] <- "#EF5600" # just a random extreme color for inf


heatmap.df3 <- heatmap.df2
heatmap.df4 <- heatmap.df2

heatmap.df3$contrib <- heatmap.df2$contrib[(apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) == 0)
                                           |(apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) == 0) ,]
heatmap.df3$sig <- heatmap.df2$sig[(apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) == 0)
                                           |(apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) == 0) ,]
heatmap.df3$sig <- heatmap.df3$sig[order(apply(heatmap.df3$contrib, 1, mean)),]
heatmap.df3$contrib <- heatmap.df3$contrib[order(apply(heatmap.df3$contrib, 1, mean)),]

heatmap.df4$contrib <- heatmap.df2$contrib[!((apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) == 0)
                                           |(apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) == 0)) ,]
heatmap.df4$sig <- heatmap.df2$sig[!((apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) == 0)
                                           |(apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) == 0)) ,]

heatmap.df4$contrib <- heatmap.df4$contrib[with(data.frame(heatmap.df4$sig == '*'), order(-Lean,-Weightloss)),]
heatmap.df4$sig <- heatmap.df4$sig[with(data.frame(heatmap.df4$sig == '*'), order(-Lean,-Weightloss)),]

heatmap.df4$change <- heatmap.df4$contrib
heatmap.df4$change[] <- 0
heatmap.df4$change[heatmap.df4$contrib > 0 & heatmap.df4$sig != ''] <- 1
heatmap.df4$change[heatmap.df4$contrib < 0 & heatmap.df4$sig != ''] <- -1
heatmap.df4$sig <- heatmap.df4$sig[order(apply(heatmap.df4$change, 1, sum), decreasing = TRUE),]
heatmap.df4$contrib <- heatmap.df4$contrib[order(apply(heatmap.df4$change, 1, sum), decreasing = TRUE),]

heatmap.df4$contrib <- heatmap.df4$contrib[order(apply(heatmap.df4$sig == '*', 1, sum)),]
heatmap.df4$sig <- heatmap.df4$sig[order(apply(heatmap.df4$sig == '*', 1, sum)),]

heatmap.df3$contrib=(heatmap.df3$contrib)*(-1)
heatmap.df4$contrib=(heatmap.df4$contrib)*(-1)

heatmap.2(t(heatmap.df3$contrib)
        , col=mypalette
        , Rowv=NULL  
        , Colv=NULL
        , dendrogram="none"
        , na.rm = TRUE
        , breaks=breaks, density.info="none", trace="none"
        , symm=F,symkey=F,symbreaks=T, scale="none"
        , margins = c(14,18)
        , lhei=c(2,4), lwid=c(1,6)
        , cellnote = t(heatmap.df3$sig)
        ,notecex=3.0
        ,notecol="black"
        ,cexRow=3
        ,cexCol=2
       )

heatmap.2(t(heatmap.df4$contrib)
        , col=mypalette
        , Rowv=NULL  
        , Colv=NULL
        , dendrogram="none"
        , na.rm = TRUE
        , breaks=breaks, density.info="none", trace="none"
        , symm=F,symkey=F,symbreaks=T, scale="none"
        , margins = c(14,18)
        , lhei=c(2,4), lwid=c(1,6)
        , cellnote = t(heatmap.df4$sig)
        ,notecex=3.0
        ,notecol="black"
        ,cexRow=3
        ,cexCol=2
       )

options(repr.plot.width = 6, repr.plot.height = 6)
control="Obese"
cellchat=cellchat_all
celltypes=levels(cellchat@idents[[control]])
genotypes <- names(cellchat@net)
heatmap.df_final <- list(contrib = data.frame(matrix(ncol = 0, nrow = 0)), pval = data.frame(matrix(ncol = 0, nrow = 0)))
used_pairs=c()
for (i in c(2:3)){
    heatmap.df <- list(contrib = data.frame(matrix(ncol = 0, nrow = 0)), pval = data.frame(matrix(ncol = 0, nrow = 0)))
    for (j in celltypes){tryCatch({for (m in celltypes){
    pair=paste0(j,"|",m)
    if (pair %in% used_pairs){}
    else{
    gg <- rankNet(cellchat, mode = "comparison", comparison=c(1,i)
            , title= genotypes[i]
            , cutoff.pvalue = 1
            , stacked = T, show.raw = T, do.stat = TRUE, return.data = FALSE,sources.use=c(j),targets.use=c(m),thresh=1)
    temp <- data.frame(gg$data[gg$data$group == levels(gg$data$group)[1],c("contribution", "pvalues")], control.contribution=gg$data[gg$data$group == levels(gg$data$group)[2],c("contribution")], row.names = paste0(rownames(gg$data)[gg$data$group == levels(gg$data$group)[2]],"__",j,"|",m))
    temp2 <- data.frame((temp[,1]/temp[,3]), row.names = rownames(temp))
    colnames(temp2) <- genotypes[i]
        
    if(sum(dim(heatmap.df$contrib) == c(0,0))){
        heatmap.df$contrib <- temp2
    } else{
        heatmap.df$contrib <- rbind(heatmap.df$contrib, temp2)

    }
    
    temp2 <- data.frame(p.adjust(temp[,2], method="bonferroni"), row.names = rownames(temp))
    colnames(temp2) <- genotypes[i]
    if(sum(dim(heatmap.df$pval) == c(0,0))){
        heatmap.df$pval <- temp2
    } else {
        heatmap.df$pval <- rbind(heatmap.df$pval, temp2)

    }

 }}}, error=function(e){})}
 if(sum(dim(heatmap.df_final$pval) == c(0,0))){
        heatmap.df_final$pval <- heatmap.df$pval 
    } else {
        heatmap.df_final$pval <- merge(heatmap.df_final$pval, heatmap.df$pval , all = TRUE, by = 0)
        rownames(heatmap.df_final$pval) <- heatmap.df_final$pval$Row.names
        heatmap.df_final$pval <- heatmap.df_final$pval[,c(-1)]
     }
 if(sum(dim(heatmap.df_final$contrib) == c(0,0))){
        heatmap.df_final$contrib <- heatmap.df$contrib 
    } else {
        heatmap.df_final$contrib <- merge(heatmap.df_final$contrib, heatmap.df$contrib , all = TRUE, by = 0)
        rownames(heatmap.df_final$contrib) <- heatmap.df_final$contrib$Row.names
        heatmap.df_final$contrib <- heatmap.df_final$contrib[,c(-1)]
     }
    }

heatmap.df_final$contrib <- log2(heatmap.df_final$contrib)

heatmap.df_final$pval <- data.frame(sapply(heatmap.df_final$pval, as.numeric), row.names=rownames(heatmap.df_final$pval))
heatmap.df_final$contrib <- data.frame(sapply(heatmap.df_final$contrib, as.numeric), row.names=rownames(heatmap.df_final$contrib))

heatmap.df_final$contrib[c('Pathway', 'Pair')] <- str_split_fixed(rownames(heatmap.df_final$contrib), '__', 2)

heatmap.df_final$pval[c('Pathway', 'Pair')] <- str_split_fixed(rownames(heatmap.df_final$pval), '__', 2)

heatmap.df<-heatmap.df_final
temp <- heatmap.df$contrib
temp$Pathway <- rownames(temp)
heatmap.df.contrib.wide <- gather(temp, Group, log2FC_contrib, Lean:Weightloss, factor_key=TRUE)
print(dim(heatmap.df.contrib.wide))
head(heatmap.df.contrib.wide)

temp <- heatmap.df$pval
temp$Pathway <- rownames(temp)
heatmap.df.pval.wide <- gather(temp, Group, adjusted_pval,Lean:Weightloss, factor_key=TRUE)
print(dim(heatmap.df.pval.wide))
head(heatmap.df.pval.wide)

heatmap.df.wide <- merge(heatmap.df.contrib.wide, heatmap.df.pval.wide, by=c("Pathway","Group"), all.x = TRUE, all.y = TRUE)
heatmap.df.wide  <- heatmap.df.wide [,!names(heatmap.df.wide) %in% c("Pair.x", "Pair.y")]

head(heatmap.df.wide)

heatmap.simple=heatmap.df

heatmap.simple$contrib=heatmap.simple$contrib[, -((ncol(heatmap.simple$contrib) - 1):ncol(heatmap.simple$contrib))]
heatmap.simple$pval=heatmap.simple$pval[, -((ncol(heatmap.simple$pval) - 1):ncol(heatmap.simple$pval))]

MAX.CONTRIB <- max(abs(heatmap.simple$contrib[!(sapply(heatmap.simple$contrib, is.infinite) | is.na(heatmap.simple$contrib))]))
MAX.CONTRIB <- round(MAX.CONTRIB)
MAX.CONTRIB

INF.VALUE <- MAX.CONTRIB+1
INF.VALUE

options(repr.plot.width = 60, repr.plot.height = 10)

heatmap.df2 <- heatmap.simple

heatmap.df2$contrib[sapply(heatmap.df2$contrib, is.infinite) & heatmap.df2$contrib > 0] <- INF.VALUE
heatmap.df2$contrib[sapply(heatmap.df2$contrib, is.infinite) & heatmap.df2$contrib < 0] <- -INF.VALUE
heatmap.final.corr=heatmap.df2

heatmap.df2$pval[is.na(heatmap.df2$pval)] <- 2


heatmap.df2$pval <- heatmap.df2$pval[(apply(is.na(heatmap.df2$contrib), 1, sum) < 2),]
heatmap.df2$contrib <- heatmap.df2$contrib[(apply(is.na(heatmap.df2$contrib), 1, sum) < 2),]

heatmap.df2$contrib <- heatmap.df2$contrib[(apply(heatmap.df2$pval > 0.01, 1, sum) < 2),]
heatmap.df2$pval <- heatmap.df2$pval[(apply(heatmap.df2$pval > 0.01, 1, sum) < 2),]

heatmap.df2$pval[heatmap.df2$pval == 2] <- NA

heatmap.df2$sig <- heatmap.df2$pval
heatmap.df2$sig[] <- ''
heatmap.df2$sig[heatmap.df2$pval <= 0.01] <- '*'

heatmap.df2$sig[is.na(heatmap.df2$pval)] <- '-'

heatmap.df2$contrib[is.na(heatmap.df2$contrib)] <- 0

heatmap.df2$sig <- heatmap.df2$sig[,c("Lean","Weightloss")]
heatmap.df2$pval <- heatmap.df2$pval[,c("Lean","Weightloss")]
heatmap.df2$contrib <- heatmap.df2$contrib[,c("Lean","Weightloss")]

scale.interval.size <- 0.05
num_breaks <- (MAX.CONTRIB * 2 / scale.interval.size)+1
breaks <- c(-(INF.VALUE+1), seq(from=-MAX.CONTRIB, to=MAX.CONTRIB, length.out=num_breaks), INF.VALUE+1)
midpoint <- 0 # the mid of the "real" scale

rampCol2 <- colorRampPalette(c("#6699FF", "white", "#FF6600"))(length(breaks)-1)
mypalette <- c(rampCol2)
mypalette[1] <- "#5689EF" # just a random extreme color for -inf
mypalette[num_breaks+1] <- "#EF5600" # just a random extreme color for inf


heatmap.df3 <- heatmap.df2
heatmap.df4 <- heatmap.df2

heatmap.df3$contrib <- heatmap.df2$contrib[(apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) == 0)
                                           |(apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) == 0) ,]
heatmap.df3$sig <- heatmap.df2$sig[(apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) == 0)
                                           |(apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) == 0) ,]
heatmap.df3$sig <- heatmap.df3$sig[order(apply(heatmap.df3$contrib, 1, mean)),]
heatmap.df3$contrib <- heatmap.df3$contrib[order(apply(heatmap.df3$contrib, 1, mean)),]

heatmap.df4$contrib <- heatmap.df2$contrib[!((apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) == 0)
                                           |(apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) == 0)) ,]
heatmap.df4$sig <- heatmap.df2$sig[!((apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) == 0)
                                           |(apply(heatmap.df2$contrib < 0 & heatmap.df2$sig != '', 1, sum) >=2 & apply(heatmap.df2$contrib > 0 & heatmap.df2$sig != '', 1, sum) == 0)) ,]

heatmap.df4$contrib <- heatmap.df4$contrib[with(data.frame(heatmap.df4$sig == '*'), order(-Lean,-Weightloss)),]
heatmap.df4$sig <- heatmap.df4$sig[with(data.frame(heatmap.df4$sig == '*'), order(-Lean,-Weightloss)),]

heatmap.df4$change <- heatmap.df4$contrib
heatmap.df4$change[] <- 0
heatmap.df4$change[heatmap.df4$contrib > 0 & heatmap.df4$sig != ''] <- 1
heatmap.df4$change[heatmap.df4$contrib < 0 & heatmap.df4$sig != ''] <- -1
heatmap.df4$sig <- heatmap.df4$sig[order(apply(heatmap.df4$change, 1, sum), decreasing = TRUE),]
heatmap.df4$contrib <- heatmap.df4$contrib[order(apply(heatmap.df4$change, 1, sum), decreasing = TRUE),]

heatmap.df4$contrib <- heatmap.df4$contrib[order(apply(heatmap.df4$sig == '*', 1, sum)),]
heatmap.df4$sig <- heatmap.df4$sig[order(apply(heatmap.df4$sig == '*', 1, sum)),]

heatmap.2(t(heatmap.df3$contrib)
        , col=mypalette
        , Rowv=NULL  
        , Colv=NULL
        , dendrogram="none"
        , na.rm = TRUE
        , breaks=breaks, density.info="none", trace="none"
        , symm=F,symkey=F,symbreaks=T, scale="none"
        , margins = c(33,35)
        , lhei=c(1,5), lwid=c(1,7)
        , cellnote = t(heatmap.df3$sig)
        ,notecex=3.0
        ,notecol="black"
        ,cexRow=3
        ,cexCol=1.5
       )

heatmap.2(t(heatmap.df4$contrib)
        , col=mypalette
        , Rowv=NULL  
        , Colv=NULL
        , dendrogram="none"
        , na.rm = TRUE
        , breaks=breaks, density.info="none", trace="none"
        , symm=F,symkey=F,symbreaks=T, scale="none"
        , margins = c(33,35)
        , lhei=c(1,5), lwid=c(1,7)
        , cellnote = t(heatmap.df4$sig)
        ,notecex=3.0
        ,notecol="black"
        ,cexRow=3
        ,cexCol=1.5
       )

lean=heatmap.final.corr$pval[,"Lean",drop=FALSE]
colnames(lean)=c("pval")
lean$contrib=heatmap.final.corr$contrib$Lean
lean$Condition="Lean"
lean$Group=rownames(lean)
rownames(lean)=NULL

head(lean)

wl=heatmap.final.corr$pval[,"Weightloss",drop=FALSE]
colnames(wl)=c("pval")
wl$contrib=heatmap.final.corr$contrib$Weightloss
wl$Condition="Weightloss"
wl$Group=rownames(wl)
rownames(wl)=NULL

head(wl)

merged.df=rbind(lean,wl)

merged.df[c('Pathway', 'Pair')] <- str_split_fixed(merged.df$Group, '__', 2)

merged.df[c('Sender', 'Receiver')] <- str_split_fixed(merged.df$Pair, '\\|', 2)

rownames(merged.df)=NULL

head(merged.df)

write.csv(merged.df,paste0(path_data,"rankNet_vs_obese_ct_fine_pairwise_scott_only_subsampled_stringent_no_lymphatic_bcells_kit.csv")) # use this data for Sankey plots



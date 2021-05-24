options(stringsAsFactors = F);library(WGCNA);enableWGCNAThreads(12)
data=read.csv("GSE43765.csv",header=T,row.names = 1)
data=as.matrix(apply(data,2,rank, ties.method= "max"))
datSummary=rownames(data)
no.samples = dim(data)[[1]];
library(preprocessCore)
datExpr=t(data)
GeneName= datSummary
ArrayName= colnames(data)
powers=c(seq(1,10,by=1),seq(12,18,by=2));
sft=pickSoftThreshold(datExpr, powerVector=powers,networkType ="signed",corFnc =cor, corOptions =list(use = 'p'))
RpowerTable=sft[[2]]
sizeGrWindow(9, 5);
pdf('choosing power.pdf');
par(mfrow = c(1,2));cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red");
dev.off()
sizeGrWindow(9, 5);
pdf('mean connectivity.pdf');
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red");
dev.off()
softPower =12
Connectivity=softConnectivity(datExpr,corFnc = "cor", corOptions ="use = 'p'",power=softPower,type="signed")
pdf("scale-free.pdf");
scaleFreePlot(Connectivity,nBreaks = 10,truncated = FALSE,removeFirst = FALSE, main = "");
dev.off()
adjacency = adjacency(datExpr,corFnc = "cor", corOptions ="use = 'p'",
                      type = "signed", power = softPower)
TOM = TOMsimilarity(adjacency,TOMType="signed");dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize =30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 0, 
                            pamRespectsDendro = FALSE,minClusterSize = minModuleSize,
                            cutHeight=0.99);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

MEList = moduleEigengenes(datExpr, colors = dynamicMods)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
MEDissThres = 0.2
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicMods, cutHeight = MEDissThres, verbose = 3);
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
pdf("DendroAndColors.pdf")
plotDendroAndColors(geneTree, cbind(dynamicMods, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, 
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleColors = mergedColors
colorOrder = c("grey", standardColors(unique(moduleColors)));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");#
pdf("METree.pdf")
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()
MEList = moduleEigengenes(datExpr, colors = dynamicMods)
nSamples=nrow(datExpr)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = cbind.data.frame(datSummary,corPvalueStudent(as.matrix(geneModuleMembership), 
                                                        nSamples));
write.table(data.frame(ArrayName,MEs),"MEs.csv",row.name=F)
kMEdat=data.frame(geneModuleMembership,MMPvalue)
write.table(data.frame(datSummary,kMEdat),"kME-MMPvalue.csv",row.names=FALSE)
k.in=intramodularConnectivity(adjacency(datExpr,corFnc = "cor", corOptions = "use ='p'",type = "signed", power = softPower),moduleColors,scaleByMax = FALSE)
datout=data.frame(datSummary, colorNEW=moduleColors, k.in)
write.table(datout, file="OutputCancerNetwork.csv", sep=",", row.names=F)
hubs    = chooseTopHubInEachModule(datExpr, moduleColors)
write.csv(data.frame(module=names(hubs),moduleColor=labels2colors(names(hubs)),hub=hubs),
          "num2color.csv",row.names=F)

gene=read.csv("OutputCancerNetwork.csv",header=T)
library(gProfileR)
for (i in unique(gene$colorNEW)[-3]){
  genes=subset(gene$datSummary,gene$colorNEW==i)
  go=gprofiler(as.vector(genes), 
               organism = "hsapiens",numeric_ns="ENTREZGENE_ACC")[,-14]
  write.table(go,"module_enrichment.csv",append =T,row.names=rep(i,nrow(go)),sep=",")}

moduleColors=gene$colorNEW #keep gene order same between gene and dat0
modules=unique(moduleColors)
n=length(modules)
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (p in (1:n)[-3]){
  inModule = is.finite(match(moduleColors, modules[p]));
  dat2=data.frame(t(datExpr))[inModule,] #make sure it is matrix data
  resamples=lapply(1:1000,function(i) a=t(sample(dat2[,1:nSamples],round(nSamples/2),replace=F)))
  K1=sapply(resamples,softConnectivity,power= softPower,type="signed") #,type="signed"?
  K=softConnectivity(t(dat2[,1:nSamples]),power= softPower,type="signed") #,type="signed"
  #outfile=paste(modules[p],"-edit.txt",sep="")
  write.table(data.frame(mean(cor(K,K1)),apply(cor(K,K1),1,sd)), file = "module-stability.csv", row.names = modules[p], append = TRUE, col.names = FALSE, sep = ", ")
  setTxtProgressBar(pb, p)}
close(pb) 

##############module preservation 
setLabels = c("Female", "Male");
datSummaryFemale=rownames(data)
datSummaryMale=rownames(data)
datExprFemale= datExpr
no.samplesFemale <- dim(datExprFemale)[[1]]
dim(datExprFemale)
datExprMale= t(read.csv("GSE43765.csv",header=T,row.names = 1)) 
colorsFemale = gene$colorNEW
colnames(datExprMale)=colnames(datExpr)
colnames(datExprFemale)=colnames(datExpr)
nSets = 2
ref = 1
test = 2
multiExpr = list(Female = list(data = datExprFemale), Male = list(data = datExprMale));
multiColor = list(Female = labels2colors(colorsFemale));
mp = modulePreservation(multiExpr, multiColor,referenceNetworks = 1,networkType="signed",nPermutations = 100,randomSeed = 1,parallelCalculation=T,quickCor = 0,verbose = 3)
save(mp, file = "modulePreservation.RData");
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];


# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
#text = modColors[plotMods];
labs = match(modColors[plotMods], standardColors(unique(modColors)-2))
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
pdf(file="FemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5,onefile=TRUE)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labs = match(modColors[plotMods], standardColors(length(unique(modColors))))
  write.table(data.frame(mod))
  #replace text to labs as number labeling: labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08); 
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], labs, cex = 1, offs = 0.08)
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
dev.off();


# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(length(unique(modColors))));  #50 should larger than module number
# Start the plot: open a suitably sized graphical window and set sectioning and margins. Alternatively,
# plot into a pdf file.
sizeGrWindow(10, 9);
pdf(file="PreservationZStatistics.pdf", w=10, h=9)
par(mfrow = c(4,4))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
for (s in 1:ncol(statsZ))
{
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  if (min > -max/12) min = -max/12
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 2.2,
       ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(30, 1200),
       cex.lab = 1.2, cex.axis = 1.2)
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 1, offs = 0.06);
  #text(moduleSizes[-1], statsZ[-c(1:2), s], labels = letter[-c(1:2)], col = "black"); #modColors[-2]);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}
# If plotting into a file, close it, otherwise it is unreadable.
dev.off();

#####other dataset projection,note match genes and moduleColors
gse9890=read.csv("GSE9890.csv",header=T,row.names = 1)
commongenes=intersect(rownames(data),rownames(gse9890)) #rownames(data) as reference
gse9890=as.matrix(apply(gse9890[commongenes,],2,rank, ties.method= "max"))
ME_gse9890 = moduleEigengenes(t(gse9890), colors = moduleColors[which(rownames(data) %in% commongenes)])
ME_gse9890=ME_gse9890$eigengenes
ArrayName_gse9890=row.names(ME_gse9890)
write.table(data.frame(ArrayName_gse9890,ME_gse9890),"ME_gse9890.csv",row.name=F)

gse50831=read.csv("GSE50831.csv",header=T,row.names = 1)
commongenes=intersect(rownames(data),rownames(gse50831)) #rownames(data) as reference
gse50831=as.matrix(apply(gse50831[commongenes,],2,rank, ties.method= "max"))
ME_gse50831 = moduleEigengenes(t(gse50831), colors = moduleColors[which(rownames(data) %in% commongenes)])
ME_gse50831=ME_gse50831$eigengenes
ArrayName_gse50831=row.names(ME_gse50831)
write.table(data.frame(ArrayName_gse50831,ME_gse50831),"ME_gse50831.csv",row.name=F)

gse26712=read.csv("GSE26712.csv",header=T,row.names = 1)
commongenes=intersect(rownames(data),rownames(gse26712)) #rownames(data) as reference
gse26712=as.matrix(apply(gse26712[commongenes,],2,rank, ties.method= "max"))
ME_gse26712 = moduleEigengenes(t(gse26712), colors = moduleColors[which(rownames(data) %in% commongenes)])
ME_gse26712=ME_gse26712$eigengenes
ArrayName_gse26712=row.names(ME_gse26712)
write.table(data.frame(ArrayName_gse26712,ME_gse26712),"ME_gse26712.csv",row.name=F)

tcga=read.csv("tcga.csv",header=T,row.names = 1)
commongenes=intersect(rownames(data),rownames(tcga)) #rownames(data) as reference
tcga=as.matrix(apply(tcga[commongenes,],2,rank, ties.method= "max"))
ME_tcga = moduleEigengenes(t(tcga), colors = moduleColors[which(rownames(data) %in% commongenes)])
ME_tcga=ME_tcga$eigengenes
ArrayName_tcga=row.names(ME_tcga)
write.table(data.frame(ArrayName_tcga,ME_tcga),"ME_tcga.csv",row.name=F)

#######survival analysis
library(tidyverse);library(tidytidbits);library(survivalAnalysis)
df <- data
map(vars(therapy_outcome, residual, neoplasm_status, indication,age,drug,stage,anatomic,
         ME1,ME2,ME3,ME4,ME5,ME6,ME7,ME8,ME9
         ,ME10,ME11,ME12,ME13,ME14,ME16,ME17,ME18,ME19,ME20,ME21,ME22,ME23,ME24), function(by)
         {
           analyse_multivariate(df,
                                vars(days_to_death, death_status),
                                covariates = list(by), # covariates expects a list
                                covariate_name_dict = covariate_names)
         }) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(os_time="OS"),
              orderer = ~order(HR),
              labels_displayed = c("endpoint", "factor", "n"),
              ggtheme = ggplot2::theme_bw(base_size = 10))

######Multivariate Cox regression analysis
library(survival);library(survminer)
res.cox <- coxph(Surv(days_to_death, death_status) ~therapy_outcome+age+ME19, data =  data)
summary(res.cox)

#####Fit survival data using the Kaplan-Meier method
surv_object <- Surv(time = data$days_to_death, event = data$death_status)
data[,11:33]=apply(data[,11:33]>0,2,as.numeric)
fit2 <- survfit(surv_object ~ data$ME19, data = data)
ggsurvplot(fit2, data =data, pval = TRUE,surv.plot.height=0.25)

#####ploting Fig. 6
pldata=read.csv("4R_gse50831.csv",header = T,row.names=1)
par(cex.lab=1.3,cex.axis=1.3)
col=standardColors(19)[c(3,5,19,2,6,7,8,11,12)]
boxplot(pldata$ME3~pldata$treatment,data=pldata,outline=F,
        ylab="Module eigengene (ME)",height=1.5,col=col[1])
boxplot(pldata$ME5~pldata$treatment,data=pldata,outline=F,
        ylab="Module eigengene (ME)",height=1.5,col=col[2])
boxplot(pldata$ME19~pldata$treatment,data=pldata,outline=F,
        ylab="Module eigengene (ME)",height=1.5,col=col[3])
boxplot(pldata$ME2~pldata$treatment,data=pldata,outline=F,
        ylab="Module eigengene (ME)",height=1.5,col=col[4])
boxplot(pldata$ME6~pldata$treatment,data=pldata,outline=F,
        ylab="Module eigengene (ME)",height=1.5,col=col[5])
boxplot(pldata$ME7~pldata$treatment,data=pldata,outline=F,
        ylab="Module eigengene (ME)",height=1.5,col="darkgrey")
boxplot(pldata$ME8~pldata$treatment,data=pldata,outline=F,
        ylab="Module eigengene (ME)",height=1.5,col=col[7])
boxplot(pldata$ME11~pldata$treatment,data=pldata,outline=F,
        ylab="Module eigengene (ME)",height=1.5,col=col[8])
boxplot(pldata$ME12~pldata$treatment,data=pldata,outline=F,
        ylab="Module eigengene (ME)",height=1.5,col=col[9])
boxplot(pldata$ME19[pldata$treatment!="ERI"]~pldata$treatment[pldata$treatment!="ERI"],data=pldata,outline=F)

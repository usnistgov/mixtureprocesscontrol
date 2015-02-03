#Required packages:  Limma,data.table,ggplot2,reshape2
require(limma)||install.packages(limma)
require(data.table)||install.packages(data.table)
require(ggplot2)||install.packages(ggplot2)
require(reshape2)||install.packages(reshape2)
require(DESeq2)||install.packages(DESeq2)
require(grid)||install.packages(grid)
require(dplyr)||install.packages(dplyr)
#The (ACTUAL) Figure code used in the manuscript for blm as of 14/09/29.  No changes to the figures will be made other than through this code!
#Except figure 1, because figure 1 is from a powerpoint.  Fig1BW.ppt

##This code all runs off of count tables - these count tables should also be added to git (and have their locations in the code suitably refactored) -
#current count table locations:
# SEQC : #ILM data as
# BLM :  #added as BLMgene.zip

#Figure2:
theme_set(theme_bw(base_size=18))
#a<-makeTargetPlot(type="BLM",df=makerefmetdfb(norm="UQN"))
#a+coord_cartesian(ylim=c(0,0.75),xlim=c(0,0.75))+scale_x_continuous(breaks=c(0,0.25,0.5,0.75))+scale_y_continuous(breaks=c(0,0.25,0.5,0.75))
#Figure3: (just call dendrogramplot)
dendrogramplot<-function(testing=FALSE){
  require(stats)
  require(RColorBrewer)
  require(DESeq2)
  SEQCDF<-makeseqcdf()
  vsndata<-varstabdata()
  test<-dist(t(assay(vsndata)))
  test<-as.matrix(test)
  rownames(test)<-colnames(test)<-with(colData(vsndata),paste(sample,site,sep=":"))
  hcr <- suppressWarnings(hclust(dist(test)))
  hcr$labels<-paste0(c(rep(c("A","B","C","D","Cm","Dm"),6)),c(rep(1,6),rep(2,6),rep(3,6),rep(4,6),rep(5,6),rep(6,6)))
  hcrd<-as.dendrogram(hcr)
  #labelColors = c("#CDB380", "#036564", "#EB6841", "#EDC951")
  labelColors<<-brewer.pal(n = 4,name = "PuOr")
  # cut dendrogram in 4 clusters
  clusMember<<- cutree(hcr, 4)
  # using dendrapply
  colLab<-function(n) {
    if (is.leaf(n)) {
      a <- attributes(n)
      labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
      attr(n, "nodePar") <- list(a$nodePar, lab.col = labCol, lab.cex=1.5)

    }
    n
  }
  clusDendro = dendrapply(hcrd, colLab)

  if(testing==TRUE){return(clusDendro)}
  # make plot
  output<-plot(clusDendro, main = "",axes=FALSE)
rm(clusMember,labelColors,pos=1)
return(output)
}
#Supplemental Figure 2b
sf1<-function(){
adf<-(subset(rbind(assignIdentity(makerefmetdfb(norm=0)),assignIdentity(makerefmetdfb(norm=1)),assignIdentity(makerefmetdfb(norm="UQN"))),identity!="unclassified"))
adf$norm<-as.factor(adf$norm);levels(adf$norm)[1]<-"None";levels(adf$norm)[2]<-"Library Size";levels(adf$norm)[3]<-"Upper Quartile" #just some sad cleaning...
adf$norm<-relevel(adf$norm,ref="None")
g<-ggplot(adf)
return(print(g+geom_point(data=subset(adf,identity!="unclassified"),aes(x=log2(a1*a2)/2,y=log2((a1)/a2),col=identity))+geom_point(data=subset(adf,identity=="external"),aes(x=log2(a1*a2)/2,y=log2((a1)/a2)),size=3)+facet_grid(. ~ norm)+geom_hline(yintercept=0)+ggtitle("Normalization")+
  theme(legend.position=c(0.9,0.9))+theme(strip.background = element_rect(fill = 'white'))+ylim(c(-2,2))+scale_colour_manual(values=c("#CC6666","black","#99CC66","#6699CC"),labels=c("Brain","ERCC","Liver","Muscle"))+xlim(c(0,18))+
  theme(legend.title=element_blank())+theme(title=element_text(size=20))+
  xlab("Mean Counts (BLM-1+BLM-2)")+ylab("Ratio of counts (BLM-1/BLM-2)")+theme(strip.text.x=element_text(size=24))+theme(axis.title=element_text(size=24))+theme(axis.text=element_text(size=20))+
  guides(colour=guide_legend(override.aes=list(size=5)))))
}#MA plots : by normalization
sf2a<-function(){
  gdfa<-(assignIdentity(makerefmetdfb(norm="UQN")))
  gdfa$plat<-"HiSeq"
  gdfa$mrat<-calcmrnafrac(gdfa)[7]/calcmrnafrac(gdfa)[1] # since the plot is comparing A1 to A2, I take the ratio of their mrna fractions as a term to use later.

  gdfb<-assignIdentity(makerefmetdfb(norm=0,platform = "5500"))
  gdfb$plat<-"SOLiD 5500"
  gdfb$mrat<-calcmrnafrac(gdfb)[8]/calcmrnafrac(gdfb)[1] #this plot is comparing A1 to B2, because i want to look at the ERCC pool effects - thus this mrna fraction ratio correction term
#  rownames(gdfa)<-gdfa$gene_id;rownames(gdfb)<-gdfb$gene_id
  gdf<-rbind(gdfa,gdfb) #it's later now.
  gdf$pme[gdf$plat=="SOLiD 5500"]<-gdfb[,9]#/sum(gdfb[,9])*sum(gdfb[,2]) #this printme variable corresponds to B2 in the data, because yeah that is what's going on in the 5500 data - but is normalized by the library size.
  gdf$pme[gdf$plat=="HiSeq"]<-gdfa[,8]#/sum(gdfa[,8])*sum(gdfa[,2]) #this "printme" variable is A2. because. A2 is what i want to look at for this dataset

  g<-ggplot(subset(gdf,identity!="unclassified"))
  g+geom_point(aes(y=log2(a1*mrat/pme),x=log2(a1/mrat*pme)/2,col=identity,size=identity))+facet_grid(. ~ plat)+scale_colour_manual(values=c("#CC6666","black","#99CC66","#6699CC"),labels=c("Brain","ERCC","Liver","Muscle"))+
    theme(legend.position=c(0.9,0.9))+labs(x="A",y=expression(paste(Log[2]," Ratio BLM-1 : BLM-2")))+scale_alpha_manual(values=c(0.7,1,0.7,0.7))+
    scale_size_manual(values=c(1.5,3,1.5,1.5),labels=c("Brain","ERCC","Liver","Muscle"))+theme(legend.title=element_blank())+theme(strip.background = element_rect(fill = 'white'))+
    scale_x_continuous(limits=c(0,20),expand=c(0,0))+scale_y_continuous(limits=c(-2,2),expand=c(0,0))+guides(colour=guide_legend(override.aes = list(size=4)))+theme(strip.text.x=element_text(size=24))+
    theme(axis.title=element_text(size=24))+theme(axis.text=element_text(size=20))+geom_point(data=subset(gdf,identity=="external"),aes(y=log2(a1*mrat/pme),x=log2(a1/mrat*pme)/2,col=identity,size=identity))
}#no longer referred to in the paper
sf4<-function(indf){
  require(reshape2)
  if(missing(indf)){
    indf<-makeseqcdf()}
  #do UQN for ILM sites:
  indfILM<-subset(indf,!site%in%c("NWU","PSU","SQW")) #subset to ILM only
  indfILM$site<-droplevels(indfILM$site) #clear out those pesky empty factor levels
  indfILM<-normalizeSEQC(indfILM,type="Both")
  indfLT<-subset(indf,site%in%c("NWU","PSU","SQW"))
  indf<-rbind(indfILM,indfLT)
  final<-NULL
  dset<-cbind(indf[,1:2],rowMeans(indf[,3:6]),rowMeans(indf[,7:10]),(indf[,11:14]),(indf[,15:18]))
  names(dset)<-c("gene_id","site","A","B","C","C","C","C","D","D","D","D")

  for(I in levels(as.factor(dset$site))){

    modeled<-SEQClm(subset(indf,site==I))
    modeled<-suppressMessages(melt(modeled));modeled$site<-I;modeled$mix<-c(rep("C",4),rep("D",4));modeled$variable<-"A"
    final<-rbind(modeled,final)
  }
  #calculates meanvalues
  final$method<-"mRNA-correcting"
  final<-cbind(subset(final,mix=="C"),subset(final,mix=="D")[,2])
  names(final)[6]<-"dval"
  final$site<-as.factor(final$site)
  levels(final$site)<-c("Lab 1 - ILM", "Lab 2 - ILM", "Lab 3 - ILM", "Lab 4 - ILM", "Lab 5 - ILM" , "Lab 6 - ILM" , "Lab 7 - LT" , "Lab 8 - LT", "Lab 9 - LT")
  require(dplyr)
  fmean<-group_by(final,method,site,variable)%>%summarise(mean(value),mean(dval),sd(value),sd(dval))
  names(fmean)<-c("method","site","variable","meanc","meand","sdc","sdd")

  g<-ggplot(final)
  g+geom_point(aes(x=value,y=dval,color=(site),pch=method),alpha=0.65,size=4)+coord_cartesian(ylim=c(0.1,0.4),xlim=c(0.6,0.9))+
    facet_wrap(~ site)+ylab("Amount of SEQC-A in SEQC-C")+xlab("Amount of SEQC-A in SEQC-D")+geom_point(aes(x=0.75,y=0.25),col="grey70")+theme(legend.position="none")+
    geom_path(data=circleFun(center=c(0.75,0.25),diameter=0.1,npoints=25),aes(x,y),col="grey")+
    geom_path(data=circleFun(center=c(0.75,0.25),diameter=0.05,npoints=25),aes(x,y),col="grey")+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(panel.margin=unit(1,"cm"))+theme(aspect.ratio=1)+
    theme(legend.text=element_text(size=rel(1.4)))+theme(axis.title=element_text(size=rel(1.6)))+theme(axis.text=element_text(size=rel(1)))+theme(strip.background = element_rect(fill = 'white'))+theme(strip.text=element_text(size=rel(1.3)))+
    geom_pointrange(data=fmean,aes(x=meanc,y=meand,ymax=meand+2*sdd,ymin=meand-2*sdd),size=1.15,shape=1)+geom_errorbarh(data=fmean,aes(x=meanc,y=meand,xmax=meanc+2*sdc,xmin=meanc-2*sdc),size=1.3)
  #if axes along each facet are needed, http://stackoverflow.com/questions/17661052/force-x-axis-text-on-for-all-facets-of-a-facet-grid-plot # it's not simple though...
}#SF4.  SEQC-C and D estimates for interlab.

#Supplemental Figure3:
sf3<-function(){figure<-makeTargetPlot(df=makerefmetdfb(norm="UQN",mapper="RSEM",method="RSEM",readtype="mostcomplex"),
                               df2 = makerefmetdfb(norm="UQN",mapper="RSEM",method="RSEM",readtype = "default"),Discriminator = "FPKM",numrings = 2)
figure+scale_alpha_manual(name="Unit",breaks=c(0,1),labels=c("Count","FPKM"),values=c(0.3,1))}#Target plot: FPKM
#Supplemental Figure4:
sf5<-function(indf){
  require(reshape2)
if(missing(indf)){
  indf<-makeseqcdf()}
#indf<-normalizeSEQC(indf) #apparently normalization of these data pre-analysis make for some issues...but only for naive/mymethod, decon is unaffected.

final<-NULL;finalMod<-NULL;finalNaive<-NULL
dset<-cbind(indf[,1:2],rowMeans(indf[,3:6]),rowMeans(indf[,7:10]),(indf[,11:14]),(indf[,15:18]))
names(dset)<-c("gene_id","site","A","B","C","C","C","C","D","D","D","D")

require("DeconRNASeq")
###SEQClm(subset(SEQCDF,site=="AGR"))  #Debugging - the data that're getting plotted are NOT these...
###value      Dval
###c1 0.7412337 0.2505380
###c2 0.7458750 0.2577192
###c3 0.7427236 0.2474030
###c4 0.7599289 0.2429792
###SEQClm(subset(SEQCDF,site=="AGR"),ignoremrna=TRUE)
###value      Dval
###c1 0.7864730 0.3006213
###c2 0.7905321 0.3086469
###c3 0.7877769 0.2971081
###c4 0.8027694 0.2921405

for(I in levels(as.factor(dset$site))){
  tmp<-DeconRNASeq(datasets=subset(dset,site==I)[c(5:12)],signatures=subset(dset,site==I)[,c(3,4)])
  tmp<-tmp$out.all
  tmp<-as.data.frame(tmp)
  tmp$mix<-c(rep("C",4),rep("D",4))
  tmp<-melt(tmp)
  tmp$site<-I
  naive<-SEQClm(subset(indf,site==I),ignoremrna=TRUE)
  naive<-melt(naive);naive$site<-I;naive$mix<-c(rep("C",4),rep("D",4));naive$variable<-"A"
  modeled<-SEQClm(subset(indf,site==I))
  modeled<-melt(modeled);modeled$site<-I;modeled$mix<-c(rep("C",4),rep("D",4));modeled$variable<-"A"
  final<-rbind(tmp,final)
  finalNaive<-rbind(naive,finalNaive)
  finalMod<-rbind(modeled,finalMod)
}
#calculates meanvalues
final$method<-"DeconRNASeq"
finalNaive$method<-"Naive"
finalMod$method<-"mRNA-correcting"
final<-rbind(final,finalNaive,finalMod)
final<-cbind(subset(final,mix=="C"),subset(final,mix=="D")[,3])
names(final)[6]<-"dval"
final$site<-as.factor(final$site)
levels(final$site)<-c("Lab 1 - ILM", "Lab 2 - ILM", "Lab 3 - ILM", "Lab 4 - ILM", "Lab 5 - ILM" , "Lab 6 - ILM" , "Lab 7 - LT" , "Lab 8 - LT", "Lab 9 - LT")
require(dplyr)
fmean<-group_by(final,method,site,variable)%>%summarise(mean(value),mean(dval),sd(value),sd(dval))
names(fmean)<-c("method","site","variable","meanc","meand","sdc","sdd")

g<-ggplot(final)
g+geom_point(aes(x=value,y=dval,color=(method),pch=method),alpha=0.65,size=4)+coord_cartesian(ylim=c(0.1,0.4),xlim=c(0.6,0.9))+
  facet_wrap(~ site)+ylab("Amount of SEQC-A in SEQC-C")+xlab("Amount of SEQC-A in SEQC-D")+geom_point(aes(x=0.75,y=0.25),col="grey70")+theme(legend.position="none")+
  geom_path(data=circleFun(center=c(0.75,0.25),diameter=0.1,npoints=25),aes(x,y),col="grey")+
  geom_path(data=circleFun(center=c(0.75,0.25),diameter=0.05,npoints=25),aes(x,y),col="grey")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(panel.margin=unit(1,"cm"))+theme(aspect.ratio=1)+
  theme(legend.text=element_text(size=rel(1.4)))+theme(axis.title=element_text(size=rel(1.6)))+theme(axis.text=element_text(size=rel(1)))+theme(strip.background = element_rect(fill = 'white'))+theme(strip.text=element_text(size=rel(1.3)))+
  geom_pointrange(data=fmean,aes(x=meanc,y=meand,ymax=meand+2*sdd,ymin=meand-2*sdd),size=1.15,shape=1)+geom_errorbarh(data=fmean,aes(x=meanc,y=meand,xmax=meanc+2*sdc,xmin=meanc-2*sdc),size=1.3)
#if axes along each facet are needed, http://stackoverflow.com/questions/17661052/force-x-axis-text-on-for-all-facets-of-a-facet-grid-plot # it's not simple though...
}#SF5.  Target plot: Naive vs Decon vs Us
#Supplemental Table1:
stab1<-function(){
Figure2<-makerefmetdfb(platform = "5500",norm=0) #we use the 5500 data here because the Illumina data has some inconsistencies with the actual spike-in mix amounts.  Those inconsistencies don't affect
#ratios, just the absolute values, and it's a pain to understand/explain, but the 5500 data illustrate it just fine.
count<-rbind(colSums(Figure2b2[match(ercc96$name[ercc96$pool=="C"],Figure2b2$gene_id),2:14]),colSums(Figure2b2[grep("ERCC-",Figure2b2$gene_id,invert=TRUE),2:14]))
#must only consider the C subpool of erccs here because pools A and B have designed differences between BLM1a and BLM1b (mixa having 2x as much pool a as mixb)
#I can't conceptually justify this distinction, but it is absolutely important.
#You can get the correct ratios from pool A or pool B after you account for the differences in designed spike-in amount,
#but for a reason that i don't understand, you can NOT just look at the entire set of 96 spike-ins as a whole and assume that the subpools even out (even though they were designed to do so!)
ercc.targ <- c(.08,.08/8,.08*8,.08,.08/8,.08*8,.08,.08,.08/8,.08*8,.003,.003,.003)
table<-rbind(count[1,]/count[2,],ercc.targ,ercc.targ*count[2,]/count[1,])
rownames(table)<-c("ERCC Count Ratio","ERCC spike proportion","Rho")

return(table)
}

sf2<-function(){ #stuff i did to compare mix models using by uqn/tmm to that achieved by calcmrnafrac
##for SEQC:
seqcmixmodeling<-function(SEQCDF=SEQCDF,nmethod="TMM",factor="ercc"){
  SEQCDFN<-NULL
  nflist<-NULL
  mflist<-NULL
  for(I in levels(SEQCDF$site)){
    tmp<-subset(SEQCDF,site==I)
    fac<-calcNormFactors(tmp[c(3:18)],method=nmethod)
    #first use the normalization factors as they were intended...
    tmp[,c(3:18)]<-tmp[,c(3:18)]*fac
    if(factor=="ercc"){    fac<-calcmrnafrac(tmp,selection=c(3:18),e5="G")}
    if(factor=="none"){  fac<-c(rep(1,8))}
    if(factor==nmethod){}
    tmp$modelC1<-(fac[1]*tmp$A1*.75)+(fac[5]*tmp$B1*.25)
    tmp$modelC1<-tmp$modelC1/sum(tmp$modelC1)*sum(tmp$C1)
    tmp$modelC2<-(fac[2]*tmp$A2*.75)+(fac[6]*tmp$B2*.25)
    tmp$modelC2<-tmp$modelC2/sum(tmp$modelC2)*sum(tmp$C2)
    tmp$modelC3<-(fac[3]*tmp$A3*.75)+(fac[7]*tmp$B3*.25)
    tmp$modelC3<-tmp$modelC3/sum(tmp$modelC3)*sum(tmp$C3)
    tmp$modelC4<-(fac[4]*tmp$A4*.75)+(fac[8]*tmp$B4*.25)
    tmp$modelC4<-tmp$modelC4/sum(tmp$modelC4)*sum(tmp$C4)
    tmp$modelD1<-fac[1]*tmp$A1*.25+fac[5]*tmp$B1*.75
    tmp$modelD1<-tmp$modelD1/sum(tmp$modelD1)*sum(tmp$D1)
    tmp$modelD2<-fac[2]*tmp$A2*.25+fac[6]*tmp$B2*.75
    tmp$modelD2<-tmp$modelD2/sum(tmp$modelD2)*sum(tmp$D2)
    tmp$modelD3<-fac[3]*tmp$A3*.25+fac[7]*tmp$B3*.75
    tmp$modelD3<-tmp$modelD3/sum(tmp$modelD3)*sum(tmp$D3)
    tmp$modelD4<-fac[4]*tmp$A4*.25+fac[8]*tmp$B4*.75
    tmp$modelD4<-tmp$modelD4/sum(tmp$modelD4)*sum(tmp$D4)
    nflist<-rbind(nflist,fac)
    mflist<-rbind(mflist,calcmrnafrac(tmp,selection=c(3:18),e5="G"))
    SEQCDFN<-rbind(SEQCDFN,tmp)

  }
  return(SEQCDFN)
}
blmmixnormmodels<-function(normtype="TMM",factor="ercc"){
  tmp<-makerefmetdfb()
  if(normtype=="TMM"){
    require(edgeR)
    fac<-c(calcNormFactors(tmp[2:14],method="TMM"))#,calcNormFactors(tmp[12:14],method="TMM"))
    #    nfac<-c(nfac,calcNormFactors(tdf[12:14]))
    for(I in 2:14){tmp[I]<-tmp[I]*fac[I-1]}
    #    return(tmp)
  }
  if(normtype=="upperquartile"|normtype=="UQN"){
    require(edgeR)
    fac<-c(calcNormFactors(tmp[2:11],method="upperquartile"),calcNormFactors(tmp[12:14],method="upperquartile"))
    #    nfac<-c(nfac,calcNormFactors(tdf[12:14]))
    for(I in 2:14){tmp[I]<-tmp[I]*fac[I-1]}
    #   return(tmp)
  }
  mixFrac1<-c(.25,.25,.5);mixFrac2<-c(.25,.5,.25)
  if(factor=="ercc"){fac<-calcmrnafrac(tmp)}


  tmp$pred1<-(tmp$bep*mixFrac1[1]*fac[11])+(tmp$lep*mixFrac1[2]*fac[12])+(tmp$mep*mixFrac1[3]*fac[13])
  tmp$pred2<-(tmp$bep*mixFrac2[1]*fac[11])+(tmp$lep*mixFrac2[2]*fac[12])+(tmp$mep*mixFrac2[3]*fac[13])
  tmp$pred1<-tmp$pred1/sum(tmp$pred1)*sum(tmp$a1)
  tmp$pred2<-tmp$pred2/sum(tmp$pred2)*sum(tmp$a2)
  return(tmp)
}
tmm<-seqcmixmodeling(makeseqcdf(),factor="none")
uqn<-seqcmixmodeling(makeseqcdf(),nmethod="upperquartile",factor="none")
tmm$norm<-"tmm"
uqn$norm<-"uqn"
fdf<-rbind(tmm,uqn)
fdf$xm<-rowMeans(fdf[,c("D1","D2","D3","D4")])
fdf$Mxm<-rowMeans(fdf[,c("modelD1","modelD2","modelD3","modelD4")])
fdf$dset<-"SEQC"
tmm<-blmmixnormmodels(factor="none")
uqn<-blmmixnormmodels(normtype="upperquartile",factor="none")
tmm$norm<-"tmm"
uqn$norm<-"uqn"
mdf<-rbind(tmm,uqn)
mdf$dset<-"BLM"
mdf$xm<-mdf$a1
mdf$Mxm<-mdf$pred1
odf<-rbind(fdf[,c("Mxm","xm","norm","dset")],mdf[,c("Mxm","xm","norm","dset")])
#none<-blmmixnormmodels()
g<-ggplot(summarise(group_by(odf,norm,dset),rmsd=sqrt(mean((Mxm-xm)^2))))
fig<-g+geom_bar(aes(x=norm,y=rmsd),stat="identity")+facet_wrap(~dset)
return(fig)
}#figure at least demonstrates the gist of what i'm going for here...

circleFun<-function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}#Generates dataframes of npoints points that form a circle centered at point center.

makerefmetdfb<-function(platform="ILLUMINA",mapper="TH",method="HTS",reference="UCSC",readtype="default",norm=0,reps=1,basedir="PAPERUSED/",...){
  #goal is to make this updated to use fread rather than fread because slow...problems:
  require(data.table)
{  string<-paste(platform,mapper,method,reference,sep="/")
   a1<-paste(basedir,string,"/BLM1a.ct",sep="")
   a1d<-paste(basedir,string,"/BLM1aD.ct",sep="")
   a1u<-paste(basedir,string,"/BLM1aU.ct",sep="")
   b1<-paste(basedir,string,"/BLM1b.ct",sep="")
   b1d<-paste(basedir,string,"/BLM1bD.ct",sep="")
   b1u<-paste(basedir,string,"/BLM1bU.ct",sep="")
   a2<-paste(basedir,string,"/BLM2a.ct",sep="")
   b2<-paste(basedir,string,"/BLM2b.ct",sep="")
   b2d<-paste(basedir,string,"/BLM2bD.ct",sep="")
   b2u<-paste(basedir,string,"/BLM2bU.ct",sep="")
   bep<-paste(basedir,string,"/BEP.ct",sep="")
   lep<-paste(basedir,string,"/LEP.ct",sep="")
   mep<-paste(basedir,string,"/MEP.ct",sep="")
   clct<-paste(basedir,string,"/CL.ct",sep="")
   bct<-paste(basedir,string,"/BCT.ct",sep="")
   lct<-paste(basedir,string,"/LCT.ct",sep="")
   mct<-paste(basedir,string,"/MCT.ct",sep="")
   if(file.exists(b1)==FALSE){if(file.exists(clct)==FALSE){return(0)}}##check to see if there is actually data to collect.  If not, don't throw an error.  How massive!
   if(length(grep("HTS",method))==1){
     int<-fread(a1,sep="\t",select=1)
     a1<-rowSums(fread(a1,sep="\t",drop=1))
     a1d<-rowSums(fread(a1d,sep="\t",drop=1))
     a1u<-rowSums(fread(a1u,sep="\t",drop=1))
     b1<-rowSums(fread(b1,sep="\t",drop=1))
     b1d<-rowSums(fread(b1d,sep="\t",drop=1))
     b1u<-rowSums(fread(b1u,sep="\t",drop=1))
     a2<-rowSums(fread(a2,sep="\t",drop=1))
     b2<-rowSums(fread(b2,sep="\t",drop=1))
     b2d<-rowSums(fread(b2d,sep="\t",drop=1))
     b2u<-rowSums(fread(b2u,sep="\t",drop=1))
     bep<-rowSums(fread(bep,sep="\t",drop=1))
     lep<-rowSums(fread(lep,sep="\t",drop=1))
     mep<-rowSums(fread(mep,sep="\t",drop=1))

   }#some specific bits to make input from hts
   else if(length(grep("LS",method))==1){
     if(platform=="5500"){
       int<-(fread(a1,sep="\t",select=1))
       a1<-fread(a1,sep="\t",select=6)#6 is count, 8 is RPKM
       a1d<-(fread(a1d,sep="\t",select=6))
       a1u<-(fread(a1u,sep="\t",select=6))
       b1<-(fread(b1,sep="\t",select=6))
       b1d<-(fread(b1d,sep="\t",select=6))
       b1u<-(fread(b1u,sep="\t",select=6))
       a2<-(fread(a2,sep="\t",select=6))
       b2<-(fread(b2,sep="\t",select=6))
       b2d<-(fread(b2d,sep="\t",select=6))
       b2u<-(fread(b2u,sep="\t",select=6))
       bep<-(fread(bep,sep="\t",select=6))
       lep<-(fread(lep,sep="\t",select=6))
       mep<-(fread(mep,sep="\t",select=6))
     }
     else if (platform=="SOLID4"){
       int<-fread(a1,select=1)
       a1<-rowSums(fread(a1)[2:3])
       a1d<-rowSums(fread(a1d)[2:3])
       a1u<-rowSums(fread(a1u)[2:3])
       b1<-rowSums(fread(b1)[2:3])
       b1d<-rowSums(fread(b1d)[2:3])
       b1u<-rowSums(fread(b1u)[2:3])
       a2<-rowSums(fread(a2)[2:3])
       b2<-rowSums(fread(b2)[2:3])
       b2d<-rowSums(fread(b2d)[2:3])
       b2u<-rowSums(fread(b2u)[2:3])
       bep<-rowSums(fread(bep)[2:3])
       lep<-rowSums(fread(lep)[2:3])
       mep<-rowSums(fread(mep)[2:3])
     }
   }# LS input
   else if (method=="CUFF"){
     if(reps==1){
       int<-(fread(a1,sep="\t",select=1))
       a1<-(fread(a1,sep="\t",header=FALSE,select=8))

       a1d<-(fread(a1d,sep="\t",header=FALSE,select=8))
       a1u<-(fread(a1u,sep="\t",header=FALSE,select=8))
       b1<-(fread(b1,sep="\t",header=FALSE,select=8))
       b1d<-(fread(b1d,sep="\t",header=FALSE,select=8))
       b1u<-(fread(b1u,sep="\t",header=FALSE,select=8))
       a2<-(fread(a2,sep="\t",header=FALSE,select=8))
       b2<-(fread(b2,sep="\t",header=FALSE,select=8))
       b2d<-(fread(b2d,sep="\t",header=FALSE,select=8)) #this one might need to switch to 8 depending on how i set it up....instead use reversecuffdiffoutput.pl
       b2u<-(fread(b2u,sep="\t",header=FALSE,select=8)) #if so, this one should be back to 7, probably...
       bep<-(fread(bep,sep="\t",header=FALSE,select=8))
       lep<-(fread(lep,sep="\t",header=FALSE,select=8))
       mep<-(fread(mep,sep="\t",header=FALSE,select=8)) #same with this one
     }
     else if(reps==0){
       int<-(fread(a1,sep="\t",select=1))
       a1<-(fread(a1,sep="\t",header=FALSE,select=8))
       a1d<-(fread(a1d,sep="\t",header=FALSE,select=8))
       a1u<-(fread(a1u,sep="\t",header=FALSE,select=8))
       b1<-(fread(b1,sep="\t",header=FALSE,select=8))
       b1d<-(fread(b1d,sep="\t",header=FALSE,select=8))
       b1u<-(fread(b1u,sep="\t",header=FALSE,select=8))
       a2<-(fread(a2,sep="\t",header=FALSE,select=8))
       b2<-(fread(b2,sep="\t",header=FALSE,select=8))
       b2d<-(fread(b2d,sep="\t",header=FALSE,select=8)) #this one might need to switch to 8 depending on how i set it up....instead use reversecuffdiffoutput.pl
       b2u<-(fread(b2u,sep="\t",header=FALSE,select=8)) #if so, this one should be back to 7, probably...
       bep<-(fread(bep,sep="\t",header=FALSE,select=8))
       lep<-(fread(lep,sep="\t",header=FALSE,select=8))
       mep<-(fread(mep,sep="\t",header=FALSE,select=8)) #same with this one
     }
   }#more specific bits for CUFF
   else if (method=="CUFFtx"|method=="CUFFTX"){
     if(reps==1){
       int<-(fread(a1,sep="\t",select=1))
       a1<-(fread(a1,sep="\t",header=FALSE,select=8))
       a1d<-(fread(a1d,sep="\t",header=FALSE,select=8))
       a1u<-(fread(a1u,sep="\t",header=FALSE,select=8))
       b1<-(fread(b1,sep="\t",header=FALSE,select=8))
       b1d<-(fread(b1d,sep="\t",header=FALSE,select=8))
       b1u<-(fread(b1u,sep="\t",header=FALSE,select=8))
       a2<-(fread(a2,sep="\t",header=FALSE,select=8))
       b2<-(fread(b2,sep="\t",header=FALSE,select=8))
       b2d<-(fread(b2d,sep="\t",header=FALSE,select=8)) #this one might need to switch to 8 depending on how i set it up....instead use reversecuffdiffoutput.pl
       b2u<-(fread(b2u,sep="\t",header=FALSE,select=8)) #if so, this one should be back to 7, probably...
       bep<-(fread(bep,sep="\t",header=FALSE,select=8))
       lep<-(fread(lep,sep="\t",header=FALSE,select=8))
       mep<-(fread(mep,sep="\t",header=FALSE,select=8)) #same with this one
     }
     else if(reps==0){
       int<-(fread(a1,sep="\t",select=1))
       a1<-(fread(a1,sep="\t",header=FALSE,select=8))
       a1d<-(fread(a1d,sep="\t",header=FALSE,select=8))
       a1u<-(fread(a1u,sep="\t",header=FALSE,select=8))
       b1<-(fread(b1,sep="\t",header=FALSE,select=8))
       b1d<-(fread(b1d,sep="\t",header=FALSE,select=8))
       b1u<-(fread(b1u,sep="\t",header=FALSE,select=8))
       a2<-(fread(a2,sep="\t",header=FALSE,select=8))
       b2<-(fread(b2,sep="\t",header=FALSE,select=8))
       b2d<-(fread(b2d,sep="\t",header=FALSE,select=8)) #this one might need to switch to 8 depending on how i set it up....instead use reversecuffdiffoutput.pl
       b2u<-(fread(b2u,sep="\t",header=FALSE,select=8)) #if so, this one should be back to 7, probably...
       bep<-(fread(bep,sep="\t",header=FALSE,select=8))
       lep<-(fread(lep,sep="\t",header=FALSE,select=8))
       mep<-(fread(mep,sep="\t",header=FALSE,select=8)) #same with this one
     }
   }#more specific bits for CUFF
   else if(method=="CLCT"){

     if(mapper=="LS"){
       a1<-(fread(clct,sep="\t",header=TRUE,select=1))
       b1<-(fread(clct,sep="\t",header=TRUE,select=16))
       a2<-(fread(clct,sep="\t",header=TRUE,select=31))
       b2<-(fread(clct,sep="\t",header=TRUE,select=36))
       b2d<-(fread(clct,sep="\t",header=TRUE,select=46))
       a1d<-(fread(clct,sep="\t",header=TRUE,select=11))
       b1d<-(fread(clct,sep="\t",header=TRUE,select=26))
       a1u<-(fread(clct,sep="\t",header=TRUE,select=6))
       b1u<-(fread(clct,sep="\t",header=TRUE,select=21))
       b2u<-(fread(clct,sep="\t",header=TRUE,select=41))
       bep<-fread(bct,sep="\t",header=TRUE,select=1)
       lep<-fread(lct,sep="\t",header=TRUE,select=6)
       mep<-fread(mct,sep="\t",header=TRUE,select=11)
       #idk what int should be here...do i have to recalibrate everything?

     }
     else if(mapper=="TH"){
       if(platform=="ILLUMINA"){
         a1<-(fread(clct,sep="\t",header=TRUE,select=1))
         b1<-(fread(clct,sep="\t",header=TRUE,select=6))
         a2<-(fread(clct,sep="\t",header=TRUE,select=11))
         b2<-(fread(clct,sep="\t",header=TRUE,select=16))
         b2d<-(fread(clct,sep="\t",header=TRUE,select=21))
         a1d<-(fread(clct,sep="\t",header=TRUE,select=26))
         b1d<-(fread(clct,sep="\t",header=TRUE,select=31))
         a1u<-(fread(clct,sep="\t",header=TRUE,select=36))
         b1u<-(fread(clct,sep="\t",header=TRUE,select=41))
         b2u<-(fread(clct,sep="\t",header=TRUE,select=46))
         bep<-fread(bct,sep="\t",header=TRUE,select=1)
         lep<-fread(lct,sep="\t",header=TRUE,select=6)
         mep<-fread(mct,sep="\t",header=TRUE,select=11)
       }
       else if(platform=="5500"){
         a1<-(fread(clct,sep="\t",header=TRUE,select=1))
         b1<-(fread(clct,sep="\t",header=TRUE,select=16))
         a2<-(fread(clct,sep="\t",header=TRUE,select=31))
         b2<-(fread(clct,sep="\t",header=TRUE,select=46))
         b2d<-(fread(clct,sep="\t",header=TRUE,select=36))
         a1d<-(fread(clct,sep="\t",header=TRUE,select=6))
         b1d<-(fread(clct,sep="\t",header=TRUE,select=26))
         a1u<-(fread(clct,sep="\t",header=TRUE,select=11))
         b1u<-(fread(clct,sep="\t",header=TRUE,select=21))
         b2u<-(fread(clct,sep="\t",header=TRUE,select=41))
         bep<-fread(bct,sep="\t",header=TRUE,select=1)
         lep<-fread(lct,sep="\t",header=TRUE,select=6)
         mep<-fread(mct,sep="\t",header=TRUE,select=11)

       }
       else if(platform=="SOLID4"){
         a1<-(fread(clct,sep="\t",header=TRUE,select=1))
         a1d<-(fread(clct,sep="\t",header=TRUE,select=6))
         a1u<-(fread(clct,sep="\t",header=TRUE,select=11))
         b1<-(fread(clct,sep="\t",header=TRUE,select=16))
         b1d<-(fread(clct,sep="\t",header=TRUE,select=21))
         b1u<-(fread(clct,sep="\t",header=TRUE,select=26))
         a2<-(fread(clct,sep="\t",header=TRUE,select=31))
         b2<-(fread(clct,sep="\t",header=TRUE,select=36))
         b2d<-(fread(clct,sep="\t",header=TRUE,select=41))
         b2u<-(fread(clct,sep="\t",header=TRUE,select=46))
         bep<-fread(bct,sep="\t",header=TRUE,select=1)
         lep<-fread(lct,sep="\t",header=TRUE,select=6)
         mep<-fread(mct,sep="\t",header=TRUE,select=11)

       }
     }

   }
   else if(method=="SAIL"){
     if(readtype=="default"){a1t<-a1
                             a1<-(fread(a1,sep="\t",select=3))###3 is TPM, 4 is FPKM
                             a1d<-(fread(a1d,sep="\t",select=3))
                             a1u<-(fread(a1u,sep="\t",select=3))
                             b1<-(fread(b1,sep="\t",select=3))
                             b1d<-(fread(b1d,sep="\t",select=3))
                             b1u<-(fread(b1u,sep="\t",select=3))
                             a2<-(fread(a2,sep="\t",select=3))
                             b2<-(fread(b2,sep="\t",select=3))
                             b2d<-(fread(b2d,sep="\t",select=3))
                             b2u<-(fread(b2u,sep="\t",select=3))
                             bep<-(fread(bep,sep="\t",select=3))
                             lep<-(fread(lep,sep="\t",select=3))
                             mep<-(fread(mep,sep="\t",select=3))
                             int<-(fread(a1t,sep="\t",select=1))
     }
     else if(readtype=="mostcomplex"){a1t<-a1
                                      a1<-(fread(a1,sep="\t",select=4))###3 is TPM, 4 is FPKM
                                      a1d<-(fread(a1d,sep="\t",select=4))
                                      a1u<-(fread(a1u,sep="\t",select=4))
                                      b1<-(fread(b1,sep="\t",select=4))
                                      b1d<-(fread(b1d,sep="\t",select=4))
                                      b1u<-(fread(b1u,sep="\t",select=4))
                                      a2<-(fread(a2,sep="\t",select=4))
                                      b2<-(fread(b2,sep="\t",select=4))
                                      b2d<-(fread(b2d,sep="\t",select=4))
                                      b2u<-(fread(b2u,sep="\t",select=4))
                                      bep<-(fread(bep,sep="\t",select=4))
                                      lep<-(fread(lep,sep="\t",select=4))
                                      mep<-(fread(mep,sep="\t",select=4))
                                      int<-(fread(a1t,sep="\t",select=1))
     }



   }
   else if (method=="RSEM"){
     if(readtype=="default"){a1t<-a1
                             a1<-(fread(a1,sep="\t",select=7))###5 is ExpectedCount, 6 is TPM, 7 is FPKM
                             a1d<-(fread(a1d,sep="\t",select=7))
                             a1u<-(fread(a1u,sep="\t",select=7))
                             b1<-(fread(b1,sep="\t",select=7))
                             b1d<-(fread(b1d,sep="\t",select=7))
                             b1u<-(fread(b1u,sep="\t",select=7))
                             a2<-(fread(a2,sep="\t",select=7))
                             b2<-(fread(b2,sep="\t",select=7))
                             b2d<-(fread(b2d,sep="\t",select=7))
                             b2u<-(fread(b2u,sep="\t",select=7))
                             bep<-(fread(bep,sep="\t",select=7))
                             lep<-(fread(lep,sep="\t",select=7))
                             mep<-(fread(mep,sep="\t",select=7))
                             int<-(fread(a1t,sep="\t",select=1))
     }
     if(readtype=="mostcomplex"){
       a1t<-a1
       a1<-(fread(a1,sep="\t",select=6))###5 is ExpectedCount, 6 is TPM, 7 is FPKM
       a1d<-(fread(a1d,sep="\t",select=6))
       a1u<-(fread(a1u,sep="\t",select=6))
       b1<-(fread(b1,sep="\t",select=6))
       b1d<-(fread(b1d,sep="\t",select=6))
       b1u<-(fread(b1u,sep="\t",select=6))
       a2<-(fread(a2,sep="\t",select=6))
       b2<-(fread(b2,sep="\t",select=6))
       b2d<-(fread(b2d,sep="\t",select=6))
       b2u<-(fread(b2u,sep="\t",select=6))
       bep<-(fread(bep,sep="\t",select=6))
       lep<-(fread(lep,sep="\t",select=6))
       mep<-(fread(mep,sep="\t",select=6))
       int<-(fread(a1t,sep="\t",select=1))
     }
     if(readtype=="mostbasic"){
       a1t<-a1
       a1<-(fread(a1,sep="\t",select=5))###5 is ExpectedCount, 6 is TPM, 7 is FPKM
       a1d<-(fread(a1d,sep="\t",select=5))
       a1u<-(fread(a1u,sep="\t",select=5))
       b1<-(fread(b1,sep="\t",select=5))
       b1d<-(fread(b1d,sep="\t",select=5))
       b1u<-(fread(b1u,sep="\t",select=5))
       a2<-(fread(a2,sep="\t",select=5))
       b2<-(fread(b2,sep="\t",select=5))
       b2d<-(fread(b2d,sep="\t",select=5))
       b2u<-(fread(b2u,sep="\t",select=5))
       bep<-(fread(bep,sep="\t",select=5))
       lep<-(fread(lep,sep="\t",select=5))
       mep<-(fread(mep,sep="\t",select=5))
       int<-(fread(a1t,sep="\t",select=1))
     }
   }#RSEM input
   else if (method=="RSEMc"){
     a1t<-a1
     a1<-(fread(a1,header=TRUE,sep="\t",select=5))###5 is ExpectedCount, 6 is TPM, 7 is FPKM
     a1d<-(fread(a1d,header=TRUE,sep="\t",select=5))
     a1u<-(fread(a1u,header=TRUE,sep="\t",select=5))
     b1<-(fread(b1,header=TRUE,sep="\t",select=5))
     b1d<-(fread(b1d,header=TRUE,sep="\t",select=5))
     b1u<-(fread(b1u,header=TRUE,sep="\t",select=5))
     a2<-(fread(a2,header=TRUE,sep="\t",select=5))
     b2<-(fread(b2,header=TRUE,sep="\t",select=5))
     b2d<-(fread(b2d,header=TRUE,sep="\t",select=5))
     b2u<-(fread(b2u,header=TRUE,sep="\t",select=5))
     bep<-(fread(bep,header=TRUE,sep="\t",select=5))
     lep<-(fread(lep,header=TRUE,sep="\t",select=5))
     mep<-(fread(mep,header=TRUE,sep="\t",select=5))
     int<-(fread(a1t,header=TRUE,sep="\t",select=1))
     int2<-(fread(a1t,header=TRUE,sep="\t",select=2))
   }
   else if (method=="RSEMt"){
     if(readtype=="default"){
       a1t<-a1
       a1<-(fread(a1,sep="\t",header=TRUE,select=5))###5 is ExpectedCount, 6 is TPM, 7 is FPKM
       a1d<-(fread(a1d,sep="\t",header=TRUE,select=5))
       a1u<-(fread(a1u,sep="\t",header=TRUE,select=5))
       b1<-(fread(b1,sep="\t",header=TRUE,select=5))
       b1d<-(fread(b1d,sep="\t",header=TRUE,select=5))
       b1u<-(fread(b1u,sep="\t",header=TRUE,select=5))
       a2<-(fread(a2,sep="\t",header=TRUE,select=5))
       b2<-(fread(b2,sep="\t",header=TRUE,select=5))
       b2d<-(fread(b2d,sep="\t",header=TRUE,select=5))
       b2u<-(fread(b2u,sep="\t",header=TRUE,select=5))
       bep<-(fread(bep,sep="\t",header=TRUE,select=5))
       lep<-(fread(lep,sep="\t",header=TRUE,select=5))
       mep<-(fread(mep,sep="\t",header=TRUE,select=5))
       int<-(fread(a1t,sep="\t",header=TRUE,select=1))
       int2<-(fread(a1t,sep="\t",header=TRUE,select=2))
     }
     else if(readtype=="mostbasic"){
       a1t<-a1
       a1<-(fread(a1,sep="\t",header=TRUE,select=6))###5 is ExpectedCount, 6 is TPM, 7 is FPKM
       a1d<-(fread(a1d,sep="\t",header=TRUE,select=6))
       a1u<-(fread(a1u,sep="\t",header=TRUE,select=6))
       b1<-(fread(b1,sep="\t",header=TRUE,select=6))
       b1d<-(fread(b1d,sep="\t",header=TRUE,select=6))
       b1u<-(fread(b1u,sep="\t",header=TRUE,select=6))
       a2<-(fread(a2,sep="\t",header=TRUE,select=6))
       b2<-(fread(b2,sep="\t",header=TRUE,select=6))
       b2d<-(fread(b2d,sep="\t",header=TRUE,select=6))
       b2u<-(fread(b2u,sep="\t",header=TRUE,select=6))
       bep<-(fread(bep,sep="\t",header=TRUE,select=6))
       lep<-(fread(lep,sep="\t",header=TRUE,select=6))
       mep<-(fread(mep,sep="\t",header=TRUE,select=6))
       int<-(fread(a1t,sep="\t",header=TRUE,select=1))
       int2<-(fread(a1t,sep="\t",header=TRUE,select=2))
     }
     if(readtype=="mostcomplex"){
       a1t<-a1
       a1<-(fread(a1,sep="\t",header=TRUE,select=7))###5 is ExpectedCount, 6 is TPM, 7 is FPKM
       a1d<-(fread(a1d,sep="\t",header=TRUE,select=7))
       a1u<-(fread(a1u,sep="\t",header=TRUE,select=7))
       b1<-(fread(b1,sep="\t",header=TRUE,select=7))
       b1d<-(fread(b1d,sep="\t",header=TRUE,select=7))
       b1u<-(fread(b1u,sep="\t",header=TRUE,select=7))
       a2<-(fread(a2,sep="\t",header=TRUE,select=7))
       b2<-(fread(b2,sep="\t",header=TRUE,select=7))
       b2d<-(fread(b2d,sep="\t",header=TRUE,select=7))
       b2u<-(fread(b2u,sep="\t",header=TRUE,select=7))
       bep<-(fread(bep,sep="\t",header=TRUE,select=7))
       lep<-(fread(lep,sep="\t",header=TRUE,select=7))
       mep<-(fread(mep,sep="\t",header=TRUE,select=7))
       int<-(fread(a1t,sep="\t",header=TRUE,select=1))
       int2<-(fread(a1t,sep="\t",header=TRUE,select=2))
     }

   }#RSEM transcriptome input
   else if (method=="xpress"|method=="xpress2"|method=="xpress3"){
     if(readtype=="default"){
       a1t<-a1
       a1<-(fread(a1,sep="\t",header=TRUE,select=8))###7 is expected counts, 8 is effective counts
       a1d<-(fread(a1d,sep="\t",header=TRUE,select=8))
       a1u<-(fread(a1u,sep="\t",header=TRUE,select=8))
       b1<-(fread(b1,sep="\t",header=TRUE,select=8))
       b1d<-(fread(b1d,sep="\t",header=TRUE,select=8))
       b1u<-(fread(b1u,sep="\t",header=TRUE,select=8))
       a2<-(fread(a2,sep="\t",header=TRUE,select=8))
       b2<-(fread(b2,sep="\t",header=TRUE,select=8))
       b2d<-(fread(b2d,sep="\t",header=TRUE,select=8))
       b2u<-(fread(b2u,sep="\t",header=TRUE,select=8))
       bep<-(fread(bep,sep="\t",header=TRUE,select=8))
       lep<-(fread(lep,sep="\t",header=TRUE,select=8))
       mep<-(fread(mep,sep="\t",header=TRUE,select=8))
       int<-(fread(a1t,sep="\t",header=TRUE,select=2)) #1 is txid, 2 is gene_id
     }
     else if(readtype=="mostbasic"){
       a1t<-a1
       a1<-(fread(a1,sep="\t",header=TRUE,select=7))###7 is expected counts, 8 is effective counts
       a1d<-(fread(a1d,sep="\t",header=TRUE,select=7))
       a1u<-(fread(a1u,sep="\t",header=TRUE,select=7))
       b1<-(fread(b1,sep="\t",header=TRUE,select=7))
       b1d<-(fread(b1d,sep="\t",header=TRUE,select=7))
       b1u<-(fread(b1u,sep="\t",header=TRUE,select=7))
       a2<-(fread(a2,sep="\t",header=TRUE,select=7))
       b2<-(fread(b2,sep="\t",header=TRUE,select=7))
       b2d<-(fread(b2d,sep="\t",header=TRUE,select=7))
       b2u<-(fread(b2u,sep="\t",header=TRUE,select=7))
       bep<-(fread(bep,sep="\t",header=TRUE,select=7))
       lep<-(fread(lep,sep="\t",header=TRUE,select=7))
       mep<-(fread(mep,sep="\t",header=TRUE,select=7))
       int<-(fread(a1t,sep="\t",header=TRUE,select=2)) #1 is txid, 2 is gene_id
     }
     else if(readtype=="mostcomplex"){
       a1t<-a1
       a1<-(fread(a1,sep="\t",header=TRUE,select=11))###7 is expected counts, 8 is effective counts
       a1d<-(fread(a1d,sep="\t",header=TRUE,select=11))
       a1u<-(fread(a1u,sep="\t",header=TRUE,select=11))
       b1<-(fread(b1,sep="\t",header=TRUE,select=11))
       b1d<-(fread(b1d,sep="\t",header=TRUE,select=11))
       b1u<-(fread(b1u,sep="\t",header=TRUE,select=11))
       a2<-(fread(a2,sep="\t",header=TRUE,select=11))
       b2<-(fread(b2,sep="\t",header=TRUE,select=11))
       b2d<-(fread(b2d,sep="\t",header=TRUE,select=11))
       b2u<-(fread(b2u,sep="\t",header=TRUE,select=11))
       bep<-(fread(bep,sep="\t",header=TRUE,select=11))
       lep<-(fread(lep,sep="\t",header=TRUE,select=11))
       mep<-(fread(mep,sep="\t",header=TRUE,select=11))
       int<-(fread(a1t,sep="\t",header=TRUE,select=2)) #1 is txid, 2 is gene_id
     }


   }#xpress transcriptome input

   if(length(names(a1))>1000){
     tdf<-data.frame(names(a1),a1,a1d,a1u,b1,b1d,b1u,a2,b2,b2d,b2u,bep,lep,mep)} #ensures naming accuracy
   else if(length(grep("RSEM|LS|RSEMt|xpress|SAIL|CUFF",method))>0){
     tdf<-data.frame(int,a1,a1d,a1u,b1,b1d,b1u,a2,b2,b2d,b2u,bep,lep,mep)} #reverts from separate RSEM naming convention to 'standard'
   else if(length(grep("RSEMt|RSEMc",method))>0){
     tdf<-data.frame(int,a1,a1d,a1u,b1,b1d,b1u,a2,b2,b2d,b2u,bep,lep,mep,int2)} #reverts from separate RSEM naming convention to 'standard'}
   else{tdf<-data.frame(int,a1,a1d,a1u,b1,b1d,b1u,a2,b2,b2d,b2u,bep,lep,mep)} #standard naming

   if(length(tdf)==14){
     colnames(tdf)<-c("gene_id","a1","a1d","a1u","b1","b1d","b1u","a2","b2","b2d","b2u","bep","lep","mep")}
   else if(length(tdf)==15){   colnames(tdf)<-c("gene_id","a1","a1d","a1u","b1","b1d","b1u","a2","b2","b2d","b2u","bep","lep","mep","geneid")}
  }
#read in and name count data from output files
{tdf<-tdf[grep("unknown|no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique",tdf$gene_id,invert=TRUE),]
 tdf$norm<-norm
 if(norm>0){tdf<-normcountdf(tdf,norm)}#do norm if norm==1
 #tdf$mean1<-rowMeans(data.frame(tdf$a1,tdf$a1d,tdf$a1u,tdf$b1,tdf$b1d,tdf$b1u))
 #tdf$mean2<-rowMeans(data.frame(tdf$a2,tdf$b2,tdf$b2d,tdf$b2u))
 #tdf$obsM<-log2((tdf$mean1)/(tdf$mean2))
 #tdf$obsM[is.na(tdf$obsM)]<-0
 #tdf$obsM[tdf$obsM==Inf]<-0
 #tdf$obsM[tdf$obsM==-Inf]<-0
 #tdf$obsA<-log2(tdf$mean1*tdf$mean2)/2
 #tdf$obsA[tdf$obsA==-Inf]<-0;tdf$obsA[is.na(tdf$obsA)]<-0
 # tdf$obsM2<-log2(round(tdf$mean1)/round(tdf$mean2))
 #tdf$expM<-log2((.25*tdf$bep+.5*tdf$mep+.25*tdf$lep)/(.25*tdf$bep+.5*tdf$lep+.25*tdf$mep))
 # tdf$expM[is.na(tdf$expM)]<-0
 # tdf$expM[tdf$expM==Inf]<-0
 # tdf$expM[tdf$expM==-Inf]<-0
 # act<-calculateexpectedmrnafractions(tdf)
 #return(act) ##comment this line out when i'm through today (06/19/13)
 # tdf$expM2<-log2((act[1]*tdf$bep+act[2]*tdf$mep+act[3]*tdf$lep)/(act[4]*tdf$bep+act[5]*tdf$lep+act[6]*tdf$mep))
 # tdf$expM2[is.na(tdf$expM2)]<-0
 # tdf$expM2[tdf$expM2==Inf]<-0
 # tdf$expM2[tdf$expM2==-Inf]<-0

 # tdf$expM3<-log2((act[1]*tdf$bep+act[2]*tdf$lep+act[3]*tdf$mep)/(act[4]*tdf$bep+act[5]*tdf$lep+act[6]*tdf$mep))
 # tdf$expM3[is.na(tdf$expM3)]<-0
 # tdf$expM3[tdf$expM3==Inf]<-0
 # tdf$expM3[tdf$expM3==-Inf]<-0

 # tdf$eM1<-((act[1]*tdf$bep+act[2]*tdf$lep+act[3]*tdf$mep))
 # tdf$eM2<-((act[4]*tdf$bep+act[5]*tdf$lep+act[6]*tdf$mep))
 # tdf$reference<-reference
 #tdf$platform<-platform
 #tdf$mapper<-mapper
 #tdf$method<-method
 tdf$uid<-paste(platform,mapper,method,reference,sep="")
 tdf[tdf==Inf]<-NA
 tdf[tdf==-Inf]<-NA
 #tdf$cutoff<-summary(tdf$obsA[tdf$obsA>0])[3]
 #trying to just do it all at once:
}#remove HTS_unknown features ; run normcountdf ; calculate means/MA, add metadata, change Inf to NA
#tdf$var1<-(((tdf$a1-tdf$mean1)^2+(tdf$a1d-tdf$mean1)^2+(tdf$a1u-tdf$mean1)^2+(tdf$b1-tdf$mean1)^2+(tdf$b1u-tdf$mean1)^2+(tdf$b1d-tdf$mean1)^2)*1/5)
#tdf$var2<-(((tdf$a2-tdf$mean2)^2+(tdf$b2-tdf$mean2)^2+(tdf$b2u-tdf$mean2)^2+(tdf$b2d-tdf$mean2)^2)*1/3)
#tdf<-predictmixfromblm(tdf)
#tdf<-bootstrapdMranges(tdf,reps)

return(tdf)



}#reads input files (Must upload those files somewhere & change the code to read those files from the remote server)
makeTargetPlot<-function(type="BLM",df,df2,Discriminator="CheeseValue",numrings=4){
  require(reshape2)
  require(ggplot2)
  theme_set(theme_bw(base_size=16))
if(type=="BLM"){ a<-blmmixfraction(df);
  a$mixnum<-c(1,1,1,1,1,1,2,2,2,2);
  a[11,]<-c(0,0,0,2);a[12,]<-c(0,0,0,2);
  a<-a[c(1:4,7:10),]#we *are* ignoring the points on the line.
  a<-as.data.frame(melt(a,id.vars=4))
  a<-data.frame(a[a$mixnum==1,],a[a$mixnum==2,])
  colnames(a)<-c("mixnum","variable","value1","blah","blah","value2")
}
  if(type=="BLM"){
    if(!missing(df2)){b<-blmmixfraction(df2);
                      b$mixnum<-c(1,1,1,1,1,1,2,2,2,2);
                      b[11,]<-c(0,0,0,2);b[12,]<-c(0,0,0,2);
                      b<-b[c(1:4,7:10),] #we *are* ignoring the points on the line.
                      b<-as.data.frame(melt(b,id.vars=4))
                      b<-data.frame(b[b$mixnum==1,],b[b$mixnum==2,])
                      colnames(b)<-c("mixnum","variable","value1","blah","blah","value2")
                      sumdf2<-data.frame(amean1=c(mean(b$value1[b$variable=="mixpropA"&b$value1>0]),mean(b$value1[b$variable=="mixpropB"&b$value1>0]),mean(b$value1[b$variable=="mixpropC"&b$value1>0])),
                                         amean2=c(mean(b$value2[b$variable=="mixpropA"&b$value2>0]),mean(b$value2[b$variable=="mixpropB"&b$value2>0]),mean(b$value2[b$variable=="mixpropC"&b$value2>0])),
                                         sd1=c(sd(b$value1[b$variable=="mixpropA"&b$value1>0]),sd(b$value1[b$variable=="mixpropB"&b$value1>0]),sd(b$value1[b$variable=="mixpropC"&b$value1>0])),
                                         sd2=c(sd(b$value2[b$variable=="mixpropA"&b$value2>0]),sd(b$value2[b$variable=="mixpropB"&b$value2>0]),sd(b$value2[b$variable=="mixpropC"&b$value2>0])))

    }
    sumdf<-data.frame(amean1=c(mean(a$value1[a$variable=="mixpropA"&a$value1>0]),mean(a$value1[a$variable=="mixpropB"&a$value1>0]),mean(a$value1[a$variable=="mixpropC"&a$value1>0])),
                      amean2=c(mean(a$value2[a$variable=="mixpropA"&a$value2>0]),mean(a$value2[a$variable=="mixpropB"&a$value2>0]),mean(a$value2[a$variable=="mixpropC"&a$value2>0])),
                      sd1=c(sd(a$value1[a$variable=="mixpropA"&a$value1>0]),sd(a$value1[a$variable=="mixpropB"&a$value1>0]),sd(a$value1[a$variable=="mixpropC"&a$value1>0])),
                      sd2=c(sd(a$value2[a$variable=="mixpropA"&a$value2>0]),sd(a$value2[a$variable=="mixpropB"&a$value2>0]),sd(a$value2[a$variable=="mixpropC"&a$value2>0])))

  }
  #do the plotting
if(type=="SEQC"){
  a<-GeneralLMest(dat,spikeID="ERCC_",components=c("A1","A2","A3","A4","B1","B2","B3","B4"),mixes=c("C1","C2","C3","C4","D1","D2","D3","D4"))

}
if(!missing(df2)&type=="SEQC"){
    g<-ggplot(subset(a,sample=="C"))
    g+geom_point(aes(x=Apct,y=a$Apct[a$sample=="D"],color=as.factor(LID)),alpha=0.7,size=6)+
      scale_y_continuous(limits=c(0,0.5),expand=c(0,0))+scale_x_continuous(limits=c(0.5,1),expand=c(0,0))+
      facet_wrap(~ LID)+ylab("Amount of SEQC-A in SEQC-C")+xlab("Amount of SEQC-A in SEQC-D")+
      geom_point(aes(x=0.75,y=0.25),col="grey70")+theme(legend.position="none")+geom_path(data=circleFun(center=c(0.75,0.25),diameter=0.1,npoints=25),aes(x,y),col="grey")+geom_path(data=circleFun(center=c(0.75,0.25),diameter=0.05,npoints=25),aes(x,y),col="grey")+geom_path(data=circleFun(center=c(0.75,0.25),diameter=0.15,npoints=25),aes(x,y),col="grey")+geom_path(data=circleFun(center=c(0.75,0.25),diameter=0.2,npoints=25),aes(x,y),col="grey")+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(panel.margin=unit(1,"cm"))+
      theme(legend.text=element_text(size=rel(1.4)))+theme(axis.title=element_text(size=rel(1.6)))+theme(axis.text=element_text(size=rel(1)))+
      theme(strip.background = element_rect(fill = 'white'))+geom_pointrange(data=sumdf,aes(x=amean1,y=amean2,ymax=amean2+sd2,ymin=amean2-sd2))+theme(aspect.ratio=1)+
      geom_errorbarh(data=sumdf,aes(x=amean1,y=amean2,xmax=amean1+sd1,xmin=amean1-sd1))+geom_point(data=subset(b,sample=="C"),aes(x=Apct,y=b$Apct[b$sample=="D"],color="red"))


  }
  if(type=="SEQC"){
    g<-ggplot(subset(a,sample=="C"))
    return(g+geom_point(aes(x=Apct,y=a$Apct[a$sample=="D"],color=as.factor(LID)),alpha=0.7,size=6)+
             scale_y_continuous(limits=c(0,0.5),expand=c(0,0))+scale_x_continuous(limits=c(0.5,1),expand=c(0,0))+
             facet_wrap(~ LID)+ylab("Amount of SEQC-A in SEQC-C")+xlab("Amount of SEQC-A in SEQC-D")+
             geom_point(aes(x=0.75,y=0.25),col="grey70")+theme(legend.position="none")+geom_path(data=circleFun(center=c(0.75,0.25),diameter=0.1,npoints=25),aes(x,y),col="grey")+
             geom_path(data=circleFun(center=c(0.75,0.25),diameter=0.05,npoints=25),aes(x,y),col="grey")+geom_path(data=circleFun(center=c(0.75,0.25),diameter=0.15,npoints=25),aes(x,y),col="grey")+
             geom_path(data=circleFun(center=c(0.75,0.25),diameter=0.2,npoints=25),aes(x,y),col="grey")+
             theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(panel.margin=unit(1,"cm"))+
             theme(legend.text=element_text(size=rel(1.4)))+theme(axis.title=element_text(size=rel(1.6)))+theme(axis.text=element_text(size=rel(1)))+
             theme(strip.background = element_rect(fill = 'white'))+geom_pointrange(data=sumdf,aes(x=amean1,y=amean2,ymax=amean2+sd2,ymin=amean2-sd2))+
             theme(aspect.ratio=1)+
             geom_errorbarh(data=sumdf,aes(x=amean1,y=amean2,xmax=amean1+sd1,xmin=amean1-sd1)))
  }

  if(!missing(df2)&type=="BLM"){
    a$mcor<-1
    b$mcor<-0
    pathdf<-NULL;pathdf2<-NULL;pathdf3<-NULL
    for(I in 1:numrings){pathdf<-rbind(pathdf,circleFun(center=c(0.25,0.25),diameter=0.05*I,npoints=25))
    pathdf2<-rbind(pathdf2,circleFun(center=c(0.25,0.5),diameter=0.05*I,npoints=25))
    pathdf3<-rbind(pathdf3,circleFun(center=c(0.5,0.25),diameter=0.05*I,npoints=25))}
    mergedm<-rbind(a,b)
    merged<-rbind(sumdf,sumdf2)
    merged$mcor<-c(1,1,1,0,0,0)
    g<-ggplot(mergedm)
    return(g+
             geom_path(data=pathdf,aes(x,y),col="grey")+geom_path(data=pathdf2,aes(x,y),col="grey")+geom_path(data=pathdf3,aes(x,y),col="grey")+
             geom_point(aes(x=value1,y=value2,col=variable,alpha=as.factor(mcor)),size=5)+
             xlab("Amount of Tissue in BLM-1")+ylab("Amount of Tissue in BLM-2")+geom_point(aes(x=c(0.25,0.25,.5),y=c(0.25,0.5,0.25)),col="grey70",size=3)+
             theme(legend.position=c(0.6,0.7))+scale_color_manual(name="Tissue",breaks=c("mixpropA","mixpropB","mixpropC","mixpropA","mixpropB","mixpropC"),
            labels=c("Brain","Liver","Muscle","Brain","Liver","Muscle"),values=rep(c("#CC6666","#99CC66","#6699CC"),2))+coord_cartesian(ylim=c(0,1),xlim=c(0,1))+
             theme(axis.text=element_text(size=rel(1.3)))+theme(axis.title=element_text(size=rel(1.5)))+theme(legend.text=element_text(size=rel(1.5)))+theme(legend.title=element_text(size=rel(1.5)))+
             geom_pointrange(data=merged,aes(x=amean1,y=amean2,ymax=amean2+sd2,ymin=amean2-sd2,alpha=as.factor(mcor)),size=1.15)+
             geom_errorbarh(data=merged,aes(x=amean1,y=amean2,xmax=amean1+sd1,xmin=amean1-sd1,height=0,alpha=as.factor(mcor)),size=1.3)+scale_alpha_manual(name=Discriminator,breaks=c(1,0),labels=c("True","False"),values=c(0.3,1))+
             theme(aspect.ratio=1)+theme(axis.ticks.margin=unit(x = 0.25,units = "cm"))+theme(axis.text=element_text(size=16)))
  }
  if(missing(df2)&type=="BLM"){
    a$mcor<-1
    mergedm<-a
    merged<-sumdf
    g<-ggplot(mergedm)
    pathdf<-NULL;pathdf2<-NULL;pathdf3<-NULL
    for(I in 1:numrings){pathdf<-rbind(pathdf,circleFun(center=c(0.25,0.25),diameter=0.05*I,npoints=25))
                         pathdf2<-rbind(pathdf2,circleFun(center=c(0.25,0.5),diameter=0.05*I,npoints=25))
                         pathdf3<-rbind(pathdf3,circleFun(center=c(0.5,0.25),diameter=0.05*I,npoints=25))
      }
    return(g+xlab("Amount of Tissue in BLM-1")+ylab("Amount of Tissue in BLM-2")+geom_path(data=pathdf,aes(x,y),col="grey")+geom_path(data=pathdf2,aes(x,y),col="grey")+geom_path(data=pathdf3,aes(x,y),col="grey")+
             geom_point(aes(x=c(0.25,0.25,.5),y=c(0.25,0.5,0.25)),col="grey70",size=3)+geom_point(aes(x=value1,y=value2,col=variable),size=5)+theme(legend.position=c(0.6,0.7))+
             scale_color_manual(name="Tissue",breaks=c("mixpropA","mixpropB","mixpropC"),labels=c("Brain","Liver","Muscle"),values=rep(c("#CC6666","#99CC66","#6699CC"),1))+
             coord_cartesian(ylim=c(0,1),xlim=c(0,1))+theme(axis.text=element_text(size=rel(1.3)))+theme(axis.title=element_text(size=rel(1.5)))+theme(legend.text=element_text(size=rel(1.5)))+
             theme(legend.title=element_text(size=rel(1.5)))+geom_pointrange(data=merged,aes(x=amean1,y=amean2,ymax=amean2+sd2,ymin=amean2-sd2),size=1.15)+
             geom_errorbarh(data=merged,aes(x=amean1,y=amean2,xmax=amean1+sd1,xmin=amean1-sd1,height=0),size=1.3)+theme(aspect.ratio=1))

  }

}#Generates a target plot discriminating between two input dataframes (from makerefmetdfb) based on discriminator value Discriminator (where df1 is assumed TRUE for this value)
blmmixfraction<-function(subdf){

  value<-unname(lm(subdf[,2]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients/sum(lm(subdf[,2]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients))
  value1<-unname(lm(subdf[,3]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients/sum(lm(subdf[,3]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients))
  value2<-unname(lm(subdf[,4]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients/sum(lm(subdf[,4]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients))
  value3<-unname(lm(subdf[,5]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients/sum(lm(subdf[,5]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients))
  value4<-unname(lm(subdf[,6]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients/sum(lm(subdf[,6]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients))
  value5<-unname(lm(subdf[,7]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients/sum(lm(subdf[,7]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients))
  value6<-unname(lm(subdf[,8]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients/sum(lm(subdf[,8]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients))
  value7<-unname(lm(subdf[,9]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients/sum(lm(subdf[,9]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients))
  value8<-unname(lm(subdf[,10]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients/sum(lm(subdf[,10]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients))
  value9<-unname(lm(subdf[,11]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients/sum(lm(subdf[,11]~subdf[,12]+subdf[,13]+subdf[,14]+0)$coefficients))
  mfrac<-calcmrnafrac(subdf,selection=12:14)
  mfrac<-mfrac/sum(mfrac)
  Amfrac<-mfrac[1]
  Bmfrac<-mfrac[2]
  Cmfrac<-mfrac[3]
  mixpropC<-(value[3]/Cmfrac)/((value[1]/Amfrac)+(value[2]/Bmfrac)+(value[3]/Cmfrac))
  mixpropB<-(value[2]/Bmfrac)/((value[1]/Amfrac)+(value[2]/Bmfrac)+(value[3]/Cmfrac))
  mixpropA<-(value[1]/Amfrac)/((value[1]/Amfrac)+(value[2]/Bmfrac)+(value[3]/Cmfrac))
  outdf<-data.frame(mixpropA,mixpropB,mixpropC)
  #mixprop<-c(mixpropA,mixpropB,mixpropC)
  mixpropC<-(value1[3]/Cmfrac)/((value1[1]/Amfrac)+(value1[2]/Bmfrac)+(value1[3]/Cmfrac))
  mixpropB<-(value1[2]/Bmfrac)/((value1[1]/Amfrac)+(value1[2]/Bmfrac)+(value1[3]/Cmfrac))
  mixpropA<-(value1[1]/Amfrac)/((value1[1]/Amfrac)+(value1[2]/Bmfrac)+(value1[3]/Cmfrac))
  outdf<-rbind(outdf,c(mixpropA,mixpropB,mixpropC))
  mixpropC<-(value2[3]/Cmfrac)/((value2[1]/Amfrac)+(value2[2]/Bmfrac)+(value2[3]/Cmfrac))
  mixpropB<-(value2[2]/Bmfrac)/((value2[1]/Amfrac)+(value2[2]/Bmfrac)+(value2[3]/Cmfrac))
  mixpropA<-(value2[1]/Amfrac)/((value2[1]/Amfrac)+(value2[2]/Bmfrac)+(value2[3]/Cmfrac))
  outdf<-rbind(outdf,c(mixpropA,mixpropB,mixpropC))
  mixpropC<-(value3[3]/Cmfrac)/((value3[1]/Amfrac)+(value3[2]/Bmfrac)+(value3[3]/Cmfrac))
  mixpropB<-(value3[2]/Bmfrac)/((value3[1]/Amfrac)+(value3[2]/Bmfrac)+(value3[3]/Cmfrac))
  mixpropA<-(value3[1]/Amfrac)/((value3[1]/Amfrac)+(value3[2]/Bmfrac)+(value3[3]/Cmfrac))
  outdf<-rbind(outdf,c(mixpropA,mixpropB,mixpropC))
  mixpropC<-(value4[3]/Cmfrac)/((value4[1]/Amfrac)+(value4[2]/Bmfrac)+(value4[3]/Cmfrac))
  mixpropB<-(value4[2]/Bmfrac)/((value4[1]/Amfrac)+(value4[2]/Bmfrac)+(value4[3]/Cmfrac))
  mixpropA<-(value4[1]/Amfrac)/((value4[1]/Amfrac)+(value4[2]/Bmfrac)+(value4[3]/Cmfrac))
  outdf<-rbind(outdf,c(mixpropA,mixpropB,mixpropC))
  mixpropC<-(value5[3]/Cmfrac)/((value5[1]/Amfrac)+(value5[2]/Bmfrac)+(value5[3]/Cmfrac))
  mixpropB<-(value5[2]/Bmfrac)/((value5[1]/Amfrac)+(value5[2]/Bmfrac)+(value5[3]/Cmfrac))
  mixpropA<-(value5[1]/Amfrac)/((value5[1]/Amfrac)+(value5[2]/Bmfrac)+(value5[3]/Cmfrac))
  outdf<-rbind(outdf,c(mixpropA,mixpropB,mixpropC))
  mixpropC<-(value6[3]/Cmfrac)/((value6[1]/Amfrac)+(value6[2]/Bmfrac)+(value6[3]/Cmfrac))
  mixpropB<-(value6[2]/Bmfrac)/((value6[1]/Amfrac)+(value6[2]/Bmfrac)+(value6[3]/Cmfrac))
  mixpropA<-(value6[1]/Amfrac)/((value6[1]/Amfrac)+(value6[2]/Bmfrac)+(value6[3]/Cmfrac))
  outdf<-rbind(outdf,c(mixpropA,mixpropB,mixpropC))
  mixpropC<-(value7[3]/Cmfrac)/((value7[1]/Amfrac)+(value7[2]/Bmfrac)+(value7[3]/Cmfrac))
  mixpropB<-(value7[2]/Bmfrac)/((value7[1]/Amfrac)+(value7[2]/Bmfrac)+(value7[3]/Cmfrac))
  mixpropA<-(value7[1]/Amfrac)/((value7[1]/Amfrac)+(value7[2]/Bmfrac)+(value7[3]/Cmfrac))
  outdf<-rbind(outdf,c(mixpropA,mixpropB,mixpropC))
  mixpropC<-(value8[3]/Cmfrac)/((value8[1]/Amfrac)+(value8[2]/Bmfrac)+(value8[3]/Cmfrac))
  mixpropB<-(value8[2]/Bmfrac)/((value8[1]/Amfrac)+(value8[2]/Bmfrac)+(value8[3]/Cmfrac))
  mixpropA<-(value8[1]/Amfrac)/((value8[1]/Amfrac)+(value8[2]/Bmfrac)+(value8[3]/Cmfrac))
  outdf<-rbind(outdf,c(mixpropA,mixpropB,mixpropC))
  mixpropC<-(value9[3]/Cmfrac)/((value9[1]/Amfrac)+(value9[2]/Bmfrac)+(value9[3]/Cmfrac))
  mixpropB<-(value9[2]/Bmfrac)/((value9[1]/Amfrac)+(value9[2]/Bmfrac)+(value9[3]/Cmfrac))
  mixpropA<-(value9[1]/Amfrac)/((value9[1]/Amfrac)+(value9[2]/Bmfrac)+(value9[3]/Cmfrac))
  outdf<-rbind(outdf,c(mixpropA,mixpropB,mixpropC))




  return(outdf)
}#uses linear models to determine the mix fraction of a dataset
varstabdata<-function(){

  library(DESeq2)

  deseqc<-NULL;deseqc$A<-rowSums(SEQCDF[,c(3,4,5,6)])
  deseqc$B<-rowSums(SEQCDF[,c(7,8,9,10)])
  deseqc$C<-rowSums(SEQCDF[,c(11,12,13,14)])
  deseqc$D<-rowSums(SEQCDF[,c(15,16,17,18)])
  SEQCDFM<-seqcmodel(SEQCDF)#this would probably be better done using a better modeling function - one that handles variable sizes, for example...
  deseqc$Cm<-round(SEQCDFM$modelC,digits=0)
  deseqc$Dm<-round(SEQCDFM$modelD,digits = 0)
  deseqc<-as.data.frame(deseqc)
  deseqctilm<-cbind(deseqc[(1:43919),],deseqc[(43919*1+1):(43919*2),],deseqc[(43919*2+1):(43919*3),],deseqc[(43919*3+1):(43919*4),],deseqc[(43919*4+1):(43919*5),],deseqc[(43919*5+1):(43919*6),])
  coldata<-data.frame(site=c("Agr","Agr","Agr","Agr","Agr","Agr","Bgi","Bgi","Bgi","Bgi","Bgi","Bgi","Cnl","Cnl","Cnl","Cnl","Cnl","Cnl","Coh","Coh","Coh","Coh","Coh","Coh","May","May","May","May","May","May","Nvs","Nvs","Nvs","Nvs","Nvs","Nvs"),sample=c("A","B","C","D","Cm","Dm"))
  rownames(coldata)<-c("AgrA","AgrB","AgrC","AgrD","AgrCm","AgrDm","BgiA","BgiB","BgiC","BgiD","BgiCm","BgiDm","CnlA","CnlB","CnlC","CnlD","CnlCm","CnlDm","CohA","CohB","CohC","CohD","CohCm","CohDm","MayA","MayB","MayC","MayD","MayCm","MayDm","NvsA","NvsB","NvsC","NvsD","NvsCm","NvsDm")
  fullctab<-DESeqDataSetFromMatrix(countData=deseqctilm,colData=coldata,design= sample ~ site)
  fullctab<-rlog(fullctab,blind=FALSE)
  return(fullctab)
} #takes the relatively ugly SEQCDF and converts it into something that DESeq can handle.  Also rlogs because why not apply a variance stabilizing transform that does great things?
makeseqcdf<-function(SEQCcountpath="SEQCtranscript/"){
  fidf<-NULL
  fdf<-NULL
  for ( I in c("AGR","BGI","CNL","COH","MAY","NVS")){
    string<-paste(SEQCcountpath,"SEQC_MAIN_ILM_",I,"_TranscriptCounts_ZSU.txt",sep="")
    tdf<-fread(string)
    tdf<-as.data.frame(tdf)
    fdf$gene_id<-tdf[,1]
    fdf$site<-c(rep(I,length(tdf[,3])))
    fdf$A1<-rowSums(tdf[,grep("_A_1_",colnames(tdf))])
    fdf$A2<-rowSums(tdf[,grep("_A_2_",colnames(tdf))])
    fdf$A3<-rowSums(tdf[,grep("_A_3_",colnames(tdf))])
    fdf$A4<-rowSums(tdf[,grep("_A_4_",colnames(tdf))])
    #fdf$A5<-rowSums(tdf[,grep("_A_5_",colnames(tdf))])
    fdf$B1<-rowSums(tdf[,grep("_B_1_",colnames(tdf))])
    fdf$B2<-rowSums(tdf[,grep("_B_2_",colnames(tdf))])
    fdf$B3<-rowSums(tdf[,grep("_B_3_",colnames(tdf))])
    fdf$B4<-rowSums(tdf[,grep("_B_4_",colnames(tdf))])
    #fdf$B5<-rowSums(tdf[,grep("_B_5_",colnames(tdf))])
    fdf$C1<-rowSums(tdf[,grep("_C_1_",colnames(tdf))])
    fdf$C2<-rowSums(tdf[,grep("_C_2_",colnames(tdf))])
    fdf$C3<-rowSums(tdf[,grep("_C_3_",colnames(tdf))])
    fdf$C4<-rowSums(tdf[,grep("_C_4_",colnames(tdf))])
    #fdf$C5<-rowSums(tdf[,grep("_C_5_",colnames(tdf))])
    fdf$D1<-rowSums(tdf[,grep("_D_1_",colnames(tdf))])
    fdf$D2<-rowSums(tdf[,grep("_D_2_",colnames(tdf))])
    fdf$D3<-rowSums(tdf[,grep("_D_3_",colnames(tdf))])
    fdf$D4<-rowSums(tdf[,grep("_D_4_",colnames(tdf))])
    #fdf$D5<-rowSums(tdf[,grep("_D_5_",colnames(tdf))])
    #fdf$E1<-rowSums(tdf[,grep("_E_1_",colnames(tdf))])
    #fdf$E2<-rowSums(tdf[,grep("_E_2_",colnames(tdf))])
    #fdf$E3<-rowSums(tdf[,grep("_E_3_",colnames(tdf))])
    #fdf$E4<-rowSums(tdf[,grep("_E_4_",colnames(tdf))])
    fdf<-as.data.frame(fdf)
    fidf<-rbind(fidf,fdf)
  }#adds in the ILM sites
  for ( I in c("NWU","PSU","SQW")){
    string<-paste(SEQCcountpath,~/Desktop/MainSEQC/SEQC_LifescopeCounts/"SEQC_LIF_",I,"_LifeScope.csv",sep="")
    fdf<-NULL
    tdf<-fread(string)
    tdf<-as.data.frame(tdf)
    fdf$gene_id<-tdf[,1]
    fdf$site<-c(rep(I,length(tdf[,3])))
    fdf$A1<-rowSums(tdf[,grep("_A_1_",colnames(tdf))])
    fdf$A2<-rowSums(tdf[,grep("_A_2_",colnames(tdf))])
    fdf$A3<-rowSums(tdf[,grep("_A_3_",colnames(tdf))])
    fdf$A4<-rowSums(tdf[,grep("_A_4_",colnames(tdf))])
    #fdf$A5<-rowSums(tdf[,grep("_A_5_",colnames(tdf))])
    fdf$B1<-rowSums(tdf[,grep("_B_1_",colnames(tdf))])
    fdf$B2<-rowSums(tdf[,grep("_B_2_",colnames(tdf))])
    fdf$B3<-rowSums(tdf[,grep("_B_3_",colnames(tdf))])
    fdf$B4<-rowSums(tdf[,grep("_B_4_",colnames(tdf))])
    #fdf$B5<-rowSums(tdf[,grep("_B_5_",colnames(tdf))])
    fdf$C1<-rowSums(tdf[,grep("_C_1_",colnames(tdf))])
    fdf$C2<-rowSums(tdf[,grep("_C_2_",colnames(tdf))])
    fdf$C3<-rowSums(tdf[,grep("_C_3_",colnames(tdf))])
    fdf$C4<-rowSums(tdf[,grep("_C_4_",colnames(tdf))])
    #fdf$C5<-rowSums(tdf[,grep("_C_5_",colnames(tdf))])
    fdf$D1<-rowSums(tdf[,grep("_D_1_",colnames(tdf))])
    fdf$D2<-rowSums(tdf[,grep("_D_2_",colnames(tdf))])
    fdf$D3<-rowSums(tdf[,grep("_D_3_",colnames(tdf))])
    fdf$D4<-rowSums(tdf[,grep("_D_4_",colnames(tdf))])
    #fdf$D5<-rowSums(tdf[,grep("_D_5_",colnames(tdf))])
    #fdf$E1<-rowSums(tdf[,grep("_E_1_",colnames(tdf))])
    #fdf$E2<-rowSums(tdf[,grep("_E_2_",colnames(tdf))])
    #fdf$E3<-rowSums(tdf[,grep("_E_3_",colnames(tdf))])
    #fdf$E4<-rowSums(tdf[,grep("_E_4_",colnames(tdf))])
    fdf<-as.data.frame(fdf)
    fidf<-rbind(fidf,fdf)
  } #adds in the SEQC sites
  return(fidf)
}#collates the many raw data files from MainSEQC into a dataframe (note that the source data needs to be moved to an online location and this function changed alongside it)
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}# function to get color labels
assignIdentity<-function(dat,index=5,forpub=1,minct=20){
  if(forpub==0){dat$identity<-"n"
                dat$identity[dat$bep>index*dat$lep&dat$bep>index*dat$mep&dat$bep>minct]<-"b"
                dat$identity[dat$lep>index*dat$bep&dat$lep>index*dat$mep&dat$lep>minct]<-"l"
                dat$identity[dat$mep>index*dat$bep&dat$mep>index*dat$lep&dat$mep>minct]<-"m"
                dat$identity[grep("ERCC-",dat$gene_id)]<-"x"
                if(sum(dat$identity=="x")==0){dat$identity[match(ercc96$V9,dat$gene_id)]<-"x"}
                #dat$varrat<-(dat$a1-dat$modeled+1)/(dat$modeled+1)
  }
  else if(forpub==1){
    dat$identity<-"unclassified"
    dat$identity[dat$bep>index*dat$lep&dat$bep>index*dat$mep&dat$bep>minct]<-"Brain"
    dat$identity[dat$lep>index*dat$bep&dat$lep>index*dat$mep&dat$lep>minct]<-"Liver"
    dat$identity[dat$mep>index*dat$bep&dat$mep>index*dat$lep&dat$mep>minct]<-"Muscle"
    dat$identity[grep("ERCC-",dat$gene_id)]<-"external"
    if(sum(dat$identity=="x")==0){dat$identity[match(ercc96$V9,dat$gene_id)]<-"external"}

  }
  return(dat)}#subsets genes within dataframe dat based on Bep Lep and Mep being > index * others (tissue-selective: >5x counts) and above a minimum threshold minct
calcmrnafrac<-function(dat,selection=2:14,e5=FALSE,type,ret=0){
    ety<-1
    if(missing(type)){type=1}
    ercc<-rownames(dat)[substr(rownames(dat),1,5)=="ERCC-"]
    if(type!=1){ercc<-rownames(dat)[match(ercc96$V1[ercc96$V4=="C"],rownames(dat))]}
    if(type=="f"){ercc<-match(ercc96$uname[ercc96$pool=="C"],dat$gene_id);ercc<-ercc[!is.na(ercc)]}
    if(length(ercc)==0){ety<-2;ercc<-match(ercc96$V9,dat$gene_id);ercc<-ercc[!is.na(ercc)]}
    if(length(ercc)==0){ety<-1;ercc<-grep("ERCC-",dat$gene_id)}
    if(length(ercc)==0){ety<-3;ercc<-grep("ERCC_",dat$gene_id)}
    if(length(ercc)==0){return("Sorry, ERROR14 LMAO YOU HAVE NO IDEA WHAT ERROR 14 IS AND NEITHER DO I")}
    non.ercc<-rownames(dat)[!rownames(dat)%in%ercc]

    if(ety<3){
      count<-rbind(colSums(dat[ercc,selection]),colSums(dat[non.ercc,selection]))
      if(length(selection)==13){ercc.targ<-c(rep(.05*.15,10),rep(.03*.1,3))}
      else{ercc.targ<-c(rep(0.08,length(selection)))}

    }
    else if(ety==3){
      count<-rbind(colSums(dat[ercc,selection]),colSums(dat[non.ercc,selection]))
      ercc.targ<-c(rep(.08,length(colnames(count))))}
    ercc.targ[grep("d",colnames(count))]<-ercc.targ[grep("d",colnames(count))]/8
    ercc.targ[grep("u",colnames(count))]<-ercc.targ[grep("u",colnames(count))]*8
    mRNA.frac<-ercc.targ*count[2,]/count[1,]
    #if(ety==3){ #this whole block is causing me problems and i don't THINK it's necessary anymore
      #wish to exclude sample5?  The flag is e5.
    #  if(e5==FALSE){
    #    mRNA.A<-mRNA.frac[grep("A",names(mRNA.frac))]
    #    mRNA.A<-mean(mRNA.A,na.rm=TRUE)
    #    mRNA.B<-mRNA.frac[grep("B",names(mRNA.frac))]
    #    mRNA.B<-mean(mRNA.B,na.rm=TRUE)
    #    mRNA.Anorm<-mRNA.A/(mRNA.A+mRNA.B)
    #    mRNA.Bnorm<-mRNA.B/(mRNA.B+mRNA.A)
    #    mRNA.frac<-c(mRNA.Anorm,mRNA.Bnorm)
    #  }
    #  if(e5==TRUE){
    #    mRNA.A<-mRNA.frac[grep("A[1-4]",names(mRNA.frac))]
    #    mRNA.B<-mRNA.frac[grep("B[1-4]",names(mRNA.frac))]
    #    mRNA.A<-mean(mRNA.A,na.rm=TRUE)
    #    mRNA.B<-mean(mRNA.B,na.rm=TRUE)
    #    mRNA.Anorm<-mRNA.A/(mRNA.A+mRNA.B)
    #    mRNA.Bnorm<-mRNA.B/(mRNA.B+mRNA.A)
    #    mRNA.frac<-c(mRNA.Anorm,mRNA.Bnorm)
    #  }
    #}
    if(ret!=0){return(count)}
    return(mRNA.frac)
  } #calculates the mRNA fraction of a makerefmetdfb-created dataframe (calcmrnafrageneral is more accomomdating of other input files)
calcmrnafracgeneral<-function(dat,spikeID="ERCC-",spikemassfraction=.1){
  #1) Identify which counts are Spike-In and which are not
#0) Identify which columns are counts and which are 'annotation':
  countcolumns<-which(unname(unlist(lapply(dat,class))=="numeric"))
  annotcolumn<-which(unname(unlist(lapply(dat,class))!="numeric"))
  ercc<-rownames(dat)[which(substr(rownames(dat),1,5)==spikeID)]           #one way to identify spikes, if row names = spikeID
  if(length(ercc)==0){ercc<-grep(spikeID,dat[,annotcolumn[1]])} #assuming that the name is in the first annotation column...
  if(length(ercc)==0){stop("I can't identify the spike-ins within your count table.   The spikeID variable should be set to something which uniquely identifies spike-ins.   Rownames are first checked for names, then if there are non-numeric columns, only the FIRST is checked for gene names. ")}
nonercc<-!(1:length(dat[,countcolumns[1]]))%in%ercc

  count<-rbind(colSums(dat[ercc,countcolumns]),colSums(dat[nonercc,countcolumns])) #determines the counts for spikes and non-spikes.
  ercc.targ<-spikemassfraction  #defines the "targeted" mass fraction for spikes : Either a vector with length = #columns,or a scalar
  mRNA.frac<-ercc.targ*count[2,]/count[1,]  #calculates an mRNA fraction based on those available data
  #this part doesn't normalize to one, but that's not exactly complicated.
  return(mRNA.frac)
}
normalizeSEQC<-function(stuff=SEQCDF,type="UQN"){
  fmol<-NULL

  if(type=="UQN"){for(I in levels(as.factor(stuff$site))){
  require(edgeR)
    tmp<-subset(stuff,site==I)
    nfac<-calcNormFactors(tmp[3:18],method="upperquartile")
    for(I in 3:18){tmp[,I]<-tmp[,I]*nfac[(I-2)]} #does the normalization

    fmol<-rbind(fmol,tmp)
  }

  return((fmol))
  }
  else if(type=="LSN"){
    for(I in levels(as.factor(stuff$site))){
      tmp<-subset(stuff,site==I)
      mfac<-max(colSums(tmp[3:18]))
      for(I in 3:18){tmp[,I]<-tmp[,I]*mfac/sum(tmp[,I])}
      fmol<-rbind(fmol,tmp)


    }
    return((fmol))

  }
  else if(type=="Both"){
    fmox<-NULL
    for(I in levels(as.factor(stuff$site))){
      tmp<-subset(stuff,site==I)
      mfac<-max(colSums(tmp[3:18]))
      for(I in 3:18){tmp[,I]<-tmp[,I]*mfac/sum(tmp[,I])}
      fmol<-rbind(fmol,tmp)
    }
    for(I in levels(as.factor(fmol$site))){
      tmp<-subset(fmol,site==I)
      nfac<-calcNormFactors(tmp[3:18],method="upperquartile")
      for(I in 3:18){tmp[,I]<-tmp[,I]*nfac[(I-2)]} #does the normalization

      fmox<-rbind(fmox,tmp)
    }
    return((fmox))
  }
}#normalizes SEQC dataframes by various methods.
SEQClm<-function(infile,e5=TRUE,e4=FALSE,ignoremrna=FALSE,filterercc=FALSE){
  if(e5==FALSE&e4==FALSE){
    infile$Amean<-rowMeans(infile[,grep("A",names(infile))])
    infile$Bmean<-rowMeans(infile[,grep("B",names(infile))])
    mfrac<-calcmrnafrac(infile,selection=c(3:14))
    #mfrac<-calcmrnafrac(infile,selection=c(3:14),type="f") #this does the calculation using "only" 1:1 subpool.  It doesn't affect anything much.
    c1<-coefficients(lm(data=infile,C1 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,C1 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    c2<-coefficients(lm(data=infile,C2 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,C2 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    c3<-coefficients(lm(data=infile,C3 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,C3 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    c4<-coefficients(lm(data=infile,C4 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,C4 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    c5<-coefficients(lm(data=infile,C5 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,C5 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    d1<-coefficients(lm(data=infile,D1 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,D1 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    d2<-coefficients(lm(data=infile,D2 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,D2 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    d3<-coefficients(lm(data=infile,D3 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,D3 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    d4<-coefficients(lm(data=infile,D4 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,D4 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    d5<-coefficients(lm(data=infile,D5 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,D5 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))

    rdf<-rbind(c1,c2,c3,c4,c5,d1,d2,d3,d4,d5)
    rdf<-data.frame(rdf,c("C1","C2","C3","C4","C5","D1","D2","D3","D4","D5"))
    names(rdf)<-c("Aest","Best","library")
  }
  else if(e5==TRUE&e4==FALSE){
    infile$Amean<-rowMeans(infile[,grep("A",names(infile))])
    infile$Bmean<-rowMeans(infile[,grep("B",names(infile))])
    mfrac<-c(1,1)
#   if(ignoremrna==FALSE){mfrac<-calcmrnafrac(infile,selection=c(3:14),type="f")}#this does the calculation using "only" 1:1 subpool.  It doesn't affect anything substantially.
    if(ignoremrna==FALSE){mfrac<-calcmrnafrac(infile,selection=c(3:14))}


    c1<-coefficients(lm(data=infile,C1 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,C1 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    c2<-coefficients(lm(data=infile,C2 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,C2 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    c3<-coefficients(lm(data=infile,C3 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,C3 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    c4<-coefficients(lm(data=infile,C4 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,C4 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    d1<-coefficients(lm(data=infile,D1 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,D1 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    d2<-coefficients(lm(data=infile,D2 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,D2 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    d3<-coefficients(lm(data=infile,D3 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,D3 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))
    d4<-coefficients(lm(data=infile,D4 ~ I(Amean*mfrac[1])+I(Bmean*mfrac[2])+0))/sum(coefficients(lm(data=infile,D4 ~ I(Amean*mfrac[1]) + I (Bmean*mfrac[2])+0)))

    rdf<-rbind(c1,c2,c3,c4,d1,d2,d3,d4)
    rdf<-data.frame(rdf,c("C1","C2","C3","C4","D1","D2","D3","D4"))
    names(rdf)<-c("Aest","Best","library")
    rdf$sample<-substr(rdf$library,1,1)
    rdf$Best<-NULL
    rdf<-data.frame(rdf[rdf$sample=="C",],AinD=rdf$Aest[rdf$sample=="D"])
    rdf$sample<-NULL
    rdf$library<-NULL
  }
  else if(e4==TRUE){
    if(ignoremrna==TRUE){
      infile$Amean<-rowMeans(infile[,grep("A",names(infile))])
      infile$Bmean<-rowMeans(infile[,grep("B",names(infile))])
      c1<-coefficients(lm(data=infile,C1 ~ I(Amean)+I(Bmean)+0))/sum(coefficients(lm(data=infile,C1 ~ I(Amean) + I (Bmean)+0)))
      c2<-coefficients(lm(data=infile,C2 ~ I(Amean)+I(Bmean)+0))/sum(coefficients(lm(data=infile,C2 ~ I(Amean) + I (Bmean)+0)))
      c3<-coefficients(lm(data=infile,C3 ~ I(Amean)+I(Bmean)+0))/sum(coefficients(lm(data=infile,C3 ~ I(Amean) + I (Bmean)+0)))

      d1<-coefficients(lm(data=infile,D1 ~ I(Amean)+I(Bmean)+0))/sum(coefficients(lm(data=infile,D1 ~ I(Amean) + I (Bmean)+0)))
      d2<-coefficients(lm(data=infile,D2 ~ I(Amean)+I(Bmean)+0))/sum(coefficients(lm(data=infile,D2 ~ I(Amean) + I (Bmean)+0)))
      d3<-coefficients(lm(data=infile,D3 ~ I(Amean)+I(Bmean)+0))/sum(coefficients(lm(data=infile,D3 ~ I(Amean) + I (Bmean)+0)))
      rdf<-rbind(c1,c2,c3,d1,d2,d3)
      rdf<-data.frame(rdf,c("C1","C2","C3","D1","D2","D3"))
      names(rdf)<-c("Aest","Best","library")
      rdf$sample<-substr(rdf$library,1,1)
      rdf$Best<-NULL
      rdf<-data.frame(rdf[rdf$sample=="C",],AinD=rdf$Aest[rdf$sample=="D"])
      rdf$sample<-NULL
      rdf$library<-NULL
      return(rdf)
    }
    infile$Amean<-rowMeans(infile[,grep("A",names(infile))])
    infile$Bmean<-rowMeans(infile[,grep("B",names(infile))])
    mfrac<-calcmrnafrac(infile,selection=c(2:15))
    c1<-coefficients(lm(data=infile,C1 ~ I(Amean*mfrac[13])+I(Bmean*mfrac[14])+0))/sum(coefficients(lm(data=infile,C1 ~ I(Amean*mfrac[13]) + I (Bmean*mfrac[14])+0)))
    c2<-coefficients(lm(data=infile,C2 ~ I(Amean*mfrac[13])+I(Bmean*mfrac[14])+0))/sum(coefficients(lm(data=infile,C2 ~ I(Amean*mfrac[13]) + I (Bmean*mfrac[14])+0)))
    c3<-coefficients(lm(data=infile,C3 ~ I(Amean*mfrac[13])+I(Bmean*mfrac[14])+0))/sum(coefficients(lm(data=infile,C3 ~ I(Amean*mfrac[13]) + I (Bmean*mfrac[14])+0)))

    d1<-coefficients(lm(data=infile,D1 ~ I(Amean*mfrac[13])+I(Bmean*mfrac[14])+0))/sum(coefficients(lm(data=infile,D1 ~ I(Amean*mfrac[13]) + I (Bmean*mfrac[14])+0)))
    d2<-coefficients(lm(data=infile,D2 ~ I(Amean*mfrac[13])+I(Bmean*mfrac[14])+0))/sum(coefficients(lm(data=infile,D2 ~ I(Amean*mfrac[13]) + I (Bmean*mfrac[14])+0)))
    d3<-coefficients(lm(data=infile,D3 ~ I(Amean*mfrac[13])+I(Bmean*mfrac[14])+0))/sum(coefficients(lm(data=infile,D3 ~ I(Amean*mfrac[13]) + I (Bmean*mfrac[14])+0)))

    rdf<-rbind(c1,c2,c3,d1,d2,d3)
    rdf<-data.frame(rdf,c("C1","C2","C3","D1","D2","D3"))
    names(rdf)<-c("Aest","Best","library")
    rdf$sample<-substr(rdf$library,1,1)
    rdf$Best<-NULL
    rdf<-data.frame(rdf[rdf$sample=="C",],AinD=rdf$Aest[rdf$sample=="D"])
    rdf$sample<-NULL
    rdf$library<-NULL
  }
  names(rdf)<-c("value","Dval")
  return(rdf)
}
deseqcdfsite<-function(){
  odf<-NULL
  for(I in levels(SEQCDFN$site)){
    tdf<-onesitedeseq(stf=I,rpoint=2)
    tdf$site<-I
    odf<-rbind(odf,as.data.frame(tdf))


  }
  return(odf)
}
onesitedeseq<-function(INDF=SEQCDF,stf="AGR",rpoint=1){

  library(DESeq2)
  #built deseqc dataframe out of SEQCDF:
    juice<-subset(INDF,site==stf)
  juice$site<-droplevels(juice$site)
  juice<-juice[rowSums(juice[,c(3:18)])>0,]
  SEQCDFM<-seqcmodel(juice)
  deseqc<-NULL
  deseqc$C1<-juice[,11]
  deseqc$C2<-juice[,12]
  deseqc$C3<-juice[,13]
  deseqc$C4<-juice[,14]
  deseqc$D1<-juice[,15]
  deseqc$D2<-juice[,16]
  deseqc$D3<-juice[,17]
  deseqc$D4<-juice[,18]
  deseqc$Cm1<-round(SEQCDFM$modelC,digits=0)
  deseqc$Cm2<-deseqc$Cm1
  deseqc$Cm3<-deseqc$Cm1
  deseqc$Cm4<-deseqc$Cm1
  deseqc$Dm1<-round(SEQCDFM$modelD,digits = 0)
  deseqc$Dm2<-deseqc$Dm1
  deseqc$Dm3<-deseqc$Dm1
  deseqc$Dm4<-deseqc$Dm1
  deseqc<-as.data.frame(deseqc)
  coldata<-data.frame(site=stf,sample=c(rep("C",4),rep("D",4),rep("Cm",4),rep("Dm",4)),replicate=c(1,2,3,4))
   fullctab<-DESeqDataSetFromMatrix(countData=deseqc,colData=coldata,design= ~ sample)
  if(rpoint==1){  return(fullctab)}
  fullctab<-DESeq(fullctab)
  fullctab<-results(fullctab,contrast=c("sample","C","Cm"))
  return(fullctab)
}
seqcmodel<-function(DF=SEQCDF){
  library(matrixStats)
  out<-NULL
  for(I in levels(DF$site)){
    tmp<-subset(DF,site==I)
    tmp$repA<-rowMeans(as.matrix(tmp[,c(3:6)]))
    tmp$sda<-rowSds(as.matrix(tmp[,c(3:6)]))
    tmp$repB<-rowMeans(as.matrix(tmp[,c(7:10)]))
    tmp$sdb<-rowSds(as.matrix(tmp[,c(7:10)]))
    tmp$repC<-rowMeans(as.matrix(tmp[,c(11:14)]))
    tmp$sdc<-rowSds(as.matrix(tmp[,c(11:14)]))
    tmp$repD<-rowMeans(as.matrix(tmp[,c(15:18)]))
    tmp$sdd<-rowSds(as.matrix(tmp[,c(15:18)]))
    mfrac<-calcmrnafrac(dat=tmp,selection=c(19,21))
    tmp$modelC<-(tmp$repA*mfrac[1]*.75)+(tmp$repB*mfrac[2]*.25)
    tmp$modelD<-(tmp$repA*mfrac[1]*.25)+(tmp$repB*mfrac[2]*.75)
    tmp$modelC<-tmp$modelC/sum(tmp$modelC)*sum(tmp$repC)
    tmp$modelD<-tmp$modelD/sum(tmp$modelD)*sum(tmp$repD)

    out<-rbind(out,tmp)
  }
  return(out)
}
GeneralLMest<-function(infile,spikeID="ERCC-",components=c("bep","lep","mep"),mixes=c("a1"),spikemassfraction=0.003,replicates=1){
  ##formatting of infile:
  # Rownames:  gene_ids, with a phrase ("spikeID") that identifies Spike-In controls
  # Columns: 1 column per component, 1 column per mix
  #must input list of components and list of mixes: For example:  the BLM mix contains components= c("bep","lep","mep") mixes=c("a1","a2")
  #spikemassfraction is the mass proportion you targeted with your spike-in controls
  #input should already be averaged across replicates unless i want to try to find a way to guess...
{replicatebehavior="basic" #I don't know what else to do with replicates right now other than average them, but i imagine i'll come up with something more interesting later.
  if(replicates==1){
  if(length(grep(paste0("^",substr(components[1],1,1)),infile))>2){message("Assuming that components with similar starting names are replicates.  Turn this behavior off with replicates=0")}
  }
  if(replicatebehavior=="basic"){
    #find all replicates and turn them into means
    RepMeans<-NULL
    for(I in 1:length(infile)){
    firstreplicates<-grep(paste0("^",substr(components[I],1,1)),components)
    assign(paste0("RepMeans$",components[I]),rowMeans(infile[,components[firstreplicates]]))

    }
  }}#testphase code trying to handle SEQCDFs and arbitrary ones with replicates.  Because I really need to spend a day making this code extensible instead of just remaking the figure in a more ad-hoc way...
  mfrac<-calcmrnafracgeneral(infile,spikeID,spikemassfraction)
  if(length(mfrac[components])!=length(components)){return("ERROR:  Number of components != Number of calculated mRNA fractions!")}
  mixval<-NULL;rdf<-NULL
  #convincing LM to handle arbitrary columns is not simple: The following block tries to do that.
{
  for(mixtext in mixes){
    modeltext<-NULL
    modeltext<-paste(modeltext,(paste("I(",mfrac[1],"*",components[1],")+0")))
    for(I in 2:length(components)){
      modeltext<-paste(modeltext,(paste("+I(",mfrac[I],"*", components[I],")+0")),sep=" ")
    }
    mixval<-coefficients(lm(data=infile, as.formula(paste(mixtext,"~",modeltext))))/sum(coefficients(lm(data=infile,as.formula(paste(mixtext,"~",modeltext)))))
    rdf<-rbind(rdf,mixval)
  }
  rownames(rdf)<-mixes
  colnames(rdf)<-components
  return(rdf)
}
}
mfdb<-function(df){
  outfrac<-NULL
  for(I in levels(as.factor(df$site))){
    tmp<-subset(df,site==I)
    for(J in 3:18){
      tmp[,J]<-tmp[,J]/sum(tmp[,J])*4e8
    }#normalization of various sites/librariesn
    tmp$repA<-rowMeans(tmp[,3:6])
    tmp$repB<-rowMeans(tmp[,7:10])
    tmp$repC<-rowMeans(tmp[,11:14])
    tmp$repD<-rowMeans(tmp[,15:18])
    #means
    mfrac<-calcmrnafrac(tmp,selection = c(19:22))
    outfrac<-rbind(outfrac,mfrac)
  }
  outfrac2<-outfrac
  outfrac2[,1]<-outfrac[,1]/(outfrac[,1]+outfrac[,2])
  outfrac2[,2]<-outfrac[,2]/(outfrac[,1]+outfrac[,2])
  outfrac2[,3]<-outfrac[,3]/(outfrac[,1]+outfrac[,2])
  outfrac2[,4]<-outfrac[,4]/(outfrac[,1]+outfrac[,2])
  #this normalizes the mRNA content - making the output effectively the A:B ratio
  return(outfrac)}

normcountdf<-function(tdf,normtype){
    if(normtype==0){
    return(tdf)}
  if(normtype==1){
    tdf$norm<-"LSN"
    maxct<-max(colSums(tdf[,2:11]))
    #    sums<-colSums(tdf[,2:11])
    #    tdf$sums<-1
    #    tdf$sums[1:10]<-sums
    for(I in 2:11){
      tdf[,I]<-tdf[,I]*maxct/sum(tdf[,I])
    }
    return(tdf)
  }
  if(normtype=="UQN"){
    library(edgeR)
    tdf$norm<-"UQN"
    nfac<-c(calcNormFactors(tdf[2:11],method="upperquartile"),calcNormFactors(tdf[12:14],method="upperquartile"))
    #    nfac<-c(nfac,calcNormFactors(tdf[12:14]))
    for(I in 2:14){tdf[I]<-tdf[I]*nfac[I-1]}
    return(tdf)
  }
  if(normtype=="TMM"){
    tdf$norm<-"TMM"
    library(edgeR)
    nfac<-c(calcNormFactors(tdf[2:11]),calcNormFactors(tdf[12:14]))
    #    nfac<-c(nfac,calcNormFactors(tdf[12:14]))
    for(I in 2:14){tdf[I]<-tdf[I]*nfac[I-1]}
    return(tdf)
  }
  if(normtype=="RLE"){
    tdf$norm<-"RLE"
    library(edgeR)
    nfac<-c(calcNormFactors(tdf[2:11],method="RLE"),calcNormFactors(tdf[12:14],method="RLE"))
    #    nfac<-c(nfac,calcNormFactors(tdf[12:14]))
    for(I in 2:14){tdf[I]<-tdf[I]*nfac[I-1]}
    return(tdf)}
  if(normtype=="RLENE"){
    tdf$norm<-"RLENE"
    library(edgeR)
    nfac<-c(calcNormFactors(tdf[grep("ERCC-",tdf$gene_id),2:11],method="RLE"),calcNormFactors(tdf[grep("ERCC-",tdf$gene_id),12:14],method="RLE"))
    #    nfac<-c(nfac,calcNormFactors(tdf[12:14]))
    for(I in 2:14){tdf[I]<-tdf[I]*nfac[I-1]}
    return(tdf)}
  return(0)}

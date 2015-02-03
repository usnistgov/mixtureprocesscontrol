#More "generalized" code to solve for unknown mixture proportions. 
#first  Get the mrna fraction
#currently works for both SEQC and BLM data!!! and arbitrary numbers of mixes (given a little massaging to get the input format right)
GeneralLMest<-function(infile,spikeID="ERCC-",components=c("bep","lep","mep"),mixes=c("a1"),spikemassfraction=c(0.003,0.003,0.003)){
##formatting of infile:
  # Rownames:  gene_ids, with a phrase ("spikeID") that identifies Spike-In controls
  # Columns: 1 column per component, 1 column per mix
  #must input list of components and list of mixes: For example:  the BLM mix contains components= c("bep","lep","mep") mixes=c("a1","a2")
  #spikemassfraction is the mass proportion you targeted with your spike-in controls
  #input should already be averaged across replicates

mfrac<-calcmrnafracgeneral(infile[,components],spikeID,spikemassfraction)  
#Mfrac goes into calcmrnafracgeneral, which calculates the mRNA fraction of each sample in the table using equation 3 - here we send it only the "components" of the mixture.
if(length(mfrac)!=length(components)){stop("ERROR:  Number of components != Number of calculated mRNA fractions!")}
{mixval<-NULL;rdf<-NULL}#set up dummy variables
#convincing LM to handle arbitrary columns is not simple: The following block tries to do that.
{
 for(mixtext in mixes){
   modeltext<-NULL
   modeltext<-paste(modeltext,(paste("I(",mfrac[1],"*",components[1],")+0")))
for(I in 2:length(components)){
  modeltext<-paste(modeltext,(paste("+I(",mfrac[I],"*", components[I],")+0")),sep=" ")
}
mixval<-coefficients(lm(data=infile, as.formula(paste(mixtext,"~",modeltext))))/sum(coefficients(lm(data=infile,as.formula(paste(mixtext,"~",modeltext)))))
#Mixval for a mixture sample in "mixes" is here determined to be the coefficients of a linear model fit to the Identity "I" of the mRNA fraction "mfrac" value of a component in "components"
#times the mixture expression of that component, plus the same value for all other components.  This is equation 2.  There is no intercept fit (+0) for this model.

rdf<-rbind(rdf,mixval)


 }
rownames(rdf)<-mixes
colnames(rdf)<-components
return(rdf)
}
}
#calculates the mRNA fraction of a given count table and returns the normalized component values.  
#thisfileneedstobeaddedtothegitrepo##
calcmrnafracgeneral<-function(dat,spikeID,spikemassfraction=0.05){
  #1) Identify which counts are Spike-In and which are not
ercc<-rownames(dat)[substr(rownames(dat),1,5)==spikeID]           #one way to identify spikes, if row names = spikeID
    if(length(ercc)==0){return("I can't identify the spike-ins within your count table.   The spikeID variable should be set to something which uniquely identifies spike-ins and rownames = gene_ids.")}
    non.ercc<-rownames(dat)[!rownames(dat)%in%ercc]
    count<-rbind(colSums(dat[ercc,]),colSums(dat[non.ercc,])) #determines the counts for spikes and non-spikes.
    ercc.targ<-spikemassfraction  #defines the "targeted" mass fraction for spikes : Either a vector with length = #columns,or a scalar
    mRNA.frac<-ercc.targ*count[2,]/count[1,]  #calculates an mRNA fraction based on those available data
    #this part doesn't normalize to one, but that's not exactly complicated.
    return(mRNA.frac)
  }
## 1.0.7 fixed bug with ComBat correction and normalize = F (line 1071).
## 1.0.6 fixed bug with variance calculation on normalized expression data (variance=apply(datExprNorm[,c((skip1+1):dim(datExprNorm)[2])],1,var,na.rm=T)
## 1.0.5 fixed bug with excludevec2 when normalize=F.
## 1.0.4 added an option to exclude genes with 0 variance in ANY group or over ALL groups of samples. Also added code before calling ModulePrinComps1 to ensure that any gene with 0 variance within the sample group was excluded prior to calculations.
## 1.0.3 adding a check to ensure there are no features with 0 variance prior to running ComBat that have been introduced by quantile normalization.
## 1.0.2 will check for missing data (NAs) and give the user the option of setting NA = 1 (recommended for RNAseq data).
## 1.0.2 will also export a list of features (probes, genes, etc.) with 0 variance as a .csv file.
## 1.0.2 can also be used locally or on cluster.
## 1.0.1 includes minor changes to ensure cluster compatibility (time stamp, .libPaths).
## This document contains R functions for constructing sample networks from genomic (or other) datasets.
## These functions are designed to perform sample outlier detection and removal, data normalization, and correction for batch effects.
## A tutorial describing the usage of these functions is available on our web site: http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/SampleNetwork.
## The SampleNetwork function was written by Mike Oldham (oldhamm@stemcell.ucsf.edu).
## DISCLAIMER: THIS CODE IS FREELY PROVIDED FOR ACADEMIC USE WITH ABSOLUTELY NO WARRANTIES.  PLEASE CONTACT MIKE OLDHAM WITH BUG REPORTS OR SUGGESTIONS.
## The ComBat function was written by WE Johnson (http://statistics.byu.edu/johnson/ComBat/Abstract.html).
## The ModulePrinComps1 function was written by Steve Horvath (shorvath@mednet.ucla.edu).

## To cite this code or the methods contained therein, please use the following references:

# 1) Oldham MC, Langfelder P, Horvath S (2011) "Sample Networks for Enhancing Cluster Analysis of Genomic Data: Application to Huntington's Disease".  Submitted.
# 2) Johnson WE, Li C, Rabinovic A (2007) "Adjusting batch effects in microarray expression data using empirical Bayes methods". Biostatistics 8: 118-127.

## The function "SampleNetwork" takes as input two files:
## 1) a feature activity file (argument "datExprT"; rows = features (e.g. probe sets), columns = feature information (e.g. probe set IDs, gene symbols, etc.) + samples);
## 2) a sample information file (argument "sampleinfo1"; rows = samples, columns = sample traits);
## Additional arguments:
## "method1" = distance measure to be used fo constructing sample networks.  choices are "correlation" and "euclidean".
## "impute1" = logical (T/F): should missing values (must be coded as "NA") be imputed (default = FALSE);
## "exclude1" = "all" or "any": should genes with 0 variance over ALL samples be excluded or should genes with 0 variance over ANY group of samples defined by indices1 be excluded (default);
## "subset1" = a binary indicator vector (0,1) that can specify a subset of genes for the analysis; by default, set to include all genes;
## "skip1" = integer describing the number of feature info columns in the expression matrix (there must be at least one);
## "indices1" = a list of vectors to subset the columns in the expression matrix (each vector in this list defines a separate group of samples for processing, e.g. "CTRL", "EXP", etc.);
## "subgroup1" = integer that points to the column in sampleinfo1 that specifies sample subgroups to be colored separately in plots (default = NULL);
## "samplelabels1" = integer that points to the column in sampleinfo1 with the sample labels that will appear in plots (note: these must be identical to the sample column headers in datExprT);
## "grouplabels1" = integer that points to the column in the sample information file that containes the group labels (note: number of groups must equal number of indices);
## "fitmodels1" = logical (T/F): should model-fitting be performed?  if FALSE (default), standard clustering by 1-ISA with accompanying Z.K plot will be produced for outlier removal;
## "whichmodel1" = should a univariate (default) or multivariate linear regression model be used?
## "whichfit1" = variable that specifies the quantity (y-variable) that will be regressed during model-fitting; choices are "pc1" (default; the 1st principal component of the entire dataset [defined by subset1]), "mean" (sample mean), or "K" (sample connectivity);
## "btrait1" = a vector of integers that identifies the columns (traits) in sampleinfo1 that should be included in model-fitting (only specified if fitmodels1 = TRUE; otherwise, NULL [default]);
## "trait1" =  a vector of integers that identifies the columns (traits) in sampleinfo1 that should be tested for the significance of individual factor levels in model-fitting; trait1 must be a subset of btrait1 and specify only categorical variables; only specified if fitmodels1 = TRUE; otherwise, NULL [default];
## "asfactors1" = a vector of integers that identifies the columns (variables) in sampleinfo1 that have been assigned to btrait1 and should be treated as factors;
## "projectname1" = a character label for the project that will appear in some output files;
## "cexlabels1" = a scaling factor for the sample labels in plots (default = 1);
## "normalize1" = logical (T/F): should quantile normalization be performed at the probe set level? default=TRUE;
## "replacenegs1" = logical (T/F): should negative expression values introduced by ComBat or already present in the matrix be replaced by the median for the corresponding probe set? default=FALSE;
## "exportfigures1" = logical (T/F): should figures generated by function be exported as pdfs? if TRUE (default), a subdirectory called "/SampleNetwork" will be created that contains each figure generated during iterative rounds of outlier removal;
## "verbose" = logical (T/F): if TRUE (default), additional network metrics will be exported as text on screen, in plots, and in "..._metrics.csv" output file.  Metrics are produced using the fundamentalNetworkConcepts function from WGCNA library;

## Note: SampleNetwork requires the following libraries: affy, cluster, impute, preprocessCore, and WGCNA;
## Note: requires WGCNA function "ModulePrinComps1" (included below);
## Note: requires "ComBat" function for batch normalization (see below).


library(affy)
library(cluster)
library(impute)
library(preprocessCore)
library(WGCNA)

allowWGCNAThreads()

## SampleNetwork function below (by Mike Oldham):

if(exists("SampleNetwork")) rm(SampleNetwork);
SampleNetwork=function(datExprT,method1="correlation",impute1=FALSE,exclude1="any",subset1=NULL,skip1,indices1,sampleinfo1,subgroup1=NULL,samplelabels1,grouplabels1,fitmodels1=FALSE,whichmodel1="univariate",whichfit1="pc1",btrait1=NULL,trait1=NULL,asfactors1=NULL,projectname1,cexlabels1=1,normalize1=TRUE,replacenegs1=FALSE,exportfigures1=TRUE,verbose=TRUE){  
  
  ## Check for numeric data:
  checknumeric=c()
  for(e in c((skip1+1):length(datExprT[1,]))){
  	checknumeric=c(checknumeric,is.numeric(datExprT[,e]))
  	}
  if(all(checknumeric)!=TRUE){
  	stop("expression matrix contains non-numeric data")
  	}
  
  ## Enforce concordance between datExprT and sampleinfo1:
  rownames(sampleinfo1)=c(1:dim(sampleinfo1)[1])
  groups1=unique(sampleinfo1[,grouplabels1])
  groups1=groups1[!is.na(groups1)]
  if(length(indices1)[1]!=length(groups1)){
    stop("number of indices does not match number of groups")
	}
  maxindex=c()
  for(i in c(1:length(indices1)[1])){
    maxindex=max(maxindex,max(indices1[[i]]))
	}
  matchlabels=data.frame(dimnames(datExprT)[[2]][(skip1+1):maxindex],sampleinfo1[,samplelabels1])
  dimnames(matchlabels)[[2]]=c("datExprT","sampleinfo1")
  Match=as.character(matchlabels[,1])==as.character(matchlabels[,2])
  if(length(Match[Match])!=length(dimnames(datExprT)[[2]][(skip1+1):maxindex])){
	stop("Sample labels in datExprT and sampleinfo1 do not match!")
	}
	  
  ## Enforce order of group indices listed in indices1 (these must be listed in the same order that groups first appear in datSample1, datExprT):
  gpordervec=c()
  for(f in c(1:length(indices1))){
    gpordervec=c(gpordervec,indices1[[f]][1])
	}
  rankorder=rank(gpordervec)
  indices1=indices1[rankorder]
  
  ## Subset expression data:
  if(is.null(subset1)){
  	subset1=rep(TRUE,length(datExprT[,1]))
  	} else {
  	  if(length(subset1)!=length(datExprT[,1])){
  	  	stop("length of subset1 must equal number of rows in datExprT")
  	  	}
  	  subset1=as.logical(subset1)
  	  }
  datExprT=datExprT[subset1,]
  
## Create SampleNetworks root directory, if necessary:
	if(length(grep(paste(projectname1,"_SampleNetworks",sep=""),getwd()))==1){
		breakdir=strsplit(getwd(),split="/")
		SNroot=c(1:length(breakdir[[1]]))[is.element(breakdir[[1]],paste(projectname1,"_SampleNetworks",sep=""))]
		setwd(paste(breakdir[[1]][1:SNroot],collapse="/"))
	} else {
		dir.create(paste(projectname1,"_SampleNetworks",sep=""))
		setwd(paste(getwd(),"/",projectname1,"_SampleNetworks",sep=""))
	}
	SNrootDir=getwd()	
	
## Check for missing data and set NA = 1 per user (recommended for RNAseq data):
	missing.values=sum(is.na(datExprT[,c((skip1+1):maxindex)]))
	total.values=nrow(datExprT)*length(c((skip1+1):maxindex))
	if(missing.values>0){
		print(paste(signif((missing.values/total.values)*100,2),"% of data are missing. Set NA = 0? (recommended for RNAseq data)",sep=""),quote=F)
		print("Alternatively, consider setting 'impute1=TRUE' for missing microarray (or other) data",quote=F)
		answer=readline(prompt="Response (y/n): ")
		while(answer==""){
			answer=readline(prompt="Response (y/n): ")
   	    }
		while(length(intersect(c("y","n"),answer))==0){
			print("Error! Response must be y/n. Please try again.")
			answer=readline(prompt="Response (y/n): ")
  	    }
		if(answer=="n"){
			break()
		} else {
			datExprT[,c((skip1+1):maxindex)][is.na(datExprT[,c((skip1+1):maxindex)])]=1
		}
	}
		
  ## Exclude probe sets with 0 variance:
  excludevec=rep(1,length(datExprT[,1]))
  if(exclude1=="all"){
	variance=apply(datExprT[,unlist(indices1)],1,var,na.rm=T)
	restvar=variance==0
	excludevec[restvar]=0
	}
  if(exclude1=="any"){
	for(a in c(1:length(indices1)[1])){
    variance=apply(datExprT[,indices1[[a]]],1,var,na.rm=T)
	restvar=variance==0
	excludevec[restvar]=0
	}
  }
  excludevec=as.logical(excludevec)
  if(length(excludevec[excludevec])>0){
      datExprOrig=datExprT
	  print(paste("Note: ",length(excludevec[!excludevec])," probe sets had 0 variance and were excluded",sep=""))
	}
  datExprT=datExprT[excludevec,]
  collectGarbage()
  
  ## Impute missing data:
  if(impute1==TRUE){
    missingtable=table(is.na(datExprT[,c((skip1+1):maxindex)]))
	if(is.na(missingtable[2])){
	  print("No missing data...no imputation")
	  } else {
	    print("Imputing missing data...")
	    datimpute=impute.knn(as.matrix(datExprT[,c((skip1+1):maxindex)]))
	    datExprT=data.frame(datExprT[,c(1:skip1)],datimpute$data)
	    }
	  }	    
  
  ## Enforce asfactors1, whichmodel1, and whichfit1:
  if(fitmodels1==TRUE){
  	if(length(asfactors1)>0){
      for(h in c(1:length(asfactors1))){
        sampleinfo1[,asfactors1[h]]=factor(sampleinfo1[,asfactors1[h]])
        }
      }
  	if(length(grep("univariate",whichmodel1))!=1&length(grep("multivariate",whichmodel1))!=1){
	  stop("whichmodel1 must equal 'univariate' or 'multivariate'")
	  }
	if(length(grep("K",whichfit1))!=1&length(grep("mean",whichfit1))!=1&length(grep("pc1",whichfit1))!=1){
  	  stop("whichfit1 must equal 'K', 'mean', or 'pc1'")
  	  }
  	}
		
  ## Initialize data frames for export:
  Outliers=data.frame()
  MetricsDF=data.frame(t(c(1:12)))
  dimnames(MetricsDF)[[2]]=c("Group","Round","Samples","Mean_IAC","Mean_Connectivity","Mean_ScaledConnectivity","Mean_ClusterCoef","Mean_MAR","Density","Decentralization","Homogeneity","PC1_VE")
  MetricsDF=MetricsDF[-1,]
  
  ## MakePlots1 function:
  MakePlots1=function(datexpr2,IAC2,adj2,sampleinfo2,meansample2,indexgp2,group2,subgroup2,btrait2,trait2,asfactors2,whichround,fitmodels2,whichmodel2,whichfit2,exportfigures2,verbose2){
	  if(whichround==1|whichround>1&is.numeric(whichround)){
	    datexprgp=datexpr2[,indexgp2]
		}
	  if(length(intersect(whichround,c("Qnorm","ComBat")))>0){
	    datexprgp=datexpr2[,c((skip1+1):length(datexpr2[1,]))]
		}
  	meansample2=apply(datexprgp,2,mean,na.rm=T)
	allsamples=c(1:length(indices1[[r]]))
	## For coloring subgroups in plots:
	subgpcolors=standardColors()[c(7,6,2,5,3,1,4,8:length(unique(sampleinfo2[indexgp2-skip1,subgroup2])))]
	colorvec=rep("black",length(datexprgp[1,]))
	if(!is.null(subgroup2)){
	  whichsubgroups=sort(unique(sampleinfo2[indexgp2-skip1,subgroup2]))
	  for(b in c(1:length(whichsubgroups))){
	    colorvec[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[b]]=subgpcolors[b]
		}
	  }
	## To set colored leaf labels (adapted from dendrapply function help):
    local({
    colLab <<- function(n,treeorder) {
      if(is.leaf(n)) {
        a <- attributes(n)
        b <<- b+1
        attr(n, "nodePar") <- c(a$nodePar, list(lab.col = colorvec[treeorder][b], lab.font= i%%3))
        }
      n
     }
    b <- 0
    })
	if(whichround==1){
	  FNC=fundamentalNetworkConcepts(adj2)
	  K2=FNC$ScaledConnectivity
	  gp.var=apply(datexprgp,1,var)
	  datexprgp=datexprgp[gp.var>0,]
      pcs1all=ModulePrinComps1(datexpr=as.matrix(t(datexprgp)),couleur=rep("gold",length(datexprgp[,1])))
	  pc2=pcs1all$PrinComps$PCgold
	  Metrics=data.frame(groups1[r],whichround,length(indexgp2),mean(IAC[upper.tri(IAC)]),mean(FNC$Connectivity),mean(FNC$ScaledConnectivity),mean(FNC$ClusterCoef),mean(FNC$MAR),FNC$Density,1-FNC$Centralization,1-FNC$Heterogeneity,pcs1all$varexplained[1,])
	  dimnames(Metrics)[[2]]=c("Group","Round","Samples","Mean_IAC","Mean_Connectivity","Mean_ScaledConnectivity","Mean_ClusterCoef","Mean_MAR","Density","Decentralization","Homogeneity","PC1_VE")
	  }
	if(whichround>1&is.numeric(whichround)){
	  samplesubset=allsamples[is.element(indices1[[r]],indexgp2)]
	  IAC=IAC2[samplesubset,samplesubset]
	  A.IAC=adj2[samplesubset,samplesubset]
	  FNC=fundamentalNetworkConcepts(A.IAC)
	  K2=FNC$ScaledConnectivity
	  newexcludevec=rep(1,length(datexprgp[,1]))
	  newvariance=apply(datexprgp,1,var)
	  newrestvar=newvariance==0
	  newexcludevec[newrestvar]=0
	  newexcludevec=as.logical(newexcludevec)
      datexprgp2=datexprgp[newexcludevec,]
	  gp.var=apply(datexprgp2,1,var)
	  datexprgp2=datexprgp2[gp.var>0,]
	  pcs1all=ModulePrinComps1(datexpr=as.matrix(t(datexprgp2)),couleur=rep("gold",length(datexprgp2[,1])))
	  pc2=pcs1all$PrinComps$PCgold
	  NewMetrics=data.frame(groups1[r],whichround,length(indexgp2),mean(IAC[upper.tri(IAC)]),mean(FNC$Connectivity),mean(FNC$ScaledConnectivity),mean(FNC$ClusterCoef),mean(FNC$MAR),FNC$Density,1-FNC$Centralization,1-FNC$Heterogeneity,pcs1all$varexplained[1,])
	  dimnames(NewMetrics)[[2]]=c("Group","Round","Samples","Mean_IAC","Mean_Connectivity","Mean_ScaledConnectivity","Mean_ClusterCoef","Mean_MAR","Density","Decentralization","Homogeneity","PC1_VE")
	  Metrics=data.frame(rbind(attr(datoutnew,"MetricsOut"),NewMetrics))
	  }
 	if(length(intersect(whichround,c("Qnorm","ComBat")))>0){
 	  print("Recalculating sample correlations...")
 	  collectGarbage()
 	  IAC=cor(datexprgp,method="p",use="p")
	  diag(IAC)=0
	  if(method1=="correlation"){
		  A.IAC=((1+IAC)/2)^2
	  }
	  if(method1=="euclidean"){
		  D.squared=as.matrix(dist(t(datexprgp),method="euclidean")^2)
		  A.IAC=1-D.squared/max(D.squared)
	  }
	  diag(A.IAC)=0
	  FNC=fundamentalNetworkConcepts(A.IAC)
	  K2=FNC$ScaledConnectivity
	  newexcludevec=rep(1,length(datexprgp[,1]))
	  newvariance=apply(datexprgp,1,var)
	  newrestvar=newvariance==0
	  newexcludevec[newrestvar]=0
	  newexcludevec=as.logical(newexcludevec)
      datexprgp2=datexprgp[newexcludevec,]
	  gp.var=apply(datexprgp2,1,var)
	  datexprgp2=datexprgp2[gp.var>0,]
	  pcs1all=ModulePrinComps1(datexpr=as.matrix(t(datexprgp2)),couleur=rep("gold",length(datexprgp2[,1])))
	  pc2=pcs1all$PrinComps$PCgold
	  NewMetrics=data.frame(groups1[r],whichround,length(indexgp2),mean(IAC[upper.tri(IAC)]),mean(FNC$Connectivity),mean(FNC$ScaledConnectivity),mean(FNC$ClusterCoef),mean(FNC$MAR),FNC$Density,1-FNC$Centralization,1-FNC$Heterogeneity,pcs1all$varexplained[1,])
	  dimnames(NewMetrics)[[2]]=c("Group","Round","Samples","Mean_IAC","Mean_Connectivity","Mean_ScaledConnectivity","Mean_ClusterCoef","Mean_MAR","Density","Decentralization","Homogeneity","PC1_VE")
	  Metrics=data.frame(rbind(attr(datoutnew,"MetricsOut"),NewMetrics))
	  }
    meanIAC=mean(IAC[upper.tri(IAC)],na.rm=T)
    cluster1=hclust(as.dist(1-A.IAC),method="average")
	cluster1order=cluster1$order
	cluster2=as.dendrogram(cluster1,hang=0.1)
	cluster2=dendrapply(cluster2,colLab,cluster1order)
    Z.K=(K2-mean(K2))/sd(K2)
	Z.C=(FNC$ClusterCoef-mean(FNC$ClusterCoef))/sd(FNC$ClusterCoef)
	Z.MAR=(FNC$MAR-mean(FNC$MAR))/sd(FNC$MAR)
	if(fitmodels2==TRUE){
	  sampleinfobtrait=data.frame(sampleinfo2[indexgp2-skip1,btrait2])
      onelevelcheckbtrait=c()
	  for(z in c(1:dim(sampleinfobtrait)[2])){
      	onelevelcheckbtrait=c(onelevelcheckbtrait,length(unique(sampleinfobtrait[,z]))>1)
      	}
	  if(!all(onelevelcheckbtrait)){
      	btrait2=btrait2[onelevelcheckbtrait]
      	}	
      if(length(trait2)>0){
	    sampleinfotrait=data.frame(sampleinfo2[indexgp2-skip1,trait2])
		onelevelchecktrait=c()
		for(y in c(1:dim(sampleinfotrait)[2])){
		  onelevelchecktrait=c(onelevelchecktrait,length(unique(sampleinfotrait[,y]))>1)
      	  }
		if(!all(onelevelchecktrait)){
		  trait2=trait2[onelevelchecktrait]
      	  }
		}
	  alltraits=union(btrait2,trait2)
      allmodelfactors=is.element(alltraits,asfactors2)
      allmodelterms=paste("sampleinfo2[indexgp2-skip1,",alltraits,"]",sep="")
      allmodelterms[allmodelfactors]=paste("factor(",allmodelterms[allmodelfactors],")",sep="")
      origmodelterms=allmodelterms
      allmodelterms=paste(allmodelterms,collapse="+")
     
	   if(whichfit2=="mean"){
        if(whichmodel2=="univariate"){
		  allsigniflist=c()
		  alltempout=c()
		  for(f in c(1:length(origmodelterms))){
		    allmodelname=paste("meansample2~",origmodelterms[f],sep="")
			alltemp=anova(lm(as.formula(allmodelname)))
			alltempout=c(alltempout,gsub(" ","",rownames(alltemp)[1]))
			allsigniflist=c(allsigniflist,alltemp[1:(dim(alltemp)[1]-1),5])
			}
		  allsigniflist[is.nan(allsigniflist)]=1
          lmoutput=lm(as.formula(allmodelname))
		  } else {
		  allmodelname=paste("meansample2~",allmodelterms,sep="")  
          alltemp=anova(lm(as.formula(allmodelname)))
		  alltempout=gsub(" ","",rownames(alltemp))
          allsigniflist=alltemp[1:(dim(alltemp)[1]-1),5]
          allsigniflist[is.nan(allsigniflist)]=1
          lmoutput=lm(as.formula(allmodelname))
          }
		bplabel="mean expression"
		}
		
      if(whichfit2=="K"){
      	if(whichmodel2=="univariate"){
		  allsigniflist=c()
		  alltempout=c()
		  for(f in c(1:length(origmodelterms))){
		    allmodelname=paste("K2~",origmodelterms[f],sep="")
			alltemp=anova(lm(as.formula(allmodelname)))
			alltempout=c(alltempout,gsub(" ","",rownames(alltemp)[1]))
			allsigniflist=c(allsigniflist,alltemp[1:(dim(alltemp)[1]-1),5])
			}
		  allsigniflist[is.nan(allsigniflist)]=1
          lmoutput=lm(as.formula(allmodelname))
		  } else {
		  allmodelname=paste("K2~",allmodelterms,sep="")  
		  alltemp=anova(lm(as.formula(allmodelname)))
		  alltempout=gsub(" ","",rownames(alltemp))
		  allsigniflist=alltemp[1:(dim(alltemp)[1]-1),5]
		  allsigniflist[is.nan(allsigniflist)]=1
		  lmoutput=lm(as.formula(allmodelname))
		  }
		bplabel="K"
		}	
		
      if(whichfit2=="pc1"){
      	if(whichmodel2=="univariate"){
		  allsigniflist=c()
		  alltempout=c()
		  for(f in c(1:length(origmodelterms))){
		    allmodelname=paste("pc2~",origmodelterms[f],sep="")
			alltemp=anova(lm(as.formula(allmodelname)))
			alltempout=c(alltempout,gsub(" ","",rownames(alltemp)[1]))
			allsigniflist=c(allsigniflist,alltemp[1:(dim(alltemp)[1]-1),5])
		  }
		  allsigniflist[is.nan(allsigniflist)]=1
          lmoutput=lm(as.formula(allmodelname))
		  } else {
		  allmodelname=paste("pc2~",allmodelterms,sep="")  
      	  alltemp=anova(lm(as.formula(allmodelname)))
		  alltempout=gsub(" ","",rownames(alltemp))
          allsigniflist=alltemp[1:(dim(alltemp)[1]-1),5]
          allsigniflist[is.nan(allsigniflist)]=1
          lmoutput=lm(as.formula(allmodelname))
          }
		bplabel="pc1"
		}
		     	    
      if(length(trait2)>0){
      varlist=vector(mode="list",length=length(trait2))
      for(q in c(1:length(trait2))){
          traitdf=data.frame(table(sampleinfo2[indexgp2-skip1,trait2[q]]))
		  varlist[[q]]=sort(unique(traitdf$Var1[traitdf$Freq>1]))
          }
      signiflist=vector(mode="list",length=length(trait2))
      for(j in c(1:length(trait2))){
        btrait3=setdiff(btrait2,trait2[j])
        whichmodelfactors=is.element(btrait3,asfactors2)
        modelterms=paste("sampleinfo2[indexgp2-skip1,",btrait3,"]",sep="")
        modelterms[whichmodelfactors]=paste("factor(",modelterms[whichmodelfactors],")",sep="")
        modelterms=paste("+",paste(modelterms,collapse="+"),sep="")
        if(length(btrait3)==0|whichmodel2=="univariate"){
        	modelterms=NULL
        	}
        whichtestfactors=is.element(trait2,asfactors2)
		if(whichfit2=="mean"){
          if(whichtestfactors[j]==TRUE){
            for(k in c(1:length(varlist[[j]]))){
              modelname=paste("meansample2~factor(sampleinfo2[indexgp2-skip1,trait2[",j,"]]==varlist[[",j,"]][",k,"])",modelterms,sep="")
              temp1=anova(lm(as.formula(modelname)))
              signiflist[[j]]=c(signiflist[[j]],temp1[1,5])
              }
            } else {
      	      for(k in c(1:length(varlist[[j]]))){
				modelname=paste("meansample2~sampleinfo2[indexgp2-skip1,trait2[",j,"]]==varlist[[",j,"]][",k,"]",modelterms,sep="")
				temp1=anova(lm(as.formula(modelname)))
				signiflist[[j]]=c(signiflist[[j]],temp1[1,5])
				}
              }
			}
        if(whichfit2=="K"){
          if(whichtestfactors[j]==TRUE){
            for(k in c(1:length(varlist[[j]]))){
              modelname=paste("K2~factor(sampleinfo2[indexgp2-skip1,trait2[",j,"]]==varlist[[",j,"]][",k,"])",modelterms,sep="")
              temp1=anova(lm(as.formula(modelname)))
              signiflist[[j]]=c(signiflist[[j]],temp1[1,5])
              }
            } else {
      	      for(k in c(1:length(varlist[[j]]))){
				modelname=paste("K2~sampleinfo2[indexgp2-skip1,trait2[",j,"]]==varlist[[",j,"]][",k,"]",modelterms,sep="")
				temp1=anova(lm(as.formula(modelname)))
				signiflist[[j]]=c(signiflist[[j]],temp1[1,5])
				}
			  }
			}
        if(whichfit2=="pc1"){
          if(whichtestfactors[j]==TRUE){
            for(k in c(1:length(varlist[[j]]))){
              modelname=paste("pc2~factor(sampleinfo2[indexgp2-skip1,trait2[",j,"]]==varlist[[",j,"]][",k,"])",modelterms,sep="")
              temp1=anova(lm(as.formula(modelname)))
              signiflist[[j]]=c(signiflist[[j]],temp1[1,5])
              }
            } else {
      	      for(k in c(1:length(varlist[[j]]))){
				modelname=paste("pc2~sampleinfo2[indexgp2-skip1,trait2[",j,"]]==varlist[[",j,"]][",k,"]",modelterms,sep="")
				temp1=anova(lm(as.formula(modelname)))
				signiflist[[j]]=c(signiflist[[j]],temp1[1,5])
				}
			  }
            }
            signiflist[[j]][is.nan(signiflist[[j]])]=1
          } ## end of for(j in c(1:length(trait2))){
        }  ## end of if(length(trait2)>0){
		
	  if(verbose2==TRUE){
	    plotrows=3
		plotcols=ceiling((length(trait2)+5)/3)
		} else {
		  plotcols=ceiling((length(trait2)+3)/3)
		  if(length(trait2)>1){
		    plotrows=3
		    } else {
	          plotrows=2
			  }
		    }
	  par(mfrow=c(plotrows,plotcols))
      par(mar=c(5,5,4,2))
	  plot(cluster2,nodePar=list(lab.cex=cexlabels1,pch=NA),main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
	  mtext(paste("distance: ",method1,sep=""),cex=0.8,line=0.2)
      par(mar=c(5,5,4,2))
      plot(Z.K,main="Connectivity", ylab="Z.K",type="n",xaxt="n",xlab="Sample",cex.main=1.8,cex.lab=1.4)
      text(Z.K,labels=sampleinfo2[indexgp2-skip1,samplelabels1],cex=cexlabels1,col=colorvec)
      abline(h=-2)
	  abline(h=-3)
	  if(verbose2==TRUE){
	    par(mar=c(5,5,4,2))
		plot(Z.C,main="ClusterCoef", ylab="Z.C",type="n",xaxt="n",xlab="Sample",cex.main=1.8,cex.lab=1.4)
	    text(Z.C,labels=sampleinfo2[indexgp2-skip1,samplelabels1],cex=cexlabels1,col=colorvec)
		abline(h=2)
		abline(h=3)
		abline(h=-2)
		abline(h=-3)		
		par(mar=c(5,5,4,2))
		plot(Z.K,Z.C,main="Connectivity vs ClusterCoef",xlab="Z.K",ylab="Z.C",col=colorvec,cex.main=1.8,cex.lab=1.4)
		  if(is.null(subgroup2)){
			  abline(lm(Z.C~Z.K),col="black",lwd=2)
			  mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
		  }
		  if(!is.null(subgroup2)){
			if(length(whichsubgroups)==2){
				for(b in c(1:length(whichsubgroups))){
					abline(lm(Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[b]]~Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[b]]),col=subgpcolors[b],lwd=2)
					}
				mtext(paste("rho = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],method="s")$estimate,2)," p = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],method="s")$p.value,2),
				"; ","rho = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],method="s")$estimate,2)," p = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],method="s")$p.value,2),
				sep=""),cex=0.8,line=0.2)
				} else {
				abline(lm(Z.C~Z.K),col="black",lwd=2)
				mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
				}
		}
		  
	  }
	  if(length(intersect(origmodelterms,alltempout))!=length(origmodelterms)){
#alltraitsbp=sort(alltraits[is.element(origmodelterms,alltempout)|is.element(alltempout,"Residuals")])
		alltraitsbp=alltraits[is.element(origmodelterms,alltempout)]
        } else {
          alltraitsbp=alltraits
          }
      allsigniflist[allsigniflist==0]=1e-300
	  par(mar=c(8,5,4,2))
	  barplot(-log(allsigniflist,10),names=as.character(dimnames(sampleinfo2)[[2]][alltraitsbp]),las=3,cex.names=0.8,ylab="-Log10 p-value",main=paste("ANOVA of ",bplabel,sep=""),cex.main=1.8,cex.lab=1.4)
	  mtext(whichmodel2,cex=0.8,line=0.2)
      abline(h=-log(.05,10),col="blue",lwd=2)
      abline(h=-log(.05/length(alltraitsbp),10),col="red",lwd=2)
      if(length(trait2)>0){
        for(m in c(1:length(signiflist)[[1]])){
          varlabel=dimnames(sampleinfo2)[[2]][trait2][m]  
          signiflist[[m]][signiflist[[m]]==0]=1e-300
		  par(mar=c(8,5,4,2))
		  plot(-log(signiflist[[m]],10),type="n",xlab=varlabel,ylab="-Log10 p-value",main=paste("ANOVA of ",bplabel,sep=""),xaxt="n",cex.main=1.8,cex.lab=1.4)
          mtext(whichmodel2,cex=0.8,line=0.2)
		  text(-log(signiflist[[m]],10),labels=varlist[[m]])
          abline(h=-log(.05,10),col="blue",lwd=2)
          abline(h=-log(.05/length(varlist[[m]]),10),col="red",lwd=2)
          } ## end of for(m in c(1:length(signiflist)[[1]])){
        } ## end of if(length(trait2)>0){
        if(verbose2==FALSE){
		  datout=data.frame(sampleinfo2[indexgp2-skip1,c(grouplabels1,samplelabels1,union(btrait2,trait2))],signif(K2,3),signif(Z.K,3))
          dimnames(datout)[[2]][length(datout[1,])-1]="K"
          dimnames(datout)[[2]][length(datout[1,])]="Z.K"
          newskip=4+length(union(btrait2,trait2))
          datout=datout[order(datout[,newskip]),]
		  } else {
		    datout=data.frame(sampleinfo2[indexgp2-skip1,c(grouplabels1,samplelabels1,union(btrait2,trait2))],signif(K2,3),signif(Z.K,3),signif(Z.C,3),signif(Z.MAR,3))
		    dimnames(datout)[[2]][length(datout[1,])-3]="K"
            dimnames(datout)[[2]][length(datout[1,])-2]="Z.K"
			dimnames(datout)[[2]][length(datout[1,])-1]="Z.C"
			dimnames(datout)[[2]][length(datout[1,])]="Z.MAR"
            newskip=4+length(union(btrait2,trait2))
            datout=datout[order(datout[,newskip]),]
		    }
      if(exportfigures2==TRUE){
		pdf(file=paste(group2,"_rd_",whichround,"_",whichfit2,".pdf",sep=""))
		par(mfrow=c(plotrows,plotcols))
        par(mar=c(5,5,4,2))
		plot(cluster2,nodePar=list(lab.cex=cexlabels1,pch=NA),main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
		mtext(paste("distance: ",method1,sep=""),cex=0.8,line=0.2)
        par(mar=c(5,5,4,2))
		plot(Z.K,main="Connectivity",ylab="Z.K",type="n",xaxt="n",xlab="Sample",cex.main=1.8,cex.lab=1.4)
        text(Z.K,labels=sampleinfo2[indexgp2-skip1,samplelabels1],cex=cexlabels1,col=colorvec)
        abline(h=-2)
		abline(h=-3)
		if(verbose2==TRUE){
	      par(mar=c(5,5,4,2))
		  plot(Z.C,main="ClusterCoef", ylab="Z.C",type="n",xaxt="n",xlab="Sample",cex.main=1.8,cex.lab=1.4)
	      text(Z.C,labels=sampleinfo2[indexgp2-skip1,samplelabels1],cex=cexlabels1,col=colorvec)
		  abline(h=2)
		  abline(h=3)
		  abline(h=-2)
		  abline(h=-3)
		  par(mar=c(5,5,4,2))
		  plot(Z.K,Z.C,main="Connectivity vs ClusterCoef",xlab="Z.K",ylab="Z.C",col=colorvec,cex.main=1.8,cex.lab=1.4)
			if(is.null(subgroup2)){
				abline(lm(Z.C~Z.K),col="black",lwd=2)
				mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
			}
			if(!is.null(subgroup2)){
			if(length(whichsubgroups)==2){
				for(b in c(1:length(whichsubgroups))){
					abline(lm(Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[b]]~Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[b]]),col=subgpcolors[b],lwd=2)
					}
				mtext(paste("rho = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],method="s")$estimate,2)," p = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],method="s")$p.value,2),
				"; ","rho = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],method="s")$estimate,2)," p = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],method="s")$p.value,2),
				sep=""),cex=0.8,line=0.2)
				} else {
				abline(lm(Z.C~Z.K),col="black",lwd=2)
				mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
				}
		  }
			}
        allsigniflist[allsigniflist==0]=1e-300
		par(mar=c(8,5,4,2))
		barplot(-log(allsigniflist,10),names=as.character(dimnames(sampleinfo2)[[2]][alltraitsbp]),las=3,cex.names=0.8,ylab="-Log10 p-value",main=paste("ANOVA of ",bplabel,sep=""),cex.main=1.8,cex.lab=1.4)
        mtext(whichmodel2,cex=0.8,line=0.2)
		abline(h=-log(.05,10),col="blue",lwd=2)
        abline(h=-log(.05/length(alltraitsbp),10),col="red",lwd=2)
        if(length(trait2)>0){
          for(m in c(1:length(signiflist)[[1]])){
            varlabel=dimnames(sampleinfo2)[[2]][trait2][m]  
            signiflist[[m]][signiflist[[m]]==0]=1e-300
			par(mar=c(8,5,4,2))
			plot(-log(signiflist[[m]],10),type="n",xlab=varlabel,ylab="-Log10 p-value",main=paste("ANOVA of ",bplabel,sep=""),xaxt="n",cex.main=1.8,cex.lab=1.4)
            mtext(whichmodel2,cex=0.8,line=0.2)
			text(-log(signiflist[[m]],10),labels=varlist[[m]])
            abline(h=-log(.05,10),col="blue",lwd=2)
            abline(h=-log(.05/length(varlist[[m]]),10),col="red",lwd=2)
            }
          }
        dev.off()
		}
	    ## end of if(exportfigures2==TRUE)
      } else {
        ## end of if(fitmodels2==TRUE){
        if(verbose2==FALSE){
		  datout=data.frame(sampleinfo2[indexgp2-skip1,c(grouplabels1,samplelabels1,union(btrait2,trait2))],signif(K2,3),signif(Z.K,3))
          dimnames(datout)[[2]][length(datout[1,])-1]="K"
          dimnames(datout)[[2]][length(datout[1,])]="Z.K"
          datout=datout[order(datout[,length(datout[1,])]),]
		  par(mfrow=c(2,1))
		  par(mar=c(5,5,4,2))
          plot(cluster2,nodePar=list(lab.cex=cexlabels1,pch=NA),main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
		  mtext(paste("distance: ",method1,sep=""),cex=0.8,line=0.2)
          par(mar=c(5,5,4,2))
		  plot(Z.K,main="Connectivity", ylab="Z.K", type="n", xaxt="n", xlab="Sample",cex.main=1.8,cex.lab=1.4)
          text(Z.K,labels=sampleinfo2[indexgp2-skip1,samplelabels1],cex=cexlabels1,col=colorvec)
          abline(h=-2)
		  abline(h=-3)
		  } else {
		    datout=data.frame(sampleinfo2[indexgp2-skip1,c(grouplabels1,samplelabels1,union(btrait2,trait2))],signif(K2,3),signif(Z.K,3),signif(Z.C,3),signif(Z.MAR,3))
            dimnames(datout)[[2]][length(datout[1,])-3]="K"
            dimnames(datout)[[2]][length(datout[1,])-2]="Z.K"
			dimnames(datout)[[2]][length(datout[1,])-1]="Z.C"
			dimnames(datout)[[2]][length(datout[1,])]="Z.MAR"
            newskip=4+length(union(btrait2,trait2))
            datout=datout[order(datout[,newskip]),]
			par(mfrow=c(2,2))
		    par(mar=c(5,5,4,2))
			plot(cluster2,nodePar=list(lab.cex=cexlabels1,pch=NA),main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
			mtext(paste("distance: ",method1,sep=""),cex=0.8,line=0.2)
            par(mar=c(5,5,4,2))
			plot(Z.K,main="Connectivity", ylab="Z.K", type="n", xaxt="n", xlab="Sample",cex.main=1.8,cex.lab=1.4)
            text(Z.K,labels=sampleinfo2[indexgp2-skip1,samplelabels1],cex=cexlabels1,col=colorvec)
            abline(h=-2)
		    abline(h=-3)
			par(mar=c(5,5,4,2))
			plot(Z.C,main="ClusterCoef", ylab="Z.C", type="n", xaxt="n", xlab="Sample",cex.main=1.8,cex.lab=1.4)
	        text(Z.C,labels=sampleinfo2[indexgp2-skip1,samplelabels1],cex=cexlabels1,col=colorvec)
		    abline(h=2)
		    abline(h=3)
			abline(h=-2)
		    abline(h=-3)
			par(mar=c(5,5,4,2))
			plot(Z.K,Z.C,main="Connectivity vs ClusterCoef",xlab="Z.K",ylab="Z.C",col=colorvec,cex.main=1.8,cex.lab=1.4)
			  if(is.null(subgroup2)){
				  abline(lm(Z.C~Z.K),col="black",lwd=2)
				  mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
			  }
			  if(!is.null(subgroup2)){
				if(length(whichsubgroups)==2){
					for(b in c(1:length(whichsubgroups))){
						abline(lm(Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[b]]~Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[b]]),col=subgpcolors[b],lwd=2)
						}
					mtext(paste("rho = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],method="s")$estimate,2)," p = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],method="s")$p.value,2),
					"; ","rho = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],method="s")$estimate,2)," p = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],method="s")$p.value,2),
					sep=""),cex=0.8,line=0.2)
					} else {
					abline(lm(Z.C~Z.K),col="black",lwd=2)
					mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
					}
			}
				}
		if(exportfigures2==TRUE){
		  pdf(file=paste(group2,"_rd_",whichround,".pdf",sep=""))
    	  if(verbose2==FALSE){
		    par(mfrow=c(2,1))
            par(mar=c(5,5,4,2))
			plot(cluster2,nodePar=list(lab.cex=cexlabels1,pch=NA),main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
			mtext(paste("distance: ",method1,sep=""),cex=0.8,line=0.2)
            par(mar=c(5,5,4,2))
			plot(Z.K,main="Connectivity", ylab="Z.K", type="n", xaxt="n", xlab="Sample",cex.main=1.8,cex.lab=1.4)
            text(Z.K,labels=sampleinfo2[indexgp2-skip1,samplelabels1],cex=cexlabels1,col=colorvec)
            abline(h=-2)
		    abline(h=-3)
			} else {
			  par(mfrow=c(2,2))
		      par(mar=c(5,5,4,2))
			  plot(cluster2,nodePar=list(lab.cex=cexlabels1,pch=NA),main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
			  mtext(paste("distance: ",method1,sep=""),cex=0.8,line=0.2)
              par(mar=c(5,5,4,2))
			  plot(Z.K,main="Connectivity", ylab="Z.K", type="n", xaxt="n", xlab="Sample",cex.main=1.8,cex.lab=1.4)
              text(Z.K,labels=sampleinfo2[indexgp2-skip1,samplelabels1],cex=cexlabels1,col=colorvec)
              abline(h=-2)
		      abline(h=-3)
			  par(mar=c(5,5,4,2))
			  plot(Z.C,main="ClusterCoef", ylab="Z.C",type="n",xaxt="n",xlab="Sample",cex.main=1.8,cex.lab=1.4)
	          text(Z.C,labels=sampleinfo2[indexgp2-skip1,samplelabels1],cex=cexlabels1,col=colorvec)
		      abline(h=2)
		      abline(h=3)
			  abline(h=-2)
		      abline(h=-3)
			  par(mar=c(5,5,4,2))
			  plot(Z.K,Z.C,main="Connectivity vs ClusterCoef",xlab="Z.K",ylab="Z.C",col=colorvec,cex.main=1.8,cex.lab=1.4)
				if(is.null(subgroup2)){
					abline(lm(Z.C~Z.K),col="black",lwd=2)
					mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
				}
				if(!is.null(subgroup2)){
				if(length(whichsubgroups)==2){
					for(b in c(1:length(whichsubgroups))){
						abline(lm(Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[b]]~Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[b]]),col=subgpcolors[b],lwd=2)
						}
					mtext(paste("rho = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],method="s")$estimate,2)," p = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[1]],method="s")$p.value,2),
					"; ","rho = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],method="s")$estimate,2)," p = ",signif(cor.test(Z.K[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],Z.C[sampleinfo2[indexgp2-skip1,subgroup2]==whichsubgroups[2]],method="s")$p.value,2),
					sep=""),cex=0.8,line=0.2)
					} else {
					abline(lm(Z.C~Z.K),col="black",lwd=2)
					mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
					}
			 }
				}
          dev.off()
          } ## end of if(exportfigures2==TRUE)    
        }
      collectGarbage()
	  attr(datout,"MetricsOut")=Metrics
	  attr(datout,"AdjMat")=A.IAC
	  print(datout)
		  cat("\n")
          if(fitmodels2==TRUE){
		    if(length(trait2)>0){
			  if(!all(onelevelchecktrait)){
			    print(paste("Warning:",dimnames(sampleinfotrait)[[2]][!onelevelchecktrait],"contains only one level and was excluded as a trait for model-fitting"))
      	        }
			  }
      	    if(!all(onelevelcheckbtrait)){
      	      print(paste("Warning:",dimnames(sampleinfobtrait)[[2]][!onelevelcheckbtrait],"contains only one level and was excluded as a btrait for model-fitting"))
      	      }	
            if(length(lmoutput$coefficients[is.na(lmoutput$coefficients)])>0){
              print(paste("Warning: the following coefficient was not defined because of a singularity:",names(lmoutput$coefficients)[is.na(lmoutput$coefficients)]))
              }
			}
          if(verbose2==TRUE){
		    cat("\n")
			print(paste("Group ",group2," mean correlation = ",signif(meanIAC,4),sep=""))
			print(paste("Group ",group2," mean connectivity = ",signif(mean(FNC$Connectivity),4),sep=""))
			print(paste("Group ",group2," mean scaled connectivity = ",signif(mean(FNC$ScaledConnectivity),4),sep=""))
			print(paste("Group ",group2," mean clustering coefficient = ",signif(mean(FNC$ClusterCoef),4),sep=""))
			print(paste("Group ",group2," mean maximum adjacency ratio = ",signif(mean(FNC$MAR),4),sep=""))
			print(paste("Group ",group2," density = ",signif(FNC$Density,4),sep=""))
			print(paste("Group ",group2," decentralization = ",signif(1-FNC$Centralization,4),sep=""))
			print(paste("Group ",group2," homogeneity = ",signif(1-FNC$Heterogeneity,4),sep=""))
			print(paste("Group ",group2," PC 1-5 var explained = ",signif(pcs1all$varexplained,4),sep=""))
			}
          cat("\n")
	  datout
	  } ## end of MakePlots1 fx
	
  for(r in c(1:length(indices1)[1])){
  
  ## Create subdirectory for group:
  tstamp1=format(Sys.time(), "%X")
  tstamp1=gsub(" AM","",tstamp1)
  tstamp1=gsub(" PM","",tstamp1)
  tstamp1=gsub(":","-",tstamp1)
  dir.create(paste(groups1[r],"_",tstamp1,sep=""))
  setwd(paste(groups1[r],"_",tstamp1,sep=""))
     
  if(exists("indexgpnew")) rm(indexgpnew);
  print("Calculating sample correlations...")
  IAC=cor(datExprT[,indices1[[r]]],method="p",use="p")
  diag(IAC)=0
  if(method1=="correlation"){
	  A.IAC=((1+IAC)/2)^2
  }
  if(method1=="euclidean"){
	  D.squared=as.matrix(dist(t(datExprT[,indices1[[r]]]),method="euclidean")^2)
	  A.IAC=1-D.squared/max(D.squared)
  }
  diag(A.IAC)=0
  Round=1
  gpoutliers=c()
  K.out=c()
  Z.Kout=c()
  Z.Cout=c()
  Z.MARout=c()
  whichround=c()
  print(paste(groups1[r]," round ",Round,sep=""))  
  datoutnew=MakePlots1(datexpr2=datExprT,IAC2=IAC,adj2=A.IAC,sampleinfo2=sampleinfo1,indexgp2=indices1[[r]],group2=groups1[r],subgroup2=subgroup1,btrait2=btrait1,trait2=trait1,asfactors2=asfactors1,whichround=Round,fitmodels2=fitmodels1,whichmodel2=whichmodel1,whichfit2=whichfit1,exportfigures2=exportfigures1,verbose2=verbose)
    
  print("Please enter the row number(s) from the output of the sample(s) you would like to remove (comma separated) or type 'None'",quote=F)
  print("Alternatively, please enter a sample connectivity threshold for removing outliers. For example, type '<-2' to remove all samples with Z.K less than -2",quote=F)
  answer=readline(prompt="Response: ")
  while(answer!="None"){
    while(answer==""){
  	  answer=readline(prompt="Response: ")
  	  }
  	if(answer[1]=="None"){
  	  break()
  	  }
    if(length(grep("<",answer))>0){
  	  answer=gsub("<","",answer)
	  answer=as.numeric(c(rownames(datoutnew)[datoutnew$Z.K<as.numeric(answer)]))
  	  } else {	
  	    if(length(grep(",",answer))>0){
  	    answer=c(as.numeric(strsplit(as.character(answer),split=",")[[1]]))
  	    }
      }
    while(length(rownames(datoutnew)[is.element(rownames(datoutnew),answer)])!=length(answer)|length(answer)==0){
      print("Error! Please enter the row number(s) from the output of the sample(s) you would like to remove (comma separated) or type 'None'",quote=F)
      print("Alternatively, please enter a sample connectivity threshold for removing outliers. For example, type '<-2' to remove all samples with Z.K less than -2",quote=F)
      answer=readline(prompt="Response: ")
  	  while(answer==""){
  	    answer=readline(prompt="Response: ")
   	    }
   	  if(answer[1]=="None"){
  	  	break()
  	  	}
   	  if(length(grep("<",answer))>0){
  	    answer=gsub("<","",answer)
		answer=as.numeric(c(rownames(datoutnew)[datoutnew$Z.K<as.numeric(answer)]))
  	    } else {	
  	      if(length(grep(",",answer))>0){
  	      answer=c(as.numeric(strsplit(as.character(answer),split=",")[[1]]))
  	      }
        }
   	  }
   	  if(answer[1]=="None"){
  	  	break()
  	  	}
   	  if(length(grep("<",answer))>0){
  	    answer=gsub("<","",answer)
		answer=as.numeric(c(rownames(datoutnew)[datoutnew$Z.K<as.numeric(answer)]))
  	    } else {	
  	      if(length(grep(",",answer))>0){
  	      answer=c(as.numeric(strsplit(as.character(answer),split=",")[[1]]))
  	      }
        }
     if(length(answer)>0){
       K.out=c(datoutnew[is.element(rownames(datoutnew),answer),is.element(dimnames(datoutnew)[[2]],"K")])
	   Z.Kout=c(datoutnew[is.element(rownames(datoutnew),answer),is.element(dimnames(datoutnew)[[2]],"Z.K")])
	   whichround=c(rep(Round,length(answer)))
	   if(verbose==FALSE){
	     Outliers=data.frame(rbind(Outliers,data.frame(sampleinfo1[as.numeric(answer),],whichround,K.out,Z.Kout,as.numeric(answer))))
	     } else {
	     Z.Cout=c(datoutnew[is.element(rownames(datoutnew),answer),is.element(dimnames(datoutnew)[[2]],"Z.C")])
	     Z.MARout=c(datoutnew[is.element(rownames(datoutnew),answer),is.element(dimnames(datoutnew)[[2]],"Z.MAR")])
	     Outliers=data.frame(rbind(Outliers,data.frame(sampleinfo1[as.numeric(answer),],whichround,K.out,Z.Kout,Z.Cout,Z.MARout,as.numeric(answer))))
		 }
	   } else {
		Outliers=Outliers
		indexgpnew=indices1[[r]]
		}
     gpoutliers=c(gpoutliers,answer)
     indexgpnew=setdiff(indices1[[r]],as.numeric(gpoutliers)+skip1)
	 Round=Round+1
     print(paste("Round ",Round,sep=""))
	 collectGarbage()
	 datoutnew=MakePlots1(datexpr2=datExprT,IAC2=IAC,adj2=A.IAC,sampleinfo2=sampleinfo1,indexgp2=indexgpnew,group2=groups1[r],subgroup2=subgroup1,btrait2=btrait1,trait2=trait1,asfactors2=asfactors1,whichround=Round,fitmodels2=fitmodels1,whichmodel2=whichmodel1,whichfit2=whichfit1,exportfigures2=exportfigures1,verbose2=verbose)
     print("Please enter the row number(s) from the output of the sample(s) you would like to remove (comma separated) or type 'None'",quote=F)
     print("Alternatively, please enter a sample connectivity threshold for removing outliers. For example, type '<-2' to remove all samples with Z.K less than -2",quote=F)
     print(paste("Note: ",length(gpoutliers)," / ",length(indices1[[r]])," (",signif(length(gpoutliers)/length(indices1[[r]]),3)*100,"%) of samples from this group have been removed",sep=""))
     answer=readline(prompt="Response: ")
     } ## end while(answer!="None"){
  
  if(!exists("indexgpnew")){
  	indexgpnew=indices1[[r]]
  	}
  datExprOut=data.frame(datExprT[,indexgpnew])
  if(normalize1==TRUE){
  	datExprOut=as.matrix(datExprOut)
  	print(paste("Quantile normalizing group ",groups1[r],sep=""))
  	datExprNorm=normalize.quantiles(datExprOut)
  	datExprNorm=data.frame(datExprT[,c(1:skip1)],datExprNorm)
  	dimnames(datExprNorm)[[2]]=c(dimnames(datExprT)[[2]][1:skip1],dimnames(datExprT)[[2]][indexgpnew])
	write.table(datExprNorm,file=paste(projectname1,"_",groups1[r],"_",length(datExprNorm[1,])-skip1,"_Qnorm.csv",sep=""),sep=",",row.names=F,col.names=T)
  	datExprNoOutliers=data.frame(datExprT[,indexgpnew])
	datExprNoOutliers=data.frame(datExprT[,c(1:skip1)],datExprNoOutliers)
	write.table(datExprNoOutliers,file=paste(projectname1,"_",groups1[r],"_",length(datExprNoOutliers[1,])-skip1,"_outliers_removed.csv",sep=""),sep=",",row.names=F,col.names=T) 
	datoutnew=MakePlots1(datexpr2=datExprNorm,IAC2=NULL,adj2=NULL,sampleinfo2=sampleinfo1,indexgp2=indexgpnew,group2=groups1[r],subgroup2=subgroup1,btrait2=btrait1,trait2=trait1,asfactors2=asfactors1,whichround="Qnorm",fitmodels2=fitmodels1,whichmodel2=whichmodel1,whichfit2=whichfit1,exportfigures2=exportfigures1,verbose2=verbose)
    if(whichfit1=="mean"){
	  print("WARNING: Regression using mean expression as outcome is not meaningful following quantile normalization!")
	  }
	} else {
	  datExprOut=data.frame(datExprT[,indexgpnew])
	  datExprOut=data.frame(datExprT[,c(1:skip1)],datExprOut)
	  write.table(datExprOut,file=paste(projectname1,"_",groups1[r],"_",length(datExprOut[1,])-skip1,"_outliers_removed.csv",sep=""),sep=",",row.names=F,col.names=T) 
	  }
  	print("Please enter the row number from the output below that you would like to use for batch normalization or type 'None'",quote=F)
  	coldf=data.frame(dimnames(datoutnew)[[2]])
  	colnames(coldf)[1]="Sample info"
  	print(coldf)
    answer=readline(prompt="Response: ")
    while(answer!="None"){
      while(answer==""){
  	    answer=readline(prompt="Response: ")
  	    }
  	  while(length(intersect(c(1:dim(datoutnew)[[2]]),as.numeric(answer)))==0&answer!="None"){
  	    print("Error! Response must match one row number in output above. Please try again or type 'None'.")
  	    answer=readline(prompt="Response: ")
  	    }
  	  if(answer[1]=="None"){
  	  	break()
  	  	}
  	  print(paste(dimnames(datoutnew)[[2]][as.numeric(answer)]," will be used. Is this correct? Type y/n.",sep=""))
  	  tempanswer=readline(prompt="Response: ")
  	  while(length(tempanswer[is.element(tempanswer,"y")])!=1&length(tempanswer[is.element(tempanswer,"n")])!=1){
  	  	tempanswer=readline(prompt="Response: ")
  	  	}
  	  if(tempanswer=="n"){
  	  	print("Please enter a new row number or type 'None'.")
  	  	answer=readline(prompt="Response: ")
  	  	while(answer==""){
  	      answer=readline(prompt="Response: ")
  	      }
  	    while(length(intersect(c(1:dim(datoutnew)[[2]]),as.numeric(answer)))==0&answer!="None"){
  	      print("Error! Response must match one row number in output above. Please try again or type 'None'.")
  	      answer=readline(prompt="Response: ")
  	      }
  	  	 if(answer[1]=="None"){
  	  	   break()
  	   	   }
  	  	print(paste(dimnames(datoutnew)[[2]][as.numeric(answer)]," will be used.",sep=""))
  	  	}
  	  if(answer[1]=="None"){
  	  	   break()
  	   	   }
  	  print("Please indicate which batch(es) you would like to correct (comma separated), or type 'All' for all batches.  If no batch correction is needed, type 'None'.",quote=F)
  	  answer2=readline(prompt="Response: ")
  	  while(answer2==""){
  	    answer2=readline(prompt="Response: ")
  	    }
  	  if(answer2[1]=="None"){
  	  	 break()
  	   	 }
	  if(answer2[1]=="All"){
		answer2=unique(factor(datoutnew[,as.numeric(answer)]))
	    }
	  if(length(grep(",",answer2))>0){
  	    answer2=factor(strsplit(answer2,split=",")[[1]])
  	    }
  	  while(length(unique(factor(datoutnew[,as.numeric(answer)]))[is.element(unique(factor(datoutnew[,as.numeric(answer)])),answer2)])!=length(answer2)){
        print("Error! Please indicate which batch(es) you would like to correct (comma separated)",quote=F)
  	    answer2=readline(prompt="Response: ")
  	    while(answer2==""){
  	      answer2=readline(prompt="Response: ")
   	      }
   	    if(length(grep(",",answer2))>0){
  	      answer2=factor(strsplit(answer2,split=",")[[1]])
  	      }
		}
        if(answer2[1]=="None"){
  	  	   break()
  	   	   }
		if(answer2[1]=="All"){
			answer2=unique(factor(datoutnew[,as.numeric(answer)]))
	       }
   	  batchlength=c()
  	  for(s in c(1:length(answer2))){
  	    batchsamples=length(datoutnew[,as.numeric(answer)][datoutnew[,as.numeric(answer)]==intersect(datoutnew[,as.numeric(answer)],answer2[s])])
  	    batchlength=c(batchlength,batchsamples)
  	    if(batchsamples<2){
  	  	  print(paste("Batch ",answer2[s]," has only one sample and cannot be corrected.",sep=""))
  	  	  }
  	    }
  	  answer2=answer2[batchlength>1]
  	  tempanswer=answer
  	  answer="None"
  	  if(length(answer2)>0){
		 if(normalize1==FALSE){
           BatchInfo=data.frame(dimnames(datExprOut)[[2]][(skip1+1):length(datExprOut[1,])])
    	   BatchInfo=data.frame(BatchInfo[,1],BatchInfo[,1],rep(1,length(BatchInfo[,1])))
  	       for(t in c(1:length(answer2))){
  	         BatchInfo[,3][datoutnew[match(BatchInfo[,1],datoutnew[,2]),as.numeric(tempanswer)]==intersect(datoutnew[,as.numeric(tempanswer)],answer2[t])]=t+1
  	         }
		   if(sum(batchlength)==dim(BatchInfo)[1]){
    	   	BatchInfo[,3]=BatchInfo[,3]-1
    	   	}
		   dimnames(BatchInfo)[[2]]=c("Array_name","Sample_name","Batch")
		   batchtable=data.frame(table(BatchInfo$Batch))
		   if(length(batchtable$Freq[batchtable$Freq==1])>0){
 		     singlebatch=batchtable$Var1[batchtable$Freq==1]
		     maxbatch=as.numeric(as.character(batchtable$Var1)[which.max(batchtable$Freq)])
			 BatchInfo[,3][BatchInfo[,3]==singlebatch]=maxbatch
			 }
		   write.table(BatchInfo,file=paste(projectname1,"_",groups1[r],"_",length(datExprOut[1,])-skip1,"_ComBat_sampleinfo.txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F)
    	   } else {
    	   BatchInfo=data.frame(dimnames(datExprNorm)[[2]][(skip1+1):length(datExprNorm[1,])])
    	   BatchInfo=data.frame(BatchInfo[,1],BatchInfo[,1],rep(1,length(BatchInfo[,1])))
    	   for(t in c(1:length(answer2))){
			     BatchInfo[,3][datoutnew[match(BatchInfo[,1],datoutnew[,2]),as.numeric(tempanswer)]==intersect(datoutnew[,as.numeric(tempanswer)],answer2[t])]=t+1
				}
  	       if(sum(batchlength)==dim(BatchInfo)[1]){
    	   	BatchInfo[,3]=BatchInfo[,3]-1
    	   	}
  	       dimnames(BatchInfo)[[2]]=c("Array_name","Sample_name","Batch")
		   ## Assign any singleton to the largest batch:
		   batchtable=data.frame(table(BatchInfo$Batch))
		   if(length(batchtable$Freq[batchtable$Freq==1])>0){
 		     singlebatch=batchtable$Var1[batchtable$Freq==1]
		     maxbatch=as.numeric(as.character(batchtable$Var1)[which.max(batchtable$Freq)])
			 BatchInfo[,3][BatchInfo[,3]==singlebatch]=maxbatch
			 }
		   write.table(BatchInfo,file=paste(projectname1,"_",groups1[r],"_",length(datExprNorm[1,])-skip1,"_ComBat_sampleinfo.txt",sep=""),sep="\t",row.names=F,quote=F)
  	       }
    
  	    print("Please enter the row number(s) from the output below that you would like to use as covariates during batch normalization or type 'None'",quote=F)
  	    print("Note: ComBat will only work with categorical covariates",quote=F)
        coldf=data.frame(dimnames(sampleinfo1)[[2]])
    	colnames(coldf)[1]="Sample info"
    	print(coldf)
        answer3=readline(prompt="Response: ")
        while(answer3!="None"){
          while(answer3==""){
  	        answer3=readline(prompt="Response: ")
  	        }
  	      if(length(grep(",",answer3))>0){
  	        answer3=c(as.numeric(strsplit(as.character(answer3),split=",")[[1]]))
  	        }
  	      while(length(intersect(c(1:dim(sampleinfo1)[[2]]),as.numeric(answer3)))==0){
  	        if(answer3[1]=="None"){
  	  	      break()
  	  	      }
  	        print("Error! Response must match at least one row number in output above. Please try again or type 'None'.",quote=F)
  	        answer3=readline(prompt="Response: ")
  	        }
  	      if(answer3[1]=="None"){
  	  	    break()
  	  	    }
  	      if(length(grep(",",answer3))>0){
  	        answer3=c(as.numeric(strsplit(as.character(answer3),split=",")[[1]]))
  	        }
  	      print(paste(paste(dimnames(sampleinfo1)[[2]][as.numeric(answer3)],collapse=" & ")," will be used. Is this correct? Type y/n.",sep=""),quote=F)
  	      tempanswer=readline(prompt="Response: ")
  	      while(length(tempanswer[is.element(tempanswer,"y")])!=1&length(tempanswer[is.element(tempanswer,"n")])!=1){
  	  	    tempanswer=readline(prompt="Response: ")
  	  	    }
  	      if(tempanswer=="n"){
  	  	    print("Please enter a new row number(s) or type 'None'.")
  	  	    answer3=readline(prompt="Response: ")
  	  	    while(answer3==""){
  	          answer3=readline(prompt="Response: ")
  	          }
  	  	    if(length(grep(",",answer3))>0){
  	          answer3=c(as.numeric(strsplit(as.character(answer3),split=",")[[1]]))
  	          }
			while(length(intersect(c(1:dim(sampleinfo1)[[2]]),as.numeric(answer3)))==0){
  	          if(answer3[1]=="None"){
  	  	        break()
  	  	        }
  	          print("Error! Response must match at least one row number in output above. Please try again or type 'None'.",quote=F)
  	          answer3=readline(prompt="Response: ")
  	          }
  	          if(answer3[1]=="None"){
  	  	        break()
  	  	        }
  	  	    print(paste(paste(dimnames(sampleinfo1)[[2]][as.numeric(answer3)],collapse=" & ")," will be used.",sep=""),quote=F)  	  	    
			}
          if(answer3[1]=="None"){
  	  	    break()
  	  	    }
          covariates=paste("Covariate_",c(1:length(answer3)),sep="")
          newdimnames=c(dimnames(BatchInfo)[[2]],covariates)
          BatchInfo=data.frame(BatchInfo,sampleinfo1[match(BatchInfo[,1],sampleinfo1[,samplelabels1]),as.numeric(answer3)])
          dimnames(BatchInfo)[[2]]=newdimnames
          if(normalize1==FALSE){
			write.table(BatchInfo,file=paste(projectname1,"_",groups1[r],"_",length(datExprOut[1,])-skip1,"_ComBat_sampleinfo.txt",sep=""),sep="\t",row.names=F,quote=F)
            } else {
			  write.table(BatchInfo,file=paste(projectname1,"_",groups1[r],"_",length(datExprNorm[1,])-skip1,"_ComBat_sampleinfo.txt",sep=""),sep="\t",row.names=F,quote=F)
              }
            answer3="None"
            } ## end of while(answer3!="None"){
  	    
  	      if(normalize1==FALSE){
  	      	datComBat=try(ComBat(paste(projectname1,"_",groups1[r],"_",length(datExprOut[1,])-skip1,"_outliers_removed.csv",sep=""),paste(projectname1,"_",groups1[r],"_",length(datExprOut[1,])-skip1,"_ComBat_sampleinfo.txt",sep=""),type="csv",write=F,par.prior=T,filter=F,skip=skip1,prior.plots=T),silent=TRUE)
  	      	datComBat=data.frame(datExprT[,c(1:skip1)],datComBat)
			} else {
## Need to ensure that quantile normalization has not introduced genes with 0 variance:
				excludevec2=rep(1,length(datExprNorm[,1]))
                variance=apply(datExprNorm[,c((skip1+1):dim(datExprNorm)[2])],1,var,na.rm=T)
				restvar=variance==0
				excludevec2[restvar]=0
				excludevec2=as.logical(excludevec2)
 			    datExprNorm=datExprNorm[excludevec2,]
				collectGarbage()

				if(length(excludevec2[excludevec2])>0){
					print(paste("Note: ",length(excludevec2[!excludevec2])," probe sets had 0 variance and were excluded",sep=""))
					write.table(datExprNorm,file=paste(projectname1,"_",groups1[r],"_",length(datExprNorm[1,])-skip1,"_Qnorm_0var_removed.csv",sep=""),sep=",",row.names=F,col.names=T)
					datComBat=try(ComBat(paste(projectname1,"_",groups1[r],"_",length(datExprNorm[1,])-skip1,"_Qnorm_0var_removed.csv",sep=""),paste(projectname1,"_",groups1[r],"_",length(datExprNorm[1,])-skip1,"_ComBat_sampleinfo.txt",sep=""),type="csv",write=F,par.prior=T,filter=F,skip=skip1,prior.plots=T),silent=TRUE)
					} else {
					datComBat=try(ComBat(paste(projectname1,"_",groups1[r],"_",length(datExprNorm[1,])-skip1,"_Qnorm.csv",sep=""),paste(projectname1,"_",groups1[r],"_",length(datExprNorm[1,])-skip1,"_ComBat_sampleinfo.txt",sep=""),type="csv",write=F,par.prior=T,filter=F,skip=skip1,prior.plots=T),silent=TRUE)
				}
			}
  	      if(is(datComBat)[1]=="try-error"){
		    print("WARNING! ComBat failed.  The following error was returned:")
			print(datComBat[1])
			} else {
			print(paste("ComBat normalized group ",groups1[r],sep=""))
			}
		  if(replacenegs1==TRUE){
		    print("Replacing negative expression values with the median for the probe set...",quote=F)
		    negatives=datComBat<0
	        totalnegs=length(negatives[negatives])
	        if(exists("MedianNegs")) rm(MedianNegs);
	          MedianNegs=function(datE){
	     	  dat2E=as.matrix(datE)
	     	  negatives=dat2E<0
	     	  for(i in c(1:length(dat2E[,1]))){
	     		if(length(dat2E[negatives[i,][negatives[i,]]])>0){
			      negatives1=dat2E[i,]<0
	   		      dat2E[i,][negatives1]=median(dat2E[i,])
	   		      }
	   		    }
	   	      dat2E
	   	      }
		datComBat=MedianNegs(datE=datComBat)
	        print(paste(signif(as.numeric(length(negatives[negatives]))/as.numeric(length(negatives))*100,3),"% were negative and replaced",sep=""))
	        }
		MeanSampleComBat=apply(datComBat,2,mean,na.rm=T)
		collectGarbage()
				if(exists("excludevec2")){
					if(length(excludevec2[excludevec2])>0){
						datComBat=data.frame(datExprT[excludevec2,c(1:skip1)],datComBat)
					} else {
						datComBat=data.frame(datExprT[,c(1:skip1)],datComBat)
					}
				}
		dimnames(datComBat)[[2]][1:skip1]=dimnames(datExprT)[[2]][1:skip1]
		datoutnew=MakePlots1(datexpr2=datComBat,IAC2=NULL,adj2=NULL,sampleinfo2=sampleinfo1,indexgp2=indexgpnew,group2=groups1[r],subgroup2=subgroup1,btrait2=btrait1,trait2=trait1,asfactors2=asfactors1,whichround="ComBat",fitmodels2=fitmodels1,whichmodel2=whichmodel1,whichfit2=whichfit1,exportfigures2=exportfigures1,verbose2=verbose)
		if(normalize1==TRUE&whichfit1=="mean"){
		  print("WARNING: Regression using mean expression as outcome is not meaningful following quantile normalization!")
		  }
		write.table(datComBat,file=paste(projectname1,"_",groups1[r],"_",length(datComBat[1,])-skip1,"_ComBat.csv",sep=""),sep=",",row.names=F)
	    } ## end if(length(answer2)>0){
      } ## end while(answer!="None"){
	write.table(attr(datoutnew,"AdjMat"),file=paste(groups1[r],"_final_adjmat.csv",sep=""),sep=",",row.names=T,col.names=T)
#tstamp=format(Sys.time(), "%X")
#tstamp=gsub(":","-",tstamp)
	datoutnew=datoutnew[order(as.character(datoutnew[,2])),]
	write.table(datoutnew,file=paste(projectname1,"_",groups1[r],"_",length(datoutnew[,1]),"_final_sample_metrics.csv",sep=""),sep=",",row.names=F,col.names=T)
	MetricsDF=data.frame(rbind(MetricsDF,attr(datoutnew,"MetricsOut")))
    TempRelMetrics=data.frame(attr(datoutnew,"MetricsOut"))
    RelFx=function(x){
      ((x-x[1])/x[1])*100
      }
    if((Round>1&is.numeric(Round))|length(intersect(Round,c("Qnorm","ComBat")))>0){
	  TempRelMetrics=data.frame(TempRelMetrics[,c(1:3)],apply(TempRelMetrics[,c(4,6:11)],2,RelFx))
      colnames(TempRelMetrics)=gsub("_"," ",colnames(TempRelMetrics))
	  pdf(file=paste(groups1[r],"_change_by_round.pdf",sep=""))
      matplot(TempRelMetrics[,4:10],type="l",lty=1,col=c("black","red","blue","green","turquoise","orange","purple"),lwd=2,ylab="Percent change",xlab="Round",main=paste(groups1[r]," sample network metrics",sep=""),xaxt="n",cex.lab=1.4,cex.main=1.8)
      axis(side=1,labels=TempRelMetrics$Round,at=c(1:length(TempRelMetrics[,1])))
      legend("topleft",legend=colnames(TempRelMetrics[4:10]),col=c("black","red","blue","green","turquoise","orange","purple"),lwd=2)
      dev.off()
	  }
	if(length(groups1)>1){
      print(paste(groups1[r]," processing complete. Hit enter to proceed to next group.",sep=""))
      answer4=readline(prompt="Press enter: ")
      }
  setwd(SNrootDir)
  } ## end for(r in c(1:length(indices1)[1])){
print("Exporting outlier summary...")
dimnames(Outliers)[[2]][is.element(dimnames(Outliers)[[2]],"K.out")]="K"
dimnames(Outliers)[[2]][is.element(dimnames(Outliers)[[2]],"Z.Kout")]="Z.K"
dimnames(Outliers)[[2]][dim(Outliers)[2]]="RowIndex"
if(verbose==TRUE){
  dimnames(Outliers)[[2]][is.element(dimnames(Outliers)[[2]],"Z.Cout")]="Z.C"
  dimnames(Outliers)[[2]][is.element(dimnames(Outliers)[[2]],"Z.MARout")]="Z.MAR"
  dimnames(Outliers)[[2]][dim(Outliers)[2]]="RowIndex"
  }
#tstamp=format(Sys.time(), "%X")
#tstamp=gsub(":","-",tstamp)
write.table(Outliers,file=paste(projectname1,"_",colnames(sampleinfo1)[grouplabels1],"_outliers","_",tstamp1,".csv",sep=""),sep=",",row.names=F,col.names=T)
#tstamp=format(Sys.time(), "%X")
#tstamp=gsub(":","-",tstamp)
write.table(MetricsDF,file=paste(projectname1,"_",colnames(sampleinfo1)[grouplabels1],"_dataset_metrics","_",tstamp1,".csv",sep=""),sep=",",row.names=F,col.names=T)
if(length(excludevec[excludevec])>0){
  write.table(datExprOrig[!excludevec,1:skip1],file=paste(projectname1,"_",colnames(sampleinfo1)[grouplabels1],"_excluded_genes_with_0_variance_",tstamp1,".csv",sep=""),sep=",",row.names=F,col.names=T)
  }
}
		  
		  
		  
## "ComBat" function below (by Evan Johnson; see http://statistics.byu.edu/johnson/ComBat/Abstract.html for more information)
## Note: 'expression_xls' is the expression index file (e.g. outputted by dChip); 
## 'sample_info_file' is a tab-delimited text file containing the colums: Array  name, sample name, Batch, and any other covariates to be included in the modeling; 
## 'type' currently supports two data file types 'txt' for a tab-delimited text file and 'csv' for an Excel .csv file (sometimes R handles the .csv file better, so use this if you have problems with a .txt file!); 
## 'write' if 'T' ComBat writes adjusted data to a file, and if 'F' and ComBat outputs the adjusted data matrix if 'F' (so assign it to an object! i.e. NewData <- ComBat('my expression.xls','Sample info file.txt', write=F));
## 'covariates=all' will use all of the columns in your sample info file in the modeling (except array/sample name), if you only want use a some of the columns in your sample info file, specify these columns here as a vector (you must include the Batch column in this list); 
## 'par.prior' if 'T' uses the parametric adjustments, if 'F' uses the nonparametric adjustments--if you are unsure what to use, try the parametric adjustments (they run faster) and check the plots to see if these priors are reasonable; 
## 'filter=value' filters the genes with absent calls in > 1-value of the samples. The defaut here (as well as in dchip) is .8. Filter if you can as the EB adjustments work better after filtering. Filter must be numeric if your expression index file contains presence/absence calls (but you can set it >1 if you don't want to filter any genes) and must be 'F' if your data doesn't have presence/absence calls; 
## 'skip' is the number of columns that contain probe names and gene information, so 'skip=5' implies the first expression values are in column 6;
## 'prior.plots' if true will give prior plots with black as a kernal estimate of the empirical batch effect density and red as the parametric estimate. 
		  
		  ComBat <- function(expression_xls, sample_info_file, type='txt', write=T, covariates='all', par.prior=T, filter=.8, skip=5, prior.plots=T){
			  
			  cat('Reading Sample Information File\n')
			  saminfo <- read.table(sample_info_file, header=T, sep='\t',comment.char='')
			  if(sum(colnames(saminfo)=="Batch")!=1){return('ERROR: Sample Information File does not have a Batch column!')}
			  
			  cat('Reading Expression Data File\n')
			  if(type=='csv'){
				  dat <- read.csv(expression_xls,header=T,comment.char='',as.is=T)
				  dat <- dat[,trim.dat(dat)]
				  colnames(dat)=scan(expression_xls,what='character',nlines=1,sep=',',quiet=T)[1:ncol(dat)]
			  }else{
				  dat <- read.table(expression_xls,header=T,comment.char='',fill=T,sep='\t', as.is=T)
				  dat <- dat[,trim.dat(dat)]
				  colnames(dat)=scan(expression_xls,what='character',nlines=1,sep='\t',quiet=T)[1:ncol(dat)]
			  }
			  if (skip>0){geneinfo <- as.matrix(dat[,1:skip]); dat <- dat[,-c(1:skip)]}else{geneinfo=NULL}
			  
			  if(filter){
				  ngenes <- nrow(dat)
				  col <- ncol(dat)/2
				  present <- apply(dat, 1, filter.absent, filter)
				  dat <- dat[present, -(2*(1:col))]
				  if (skip>0){geneinfo <- geneinfo[present,]}
				  cat('Filtered genes absent in more than',filter,'of samples. Genes remaining:',nrow(dat),'; Genes filtered:',ngenes-nrow(dat),'\n')
			  }
			  
			  if(any(apply(dat,2,mode)!='numeric')){return('ERROR: Array expression columns contain non-numeric values! (Check your .xls file for non-numeric values and if this is not the problem, make a .csv file and use the type=csv option)')}
			  
			  tmp <- match(colnames(dat),saminfo[,1])
			  if(any(is.na(tmp))){return('ERROR: Sample Information File and Data Array Names are not the same!')}
			  tmp1 <- match(saminfo[,1],colnames(dat))
			  saminfo <- saminfo[tmp1[!is.na(tmp1)],]		
			  
			  if(any(covariates != 'all')){saminfo <- saminfo[,c(1:2,covariates)]}
			  design <- design.mat(saminfo)	
			  
			  batches <- list.batch(saminfo)
			  n.batch <- length(batches)
			  n.batches <- sapply(batches, length)
			  n.array <- sum(n.batches)
			  
##Standardize Data across genes
			  cat('Standardizing Data across genes\n')
			  B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat)) #Standarization Model
			  grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
			  var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)
			  
			  stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
			  if(!is.null(design)){tmp <- design;tmp[,c(1:n.batch)] <- 0;stand.mean <- stand.mean+t(tmp%*%B.hat)}	
			  s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))
			  
##Get regression batch effect parameters
			  cat("Fitting L/S model and finding priors\n")
			  batch.design <- design[,1:n.batch]
			  gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))
			  delta.hat <- NULL
			  for (i in batches){
				  delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var))
			  }
			  
##Find Priors
			  gamma.bar <- apply(gamma.hat, 1, mean)
			  t2 <- apply(gamma.hat, 1, var)
			  a.prior <- apply(delta.hat, 1, aprior)
			  b.prior <- apply(delta.hat, 1, bprior)
			  
			  
##Plot empirical and parametric priors
			  
			  if (prior.plots){
				  par(mfrow=c(2,2))
				  tmp <- density(gamma.hat[1,])
				  plot(tmp,  type='l', main="Density Plot")
				  xx <- seq(min(tmp$x), max(tmp$x), length=100)
				  lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
				  qqnorm(gamma.hat[1,])	
				  qqline(gamma.hat[1,], col=2)	
				  
				  tmp <- density(delta.hat[1,])
				  invgam <- 1/rgamma(ncol(delta.hat),a.prior[1],b.prior[1])
				  tmp1 <- density(invgam)
				  plot(tmp,  typ='l', main="Density Plot", ylim=c(0,max(tmp$y,tmp1$y)))
				  lines(tmp1, col=2)
				  qqplot(delta.hat[1,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')	
				  lines(c(0,max(invgam)),c(0,max(invgam)),col=2)	
				  title('Q-Q Plot')
			  }
			  
##Find EB batch adjustments
			  
			  gamma.star <- delta.star <- NULL
			  if(par.prior){
				  cat("Finding parametric adjustments\n")
				  for (i in 1:n.batch){
					  temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
					  gamma.star <- rbind(gamma.star,temp[1,])
					  delta.star <- rbind(delta.star,temp[2,])
				  }
			  }else{
				  cat("Finding nonparametric adjustments\n")
				  for (i in 1:n.batch){
					  temp <- int.eprior(as.matrix(s.data[,batches[[i]]]),gamma.hat[i,],delta.hat[i,])
					  gamma.star <- rbind(gamma.star,temp[1,])
					  delta.star <- rbind(delta.star,temp[2,])
				  }
			  }
			  
			  
### Normalize the Data ###
			  cat("Adjusting the Data\n")
			  
			  bayesdata <- s.data
			  j <- 1
			  for (i in batches){
				  bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
				  j <- j+1
			  }
			  
			  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean
#if(skip>0){bayesdata <- cbind(geneinfo,bayesdata)}
			  
			  
			  if(write){
				  output_file <- paste('Adjusted',expression_xls,'.xls',sep='_')
				  cat(c(colnames(geneinfo),colnames(dat),'\n'),file=output_file,sep='\t')
				  suppressWarnings(write.table(cbind(geneinfo,formatC(as.matrix(bayesdata), format = "f")), file=output_file, sep="\t", quote=F,row.names=F,col.names=F,append=T))
				  cat("Adjusted data saved in file:",output_file,"\n")
			  }else{return(bayesdata)}
		  }
		  
# filters data based on presence/absence call
		  filter.absent <- function(x,pct){
			  present <- T
			  col <- length(x)/2
			  pct.absent <- (sum(x[2*(1:col)]=="A") + sum(x[2*(1:col)]=="M"))/col
			  if(pct.absent > pct){present <- F}
			  present
		  }
		  
# Next two functions make the design matrix (X) from the sample info file 
		  build.design <- function(vec, des=NULL, start=2){
			  tmp <- matrix(0,length(vec),nlevels(vec)-start+1)
			  for (i in 1:ncol(tmp)){tmp[,i] <- vec==levels(vec)[i+start-1]}
			  cbind(des,tmp)
		  }
		  
		  design.mat <- function(saminfo){
			  tmp <- which(colnames(saminfo) == 'Batch')
			  tmp1 <- as.factor(saminfo[,tmp])
			  cat("Found",nlevels(tmp1),'batches\n')
			  design <- build.design(tmp1,start=1)
			  ncov <- ncol(as.matrix(saminfo[,-c(1:2,tmp)]))
			  cat("Found",ncov,'covariate(s)\n')
			  if(ncov>0){
				  for (j in 1:ncov){
					  tmp1 <- as.factor(as.matrix(saminfo[,-c(1:2,tmp)])[,j])
					  design <- build.design(tmp1,des=design)
				  }
			  }
			  design
		  }
		  
# Makes a list with elements pointing to which array belongs to which batch
		  list.batch <- function(saminfo){
			  tmp1 <- as.factor(saminfo[,which(colnames(saminfo) == 'Batch')])
			  batches <- NULL
			  for (i in 1:nlevels(tmp1)){batches <- append(batches, list((1:length(tmp1))[tmp1==levels(tmp1)[i]]))}
			  batches
		  }
		  
# Trims the data of extra columns, note your array names cannot be named 'X' or start with 'X.'
		  trim.dat <- function(dat){
			  tmp <- strsplit(colnames(dat),'\\.')
			  tr <- NULL
			  for (i in 1:length(tmp)){tr <- c(tr,tmp[[i]][1]!='X')}
			  tr
		  }
		  
# Following four find empirical hyper-prior values
		  aprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (2*s2+m^2)/s2}
		  bprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (m*s2+m^3)/s2}
		  postmean <- function(g.hat,g.bar,n,d.star,t2){(n*t2*g.hat+d.star*g.bar)/(n*t2+d.star)}
		  postvar <- function(sum2,n,a,b){(.5*sum2+b)/(n/2+a-1)}
		  
		  
# Pass in entire data set, the design matrix for the entire data, the batch means, the batch variances, priors (m, t2, a, b), columns of the data  matrix for the batch. Uses the EM to find the parametric batch adjustments
		  
		  it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
			  n <- ncol(sdat)
			  g.old <- g.hat
			  d.old <- d.hat
			  change <- 1
			  count <- 0
			  while(change>conv){
				  g.new <- postmean(g.hat,g.bar,n,d.old,t2)
				  sum2 <- apply((sdat-g.new%*%t(rep(1,n)))^2, 1, sum)
				  d.new <- postvar(sum2,n,a,b)
				  change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
				  g.old <- g.new
				  d.old <- d.new
				  count <- count+1
			  }
#cat("This batch took", count, "iterations until convergence\n")
			  adjust <- rbind(g.new, d.new)
			  rownames(adjust) <- c("g.star","d.star")
			  adjust
		  }
		  
#likelihood function used below
		  L <- function(x,g.hat,d.hat){prod(dnorm(x,g.hat,sqrt(d.hat)))}
		  
# Monte Carlo integration function to find the nonparametric adjustments
		  int.eprior <- function(sdat,g.hat,d.hat){
			  g.star <- d.star <- NULL
			  r <- nrow(sdat)
			  for(i in 1:r){
				  g <- g.hat[-i]
				  d <- d.hat[-i]		
				  x <- sdat[i,]
				  n <- length(x)
				  j <- numeric(n)+1
				  dat <- matrix(x,length(g),n,byrow=T)
				  resid2 <- (dat-g)^2
				  sum2 <- resid2%*%j
				  LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
				  g.star <- c(g.star,sum(g*LH)/sum(LH))
				  d.star <- c(d.star,sum(d*LH)/sum(LH))
#if(i%%1000==0){cat(i,'\n')}
			  }
			  adjust <- rbind(g.star,d.star)
			  rownames(adjust) <- c("g.star","d.star")
			  adjust	
		  } 
		  
		  
		  
# The function ModulePrinComps1 (by Steve Horvath) finds the first principal component (eigengene) in each 
# module defined by the colors of the input vector "couleur" (Pardon my French).
# It also reports the variances explained by the first 5 principal components.
# This requires the R library impute
# The output is a list with 2 components: 
# 1) a data frame of module eigengenes (PCs), 
# 2) a data frame that lists the percent variance explained by the first 5 PCs of a module
# Options: if removeGrey=T, then no output is generated for the grey module.
# Recall that grey often denotes genes outside proper modules. 
		  if (exists("ModulePrinComps1" ) ) rm(ModulePrinComps1);
		  ModulePrinComps1=function(datexpr,couleur,removeGrey=F, FiveComponents=F) {
			  modlevels=levels(factor(couleur))
			  if ( removeGrey ) modlevels=setdiff(modlevels, c("grey") ); 
			  if (FiveComponents ) {print("To speed up the calculation, we only compute the five principal components of each module.  Therefore, the estimate of the proportion of variance explained is no longer accurate. If you want an accurate estimate of the proportion of var explained, please choose  the option FiveComponents=F ")  ;} 
			  PrinComps=data.frame(matrix(666,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
			  varexplained= data.frame(matrix(666,nrow= 5,ncol= length(modlevels)))
			  names(PrinComps)=paste("PC",modlevels,sep="")
			  for(i in c(1:length(modlevels)) ){
				  print(i)   
				  modulename    = modlevels[i]
				  restrict1= as.character(couleur)== modulename
# in the following, rows are genes and columns are samples
				  datModule=t(datexpr[, restrict1])
				  is.saved = FALSE;
				  if (exists(".Random.seed"))
				  {
					  saved.seed = .Random.seed;
					  is.saved = TRUE;
				  }
				  datModule=impute.knn(as.matrix(datModule))
				  if ( (length(datModule)==3) && 
					  (!is.null(names(datModule))) &&
					  (names(datModule)[1]=="data") )
				  datModule = datModule$data;
				  if (is.saved) .Random.seed = saved.seed;
				  datModule=t(scale(t(datModule)))
				  if (FiveComponents ) { svd1 =svd(datModule, nu = 5, nv = 5) } else {svd1=svd(datModule)}
				  mtitle=paste("PCs of ", modulename," module", sep="")
				  varexplained[,i]= (svd1$d[1:5])^2/sum(svd1$d^2)
# this is the first principal component
				  pc1=svd1$v[,1]
				  signh1=sign(sum(cor(pc1,  t(datModule))))
				  if (signh1 != 0)  pc1=signh1* pc1
				  PrinComps[,i]= pc1
			  }
			  list(PrinComps=PrinComps, varexplained=varexplained)
		  }
		  
		  
## END

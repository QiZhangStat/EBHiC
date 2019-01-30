###' The goal of the demo is to analyze a simulated sample dataset using Fit-Hi-C and EBHiC, 
###' and compare their accuracy of the top candidates and AUPRC.
###'  	
rm(list=ls())
library(parallel)
library(data.table)
library(FitHiC) ### used only for comparison		
library(PRROC) ### used only for calculating AUPRC
source('./nphic.r')

ncores=1 ### number of cores EBHiC uses
options(mc.cores=ncores,mc.preschedule=F)

### run FitHiC
outdir = './sample_data/'

	fragsfile  = paste0(outdir,"fragmentLists.gz",sep="")
	intersfile = paste(outdir,"contactCounts.gz",sep="")
    ## FitHiC method will generate output into the "outdir" directory	
	FitHiC(fragsfile, intersfile, outdir, libname="Simulation", distUpThres=100, distLowThres=0)
	
### run NebHiC
## number of breaks for spline basis
nbreaks    = 101
## the maximum number of iterations
iter.max   = 100
## initial values
init       = NULL
## the maximum error control for convergence
error.max  = 1e-5

	## Set up output file names. In this demo, we use the FitHi-C output file as the input file
	input.full.name = paste0(outdir, "Simulation.spline_pass2.significances.txt.gz",sep="")
	chr.data  = read.table(gzfile(input.full.name), header=T)
	## The wrapper for EBHiC method 
	out.ebhic = ebhicWrapper(dt=chr.data, nbreaks=nbreaks, ncores=ncores, iter.max=iter.max, init=init, error.max=error.max)  
	rm(chr.data)
	save(out.ebhic, file=paste0(outdir,'ebhic.rdata'), version = NULL)
	#data.csv = out.ebhic[['output']]
	write.table(out.ebhic[['output']],file=paste0(outdir,'ebhic.csv'),sep=",",col.names = TRUE,row.names=FALSE)
    rm(out.ebhic)

 


### ACCURACY plot shows that NebHiC is superior compared with FitHiC as the proportion of number of pairs selected is small 

	## read the outputs from NebHiC
	dthic = read.csv(paste0(outdir,'ebhic.csv'),header=T)
	## delta is the true peak status of the simulated data.
	delta = read.table(paste0(outdir,'truth.txt'),header=F)$V1
	## rank the p-value of FitHiC by ascending order and get the positions of each value
    dthic$rankfit = rank(dthic$pval_fithic,ties.method='min')
	## rank the p-value of NebHiC by ascending order and get the positions of each value
    dthic$ranknp  = rank(dthic$locfdr,ties.method='min')
    np = nrow(dthic)
	## Set up different thresholds of proportion for the p-values
    thr.prop = c(0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1)
    dtplot = t(sapply(thr.prop, function(thr) c(mean(delta[dthic$ranknp<=thr*np]),mean(delta[dthic$rankfit<=thr*np]))))
    colnames(dtplot) = c('np','fit')
	
	## Accuracy Plot
    legd = c('NebHiC','FitHiC')
pdf(paste0(outdir,"acc.pdf"))
	plot(0.01,0.01,type='n',log='x',xlim=c(min(thr.prop),max(thr.prop)),ylim=c(0,1),main='Accuracy plot',xlab="Proportion Selected",ylab="Accuracy")
    legend(median(thr.prop),0.8,legd,col=1:2,lty=1:2,lwd=5,cex=2)
    lines(thr.prop,dtplot[,'np'], col=1,lty=1)
    lines(thr.prop,dtplot[,'fit'],col=2,lty=2)
dev.off()


### Calculate the area under precision recall curve (AUPRC)
	PRC_input = data.frame(pval_fithic=dthic$pval_fithic, locfdr=dthic$locfdr, delta=delta)
	## rank the p-value of FitHiC by ascending order and get the positions of each value
	PRC_curve_fithic = pr.curve(scores.class0=rank(-PRC_input$pval_fithic), weights.class0=PRC_input$delta, curve=T)
	## rank the p-value of FitHiC by ascending order and get the positions of each value	
	PRC_curve_locfdr = pr.curve(scores.class0=rank(-PRC_input$locfdr), weights.class0=PRC_input$delta, curve=T)
  

### print the AUPRC 
print(c(PRC_curve_fithic$auc.integral,PRC_curve_locfdr$auc.integral))


	  

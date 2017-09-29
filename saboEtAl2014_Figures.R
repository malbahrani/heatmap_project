source('analysisEnvironment.R')

# Myc promoter binding events in EuMyc
	mycEuMycPromPeaks= peaksRefSimple$EuMyc[1:5]
	for(i in 1:5) {
		res= countOverlaps(mycEuMycPromPeaks[[i]], promoterGR)
		mycEuMycPromPeaks[[i]]= mycEuMycPromPeaks[[i]][res>0]
	}

# Myc distal binding events in EuMyc
	mycEuMycDistalPeaks= peaksRefSimple$EuMyc[1:5]
	for(i in 1:5) {
		res= countOverlaps(mycEuMycDistalPeaks[[i]], promoterGR)
		mycEuMycDistalPeaks[[i]]= mycEuMycDistalPeaks[[i]][res==0]
	}

# Fig 1a: heatmap of promoter Myc peaks and additional marks
	# myc promoter peaks in either C, P or one of the T
	mycEuMycPromUnion= mycEuMycPromPeaks[[1]]
	for(i in 2:length(mycEuMycPromPeaks)) mycEuMycPromUnion= union(mycEuMycPromUnion, mycEuMycPromPeaks[[i]]) # 13628
	mycEuMycPromUnion6Kb= GRsetwidth(mycEuMycPromUnion, newWidth=6000) # 3Kb upstream and downstream the mid point of each region (negative starts are corrected)
	mycEuMycPromUnion6KbChr1= mycEuMycPromUnion6Kb[seqnames(mycEuMycPromUnion6Kb)=='chr1']
	grl= c(peaksRefSimple$EuMyc[c(1:10, 16:20, 11:15, 21:25)], CGI=CGIgr) # Myc, Pol2, K4me3, K4me1, K27ac
	EuMycPromData= GRheatmapData(grl= grl, refgr = mycEuMycPromUnion6KbChr1, type= rep('gr', length(grl)), nbins= 10, txdb = txdb)
	png('./figures/Fig1a.png', 3000, 2500, res=400)
		GRheatmapPlot(matList= EuMycPromData, rowLab = FALSE, colLab = TRUE, clusterInds = 1:length(EuMycPromData))
	dev.off()

# Fig 1b: heatmap of distal Myc peaks and additional marks
	# myc distal peaks in either C, P or one of the T
	mycEuMycDistalUnion= mycEuMycDistalPeaks[[1]]
	for(i in 2:length(mycEuMycDistalPeaks)) mycEuMycDistalUnion= union(mycEuMycDistalUnion, mycEuMycDistalPeaks[[i]]) # 36048
	mycEuMycDistalUnion6Kb= GRsetwidth(mycEuMycDistalUnion, newWidth=6000) # 3Kb upstream and downstream the mid point of each region (negative starts are corrected)
	mycEuMycDistalUnion6KbChr1= mycEuMycDistalUnion6Kb[seqnames(mycEuMycDistalUnion6Kb)=='chr1']
	EuMycDistalData= GRheatmapData(grl= grl, refgr = mycEuMycDistalUnion6KbChr1, type= rep('gr', length(grl)), nbins= 10, txdb = txdb)
	png('./figures/Fig1b.png', 3000, 2500, res=400)
		GRheatmapPlot(matList= EuMycDistalData, rowLab = FALSE, colLab = TRUE, clusterInds = 1:length(EuMycDistalData))
	dev.off()

# Fig 1c: heatmap of promoter not Myc bound and additional marks
	# myc unbound promoters for each condition
	EuMycUnboundProm= list()
	for(cond in names(peaksRefSimple$EuMyc)[1:5]) {
		res= countOverlaps(promoterGR, peaksRefSimple$EuMyc[[cond]])
		EuMycUnboundProm[[cond]]= promoterGR[res==0]
	}
	# promoter not Myc bound in any condition
	EuMycUnboundPint= EuMycUnboundProm[[1]]
	for(i in 2:length(EuMycUnboundProm)) EuMycUnboundPint= intersect(EuMycUnboundPint, EuMycUnboundProm[[i]]) # 16753
	EuMycUnboundPint6Kb= GRsetwidth(EuMycUnboundPint, newWidth=6000) # 3Kb upstream and downstream the mid point of each region (negative starts are corrected)
	EuMycUnboundPint6KbChr1= EuMycUnboundPint6Kb[seqnames(EuMycUnboundPint6Kb)=='chr1'] # 1512
	grl= c(peaksRefSimple$EuMyc[c(1:10, 16:20, 11:15, 21:25)], CGI=CGIgr) # Myc, Pol2, K4me3, K4me1, K27ac
	EuMycUnboundPromData= GRheatmapData(grl= grl, refgr = EuMycUnboundPint6KbChr1, type= rep('gr', length(grl)), nbins= 10, txdb = txdb)
	png('./figures/Fig2c.png', 3000, 2500, res=400)
		GRheatmapPlot(matList= EuMycUnboundPromData, rowLab = FALSE, colLab = TRUE, clusterInds = 1:length(EuMycUnboundPromData))
	dev.off()


# Fig 1d: heatmap of enhancer sites not Myc bound and additional marks
	# enhancers in at least on EuMyc condition
	EuMycEnhUnion= enhancerList$EuMyc[[1]]
	for(i in 2:length(enhancerList$EuMyc)) EuMycEnhUnion= union(EuMycEnhUnion, enhancerList$EuMyc[[i]]) # 65649	
	# regions boud by myc in at least one EuMyc condition
	EuMycMycUnion= peaksRefSimple$EuMyc[[1]]
	for(i in 2:5) EuMycMycUnion= union(EuMycMycUnion, peaksRefSimple$EuMyc[[i]]) # 48649	
	# enhancer not bound by myc in any EuMyc condition
	res= countOverlaps(EuMycEnhUnion, EuMycMycUnion)
	EuMycUnbEnhUnion= EuMycEnhUnion[res==0] # 43334
	# 6Kb around the enhancers and seleection of Chr1 regions for the heatmap
	EuMycUnbEnhUnion6Kb= GRsetwidth(EuMycUnbEnhUnion, newWidth=6000)
	EuMycUnbEnhUnion6KbChr1= EuMycUnbEnhUnion6Kb[seqnames(EuMycUnbEnhUnion6Kb)=='chr1'] # 3276
	grl= c(peaksRefSimple$EuMyc[c(1:10, 16:20, 11:15, 21:25)], CGI=CGIgr) # Myc, Pol2, K4me3, K4me1, K27ac
	EuMycUnbEnhData= GRheatmapData(grl= grl, refgr = EuMycUnbEnhUnion6KbChr1, type= rep('gr', length(grl)), nbins= 10, txdb = txdb)
	png('./figures/ExtDataFig3c.png', 3000, 2500, res=400)
		GRheatmapPlot(matList= EuMycUnbEnhData, rowLab = FALSE, colLab = TRUE, clusterInds = 1:length(EuMycUnbEnhData))
	dev.off()


# assigning Myc promoter peaks to genes
	mycBoundGSref= list()
	for(i in 1:4) {
		message(i)
		mycBoundGSref[[i]]= list()
		inds= grep('Myc', names(peaksRef[[i]]))
		for(cond in names(peaksRef[[i]])[inds]) {
			cat(cond, ' ')
			res = findOverlaps(peaksRef[[i]][[cond]], promoterGR, maxgap=0L, minoverlap=1L, type='any', select='all')
			resGS= peaksRef[[i]][[cond]]$enrichment[queryHits(res)]
			names(resGS)= promoterGR[subjectHits(res)]$gene_symbol
			mycBoundGSref[[i]][[cond]]= resGS
		}
	}
	names(mycBoundGSref)= names(peaksRef)

# building gene expression and Myc binding matrix
	# EuMyc
		EuMycExprBindMat= deseqRes[,-c(grep('vsn', colnames(deseqRes)))]
		for(cond in names(mycBoundGSref$EuMyc)) {
			message(cond)
			ints= rep(NA, nrow(EuMycExprBindMat))
			names(ints)= rownames(EuMycExprBindMat)
			int= mycBoundGSref$EuMyc[[cond]]
			gs= names(int)
			maxInt= tapply(int, INDEX=as.factor(gs), FUN=max)
			ints[names(maxInt)]= maxInt
			EuMycExprBindMat= cbind(EuMycExprBindMat, ints= ints[rownames(EuMycExprBindMat)])
		}
		colnames(EuMycExprBindMat)[20:24]= paste0(names(mycBoundGSref$EuMyc), '.int')
		save(EuMycExprBindMat, file='./data/EuMycExprBindMat.rda', compress=TRUE)
	# MycER
		mycERexprBindMat= mycERrnaseq
		for(cond in names(mycBoundGSref$mycER)) {
			message(cond)
			ints= rep(NA, nrow(mycERexprBindMat))
			names(ints)= rownames(mycERexprBindMat)
			int= mycBoundGSref$mycER[[cond]]
			gs= names(int)
			maxInt= tapply(int, INDEX=as.factor(gs), FUN=max)
			ints[names(maxInt)]= maxInt
			mycERexprBindMat= cbind(mycERexprBindMat, ints= ints[rownames(mycERexprBindMat)])
		}
		colnames(mycERexprBindMat)[1:2]= sub('X', 't', colnames(mycERexprBindMat)[1:2])
		colnames(mycERexprBindMat)[11:12]= paste0(names(mycBoundGSref$mycER), '.int')
		save(mycERexprBindMat, file='./data/mycERexprBindMat.rda', compress=TRUE)

	

# Differentially expressed genes (DEG)
	DEGlist= list(EuMyc=list(), MycER=list())
	# EuMyc (containing numbers for Fig2b Venn diagrams, marked with "Fig2b")
		over3inds= which(apply(EuMycExprBindMat[,1:11], 1, max)>=3) # 14840 over 21819
		mat= EuMycExprBindMat[over3inds,]
		# for each condition
			DEGlist$EuMyc$P= rownames(mat)[which(mat$qvalue.CP<0.05)] # 4352
			DEGlist$EuMyc$T1= rownames(mat)[which(mat$qvalue.CT1<0.05)] # 3968 "Fig2b"
			DEGlist$EuMyc$T2= rownames(mat)[which(mat$qvalue.CT2<0.05)] # 3947 "Fig2b"
			DEGlist$EuMyc$T3= rownames(mat)[which(mat$qvalue.CT3<0.05)] # 4592 "Fig2b"
			length(intersect(DEGlist$EuMyc$T1, intersect(DEGlist$EuMyc$T2, DEGlist$EuMyc$T3))) # 1914 "Fig2b"
		# for each group
			DEGlist$EuMyc$T= unique(c(DEGlist$EuMyc$T1, DEGlist$EuMyc$T2, DEGlist$EuMyc$T3)) # 6866
			DEGlist$EuMyc$Psp= DEGlist$EuMyc$P[!(DEGlist$EuMyc$P %in% DEGlist$EuMyc$T)] # 1018 "Fig2b"
			DEGlist$EuMyc$PT= intersect(DEGlist$EuMyc$P, DEGlist$EuMyc$T) # 3334 "Fig2b"
			DEGlist$EuMyc$Tsp= DEGlist$EuMyc$T[!(DEGlist$EuMyc$T %in% DEGlist$EuMyc$P)] # 3532 "Fig2b"
			DEGlist$EuMyc$allDEG= unique(c(DEGlist$EuMyc$P, DEGlist$EuMyc$T)) # 7884
		# breaking down in Up and Down
			DEGlist$EuMyc$Pup= rownames(mat)[which(mat$qvalue.CP<0.05 & mat$log2FC.CP>0)] # 2514
			DEGlist$EuMyc$Pdown= rownames(mat)[which(mat$qvalue.CP<0.05 & mat$log2FC.CP<0)] # 1838
			DEGlist$EuMyc$T1up= rownames(mat)[which(mat$qvalue.CT1<0.05 & mat$log2FC.CT1>0)] # 2600
			DEGlist$EuMyc$T1down= rownames(mat)[which(mat$qvalue.CT1<0.05 & mat$log2FC.CT1<0)] # 1368
			DEGlist$EuMyc$T2up= rownames(mat)[which(mat$qvalue.CT2<0.05 & mat$log2FC.CT2>0)] # 2648
			DEGlist$EuMyc$T2down= rownames(mat)[which(mat$qvalue.CT2<0.05 & mat$log2FC.CT2<0)] # 1299
			DEGlist$EuMyc$T3up= rownames(mat)[which(mat$qvalue.CT3<0.05 & mat$log2FC.CT3>0)] # 3110
			DEGlist$EuMyc$T3down= rownames(mat)[which(mat$qvalue.CT3<0.05 & mat$log2FC.CT3<0)] # 1482
	# MycER
		over3inds= which(apply(mycERexprBindMat[,1:11], 1, max)>=3) # 8584 over 21676
		mat= mycERexprBindMat[over3inds,]
		# for each time point
			DEGlist$mycER$t4h= rownames(mat)[which(mat$qvalue.4hvs0h<0.05)] # 1869
			DEGlist$mycER$t8h= rownames(mat)[which(mat$qvalue.8hvs0h<0.05)] # 1435
			DEGlist$mycER$t16h= rownames(mat)[which(mat$qvalue.16hvs0h<0.05)] # 2333
		# for each group
			DEGlist$mycER$allDEG= unique(c(DEGlist$mycER$t4h, DEGlist$mycER$t8h, DEGlist$mycER$t16h)) # 2587
		# breaking down in Up and Down
			DEGlist$mycER$t4hUp= rownames(mat)[which(mat$qvalue.4hvs0h<0.05 & mat$log2FC.4hvs0h>0)] # 1008
			DEGlist$mycER$t4hDown= rownames(mat)[which(mat$qvalue.4hvs0h<0.05 & mat$log2FC.4hvs0h<0)] # 861
			DEGlist$mycER$t8hUp= rownames(mat)[which(mat$qvalue.8hvs0h<0.05 & mat$log2FC.8hvs0h>0)] # 932
			DEGlist$mycER$t8hDown= rownames(mat)[which(mat$qvalue.8hvs0h<0.05 & mat$log2FC.8hvs0h<0)] # 503
			DEGlist$mycER$t16hUp= rownames(mat)[which(mat$qvalue.16hvs0h<0.05 & mat$log2FC.16hvs0h>0)] # 1337
			DEGlist$mycER$t16hDown= rownames(mat)[which(mat$qvalue.16hvs0h<0.05 & mat$log2FC.16hvs0h<0)] # 996
	save(DEGlist, file='./data/DEGlist.rda', compress=TRUE)

# heatmap, hierarchical clustering and scatter plots of DEGs
	# colorbar for DEG specific for P or T, or shared.
	sidecols= rep('black', length(DEGlist$EuMyc$allDEG))
	names(sidecols)=DEGlist$EuMyc$allDEG
	sidecols[DEGlist$EuMyc$Tsp]='brown' # specific for T
	sidecols[DEGlist$EuMyc$PT]='gold' # shared
	DEGmat= EuMycExprBindMat[DEGlist$EuMyc$allDEG,1:11]
	DEGmat[DEGmat==0]= min(DEGmat[DEGmat>0]) # replacing 0 with the minimal positive number
	DEGmat= log2(DEGmat/rowMeans(DEGmat[,1:4])) # normalizing to the baseline C
	# ExtDataFig4a: Hierarchical clustering of DEG
		DEGmat[DEGmat>5]=5
		DEGmat[DEGmat<(-5)]=-5
		pdf('./figures/ExtDataFig4a.pdf', 6, 6)
			plot(hclust(dist(t(DEGmat))), xlab='', ylab='Euclidean distance')
		dev.off()
	# Fig 2a: appending info on Myc binding
	bindMat= EuMycExprBindMat[DEGlist$EuMyc$allDEG,20:24]
	bindMat= bindMat/5
	bindMat[bindMat>2]= 2 # binding up to 10 to a scale of 2 in the heatmap
	bindMat[is.na(bindMat)]= 0
	# Fig2a: heatmap of gene expression and myc binding for EuMyc DEG
	DEGmat[DEGmat>2]=2
	DEGmat[DEGmat<(-2)]=-2	
	mat= cbind(bindMat, DEGmat)
	cpBB= colorRampPalette(c('blue','beige'))
	cpBR= colorRampPalette(c('beige','red'))
	cols= c(cpBB(64), cpBR(65)[-1])
	pdf('./figures/Fig2a.pdf', 6, 9)
		heatmap.2(as.matrix(mat), breaks=seq(-2, 2, length.out=129), col=cols, Colv=NULL, dendrogram='row',
			trace='none', density.info='none', RowSideColors=sidecols, labRow=FALSE, margins=c(8,1))
	dev.off()

# Fig2c: barplot of EuMyc DEG numbers
		degmatrix= matrix(NA, 4,3)
		rownames(degmatrix)= c('P','T1','T2','T3')
		colnames(degmatrix)= c('total','up','down')
		degmatrix[1,]=c(length(DEGlist$EuMyc$P), length(DEGlist$EuMyc$Pup), length(DEGlist$EuMyc$Pdown))
		degmatrix[2,]= c(length(DEGlist$EuMyc$T1), length(DEGlist$EuMyc$T1up), length(DEGlist$EuMyc$T1down))
		degmatrix[3,]= c(length(DEGlist$EuMyc$T2), length(DEGlist$EuMyc$T2up), length(DEGlist$EuMyc$T2down))
		degmatrix[4,]= c(length(DEGlist$EuMyc$T3), length(DEGlist$EuMyc$T3up), length(DEGlist$EuMyc$T3down))
		degmatrixBound= matrix(NA,4,3)
		Pbound= rownames(EuMycExprBindMat)[!is.na(EuMycExprBindMat$Myc.P.int)]
		degmatrixBound[1,]= c(length(intersect(DEGlist$EuMyc$P, Pbound)), length(intersect(DEGlist$EuMyc$Pup, Pbound)), length(intersect(DEGlist$EuMyc$Pdown, Pbound)))
		T1bound= rownames(EuMycExprBindMat)[!is.na(EuMycExprBindMat$Myc.T27.int)]
		degmatrixBound[2,]= c(length(intersect(DEGlist$EuMyc$T1, T1bound)), length(intersect(DEGlist$EuMyc$T1up, T1bound)), length(intersect(DEGlist$EuMyc$T1down, T1bound)))
		T2bound= rownames(EuMycExprBindMat)[!is.na(EuMycExprBindMat$Myc.T29.int)]
		degmatrixBound[3,]= c(length(intersect(DEGlist$EuMyc$T2, T2bound)), length(intersect(DEGlist$EuMyc$T2up, T2bound)), length(intersect(DEGlist$EuMyc$T2down, T2bound)))
		T3bound= rownames(EuMycExprBindMat)[!is.na(EuMycExprBindMat$Myc.T31.int)]
		degmatrixBound[4,]= c(length(intersect(DEGlist$EuMyc$T3, T3bound)), length(intersect(DEGlist$EuMyc$T3up, T3bound)), length(intersect(DEGlist$EuMyc$T3down, T3bound)))
		pdf('./figures/Fig2c.pdf', 3, 4)
			barplot(t(degmatrix), beside=TRUE, las=3)
			barplot(t(degmatrixBound), beside=TRUE, add=TRUE, names.arg=rep('', 4), density=15, col='gray60')
		dev.off()

# Fig2f: comparing nanostring and RNAseq results
	Nnorm= read.csv('./data/nanostring_normalized.csv', row.names=1) # Nanostring normalized
    sf=as.numeric(Nnorm['Tbp',2:10]/mean(as.numeric(Nnorm['Tbp',2:10])))
    for(i in 2:10) Nnorm[,i]= Nnorm[,i]*sf[i-1] # nornalizing by Tbp housekeeper
	NnormCE= read.csv('./data/nanostringCE_250000cells.csv', row.names=1) # Nanostring normalized and adjusted for cell equivalent
	Nnorm= cbind(Nnorm, Cavg= rowMeans(Nnorm[,c('C1','C2','C3')]), Pavg= rowMeans(Nnorm[,c('P1','P2','P3')]))
	NnormCE= cbind(NnormCE, Cavg= rowMeans(NnormCE[,c('C1.cell.eq','C2.cell.eq','C3.cell.eq')]),
	Pavg= rowMeans(NnormCE[,c('P1.cell.eq','P2.cell.eq','P3.cell.eq')]))
	NnormCE= NnormCE[,-(1:5)]
	bpfunComb= function(baseC, cond, Nmat) {
        logRs= log2(Nmat[,cond]/Nmat[,baseC])
        return(logRs)
    }
    resP= bpfunComb(baseC='Cavg', cond='Pavg', Nmat=NnormCE)
    resT1= bpfunComb(baseC='Cavg', cond='T1.cell.eq', Nmat=NnormCE)
    resT2= bpfunComb(baseC='Cavg', cond='T2.cell.eq', Nmat=NnormCE)
    resT3= bpfunComb(baseC='Cavg', cond='T3.cell.eq', Nmat=NnormCE)
    boxplotList= list(P=resP, T1=resT1, T2=resT2, T3=resT3)
    pdf('./figures/Fig2f.pdf', 3, 3.5)
        boxplot(boxplotList, outline=F, las=3, ylab='log2(Ratio vs C)')
        abline(h=0, lwd=2.5, col='gray50')
    dev.off()

# Fig3h_promoter: heatmap of promoter 3T9 Serum and MycER Myc peaks
	peaks3T9= c(peaksRefSimple$Serum[1:3], peaksRefSimple$mycER[1:2])
	names(peaks3T9)= paste(c(rep('Serum',3), rep('mycER',2)), names(peaks3T9), sep='.')
	peaks3T9prom= peaks3T9
	for(i in 1:length(peaks3T9prom)) {
		res= countOverlaps(peaks3T9prom[[i]], promoterGR)
		peaks3T9prom[[i]]= peaks3T9prom[[i]][res==1]
	}
	peaks3T9promUnion= peaks3T9prom[[1]]
	for(i in 2:length(peaks3T9prom)) peaks3T9promUnion= union(peaks3T9promUnion, peaks3T9prom[[i]]) # 4182
	peaks3T9promUnion6Kb= GRsetwidth(peaks3T9promUnion, newWidth=6000) # 3Kb upstream and downstream the mid point of each region (negative starts are corrected)
 	data3T9prom= GRheatmapData(grl= peaks3T9, refgr = peaks3T9promUnion6Kb, type= rep('gr', length(peaks3T9)), nbins= 10, txdb = txdb)
	png('./figures/Fig3h_promoter.png', 2500, 2500, res=400)
		GRheatmapPlot(matList= data3T9prom, rowLab = FALSE, colLab = TRUE, clusterInds = 1:length(data3T9prom))
	dev.off()

# Fig3h_distal: heatmap of distal 3T9 Serum and MycER Myc peaks
	peaks3T9distal= peaks3T9
	for(i in 1:length(peaks3T9distal)) {
		res= countOverlaps(peaks3T9distal[[i]], promoterGR)
		peaks3T9distal[[i]]= peaks3T9distal[[i]][res==0]
	}
	peaks3T9distalUnion= peaks3T9distal[[1]]
	for(i in 2:length(peaks3T9distal)) peaks3T9distalUnion= union(peaks3T9distalUnion, peaks3T9distal[[i]]) # 19305
	peaks3T9distalUnion6Kb= GRsetwidth(peaks3T9distalUnion, newWidth=6000) # 3Kb upstream and downstream the mid point of each region (negative starts are corrected)
	peaks3T9distalUnion6KbChr1= peaks3T9distalUnion6Kb[seqnames(peaks3T9distalUnion6Kb)=='chr1'] # 1174
 	data3T9distal= GRheatmapData(grl= peaks3T9, refgr = peaks3T9distalUnion6KbChr1, type= rep('gr', length(peaks3T9)), nbins= 10, txdb = txdb)
	png('./figures/Fig3h_distal.png', 2500, 2500, res=400)
		GRheatmapPlot(matList= data3T9distal, rowLab = FALSE, colLab = TRUE, clusterInds = 1:length(data3T9distal))
	dev.off()

# measuring Myc invasion of active promoters in EuMyc (proportion of H3K4me3 marked promoters that are Myc bound)
		EuMycInvasionMat= matrix(NA,5,5)
		rownames(EuMycInvasionMat)= c('C','P','T1','T2','T3')
		colnames(EuMycInvasionMat)= c('K4me3.total','%','Myc.total','%','Myc/K4me3 %')
		for(i in 1:5) {
			message(i)
			res= countOverlaps(promoterGR, peaksRefSimple$EuMyc[[i+15]]) # k4me3
			k4me3prom= promoterGR[res>0]
			res= countOverlaps(promoterGR, peaksRefSimple$EuMyc[[i]]) #myc
			mycProm= promoterGR[res>0]
			mycK4me3prom= countOverlaps(mycProm, k4me3prom) # myc & K4me3
			EuMycInvasionMat[i, 1]= length(k4me3prom)
			EuMycInvasionMat[i, 2]= round(100*length(k4me3prom) / length(promoterGR))
			EuMycInvasionMat[i, 3]= length(mycProm)
			EuMycInvasionMat[i, 4]= round(100*length(mycProm) / length(promoterGR))
			EuMycInvasionMat[i, 5]= round(100*length(mycK4me3prom) / length(k4me3prom))
		}
		show(EuMycInvasionMat)
			#    K4me3.total  % Myc.total  % Myc/K4me3 %
			# C        30811 56     10465 19          34
			# P        30462 55     20036 36          66
			# T1       27787 50     26070 47          94
			# T2       29097 53     25686 46          88
			# T3       29524 53     25820 47          87

# measuring Myc invasion of active enhancers (K27ac distal peaks) in EuMyc and mycER
	# EuMyc
		# C
			k27=peaksRef$EuMyc$K27ac.C1
			distalk27=k27[countOverlaps(k27, promoterGR)==0]
			length(which((countOverlaps(peaksRef$EuMyc$Myc.C1, distalk27))>0))/length(distalk27) # 3%
		# T1
			k27=peaksRef$EuMyc$K27ac.T27
			distalk27=k27[countOverlaps(k27, promoterGR)==0]
			length(which((countOverlaps(peaksRef$EuMyc$Myc.T27, distalk27))>0))/length(distalk27) # 60%
	# 3T9 MycER
		# 0h
			k27=peaksRef$mycER$K27ac.0h
			distalk27=k27[countOverlaps(k27, promoterGR)==0]
			length(which((countOverlaps(peaksRef$mycER$Myc.0h, distalk27))>0))/length(distalk27) # 33%
		# T1
			k27=peaksRef$mycER$K27ac.4h
			distalk27=k27[countOverlaps(k27, promoterGR)==0]
			length(which((countOverlaps(peaksRef$mycER$Myc.4h, distalk27))>0))/length(distalk27) # 59%


#####################################
######### Session info ############
#####################################
sessionInfo()
	# R version 3.0.2 (2013-09-25)
	# Platform: x86_64-unknown-linux-gnu (64-bit)

	# locale:
	# [1] C

	# attached base packages:
	# [1] parallel  stats     graphics  grDevices utils     datasets  methods  
	# [8] base     

	# other attached packages:
	#  [1] BSgenome.Mmusculus.UCSC.mm9_1.3.19      
	#  [2] BSgenome_1.30.0                         
	#  [3] org.Hs.eg.db_2.10.1                     
	#  [4] TxDb.Hsapiens.UCSC.hg19.knownGene_2.10.1
	#  [5] compEpiTools_0.1                        
	#  [6] org.Mm.eg.db_2.10.1                     
	#  [7] TxDb.Mmusculus.UCSC.mm9.knownGene_2.10.1
	#  [8] GenomicFeatures_1.14.5                  
	#  [9] gplots_2.13.0                           
	# [10] Rsamtools_1.14.3                        
	# [11] Biostrings_2.30.1                       
	# [12] GenomicRanges_1.14.4                    
	# [13] XVector_0.2.0                           
	# [14] IRanges_1.20.7                          
	# [15] topGO_2.14.0                            
	# [16] SparseM_1.03                            
	# [17] GO.db_2.10.1                            
	# [18] RSQLite_0.11.4                          
	# [19] DBI_0.2-7                               
	# [20] AnnotationDbi_1.24.0                    
	# [21] Biobase_2.22.0                          
	# [22] BiocGenerics_0.8.0                      
	# [23] graph_1.40.1                            

	# loaded via a namespace (and not attached):
	#  [1] KernSmooth_2.23-12 RCurl_1.95-4.1     XML_3.98-1.1       biomaRt_2.18.0    
	#  [5] bitops_1.0-6       caTools_1.16       gdata_2.13.3       grid_3.0.2        
	#  [9] gtools_3.3.1       lattice_0.20-29    rtracklayer_1.22.7 stats4_3.0.2      
	# [13] tools_3.0.2        zlibbioc_1.8.0  

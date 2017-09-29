### Mustafa: import needed packages ###
require(compEpiTools)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(org.Hs.eg.db)
require(BSgenome.Mmusculus.UCSC.mm9)

### Mustafa: loading data from local folder ###
load('/Users/mustafa.albahrani/workspace/clustering_tool/saboEtAl2014_sourceCode/data/peaksRef.rda')
load('/Users/mustafa.albahrani/workspace/clustering_tool/saboEtAl2014_sourceCode/data/mycSummitsRef.rda')
load('/Users/mustafa.albahrani/workspace/clustering_tool/saboEtAl2014_sourceCode/data/CGIgr.rda')
load('/Users/mustafa.albahrani/workspace/clustering_tool/saboEtAl2014_sourceCode/data/deseqRes.rda')
### Mustafa: NOTE: I could not find a file in the signified folder with this name! [I might send the papers authors about this dataset] ###
load('./data/ChIPfiles.rda')

########################################
###### compEpiTools additional functions
########################################

# - topGOres: species not recognized as method of OrgDb in topGOres, replace code with:
	topGOres= function(ids, Pthr=1e-5, maxN= 5000, minN=5, orgdb) {
		if(class(ids)!='character') stop('topGOres: ids has to be of class character ...')
		if(class(Pthr)!='numeric') stop('topGOres: Pthr has to be of class character ...')
		if(class(maxN)!='numeric') stop('topGOres: maxN has to be of class character ...')
		if(class(minN)!='numeric') stop('topGOres: minN has to be of class character ...')
		if(class(orgdb)!='OrgDb') stop('topGOres: orgdb has to be of class OrgDb ...')
		if(AnnotationDbi::species(orgdb)!='Homo sapiens' && AnnotationDbi::species(orgdb)!='Mus musculus')
			stop('topGOres: only org.Mm.eg.db and org.Hs.eg.db are supported ...')

		allEG= keys(orgdb)
		EGvec= rep(0, length(allEG))
		names(EGvec)= allEG
		EGvec[ids]= 1
		EGvec= factor(EGvec)
		if(AnnotationDbi::species(orgdb)=='Homo sapiens') dbname= 'org.Hs.eg.db'
		if(AnnotationDbi::species(orgdb)=='Mus musculus') dbname= 'org.Mm.eg.db'		
		GOdata= new("topGOdata", ontology = "BP", allGenes = EGvec, annot = annFUN.org, mapping= dbname, nodeSize=minN)
		res=runTest(GOdata, algorithm = "classic", statistic = "fisher")
		resdf= GenTable(GOdata, classic=res, topNodes = 50)
		pvals= resdf$classic
		pvals= as.numeric(sub('<', '', pvals))
		inds= which(resdf$Annotated<maxN & pvals< Pthr)
		if(length(inds)==0) return(NULL)
		resdf= resdf[inds,]
		return(resdf)
	}

# - simplifyGOterms: a function to eliminate GOterms having more than X% (maxOv)erlap with their children terms
	simplifyGOterms= function(goterms, maxOverlap= 0.8, ontology) {
		if(class(goterms)!='character') stop('simplifyGOterms: goterms has to be of class character ...')
		if(class(maxOverlap)!='numeric') stop('simplifyGOterms: maxOverlap has to be of class numeric ...')
		if(maxOverlap<0 || maxOverlap>1) stop('simplifyGOterms: maxOverlap is a percentage and has to range in [0,1] ...')
		if(!all(ontology %in% c('BP','CC','MF'))) stop('simplifyGOterms: ontology has to be one of: CC, BP, MF ...')
		
		### Mustafa: org.Mm.egGO2ALLEGS bellow is an R object that provides mappings between a given GO identifier and all of the Entrez Gene identifiers annotated at that GO term OR TO ONE OF IT’S CHILD NODES in the GO ontology ###
		go2allEGs= as.list(org.Mm.egGO2ALLEGS)
		if(ontology=='CC') go2parents= as.list(GOCCPARENTS) # Mustafa: this line means that if the `ontology` argument is set to `CC` then access the  cellular component data set (`GOCCPARENTS`) descriping the association of this GO term with their direct parent CC terms. ###
		if(ontology=='BP') go2parents= as.list(GOBPPARENTS) # Mustafa: this line means that if the `ontology` argument is set to `BP` then access the  biological process data set descriping the association of this GO term with their direct parent BP terms. ###
		if(ontology=='MF') go2parents= as.list(GOMFPARENTS) # Mustafa: this line means that if the `ontology` argument is set to `MF` then access the  molecular function data set descriping the association of this GO term with their direct parent MF terms. ###
		go2discard= NULL
		for(goterm in goterms) { ### Mustafa: this line iterates through each given GO term in `goterms` argument part of `simplifyGOterms` function ###
			parents= go2parents[[goterm]] ### Mustafa: this line reference (i.e. copy) the direct parents from of a `goterm` from the `go2parents` list ###
			parents= intersect(parents, goterms) ### Mustafa: caluclate the interesected GO terms between `parents` and `goterms` objects and then store them back into `parents` object ###
			if(length(parents)==0) next # no parents are found for a given GO term, check the others
			gotermEGs= go2allEGs[[goterm]] # EGs associated to a given GO term
			for(parent in parents) {
				parentEGs= go2allEGs[[parent]] # EGs associated to its parent
				commonEGs= intersect(gotermEGs, parentEGs)
				if(length(commonEGs)/length(parentEGs)>maxOverlap) go2discard= c(go2discard, parent)
			}
		}
		if(length(go2discard)>0) # discard redundant parents
		goterms= goterms[-which(goterms %in% go2discard)]
		return(goterms)
	}
	# example:
	# simplifyGOterms(goterms=c('GO:0002320','GO:0002244'), maxOverlap= 0.4, ontology='BP')

# - GR2fasta: a new GRanges method to extract and write to the disk a fasta file containing genomic seqs for the GRanges regions in a genome
''' I have stopped here !!!!!! '''
	setGeneric('GR2fasta', function(gr, bsgenome, fileout=NULL) standardGeneric('GR2fasta'))
	setMethod('GR2fasta','GRanges', function(gr, bsgenome, fileout=NULL) {
		if(class(bsgenome)!='BSgenome') stop('GR2fasta: org has to be of class BSgenome ...')
		if(!is.null(fileout) && class(fileout)!='character') stop('GR2fasta: fileout has to be either NULL or of class character ...')
		if(!is.null(fileout) && file.exists(fileout)) stop('GR2fasta: fileout already exists! ...')

		chrs= as.character(seqnames(gr))
		allseqs= NULL
		for(chr in unique(chrs)) {
			message(chr)
			chrinds= which(chrs==chr)
			chrseq= bsgenome[[chr]]
			starts= start(gr[chrinds])
			ends= end(gr[chrinds])
			seqs= Views(chrseq, starts, ends)
			seqnames= paste0(chr,':', starts, '-', ends)
			names(seqs)= seqnames
			seqs= DNAStringSet(as.character(seqs), use.names=T)
			allseqs= c(allseqs, seqs)
			if(!is.null(fileout)) writeXStringSet(seqs, filepath=fileout, append=TRUE, format="fasta")
		}
		invisible(allseqs)	
	})
	# example:
	# require(BSgenome.Mmusculus.UCSC.mm9)
	# gr= GRanges(Rle(c('chr1','chr2')), ranges=IRanges(start=c(1e7, 2e7), end=c(1e7+19, 2e7+19)))
	# show(GR2fasta(gr= gr, bsgenome= Mmusculus, fileout=NULL))



# - unionMaxPvalue: new GRanges method to perform union of peaks keeping the pvalue of the most significant peak:

	setGeneric('unionMaxPvalue', function(gr1, gr2, pvalue1= gr1$pvalue, pvalue2= gr2$pvalue) standardGeneric('unionMaxPvalue'))
	setMethod('unionMaxPvalue','GRanges', function(gr1, gr2, pvalue1= gr1$pvalue, pvalue2= gr2$pvalue) {
		if(class(gr2)!='GRanges') stop('unionMaxPvalue: gr2 has to be of class GRanges ...')
		if(class(pvalue1)!='numeric') stop('unionMaxPvalue: pvalue1 has to be of class numeric ...')
		if(class(pvalue2)!='numeric') stop('unionMaxPvalue: pvalue2 has to be of class numeric ...')

		grU= union(gr1, gr2)
		pvals= matrix(NA, length(grU), 2)
		# determining max pvalue over gr1
		ov1= findOverlaps(query= grU, subject= gr1, maxgap=0L, minoverlap=1L, type='any', select='all')
		ov1= tapply(pvalue1[subjectHits(ov1)], INDEX= as.factor(queryHits(ov1)), FUN= max)
		pvals[as.numeric(names(ov1)), 1]= ov1
		# determining max pvalue over gr2
		ov2= findOverlaps(query= grU, subject= gr2, maxgap=0L, minoverlap=1L, type='any', select='all')
		ov2= tapply(pvalue2[subjectHits(ov2)], INDEX= as.factor(queryHits(ov2)), FUN= max)
		pvals[as.numeric(names(ov2)), 2]= ov2
		# determine max pvalue over gr1 and gr2
		pvals= apply(pvals, 1, function(x) max(x, na.rm=TRUE))
		mcols(grU)= data.frame(pvalue= pvals)
		return(grU)
	})
	# example:
	# gr1= GRanges(Rle(c('chr1','chr1')), ranges= IRanges(start=c(100,200), end=c(150,250)), pvalue= c(200, 100))
	# gr2= GRanges(Rle(c('chr1','chr1')), ranges= IRanges(start=c(80,400), end=c(300,500)), pvalue= c(50, 150))
	# unionMaxPvalue(gr1, gr2)

# - GRoverlap: new function to determine the overlap between N GRanges and display an heatmap

	GRoverlap= function(GRlist, plot=TRUE) {
		if(class(GRlist)!='list') stop('GRoverlapFun: GRlist has to be a list ...')
		if(!all(sapply(GRlist, function(x) class(x)=='GRanges'))) stop('GRoverlapFun: GRlist has to be a list of GRanges ...')
		if(class(plot)!='logical') stop('GRoverlapFun: plot has to be of class logical ...')

		grL= length(GRlist)
		overlapMatrix=matrix(100, grL, grL)
		rownames(overlapMatrix)= names(GRlist)
		colnames(overlapMatrix)= names(GRlist)
		for(i in 1:grL) {
			message(i)
			for(j in (1:grL)[-i]) {
				message(j, ' ', appendLF = FALSE)
				res= countOverlaps(query= GRlist[[i]], subject= GRlist[[j]], maxgap=0L, minoverlap=1L, type='any')
				count= length(which(res>0))
				overlapMatrix[i,j]= 100*count/length(GRlist[[i]])
			}
		}
		if(plot) {
			cpWB= colorRampPalette(c('white','beige'))
			cpBR= colorRampPalette(c('beige','red'))
			cols= c(cpWB(50), cpBR(51)[-1])
			heatmap.2(overlapMatrix, cexRow=.7, cexCol=.7, col=cols, Rowv=NULL, Colv=NULL, cellnote=signif(overlapMatrix, 2),
				dendrogram='none', trace='none', density.info='none', notecol='black')
		}
		invisible(overlapMatrix)
	}
	# example:
	# starts= seq(100, 500, length.out=5)
	# gr1= GRanges(Rle('chr1'), ranges= IRanges(start=starts, end=starts+100))
	# starts= seq(300, 700, length.out=5)
	# gr2= GRanges(Rle('chr1'), ranges= IRanges(start=starts+50, end=starts+120))
	# GRoverlap(GRlist= list(gr1, gr2), plot=FALSE)


# - GRsetwidth: new GRanges method to set the width of a GRanges based on the mid point of each region:

	setGeneric('GRsetwidth', function(gr, newWidth) standardGeneric('GRsetwidth'))
	setMethod('GRsetwidth','GRanges', function(gr, newWidth) {
		if(class(newWidth)!='numeric') stop('newWidth: extension has to be of class numeric ...')
		mp =GRmidpoint(gr)
		newS= start(mp)- round(newWidth)/2
		newS[newS<0]= 1 # width is not guarantee to be exactly newWidth
		newE= start(mp)+ round(newWidth)/2 -1 # this might go beyond the chromosome length!
		start(gr)= newS
		end(gr)= newE
		return(gr)
	})
	# example:
	# gr= GRanges(Rle(c('chr1','chr1')), ranges= IRanges(start=c(100,200), end=c(150,250)))
	# GRsetwidth(gr, 1000)

# - update of the TSS function with the following for a complete TSS DB:

	TSS= function(txdb) {
		TSSr= promoters(txdb, upstream=0, downstream=0)
		ann=select(txdb, keys=transcripts(txdb)$tx_name, keytype='TXNAME', columns = c('TXNAME','GENEID'))
		rownames(ann)= ann[,1]
		names(TSSr)= ann[TSSr$tx_name,'GENEID']
		end(TSSr[strand(TSSr)=='+'])= start(TSSr[strand(TSSr)=='+']) # take min pos for - and max for + strand genes
		start(TSSr[strand(TSSr)=='-'])= end(TSSr[strand(TSSr)=='-'])
		return(TSSr)
	}


# - GRannotateSimple: a GRanges method to split a GRanges in three GRanges (promoter, intragenic and intergenic)
	setGeneric('GRannotateSimple', function(gr, txdb, upstream= 2000, downstream= 1000) standardGeneric('GRannotateSimple'))
	setMethod('GRannotateSimple','GRanges', function(gr, txdb, upstream= 2000, downstream= 1000) {
		if(class(txdb)!='TranscriptDb') stop('GRannotateSimple: txdb has to be an object of class TranscriptDb ...')
		if(class(upstream)!='numeric') stop('GRannotateSimple: upstream has to be an object of class numeric ...')
		if(class(downstream)!='numeric') stop('GRannotateSimple: downstream has to be an object of class numeric ...')

		grList= list()

		# promoter
		promoterGR= promoters(txdb, upstream= upstream, downstream= downstream)
		promoterCounts= countOverlaps(gr, promoterGR)
		promoterInds= which(promoterCounts>0) # gr ranges in promoters
		if(length(promoterInds)>0) {
			grList$promoter= gr[promoterInds]
			gr= gr[-promoterInds]
		}
		else grList$promoter= NA

		# intragenic
		trGR= transcripts(txdb)
		intragenicCounts= countOverlaps(gr, trGR)
		intragenicInds= which(intragenicCounts>0) # gr ranges in genebodies AND not in promoters
		if(length(intragenicInds)>0) {
			grList$intragenic= gr[intragenicInds]
			gr= gr[-intragenicInds]
		}
		else grList$intragenic= NA

		# intergenic
		if(length(gr)>0) {
			grList$intergenic= gr
		}
		else grList$intergenic= NA

		return(grList)
	})
	# example:
	# require(TxDb.Mmusculus.UCSC.mm9.knownGene)
	# txdb= TxDb.Mmusculus.UCSC.mm9.knownGene
	# gr= GRanges(Rle(c('chr1','chr1')), ranges= IRanges(start=c(100,200), end=c(150,250)))
	# GRannotateSimple(gr, txdb)

	
# - enhancers: a GRanges method to define enhancers as distal K4me1 peaks (GRanges) not associated to promoters
	setGeneric('enhancers', function(gr, txdb, upstream= 2000, downstream= 1000, CGIgr=NULL) standardGeneric('enhancers'))
	setMethod('enhancers','GRanges', function(gr, txdb, upstream= 2000, downstream= 1000, CGIgr=NULL) {
		if(class(txdb)!='TranscriptDb') stop('enhancers: txdb has to be an object of class TranscriptDb ...')
		if(class(upstream)!='numeric') stop('enhancers: upstream has to be an object of class numeric ...')
		if(class(downstream)!='numeric') stop('enhancers: downstream has to be an object of class numeric ...')
		if(!is.null(CGIgr) && class(CGIgr)!='GRanges') stop('enhancers: CGIgr has to be either NULL or an object of class GRanges ...') 

		promoterGR= promoters(txdb, upstream= upstream, downstream= downstream)
		# identification of gr GRanges not overlapping with promoters
		res= countOverlaps(gr, promoterGR)
		gr= gr[res==0]
		# optional identification of gr GRanges not overlapoing with CpG islands
		# might be useful to discard possible unannotated TSS
		if(!is.null(CGIgr)) {
			res= countOverlaps(gr, CGIgr)
			gr= gr[res==0]
		}

		return(gr)
	})
	# example:
	# require(TxDb.Mmusculus.UCSC.mm9.knownGene)
	# txdb= TxDb.Mmusculus.UCSC.mm9.knownGene
	# gr= GRanges(Rle(c('chr1','chr1')), ranges= IRanges(start=c(100,200), end=c(150,250)))
	# enhancers(gr, txdb)



#####################################
######### pre-processing ############
#####################################

	peaksRefSimple= peaksRef
	peaksRefSimple$EuMyc$Myc.P= union(peaksRef$EuMyc$Myc.P, peaksRef$EuMyc$Myc.P) # eliminating the redundancy between 3 replicates of Myc P
	txdb= TxDb.Mmusculus.UCSC.mm9.knownGene
	EG2GSlist=as.list(org.Mm.egSYMBOL)
	EG2GSlist.hsa=as.list(org.Hs.egSYMBOL)
	# Mmu promoters
		promoterGR= promoters(txdb, upstream= 2000, downstream= 1000)
		# adding gene ids and gene symbol to the promoter GRanges
		ann=select(txdb, keys=transcripts(txdb)$tx_name, keytype='TXNAME', columns = c('TXNAME','GENEID'))
		rownames(ann)= ann[,1]
		promoterGR$gene_id=ann[promoterGR$tx_name,'GENEID']
		promoterGR$gene_symbol=as.character(EG2GSlist[promoterGR$gene_id])
	# Hsa promoters
		txdb.hsa= TxDb.Hsapiens.UCSC.hg19.knownGene
		promoterGR.hsa= promoters(txdb.hsa, upstream= 2000, downstream= 1000)
		# adding gene ids and gene symbol to the promoter GRanges
		ann=select(txdb.hsa, keys=transcripts(txdb.hsa)$tx_name, keytype='TXNAME', columns = c('TXNAME','GENEID'))
		rownames(ann)= ann[,1]
		promoterGR.hsa$gene_id=ann[promoterGR.hsa$tx_name,'GENEID']
		promoterGR.hsa$gene_symbol=as.character(EG2GSlist.hsa[promoterGR.hsa$gene_id])




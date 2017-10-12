# BAM and BED files listed below can be retrieved in the GEO DB series id GSE51011 with the specified GSM sample id

# set of peaks were converted into GRanges objects based on this file (see below), saved in the ./data/ directory and will be automatically loaded with the R command line:
# source('analysisEnvironment.R')
# which is the first line of the saboEtAl2014_Figures.R, saboEtAl2014_ExtData.R and saboEtAl2014_ExtData10.R files

library(compEpiTools)

####### DNAseI-seq enriched regions
DNAseIpeaks= list()
DNAseIpeaks$t0h.A='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1230377_3T9.mycER.DNAseI.0hOHT.bed' 
DNAseIpeaks$t0h.B='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1395176_3T9.mycER.DNAseI.0hOHT.B.bed'
DNAseIpeaks$t4h.A='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1230378_3T9.mycER.DNAseI.4hOHT.bed'
DNAseIpeaks$t4h.B='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1395177_3T9.mycER.DNAseI.4hOHT.B.bed'

save(DNAseIpeaks, file='./data/DNAseIpeaks.rda', compress=TRUE)

##############################
####### listing BAM ChIP files
##############################
ChIPfiles= list(EuMyc=list(), mycER= list(), P493= list(), Serum= list())

####### EuMyc B-cell lymphoma
# Myc
	ChIPfiles$EuMyc$Myc.C1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234471_Eu-myc.Myc.C.1.bed' 	
	ChIPfiles$EuMyc$Myc.P1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234472_Eu-myc.Myc.P.1.bed' 	
	ChIPfiles$EuMyc$Myc.P2='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234473_Eu-myc.Myc.P.2.bed' 	
	ChIPfiles$EuMyc$Myc.P3='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234474_Eu-myc.Myc.P.3.bed'
	ChIPfiles$EuMyc$Myc.T27='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234475_Eu-myc.Myc.T.1.bed' 	
	ChIPfiles$EuMyc$Myc.T29='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234476_Eu-myc.Myc.T.2.bed' 	
	ChIPfiles$EuMyc$Myc.T31='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234477_Eu-myc.Myc.T.3.bed' 	
# Pol2
	ChIPfiles$EuMyc$Pol2.C1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234478_Eu-myc.Pol2.C.1.bed' 	
	ChIPfiles$EuMyc$Pol2.P1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234479_Eu-myc.Pol2.P.1.bed' 	
	ChIPfiles$EuMyc$Pol2.T27='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234480_Eu-myc.Pol2.T.1.bed'
	ChIPfiles$EuMyc$Pol2.T29='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234481_Eu-myc.Pol2.T.2.bed' 	
	ChIPfiles$EuMyc$Pol2.T31='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234482_Eu-myc.Pol2.T.3.bed' 	
# H3K4me1
	ChIPfiles$EuMyc$K4me1.C1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234488_Eu-myc.H3K4me1.C1.bed' 	
	ChIPfiles$EuMyc$K4me1.C2='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234519_Eu-myc.H3K4me1.C2.bed' 	
	ChIPfiles$EuMyc$K4me1.P1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234489_Eu-myc.H3K4me1.P1.bed' 	
	ChIPfiles$EuMyc$K4me1.P2='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234520_Eu-myc.H3K4me1.P2.bed' 
	ChIPfiles$EuMyc$K4me1.T27='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234490_Eu-myc.H3K4me1.T1.bed' 	
	ChIPfiles$EuMyc$K4me1.T29='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234491_Eu-myc.H3K4me1.T2.bed'
	ChIPfiles$EuMyc$K4me1.T31='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234492_Eu-myc.H3K4me1.T3.bed' 	
# H3K4me3
	ChIPfiles$EuMyc$K4me3.C1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234483_Eu-myc.H3K4me3.C1.bed' 	
	ChIPfiles$EuMyc$K4me3.C2='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234517_Eu-myc.H3K4me3.C2.bed' 	
	ChIPfiles$EuMyc$K4me3.P1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234484_Eu-myc.H3K4me3.P1.bed' 	
	ChIPfiles$EuMyc$K4me3.P2='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234518_Eu-myc.H3K4me3.P2.bed' 	
	ChIPfiles$EuMyc$K4me3.T27='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234485_Eu-myc.H3K4me3.T1.bed' 	
	ChIPfiles$EuMyc$K4me3.T29='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234486_Eu-myc.H3K4me3.T2.bed' 	
	ChIPfiles$EuMyc$K4me3.T31='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234487_Eu-myc.H3K4me3.T3.bed'
# H3K27ac
	ChIPfiles$EuMyc$K27ac.C1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234493_Eu-myc.H3K27ac.C.bed' 	
	ChIPfiles$EuMyc$K27ac.P1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234494_Eu-myc.H3K27ac.P.bed' 	
	ChIPfiles$EuMyc$K27ac.T27='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234495_Eu-myc.H3K27ac.T1.bed' 	
	ChIPfiles$EuMyc$K27ac.T29='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234496_Eu-myc.H3K27ac.T2.bed' 	
	ChIPfiles$EuMyc$K27ac.T31='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234497_Eu-myc.H3K27ac.T3.bed' 	

####### mycER 3T9 fibroblasts
# Myc
	ChIPfiles$mycER$Myc.0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234508_3T9.mycER.Myc.0hOHT.bed'
	ChIPfiles$mycER$Myc.4h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234509_3T9.mycER.Myc.4hOHT.bed' 
# Pol2
	ChIPfiles$mycER$Pol2.0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234505_3T9.mycER.Pol2.0hOHT.bed' 	
	ChIPfiles$mycER$Pol2.4h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234506_3T9.mycER.Pol2.4hOHT.bed' 	
# H3K4me1
	ChIPfiles$mycER$K4me1.0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234512_3T9.mycER.H3K4me1.0hOHT.bed'
	ChIPfiles$mycER$K4me1.4h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234513_3T9.mycER.H3K4me1.4hOHT.bed' 	
# H3K4me3
	ChIPfiles$mycER$K4me3.0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234510_3T9.mycER.H3K4me3.0hOHT.bed' 	
	ChIPfiles$mycER$K4me3.4h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234511_3T9.mycER.H3K4me3.4hOHT.bed' 	
# H3K27ac
	ChIPfiles$mycER$K27ac.0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234514_3T9.mycER.H3K27ac.0hOHT.bed' 	
	ChIPfiles$mycER$K27ac.4h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234515_3T9.mycER.H3K27ac.4hOHT.bed' 

####### P493
# Myc
	ChIPfiles$P493$Myc.t0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234499_P493.Myc.t0h.bed' 	
	ChIPfiles$P493$Myc.t1h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386342_P493.Myc.t1h.bed'
	ChIPfiles$P493$Myc.t24h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386343_P493.Myc.t24h.bed'
	ChIPfiles$P493$Myc.growing='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234501_P493.Myc.HighMyc.bed' 	
	ChIPfiles$P493$Myc.low='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234500_P493.Myc.LowMyc.bed'
# Pol2
	ChIPfiles$P493$Pol2.t0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234502_P493.Pol2.t0h.bed' 	
	ChIPfiles$P493$Pol2.t24h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386344_P493.Pol2.t24h.bed' 	
# H3K4me3
	#ChIPfiles$P493$K4me3.t0h= TBD
	ChIPfiles$P493$K4me3.t24h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386345_P493.H3K4me3.t24h.bed'
# H3K27ac
	#ChIPfiles$P493$K27ac.t0h= TBD 
	ChIPfiles$P493$K27ac.t24h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386346_P493.H3K27ac.t24h.bed'

####### Serum
# Myc
	ChIPfiles$Serum$Myc.0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386348_3T9.Serum.Myc.t0h.bed' 
	ChIPfiles$Serum$Myc.1h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386349_3T9.Serum.Myc.t1h.bed' 
	ChIPfiles$Serum$Myc.2h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386350_3T9.Serum.Myc.t2h.bed' 

# the following file is not provided as it contained private links. It should be created based on GEO samples files.
##save(ChIPfiles, file='./data/ChIPfiles.rda', compress=TRUE)





##############################
####### listing BAM input files
##############################
InputFiles= list(EuMyc=list(), mycER= list(), P493= list(), Serum= list())

####### EuMyc
# Myc
	InputFiles$EuMyc$Myc.C1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234471_Eu-myc.Myc.C.1.bed' 
	InputFiles$EuMyc$Myc.P1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234472_Eu-myc.Myc.P.1.bed' 
	InputFiles$EuMyc$Myc.P2='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234473_Eu-myc.Myc.P.2.bed'
	InputFiles$EuMyc$Myc.P3='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234474_Eu-myc.Myc.P.3.bed'
	InputFiles$EuMyc$Myc.T27='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234475_Eu-myc.Myc.T.1.bed'
	InputFiles$EuMyc$Myc.T29='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234476_Eu-myc.Myc.T.2.bed'
	InputFiles$EuMyc$Myc.T31='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234477_Eu-myc.Myc.T.3.bed'

####### mycER
# Myc
	InputFiles$mycER$Myc.0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234508_3T9.mycER.Myc.0hOHT.bed' 
	InputFiles$mycER$Myc.4h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234509_3T9.mycER.Myc.4hOHT.bed'

####### P493
# Myc
	#InputFiles$P493$Myc.no='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM894093
	InputFiles$P493$Myc.low='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234500_P493.Myc.LowMyc.bed'
	InputFiles$P493$Myc.high='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234501_P493.Myc.HighMyc.bed'
	#InputFiles$P493$Myc.growing='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM894093
	#InputFiles$P493$Myc.low='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM894093

####### Serum
# Myc
	InputFiles$Serum$Myc.0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386348_3T9.Serum.Myc.t0h.bed'
	InputFiles$Serum$Myc.1h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386349_3T9.Serum.Myc.t1h.bed' 
	InputFiles$Serum$Myc.2h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386350_3T9.Serum.Myc.t2h.bed' 

# the following file is not provided as it contained private links. It should be created based on GEO samples files.
##save(InputFiles, file='./data/InputFiles.rda', compress=TRUE)







##############################
####### listing MACS peaks BED files
##############################
macsfiles= list(EuMyc=list(), mycER= list(), P493= list(), Serum= list())

####### EuMyc
# Myc
	macsfiles$EuMyc$Myc.C1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234471_Eu-myc.Myc.C.1.bed' 
	macsfiles$EuMyc$Myc.P1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234472_Eu-myc.Myc.P.1.bed' 
	macsfiles$EuMyc$Myc.P2='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234473_Eu-myc.Myc.P.2.bed' 
	macsfiles$EuMyc$Myc.P3='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234474_Eu-myc.Myc.P.3.bed' 
	macsfiles$EuMyc$Myc.T27='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234475_Eu-myc.Myc.T.1.bed' 
	macsfiles$EuMyc$Myc.T29='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234476_Eu-myc.Myc.T.2.bed' 
	macsfiles$EuMyc$Myc.T31='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234477_Eu-myc.Myc.T.3.bed' 
# Pol2
	macsfiles$EuMyc$Pol2.C1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234478_Eu-myc.Pol2.C.1.bed'
	macsfiles$EuMyc$Pol2.P1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234479_Eu-myc.Pol2.P.1.bed'
	macsfiles$EuMyc$Pol2.T27='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234480_Eu-myc.Pol2.T.1.bed'
	macsfiles$EuMyc$Pol2.T29='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234481_Eu-myc.Pol2.T.2.bed'
	macsfiles$EuMyc$Pol2.T31='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234482_Eu-myc.Pol2.T.3.bed'
# H3K4me1
	macsfiles$EuMyc$K4me1.C1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234488_Eu-myc.H3K4me1.C1.bed'
	macsfiles$EuMyc$K4me1.C2='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234519_Eu-myc.H3K4me1.C2.bed'
	macsfiles$EuMyc$K4me1.P1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234489_Eu-myc.H3K4me1.P1.bed'
	macsfiles$EuMyc$K4me1.P2='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234520_Eu-myc.H3K4me1.P2.bed'
	macsfiles$EuMyc$K4me1.T27='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234490_Eu-myc.H3K4me1.T1.bed'
	macsfiles$EuMyc$K4me1.T29='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234491_Eu-myc.H3K4me1.T2.bed'
	macsfiles$EuMyc$K4me1.T31='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234492_Eu-myc.H3K4me1.T3.bed'
# H3K4me3
	macsfiles$EuMyc$K4me3.C1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234483_Eu-myc.H3K4me3.C1.bed'
	macsfiles$EuMyc$K4me3.C2='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234517_Eu-myc.H3K4me3.C2.bed'
	macsfiles$EuMyc$K4me3.P1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234484_Eu-myc.H3K4me3.P1.bed'
	macsfiles$EuMyc$K4me3.P2='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234518_Eu-myc.H3K4me3.P2.bed'
	macsfiles$EuMyc$K4me3.T27='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234485_Eu-myc.H3K4me3.T1.bed'
	macsfiles$EuMyc$K4me3.T29='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234486_Eu-myc.H3K4me3.T2.bed'
	macsfiles$EuMyc$K4me3.T31='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234487_Eu-myc.H3K4me3.T3.bed'
# H3K27ac
	macsfiles$EuMyc$K27ac.C1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234493_Eu-myc.H3K27ac.C.bed'
	macsfiles$EuMyc$K27ac.P1='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234494_Eu-myc.H3K27ac.P.bed'
	macsfiles$EuMyc$K27ac.T27='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234495_Eu-myc.H3K27ac.T1.bed'
	macsfiles$EuMyc$K27ac.T29='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234496_Eu-myc.H3K27ac.T2.bed'
	macsfiles$EuMyc$K27ac.T31='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234497_Eu-myc.H3K27ac.T3.bed'

####### mycER
# Myc
	macsfiles$mycER$Myc.0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234508_3T9.mycER.Myc.0hOHT.bed'
	macsfiles$mycER$Myc.4h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234509_3T9.mycER.Myc.4hOHT.bed'
# Pol2
	macsfiles$mycER$Pol2.0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234505_3T9.mycER.Pol2.0hOHT.bed'
	macsfiles$mycER$Pol2.4h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234506_3T9.mycER.Pol2.4hOHT.bed'
# H3K4me1
	macsfiles$mycER$K4me1.0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234512_3T9.mycER.H3K4me1.0hOHT.bed'
	macsfiles$mycER$K4me1.4h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234513_3T9.mycER.H3K4me1.4hOHT.bed'
# H3K4me3
	macsfiles$mycER$K4me3.0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234510_3T9.mycER.H3K4me3.0hOHT.bed'
	macsfiles$mycER$K4me3.4h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234511_3T9.mycER.H3K4me3.4hOHT.bed'
# H3K27ac
	macsfiles$mycER$K27ac.0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234514_3T9.mycER.H3K27ac.0hOHT.bed'
	macsfiles$mycER$K27ac.4h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234515_3T9.mycER.H3K27ac.4hOHT.bed'

####### P493
# Myc
	macsfiles$P493$Myc.t0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234499_P493.Myc.t0h.bed'
	macsfiles$P493$Myc.t1h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386342_P493.Myc.t1h.bed'
	macsfiles$P493$Myc.t24h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386343_P493.Myc.t24h.bed'
	macsfiles$P493$Myc.growing='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234501_P493.Myc.HighMyc.bed'
	macsfiles$P493$Myc.low='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234500_P493.Myc.LowMyc.bed'
# Pol2
	macsfiles$P493$Pol2.t0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1234502_P493.Pol2.t0h.bed'
	macsfiles$P493$Pol2.t24h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386344_P493.Pol2.t24h.bed'
# H3K4me3
	#macsfiles$P493$K4me3.t0h= TBD
	macsfiles$P493$K4me3.t24h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386345_P493.H3K4me3.t24h.bed'
# H3K27ac
	#macsfiles$P493$K27ac.t0h= TBD
	macsfiles$P493$K27ac.t24h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386346_P493.H3K27ac.t24h.bed'

####### Serum
# Myc
	macsfiles$Serum$Myc.0h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386348_3T9.Serum.Myc.t0h.bed'
	macsfiles$Serum$Myc.1h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386349_3T9.Serum.Myc.t1h.bed'
	macsfiles$Serum$Myc.2h='/Users/mustafa.albahrani/Desktop/heatmap_project/data/GSE51011_RAW/GSM1386350_3T9.Serum.Myc.t2h.bed'

# the following file is not provided as it contained private links. It should be created based on GEO samples files.
## save(macsfiles, file='./data/macsfiles.rda', compress=TRUE)



# loading MACS peaks (pvalue<1e-8) and transforming in GRanges
	peaksRef= macsfiles
	for(cell in names(peaksRef)) {
		message(cell, '\n')
		for(cond in names(macsfiles[[cell]])) {
			cat(cond, ' ')
			data= read.table(macsfiles[[cell]][[cond]], header=T, stringsAsFactors=FALSE)
			inds= which(data[,'X.10.log10.pvalue.']>80) # header of the MACS xls file, corresponding to the MACS BED file uploaded on GEO
			if(length(inds)>0) data= data[inds,]
			gr= GRanges(seqnames=Rle(data[,1]), ranges= IRanges(start=data[,2], end=data[,3]), pvalue= data[,'X.10.log10.pvalue.'])
			peaksRef[[cell]][[cond]]= gr
		}
	}

# determining enrichment for Myc peaks
	for(i in 1:4) {
		message(i)
		mycinds= grep('Myc', names(peaksRef[[i]]))
		for(j in mycinds) {
			cat(j, ' ')
			chipCov= GRcoverage(Object= peaksRef[[i]][[j]], bam= ChIPfiles[[i]][[j]], Nnorm=TRUE, Snorm=FALSE)
			inputCov= GRcoverage(Object= peaksRef[[i]][[j]], bam= InputFiles[[i]][[j]], Nnorm=TRUE, Snorm=FALSE)
			enrichment= signif(log2(chipCov - inputCov),3)
			mcols(peaksRef[[i]][[j]])=  data.frame(pvalue= peaksRef[[i]][[j]]$pvalue, enrichment= enrichment) 
		}
	}
	save(peaksRef, file='./data/peaksRef.rda', compress=TRUE)


# determining summit for Myc peaks
	mycSummits= list()
	for(i in 1:4) {
		message(i)
		mycSummits[[i]]= list()
		mycinds= grep('Myc', names(peaksRef[[i]]))
		for(j in names(peaksRef[[i]])[mycinds]) {
			cat(j, ' ')
			mycSummits[[i]][[j]]= GRcoverageSummit(Object= peaksRef[[i]][[j]], bam= ChIPfiles[[i]][[j]])
			pvalue= peaksRef[[i]][[j]]$pvalue
			enrichment= peaksRef[[i]][[j]]$enrichment
			mcols(mycSummits[[i]][[j]])=  data.frame(pvalue= pvalue, enrichment= enrichment) 
		}
	}
	names(mycSummits)= names(peaksRef)
	save(mycSummits, file='./data/mycSummits.rda', compress=TRUE)



# refining peaks and summits
	# concatenating myc peaks over P replicates
	temp= c(peaksRef$EuMyc$Myc.P1, peaksRef$EuMyc$Myc.P2, peaksRef$EuMyc$Myc.P3)
	peaksRefRef.EuMyc= c(peaksRef$EuMyc[1], list(Myc.P= temp), peaksRef$EuMyc[5:12]) # 31187 pooled peaks, 17340 peaks after union
	# union and max pvalue over replicates for histone marks
	temp= unionMaxPvalue(peaksRef$EuMyc$K4me1.C1, peaksRef$EuMyc$K4me1.C2)
	peaksRefRef.EuMyc= c(peaksRefRef.EuMyc, list(K4me1.C= temp))
	temp= unionMaxPvalue(peaksRef$EuMyc$K4me1.P1, peaksRef$EuMyc$K4me1.P2)
	peaksRefRef.EuMyc= c(peaksRefRef.EuMyc, list(K4me1.P= temp), peaksRef$EuMyc[17:19])
	temp= unionMaxPvalue(peaksRef$EuMyc$K4me3.C1, peaksRef$EuMyc$K4me3.C2)
	peaksRefRef.EuMyc= c(peaksRefRef.EuMyc, list(K4me3.C= temp))
	temp= unionMaxPvalue(peaksRef$EuMyc$K4me3.P1, peaksRef$EuMyc$K4me3.P2)
	peaksRefRef.EuMyc= c(peaksRefRef.EuMyc, list(K4me3.P= temp), peaksRef$EuMyc[24:31])
	peaksRef= peaksRef
	peaksRef$EuMyc= peaksRefRef.EuMyc
	save(peaksRef, file='./data/peaksRef.rda', compress=TRUE)

	temp= c(mycSummits$EuMyc$Myc.P1, mycSummits$EuMyc$Myc.P2, mycSummits$EuMyc$Myc.P3)
	mycSummitsRef.EuMyc= c(mycSummits$EuMyc[1], list(Myc.P= temp), mycSummits$EuMyc[5:7])
	mycSummitsRef= mycSummits
	mycSummitsRef$EuMyc= mycSummitsRef.EuMyc
	save(mycSummitsRef, file='./data/mycSummitsRef.rda', compress=TRUE)
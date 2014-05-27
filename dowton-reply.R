require(ape)
require(spider)#?spider
require(adhoc)#?adhoc
require(phangorn)
#clear memory
rm(list = ls())
save(list=c("mat", "spv", "fp", "sp"), file="/media/1TB/auto_backed_up/R/dowton_reply/dowton_analysis.RData", compress="bzip2", compression_level=9)
load(file="/media/1TB/auto_backed_up/R/dowton_reply/dowton_analysis.RData")
ls()
#load
dat <- read.dna(file="/home/rupert/aves_sp.fas", format="fasta")
#delete gaps
dat2 <- del.gaps(dat)

#trim up names
spl <- strsplit(names(dat2), "\\|")
nnam <- sapply(spl, function(xx) paste(gsub(" ", "_", xx[2]), xx[1], sep="_"))

#rename
names(dat2) <- nnam

#filter the less than 500 bp seqs
lx <- sapply(dat2, function(x) length(x))
dat3 <- dat2[which(lx>=500)]#to combine... which((lx>=500)&(lx<=1000))
dat4 <- dat3[sort(names(dat3))]


#remove shitty sequences - need to check alignment first 
rem <- c(
"Pterodroma_lessonii_PTERO001-14",
"Pterodroma_lessonii_PTERO002-14",
"Pterodroma_lessonii_PTERO003-14",
"Accipiter_cooperii_KBNA489-04",
"Pterodroma_macgillivrayi_GBIR2920-12",
"Pterodroma_macgillivrayi_GBIR2921-12",
"Gallus_domesticus_RSMS023-11",
"Tyto_alba_GBIR0313-06",
"Dendroica_coronata_KBNA033-04",
"Dendroica_coronata_KBNA034-04",
"Regulus_calendula_KBNA028-04",
"Vermivora_ruficapilla_KBNA045-04",
"Buceros_hydrocorax_CMCPB055-10",
"Egretta_garzetta_SIBIQ131-12",
"Accipiter_striatus_KBNA491-04",
"Sterna_paradisaea_KBNA506-04",
"Aethopyga_primigenia_CMCPB040-10",
"Agelaius_phoeniceus_KBNA038-04",
"Clamator_glandarius_SIBIQ106-12",
"Dicrurus_balicassius_CMCPB009-10",
"Euphagus_cyanocephalus_CDLSU035-05",
"Lonchura_leucogastra_CMCPB048-10",
"Molothrus_ater_BOTW211-04",
"Nothoprocta_ornata_USNMB012-10",
"Picoides_borealis_GBIR595-07",
"Vireo_flavifrons_BOTW479-06",
"Alcedo_cyanopectus_CMCPB064-10"
)

#remove shitty sequences
dat5 <- dat4[which(!names(dat4) %in% rem)]
length(dat4)
length(dat5)


#cleanup
ch <- checkDNA(dat5)#?checkDNA
dat6 <- dat5[which(!ch >=10)]#651*0.02
length(dat5)
length(dat6)
write.dna(dat6, file="aves_unaligned.fas", colw=100000, format="fasta")


#start alignment on "LYLIFG"
#end alignment on "LYQHL"


##anlaysis
dat <- read.dna(file="/media/1TB/auto_backed_up/R/dowton_reply/aves_trimmed.fas", format="fasta")


#get stats
uspl <- strsplit(dimnames(dat)[[1]], "_")
spv <- sapply(uspl, function(yy) paste(yy[1], yy[2], sep="_"))
gnv <- sapply(uspl, function(nn) nn[1])

length(dat)
dataStat(spv, gnv, thresh=2)#?dataStat

mat <- dist.dna(dat, model="raw", pairwise.deletion=TRUE)

#local minima good to remove high values first so more smoothed
mk <- mat[which(mat < 0.05)]
mi <- localMinima(mk)
plot(mi)


##threshopt
fp <- threshOpt(mat, spv, threshold=0.04)
sp <- threshOpt(mat, spv, threshold=0.006)

#cum error props
fp["Cumulative error"]/(fp["Cumulative error"]+fp["True pos"]+fp["True neg"])
sp["Cumulative error"]/(sp["Cumulative error"]+sp["True pos"]+sp["True neg"])

#false pos props
fp["False pos"]/(fp["Cumulative error"]+fp["True pos"]+fp["True neg"])
sp["False pos"]/(sp["Cumulative error"]+sp["True pos"]+sp["True neg"])

#false neg props
fp["False neg"]/(fp["Cumulative error"]+fp["True pos"]+fp["True neg"])
sp["False neg"]/(sp["Cumulative error"]+sp["True pos"]+sp["True neg"])

#dists
nnd <- nonConDist(mat, spv, propZero = FALSE, rmNA = FALSE)
summary(nnd)
quantile(nnd, probs=0.06)#?quantile
plot(sort(nnd))


threshVal <- seq(0.001,0.015, by = 0.001)
opt <- lapply(threshVal, function(x) threshOpt(mat, spv, thresh = x))
optMat <- do.call(rbind, opt)
barplot(t(optMat)[4:5,], names.arg=optMat[,1], xlab="Threshold values", ylab="Cumulative error")
legend(x = 2.5, y = 29, legend = c("False positives", "False negatives"), fill = c("grey75", "grey25"))

#best thresh=0.006
bc <- bestCloseMatch(mat, spv, threshold=0.006)
ti <- table(bc)
ti["correct"] / (ti["ambiguous"]+ti["correct"]+ti["incorrect"]+ti["no id"])#=0.8407376


#with singletons removed
sing <- rmSingletons(spv, exclude = TRUE)
ta <- table(bc[sing])
ta["correct"] / (ta["ambiguous"]+ta["correct"]+ta["incorrect"]+ta["no id"])#=0.8806994

# best thresh=0.006 including the no id singletons id'd as correct
singF <- rmSingletons(spv, exclude=FALSE)
ns <- table(bc[singF])
(ti["correct"]+ns["no id"]) / length(bc)# =0.8819149

#with 4%
fo <- bestCloseMatch(mat, spv, threshold=0.04)
fot <- table(fo[sing])
fot["correct"] / (fot["ambiguous"]+fot["correct"]+fot["incorrect"]+fot["no id"])#=0.915406

# NN with singletons removed
nn <- nearNeighbour(mat, spv, names = FALSE)
nnt <- table(nn[sing])
nnt["TRUE"]/(nnt["FALSE"]+nnt["TRUE"])#=0.9654505



#free up some memory
gc()

#Categorisation of these error rates follows Meyer and Paulay [82]: 
#"False positives are the identification of spurious novel taxa (splitting) within a species whose intraspecific variation extends deeper than the threshold value; 
#false negatives are inaccurate identification (lumping) within a cluster of taxa whose interspecific divergences are shallower than the proposed value” (p. 2230). 
#The optimum threshold is found where cumulative errors are minimised. 
#True positive identifications were recorded when only conspecific matches were delivered within the threshold percent of the query.
#False negative identifications occurred when more than one species was recorded within the threshold, and a 
#False positive was returned when there were no matches within the threshold value although conspecific species were available in the dataset.


##using adhoc

out1 <- checkDNAbcd(dat, DistModel="raw")
out2 <- adhocTHR(out1)

newteph <- tephdata
dimnames(tephdata["Bactrocera_cucurbitae_Benin_AB33598852_JEMU",])[[1]][1] <- "Newgenus_mesomelas_Cameroon_718_JEMU"


dimnames(newteph["Bactrocera_cucurbitae_Benin_AB33598852_JEMU",])[[1]][1] <- "Differentgenus_cucurbitae_Benin_AB33598852_JEMU"

attr(newteph, "dimnames")[[1]][3] <- "Differentgenus_cucurbitae_Benin_AB33598852_JEMU"
labels(newteph)
out1 <- checkDNAbcd(newteph, DistModel="raw")
out2 <- adhocTHR(out1)
  


data(tephdata);
out1<-checkDNAbcd(tephdata, DistModel="raw");
out2<-adhocTHR(out1);

layout(matrix(1,1,1));
par(font.sub=8);
plot(RE~thres,out2$IDcheck,ylim=c(0,max(c(out2$IDcheck$RE,out2$ErrProb))),xlab=NA,ylab=NA);
title(main="Ad hoc threshold",xlab="Distance", ylab="Relative identification error (RE)")
title(sub=paste("For a RE of",round(out2$ErrProb,4), "use a threshold of", round(out2$THR,4)));
regcoef<-out2$reg$coefficient;
curve(regcoef[1] + regcoef[2]*x + regcoef[3]*x^2 + regcoef[4]*x^3, add=TRUE);
segments(-1,out2$ErrProb,out2$THR,out2$ErrProb);segments(out2$THR,-1,out2$THR,out2$ErrProb);

#check if species names are unique
unique(spv)
unique(fghg)
fghg <- sapply(uspl, function(nn) nn[2])

str(out1)

## fleshfly analysis
tab <- read.table("/home/rupert/Dropbox/Dowton/genbank.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

seqs <- read.GB(tab$accession)
names(seqs) <- paste(attr(seqs, "species"), attr(seqs, "accession_num"), sep="_")
write.dna(seqs, file="/home/rupert/Dropbox/Dowton/sarcophaga_unaligned.fas", colw=100000, format="fasta")


sarc <- read.dna(file="/home/rupert/Dropbox/Dowton/sarcophaga_aligned.fas", format="fasta")

#get stats
sasp <- strsplit(dimnames(sarc)[[1]], "_")
sspv <- sapply(sasp, function(yy) paste(yy[1], yy[2], sep="_"))
sgnv <- sapply(sasp, function(nn) nn[1])

length(labels(sarc))
dataStat(sspv, sgnv, thresh=2)#?dataStat

smat <- dist.dna(sarc, model="raw", pairwise.deletion=TRUE)

threshVal <- seq(0.001,0.10, by = 0.001)
opt <- lapply(threshVal, function(x) threshOpt(smat, sspv, thresh = x))
optMat <- do.call(rbind, opt)
barplot(t(optMat)[4:5,], names.arg=optMat[,1], xlab="Threshold values", ylab="Cumulative error")
legend(x = 2.5, y = 29, legend = c("False positives", "False negatives"), fill = c("grey75", "grey25"))


#          Threshold True neg True pos False neg False pos Cumulative error
#[23,]     0.023       13      414         0         1                1
#[40,]     0.040       11      388        29         0               29

smi <- localMinima(smat)#local min says 2.3 best

to <- threshOpt(smat, sspv, thresh = 0.023)
fo <- threshOpt(smat, sspv, thresh = 0.04)

round((to["Cumulative error"]/(to["Cumulative error"]+to["True pos"]+to["True neg"])), digits=4)
round((fo["Cumulative error"]/(fo["Cumulative error"]+fo["True pos"]+fo["True neg"])), digits=4)

#best thresh=0.023
bc <- bestCloseMatch(smat, sspv, threshold=0.023)
ti <- table(bc)
ti["correct"] / (ti["correct"]+ti["no id"])# =0.9672897

# best thresh=0.023 with singletons removed
sing <- rmSingletons(sspv, exclude = TRUE)
ta <- table(bc[sing])
ta["correct"] / (ta["correct"]+ta["no id"])# =0.9975904

# best thresh=0.023 including the no id singletons id'd as correct
singF <- rmSingletons(sspv, exclude=FALSE)
ns <- table(bc[singF])
(ti["correct"]+ns["no id"]) / (ti["correct"]+ti["no id"])


#poor thresh=0.04
bn <- bestCloseMatch(smat, sspv, threshold=0.04)
ty <- table(bn)
ty["correct"] / (ty["correct"]+ty["incorrect"]+ty["no id"])# = 0.9696262

# poor thresh=0.04 with singletons removed
tl <- table(bn[sing])
tl["correct"] / (tl["correct"])# =1


#assess which specimens cause problems 0.023
sspv[which(bc=="incorrect")]
sspv[which(bc=="ambiguous")]
singlist <- sspv[rmSingletons(sspv, exclude=FALSE)]
sspv[which(bc=="no id")]
sspv[which(bc=="no id")][which(!sspv[which(bc=="no id")] %in% singlist)]

#assess which specimens cause problems 0.04
sspv[which(bn=="incorrect")]
sspv[which(bn=="ambiguous")]
sspv[which(bn=="no id")]
sspv[which(bn=="no id")][which(!sspv[which(bn=="no id")] %in% singlist)]


cbind(sspv,bc)
cbind(sspv,bn)



# problems with 0.023
# hap sharing between two S. zeta (JN965170, JN965167) and the only S. beta (GQ254434) # the S. beta specimen (as per GenBank) is identified as S. zeta in Meiklejohn (2012) suppl
# S. hardyi (GQ254450 at 6.11% div)  could not be ID but conspecifics present. the specimen is S. furcata according to Meiklejohn (2012) suppl
# and S. australis (JN964738 KM673 at 3.81%) could not be ID but conspecifics present.

# problems with 0.04
# ID S. megafilosia (GQ254479) as S. meiofilosia (GQ254481) and vice versa at 2.58%
# S. meiofilosia bin BOLD:AAI0956. S. megafilosia bin BOLD:AAI0955
# same zeta/beta problem
# ID S. australis correctly, but still not hardyi


#dists
snnd <- nonConDist(smat, sspv, propZero = FALSE, rmNA = FALSE)
summary(snnd)
cbind(labels(sarc), snnd)

smid <- maxInDist(smat, sspv, propZero = FALSE, rmNA = FALSE)#?maxInDist
summary(smid)
cbind(labels(sarc), smid)


unique(sspv[which(sgnv == "Sarcophaga")])

#check nj
pdf(file="pdf.pdf", width=10, height=60)
plot(midpoint(ladderize(nj(smat))))
dev.off()

#try adhoc
out1 <- checkDNAbcd(sarc, DistModel="raw")
out2 <- adhocTHR(out1, Reg="linear")#?adhocTHR


####################################################### USING THE REAL DATA #################################################################################3

dat <- read.nexus.data(file="/home/rupert/Dropbox/Dowton/COI_CAD_SystBiol.nex")#?read.nexus.data
#remove the non-Sarcophaga specimens and the unknown specimens
#nonsarc <- c(
#grep("mil", labels(dat)),
#grep("tach", labels(dat)),
#grep("dome", labels(dat)),
#grep("bre", labels(dat)),
#grep("var", labels(dat)),
#grep("blA", labels(dat)),
#grep("blA", labels(dat)))
#
#datred <- dat[-nonsarc]
#labels(dat)

#split up the data into COI and CAD
datbin <- as.matrix(as.DNAbin(dat))#?ape
coi <- datbin[1:length(labels(datbin)),1:657]
cad <- datbin[1:length(labels(datbin)),658:1434]



tr <- midpoint(ladderize(nj(dist.dna(coi, model="raw", pairwise.deletion=TRUE))))
pdf(file="dow_coi_red.pdf", width=6, height=40)
plot(tr, no.margin=TRUE, cex=0.6,  edge.width=1, label.offset=0.0005, font=1)
#nodelabels(frame="n")
dev.off()

#relabel with species vector and rename outgroups
sspv <- sapply(strsplit(labels(coi), "KM[0-9][0-9][0-9]|JW[0-9][0-9][0-9]|JW[0-9][0-9][0-9]v1|JW[0-9][0-9][0-9]v2|JW[0-9][0-9]v2"), function(xx) paste(xx[2], sep=""))
sspv[434] <- "tach"
sspv[435] <- "M_dome"


#get number of species
length(unique(labels(sspv)))

smat <- dist.dna(coi, model="raw", pairwise.deletion=TRUE)

threshVal <- seq(0.001,0.10, by = 0.001)
opt <- lapply(threshVal, function(x) threshOpt(smat, sspv, thresh = x))
optMat <- do.call(rbind, opt)
barplot(t(optMat)[4:5,], names.arg=optMat[,1], xlab="Threshold values", ylab="Cumulative error")
legend(x = 2.5, y = 400, legend = c("False positives", "False negatives"), fill = c("grey75", "grey25"))

#local minima estimation
smi <- localMinima(smat)#local min says 2.13% best

to <- threshOpt(smat, sspv, thresh = 0.0213)
fo <- threshOpt(smat, sspv, thresh = 0.04)

round((to["Cumulative error"]/(to["Cumulative error"]+to["True pos"]+to["True neg"])), digits=4)
round((fo["Cumulative error"]/(fo["Cumulative error"]+fo["True pos"]+fo["True neg"])), digits=4)


#best thresh=0.0213
bc <- bestCloseMatch(smat, sspv, threshold=0.0213)
ti <- table(bc)
ti["correct"] / length(bc)# = 0.954023


# best thresh=0.0213 with singletons removed
sing <- rmSingletons(sspv, exclude = TRUE)
ta <- table(bc[sing])
ta["correct"] / length(sing)# = 0.9857482

# best thresh=0.0213 including the no id singletons id'd as correct
singF <- rmSingletons(sspv, exclude=FALSE)
ns <- table(bc[singF])
(ti["correct"]+ns["no id"]) / length(bc)# = 0.9862069
 
#assess which specimens cause problems 0.0213
labels(coi)[which(bc=="incorrect")]
labels(coi)[which(bc=="ambiguous")]
labels(coi)[which(bc=="no id")]


#get indexes of only those in the Dowton paper
sarc <- c(
grep("afr", labels(coi)),
grep("alc", labels(coi)),
grep("alp", labels(coi)),
grep("iot", labels(coi)),#aurifrons
grep("aus", labels(coi)),
grep("koh", labels(coi)),#bancroftorum clade 1
grep("per", labels(coi)),#bancroftorum clade 2
grep("bid", labels(coi)),
grep("bif", labels(coi)),
grep("cra", labels(coi)),
grep("dux", labels(coi)),
grep("fro", labels(coi)),
grep("fur", labels(coi)),
grep("imp", labels(coi)),
grep("kap", labels(coi)),
grep("unE", labels(coi)),#kohla
grep("blC", labels(coi)),#littoralis
grep("meg", labels(coi)),
grep("mei", labels(coi)),
grep("unF", labels(coi)),#misera
grep("mul", labels(coi)),
grep("omi", labels(coi)),
grep("unC", labels(coi)),#peregrina
grep("sen", labels(coi)),#piva
grep("pra", labels(coi)),
grep("ruf", labels(coi)),
grep("aur", labels(coi)),#sigma
grep("spf", labels(coi)),
grep("spg", labels(coi)),
grep("mis", labels(coi)),#taenionota
grep("tor", labels(coi)),
grep("lon", labels(coi)),#villisterna
grep("zet", labels(coi)))

#with singletons removed
singlist <- sspv[rmSingletons(sspv, exclude=FALSE)]
labels(coi)[which(bc=="incorrect")][which(!sspv[which(bc=="incorrect")] %in% singlist)]
labels(coi)[which(bc=="ambiguous")][which(!sspv[which(bc=="ambiguous")] %in% singlist)]
labels(coi)[which(bc=="no id")][which(!sspv[which(bc=="no id")] %in% singlist)]
#only the non-singletons in sarc (queries from the paper)
jh <- labels(coi)[which(bc=="no id")][which(!sspv[which(bc=="no id")] %in% singlist)]
jh[jh %in% labels(coi)[sarc]]

#no id 3
#incorrect 1
####

labels(coi)
#just the 32 species in the paper
yu <- table(bc[sarc])
yu["correct"] / length(sarc)# = 0.9703704
#check the names
labels(coi)[sarc[which(bc[sarc] == "no id")]]
#gf <- labels(coi)[sarc[which(bc[sarc] == "no id")]]





#excluding the singletons in the species list from the paper
fg <- bc[sarc[which(labels(coi)[sarc] %in% labels(coi)[sing])]]
kl <- table(fg)
kl["correct"] / length(fg)# = 0.9899244



#poor thresh=0.04
bn <- bestCloseMatch(smat, sspv, threshold=0.04)
ty <- table(bn)
ty["correct"] / length(bn)# = 0.9586207

# poor thresh=0.04 with singletons removed
tl <- table(bn[sing])
tl["correct"] / length(bn[sing])# = 0.9904988




#singlist <- sspv[rmSingletons(sspv, exclude=FALSE)]
lk <- rmSingletons(sspv, exclude=FALSE)
io <- sspv[sarc[sarc %in% lk]]
labels(datbin)[which(bc=="no id")][which(!sspv[which(bc=="no id")] %in% )]

sspv[sarc]

which(!sspv[which(bn=="no id")] %in% singlist)

which(!sspv[which(bc=="no id")] %in% sarc)

sarc[which(bc=="no id") %in% sarc]

# KM592koh is incorrect. Clusters with "mis" in the tree. labelled as "taenionota" in MPE and GenBank. 



#assess which specimens cause problems 0.04
sspv[which(bn=="incorrect")]
sspv[which(bn=="ambiguous")]
sspv[which(bn=="no id")]
sspv[which(bn=="no id")][which(!sspv[which(bn=="no id")] %in% singlist)]


cbind(sspv,bc)
cbind(sspv,bn)



#nonsarc <- c(
#grep("mil", labels(dat)),
#grep("tach", labels(dat)),
#grep("dome", labels(dat)),
#grep("bre", labels(dat)),
#grep("var", labels(dat)),
#grep("blA", labels(dat)),
#grep("blA", labels(dat)))
#
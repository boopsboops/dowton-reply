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

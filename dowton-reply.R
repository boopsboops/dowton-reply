## R script to repeat analyses from Collins & Cruickshank: "Known knowns, known unknowns, unknown unknowns and unknown knowns in DNA barcoding: A reply to Dowton et al."
## Date 27/05/14

# load libraries
require("ape")# v3.1-1
require("spider")# v1.3-0
require("phangorn")# v1.99-7

# read the data is nexus format from Dryad
dat <- read.nexus.data(file="http://datadryad.org/bitstream/handle/10255/dryad.53505/COI_CAD_SystBiol.nex?sequence=1")

#split up the data into COI and CAD
datbin <- as.matrix(as.DNAbin(dat))#?ape
coi <- datbin[1:length(labels(datbin)),1:657]
cad <- datbin[1:length(labels(datbin)),658:1434]

# to make a quick NJ tree
#tr <- midpoint(ladderize(nj(dist.dna(coi, model="raw", pairwise.deletion=TRUE))))
#pdf(file="coi.pdf", width=6, height=40)
#plot(tr, no.margin=TRUE, cex=0.6,  edge.width=1, label.offset=0.0005, font=1)
#dev.off()

# re-label with species vector and rename outgroups
sspv <- sapply(strsplit(labels(coi), "KM[0-9][0-9][0-9]|JW[0-9][0-9][0-9]|JW[0-9][0-9][0-9]v1|JW[0-9][0-9][0-9]v2|JW[0-9][0-9]v2"), function(xx) paste(xx[2], sep=""))
sspv[434] <- "tach"
sspv[435] <- "M_dome"

# get number of species
length(unique((sspv)))

# make distance matrix using p distances
smat <- dist.dna(coi, model="raw", pairwise.deletion=TRUE)

# get false pos and false neg errors for varying thresholds
threshVal <- seq(0.001,0.10, by = 0.001)
opt <- lapply(threshVal, function(x) threshOpt(smat, sspv, thresh = x))
optMat <- do.call(rbind, opt)
print(optMat)
# plot
barplot(t(optMat)[4:5,], names.arg=optMat[,1], xlab="Threshold values", ylab="Cumulative error")
legend(x = 2.5, y = 400, legend = c("False positives", "False negatives"), fill = c("grey75", "grey25"))
# optimum threshold is between 0.017 and 0.025

# Categorisation of these error rates follows Meyer and Paulay [82]: 
# "False positives are the identification of spurious novel taxa (splitting) within a species whose intraspecific variation extends deeper than the threshold value; 
# false negatives are inaccurate identification (lumping) within a cluster of taxa whose interspecific divergences are shallower than the proposed value” (p. 2230). 
# The optimum threshold is found where cumulative errors are minimised. 
# True positive identifications were recorded when only conspecific matches were delivered within the threshold percent of the query.
# False negative identifications occurred when more than one species was recorded within the threshold, and a 
# False positive was returned when there were no matches within the threshold value although conspecific species were available in the dataset.


# accurate estimation of local minima threshold
smi <- localMinima(smat)# local minima says 0.0213 (= 2.13%) is optimum. This lies within the range estimated by 'threshOpt' above.

# compare the error rates for the optimum threshold (2.13%) and the arbitrary threshold (4%)
to <- threshOpt(smat, sspv, thresh = 0.0213)
fo <- threshOpt(smat, sspv, thresh = 0.04)
round((to["Cumulative error"]/(to["Cumulative error"]+to["True pos"]+to["True neg"])), digits=4)
round((fo["Cumulative error"]/(fo["Cumulative error"]+fo["True pos"]+fo["True neg"])), digits=4)


# run the best close match on the optimum thresholf (= 0.0213)
bc <- bestCloseMatch(smat, sspv, threshold=0.0213)

# index the singletons
sing <- rmSingletons(sspv, exclude = TRUE)

#get indexes of only those 33 species queried in the Dowton paper
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


# comparison and tabulate results
# for best close match (BCM) at 0.0213 (= 2.13% p-dist)
vc <- data.frame(cbind(labels(coi),sspv,bc))
vc$new1 <- NA
vc$new2 <- NA
vc$new1[match(labels(coi)[sarc], vc$V1)] <- "query"
vc$new1[match(labels(coi)[-sarc], vc$V1)] <- "non-query"
vc$new2[match(labels(coi)[sing], vc$V1)] <- "non-singleton"
vc$new2[match(labels(coi)[-sing], vc$V1)] <- "singleton"
vc$new3 <- nonConDist(smat, sspv, propZero = FALSE, rmNA = FALSE)# distance to closest non-conspecific
vc$new4 <- maxInDist(smat, sspv, propZero = FALSE, rmNA = FALSE)# maximum intraspecific distance
vc$new4[vc$new2 == "singleton"] <- NA # insert NAs for singletons
colnames(vc) <- c("taxon_label", "species_vector", "bcm_id_result", "query_status", "singleton_status", "non_con_dist", "max_in_dist")
print(vc)

# all results 
t1 <- table(vc$bcm_id_result)
t1["correct"] / length(vc$bcm_id_result)# = 0.954023

# all excluding singletons
t2 <- table(vc$bcm_id_result[vc$singleton_status == "non-singleton"])
t2["correct"] / length(vc$bcm_id_result[vc$singleton_status == "non-singleton"])# = 0.9857482

# just those used as queries by Dowton et al. (i.e. sarc)
t3 <- table(vc$bcm_id_result[vc$query_status == "query"])
t3["correct"] / length(vc$bcm_id_result[vc$query_status == "query"])# = 0.9703704

# just those used as queries by Dowton et al., but excluding singletons
t4 <- table(vc$bcm_id_result[vc$query_status == "query" & vc$singleton_status == "non-singleton"])
t4["correct"] / length(vc$bcm_id_result[vc$query_status == "query" & vc$singleton_status == "non-singleton"])# = 0.9899244

# notes 
# KM592koh is incorrect. Clusters with "mis" in the tree. labeled as "taenionota" in MPE and GenBank. 
# KM673aus is "no id". 
# KM311omi is "no id". 
# KM260spg is "no id". 



## ANALYSIS of CORRECTED COI sequences
ccoi <- read.dna(file="/home/rupert/LaTeX/dowton-reply/COI_corrected.fasta", format="fasta")

# re-label with species vector and rename outgroups
csspv <- sapply(strsplit(labels(ccoi), "KM[0-9][0-9][0-9]|JW[0-9][0-9][0-9]|JW[0-9][0-9][0-9]v1|JW[0-9][0-9][0-9]v2|JW[0-9][0-9]v2"), function(xx) paste(xx[2], sep=""))
csspv[434] <- "tach"
csspv[435] <- "M_dome"

# make distance matrix using p distances
csmat <- dist.dna(ccoi, model="raw", pairwise.deletion=TRUE)

# run the best close match
cbc <- bestCloseMatch(csmat, csspv, threshold=0.0213)

#check new result for KM311omi
cbc[which(labels(ccoi) == "KM311omi")]

#check new result for KM260spg
cbc[which(labels(ccoi) == "KM260spg")]



## ANALYSIS of CAD

o <- (1-(1/405))*100
round(o, digits=1)

vc$max_in_dist[vc$taxon_label == "KM673aus"] *100

?spider
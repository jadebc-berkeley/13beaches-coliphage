
# ------------------------------------
# 10-aim1-wq-summary-tables.R
#
# Create summary statistics for
# water quality data by 
# beach and overall
# for Enterococcus EPA 1600 and qPCR 1611
#
# ------------------------------------

# ------------------------------------
# preamble
# ------------------------------------
rm(list=ls())



# ------------------------------------
# load the water quality sample data
# ------------------------------------
d <- read.csv("~/dropbox/13beaches/data/final/13beaches-wq-samples.csv",stringsAsFactors=F)

d$beachcode[d$beachcode==""] <- d$beach[d$beachcode==""]

# ------------------------------------
# At Mission Bay, recode the 
# EPA 1600 values with Enterolert,
# which is the assay used in the analysis
# at that beach
# ------------------------------------
d$entero1600cfu[d$beach=="Mission Bay"] <- d$enteroELTmpn[d$beach=="Mission Bay"]
d$entero1600cfu_nd[d$beach=="Mission Bay"] <- d$enteroELTmpn_nd[d$beach=="Mission Bay"]


# ------------------------------------
# Summarize non-detects by beach
# and sampling location
# ------------------------------------
nd.qpcr <- table(d$beachcode,d$enteroQPCRcce_nd)[,2:3]
nd.1600 <- table(d$beachcode,d$entero1600cfu_nd)[,2:3]

nd.qpcr
nd.1600


tot.nd.qpcr <- table(d$enteroQPCRcce_nd)[2:3]
tot.nd.1600 <- table(d$entero1600cfu_nd)[2:3]


# ------------------------------------
# impute non-detects at 0.1 before
# calculating summary statistics
# ------------------------------------

d$enteroQPCRcce[d$enteroQPCRcce_nd=="Below detection"] <- 0.1
d$entero1600cfu[d$entero1600cfu_nd=="Below detection"] <- 0.1

# ------------------------------------
# Summarize min and max by beach
# and sampling location
# ------------------------------------

min.qpcr <- tapply(d$enteroQPCRcce,d$beachcode,function(x) min(x,na.rm=T))
max.qpcr <- tapply(d$enteroQPCRcce,d$beachcode,function(x) max(x,na.rm=T))

min.1600 <- tapply(d$entero1600cfu,d$beachcode,function(x) min(x,na.rm=T))
max.1600 <- tapply(d$entero1600cfu,d$beachcode,function(x) max(x,na.rm=T))

round(cbind(min.qpcr, max.qpcr),1)
round(cbind(min.1600, max.1600),1)

# overall min/max
tot.min.qpcr <- min(d$enteroQPCRcce,na.rm=T)
tot.max.qpcr <- max(d$enteroQPCRcce,na.rm=T)

tot.min.1600 <- min(d$entero1600cfu,na.rm=T)
tot.max.1600 <- max(d$entero1600cfu,na.rm=T)

# ------------------------------------
# Summarize geometric mean by beach
# and sampling location
# ------------------------------------
geomean <- function(x) {
	logx   <- log(x)
	mulogx <- mean(logx,na.rm=T)
	return( exp(mulogx) )
}
mu.qpcr <- tapply(d$enteroQPCRcce,d$beachcode,geomean)
mu.1600 <- tapply(d$entero1600cfu,d$beachcode,geomean)

# overall means
tot.mu.qpcr <- geomean(d$enteroQPCRcce)
tot.mu.1600 <- geomean(d$entero1600cfu)


# ------------------------------------
# combine results into a single table
# for each of the 1600 and qPCR assays
# ------------------------------------

total.1600  <- c(sum(tot.nd.1600),tot.nd.1600[1],tot.min.1600,tot.max.1600,tot.mu.1600)
total.qpcr  <- c(sum(tot.nd.qpcr),tot.nd.qpcr[1],tot.min.qpcr,tot.max.qpcr,tot.mu.qpcr)


wqtab.1600 <- cbind(
	rowSums(nd.1600),
	nd.1600[,1],
	min.1600,
	max.1600,
	mu.1600
)

wqtab.qpcr <- cbind(
	rowSums(nd.qpcr),
	nd.qpcr[,1],
	min.qpcr,
	max.qpcr,
	mu.qpcr
)

colnames(wqtab.1600) <- colnames(wqtab.qpcr) <- c("Nsamples","NonDetects","Min","Max","Geomean")
names(total.1600) <- names(total.qpcr) <- c("Nsamples","NonDetects","Min","Max","Geomean")

round(wqtab.1600,1)
round(wqtab.qpcr,1)

round(total.1600,1)
round(total.qpcr,1)
# ------------------------------------
# write summary tables to a file
# to be formatted in TeX
# ------------------------------------

save(wqtab.1600,total.1600,wqtab.qpcr,total.qpcr,file="~/dropbox/13beaches/aim1-results/rawoutput/aim1-wq-summary-tables.RData")






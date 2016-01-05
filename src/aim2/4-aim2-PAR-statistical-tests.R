



# --------------------------------------
# 4-aim2-PAR-statistical-tests.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# conduct formal statistical tests for
# effect modification of the population
# attributable risk (PAR) by age
# using bootstrap replicates
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------
rm(list=ls())


# --------------------------------------
# calculate differences in the bootstrap 
# replicates by different ages
# to compare 0-4 y vs 5-10 y and vs >10 y
# --------------------------------------

PARcomp <- function(x1,x2) {
	diff <- x1$stats[,1] - x2$stats[,1]
	mu.diff <- mean(diff)
	se.diff <- sd(diff)
	Z.diff <- mu.diff/se.diff
	P.diff <- 2*pnorm(-abs(Z.diff))
	res <- c(mu.diff,se.diff,Z.diff,P.diff)
	names(res) <- c("diff","SEdiff","Z","P")
	return(res)
}

# --------------------------------------
# Diarrhea
# --------------------------------------
load("~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARswimex-diar.RData")
PARcomp(PARswimex.diar.0to4,PARswimex.diar.5to10)
PARcomp(PARswimex.diar.0to4, PARswimex.diar.11plus)
PARcomp(PARswimex.diar.5to10, PARswimex.diar.11plus)

# --------------------------------------
# GI illness
# --------------------------------------
load("~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARswimex-gi.RData")
PARcomp(PARswimex.gi.0to4,PARswimex.gi.5to10)
PARcomp(PARswimex.gi.0to4, PARswimex.gi.11plus)
PARcomp(PARswimex.gi.5to10, PARswimex.gi.11plus)




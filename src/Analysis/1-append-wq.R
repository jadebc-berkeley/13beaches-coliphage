# ------------------------------------
#
# Create dataset with water quality
# sample-specific information
#
# version 1 
# ------------------------------------

# ------------------------------------
# preamble
# ------------------------------------
rm(list=ls())


# ------------------------------------
# Load the datasets
# ------------------------------------

neear <- read.csv("~/Documents/CRG/coliphage/13beaches-data/final/neear-wq-samples.csv",stringsAsFactors=F)
adm   <- read.csv("~/Documents/CRG/coliphage/13beaches-data/final/adm-wq-samples.csv",stringsAsFactors=F)
mb    <- read.csv("~/Documents/CRG/coliphage/13beaches-data/final/mb-wq-samples.csv",stringsAsFactors=F)

# ------------------------------------
# Import berm and groundwater flow 
# variables into the adm dataset
# ------------------------------------
epi=read.csv("~/Documents/CRG/coliphage/13beaches-data/final/13beaches-epi.csv")
epi$bermyn=ifelse(epi$berm=="Open",1,0)
epi$bermyn[epi$berm==""]=NA
epi$groundyn=ifelse(epi$groundwater=="Above median flow",1,0)
epi$groundyn[epi$groundwater==""]=NA
colnames(epi)[grep("intdate",colnames(epi))]="coldate"

# collapse rows to single value by beach and date
berm=aggregate(bermyn~beach+beachcode+coldate,data=epi,mean,na.rm=TRUE)
gw=aggregate(groundyn~beach+beachcode+coldate,data=epi,mean,na.rm=TRUE)

adm=merge(adm,berm,by=c("beach","beachcode","coldate"),all.x=TRUE)
adm=merge(adm,gw,by=c("beach","beachcode","coldate"),all.x=TRUE)

adm=adm[order(adm$beach,adm$beachcode,adm$coldate),]
colnames(adm)[grep("berm",colnames(adm))]="berm"
colnames(adm)[grep("ground",colnames(adm))]="groundwater"


# ------------------------------------
# Subset the data to common variables
# and append
# ------------------------------------

neear$stationid <- paste(neear$stationid)
neear$beachcode <- neear$beach
neear$sampleid <- neear$sampleid_neear1600
neear$coldate <- as.Date(neear$coldate,"%d%b%Y")
adm$coldate   <- as.Date(adm$coldate,"%d%b%Y")
mb$coldate    <- as.Date(mb$coldate,"%d%b%Y")


varselect.adm <- c("beach","beachcode","stationid","sampleid","coldate","fmc1601mpn",
               "fmc1601mpn_nd","fmc1602mpn","fmc1602mpn_nd","fpc1601mpn",
               "fpc1601mpn_nd","fpc1602mpn","fpc1602mpn_nd","entero1600cfu","entero1600cfu_nd",
               "berm","groundwater")

varselect.mb <- c("beach","beachcode","stationid","coldate","fmc1601mpn",
                   "fmc1601mpn_nd","fpc1601mpn","fpc1601mpn_nd","enteroELTmpn","enteroELTmpn_nd")

varselect.neear <- c("beach","beachcode","stationid","sampleid","coldate","fpc1601mpn",
                   "fpc1601mpn_nd","entero1600cfu","entero1600cfu_nd")


nd <- subset(neear,select=varselect.neear)
nd$fmc1601mpn=NA
nd$fmc1601mpn_nd=NA
nd$fmc1602mpn=NA
nd$fmc1602mpn_nd=NA
nd$fpc1602mpn=NA
nd$fpc1602mpn_nd=NA
nd$berm=NA
nd$groundwater=NA
nd$entero=nd$entero1600cfu
nd$entero_nd=nd$entero1600cfu_nd
nd$entero1600cfu=NULL
nd$entero1600cfu_nd=NULL

ad <- subset(adm,select=varselect.adm)
ad$entero=ad$entero1600cfu
ad$entero_nd=ad$entero1600cfu_nd
ad$entero1600cfu=NULL
ad$entero1600cfu_nd=NULL

mb <- subset(mb,select=varselect.mb)
mb$fmc1602mpn=NA
mb$fmc1602mpn_nd=NA
mb$fpc1602mpn=NA
mb$fpc1602mpn_nd=NA
mb$groundwater=NA
mb$berm=NA
mb$entero=mb$enteroELTmpn
mb$entero_nd=mb$enteroELTmpn_nd
mb$enteroELTmpn=NULL
mb$enteroELTmpn_nd=NULL
mb$sampleid=NA

wq <- rbind(nd,ad,mb)

# ------------------------------------
# Subset to beaches with coliphage assays
# ------------------------------------
beach.list=c("Avalon","Doheny","Malibu","Mission Bay","Goddard","Fairhope")
wq=subset(wq,wq$beach %in% beach.list)


# ------------------------------------
# impute non-detects at 0.1 before
# calculating summary statistics
# ------------------------------------

wq$entero[wq$entero_nd=="Below detection"] <- 0.1
wq$fmc1601mpn[wq$fmc1601mpn_nd=="Below detection"] <- 0.1
wq$fmc1602mpn[wq$fmc1602mpn_nd=="Below detection"] <- 0.1
wq$fpc1601mpn[wq$fpc1601mpn_nd=="Below detection"] <- 0.1
wq$fpc1602mpn[wq$fpc1602mpn_nd=="Below detection"] <- 0.1

# ------------------------------------
# rename variables for ease of coding
# ------------------------------------
wq$entero1600cfu=wq$entero
wq$fmc1601=wq$fmc1601mpn
wq$fmc1602=wq$fmc1602mpn
wq$fpc1601=wq$fpc1601mpn
wq$fpc1602=wq$fpc1602mpn

wq$fmc1601_nd=wq$fmc1601mpn_nd
wq$fmc1602_nd=wq$fmc1602mpn_nd
wq$fpc1601_nd=wq$fpc1601mpn_nd
wq$fpc1602_nd=wq$fpc1602mpn_nd

drops=c("entero1600cfu_nd","fmc1601mpn","fmc1602mpn","fpc1601mpn","fpc1602mpn")

wq=wq[,!names(wq)%in% drops]

# ------------------------------------
# make risk variable
# ------------------------------------
wq$risk = NA
wq$risk[wq$beach=="Fairhope"| wq$beach=="Goddard"|
          (wq$beach=="Doheny" & wq$berm==1)|
          (wq$beach=="Avalon" & wq$groundwater==1)]=1
wq$risk[wq$beach=="Malibu"|wq$beach=="Mission Bay"|
          (wq$beach=="Doheny" & wq$berm==0)|
          (wq$beach=="Avalon" & wq$groundwater==0)]=0
wq$risk=factor(wq$risk)
levels(wq$risk)=c("Low","High")

wq$berm=NULL
wq$groundwater=NULL


# ------------------------------------
# drop 1602 values for Malibu
# ------------------------------------
wq$fmc1601[wq$beach=="Malibu"]=NA
wq$fmc1602[wq$beach=="Malibu"]=NA
wq$fmc1601_nd[wq$beach=="Malibu"]=NA
wq$fmc1602_nd[wq$beach=="Malibu"]=NA

save(wq,file="~/Documents/CRG/coliphage/13beaches-data/temp/beaches-coli-ent-wq.RData")





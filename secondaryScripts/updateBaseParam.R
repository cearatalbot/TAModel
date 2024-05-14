##Base parameter changes from watershed spot calibration, adding params for dynamic phenology 
#By: Ceara J. Talbot


setwd("")
PFTtable<-read.csv("parameters/paramTable.csv",stringsAsFactors = F)

PFTtable[37,2:8]<-c(0, 0.0027, 0.0027, 0,0.0027,0,0.0027)
PFTtable[35,2:8]<-c(0.1, 0.1, 0.1, 0.1, 0, 0.1, 0)
PFTtable[54,2:8]<-0.001
PFTtable[12,2:8]<-c(5,5,4,4,4,2,4)
PFTtable[17,2:8]<-0.4
PFTtable[15,2:8]<-c(rep(0.07, times= 7))  #old val is .004
PFTtable[18,2:8]<-c(7, 7, 12.5, 12.5, 12, 8, 12)#c(10.9, 10.9, 12.5, 12.5, 12,8,15) #from Zhou et al 2016... NF modified to account for maritime -- 37% higher than temperate (so, scaled up by 0.5*.37=0.185)
PFTtable[55,1]<-"Lmin"
PFTtable[55,2:8]<-c(0.2, 4, 3, 0.2, 0.1, 0.2, 0.1)
PFTtable[20,2:8]<-c(70, 70, 80, 80, 30, 180, 30) ####reduce vpdlim
PFTtable[59,1]<-"Kvpd2"
PFTtable[59,2:8]<-2
PFTtable[3,2:8]<-0.08
PFTtable[9,4:8]<-c(0.35, 0.35, 0.62, 0.68, 0.62)

PFTtable[60,1]<-"Bp1"
PFTtable[60,2:8]<-0.55 #best is 0.55
PFTtable[61,1]<-"Bp2"
PFTtable[61,2:8]<-0.30 #best is 0.40

PFTtable[62,1]<-"W30"
PFTtable[62,2:8]<-40 #best was 60

PFTtable[33,2:8]<-10 #W20, best was 30

PFTtable[32,2:8]<-100 #Tstar
PFTtable[63,1]<-"Tstar2"
PFTtable[63,2:8]<-15 #used to be 50, best was 200

PFTtable[52,2:8]<-c(0.23, 0.23, 0.11, 0.11, 0.06, 0.04, 0.06)

#chilling period params: 
PFTtable[56,1]<-"pk"
PFTtable[56,2:8]<-c(0.05, 0.05, 0.05, 0.05, 0.65, 0.65, 0.65) #GR, SH, CR guessed based on debr
PFTtable[57,1]<-"pa"
PFTtable[57,2:8]<-c(0.65, 0.70, 0.70, 0.65, 0.65, 0.65, 0.65)
PFTtable[58,1]<-"pb"
PFTtable[58,2:8]<-c(200, 100, 100, 200, 100, 100, 100)

write.csv(PFTtable, "parameters/paramTable_2024.csv", row.names=F)

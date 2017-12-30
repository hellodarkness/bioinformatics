source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("statmod")
library(edgeR)
library(statmod)

source("https://bioconductor.org/biocLite.R")
biocLite("limma")
biocLite("methylMnM")
library(methylMnM)

setwd("/Users/lucydiamondsky/mem")
datafile<-("/Users/lucydiamondsky/mem")
filepath<-datafile[1]
dirwrite<-paste(setwd(getwd()),"/",sep="")

binlength<-500
file.cpgsite<-paste(filepath,"/hg19_CpGsite.bed",sep="")
writefile<-paste(dirwrite,"numallcpg.bed",sep="")
reportfile<-paste(dirwrite,"report_numallcpg.txt",sep="")
countcpgbin(file.cpgsite,file.blacklist=NULL,file.bin=NULL, writefile=writefile, reportfile = reportfile)

file<-paste(filepath,"/hg19_MRE_CpGsite.bed",sep="")
file1<-paste(filepath,"/hg19_CpGsite.bed",sep="")
allcpgfile<-paste(dirwrite,"numallcpg.bed",sep="")
five_Mre_CpGsite<-read.table(file, header=FALSE, as.is=TRUE)
four_Mre_CpGsite<-five_Mre_CpGsite[five_Mre_CpGsite[,4]!="ACGT",]
mrecpg.site<-four_Mre_CpGsite[four_Mre_CpGsite[,4]!="CGCG",]
writefile<-paste(dirwrite,"three_mre_cpg.bed",sep="")
countMREcpgbin(mrecpg.site,file.allcpgsite=file1,file.bin=allcpgfile,writefile=writefile)

file3<-paste(filepath,"/H_2034_Medip.extended.bed",sep="")
allcpgfile<-paste(dirwrite,"numallcpg.bed",sep="")
writefile<-paste(dirwrite,"H_2034_Medip_num500.bed",sep="")
reportfile<-paste(dirwrite,"H_2034_Medip_num500_report.txt",sep="")
countMeDIPbin(file.Medipsite=file3,file.blacklist=NULL,file.bin=allcpgfile,file.CNV=NULL,writefile=writefile, reportfile = reportfile)

file4<-paste(filepath,"/human2034_Mre.bed",sep="")
writefile<-paste(dirwrite,"human2034_MRE_num500.bed",sep="")
reportfile<-paste(dirwrite,"human2034_MRE_num500_report.bed",sep="")
countMREbin(file.MREsite=file4,file.blacklist=NULL, file.bin=allcpgfile,file.CNV=NULL,writefile=writefile, reportfile = reportfile

file5<-paste(filepath,"/H_2080_Medip.extended.bed",sep="")
allcpgfile<-paste(dirwrite,"numallcpg.bed",sep="")
writefile<-paste(dirwrite,"H_2080_MeDIP_num500.bed",sep="")
reportfile<-paste(dirwrite,"H_2080_MeDIP_num500_report.txt",sep="")
countMeDIPbin(file.Medipsite=file5,file.blacklist=NULL,file.bin=allcpgfile,file.CNV=NULL,writefile=writefile, reportfile = reportfile)

file6<-paste(filepath,"/human2080_mre.bed",sep="")
writefile<-paste(dirwrite,"human2080_MRE_num500.bed",sep="")
reportfile<-paste(dirwrite,"human2080_MRE_num500_report.txt",sep="")
countMREbin(file.MREsite=file6,file.blacklist=NULL, file.bin=allcpgfile,file.CNV=NULL,writefile=writefile, reportfile = reportfile)

datafile1<-paste(dirwrite,"H_2034_Medip_num500.bed",sep="")
datafile2<-paste(dirwrite,"H_2080_MeDIP_num500.bed",sep="")
datafile3<-paste(dirwrite,"human2034_MRE_num500.bed",sep="")
datafile4<-paste(dirwrite,"human2080_MRE_num500.bed",sep="")
datafile<-c(datafile1,datafile2,datafile3,datafile4)
chrstring<-NULL
cpgfile<-paste(dirwrite,"numallcpg.bed",sep="")
mrecpgfile<-paste(dirwrite,"three_mre_cpg.bed",sep="")
writefile<-paste(dirwrite,"pval_H2034_H2080.bed",sep="")
reportfile<-paste(dirwrite,"report_pval_H2034_H2080.txt",sep="")
MnM.test(file.dataset=datafile,chrstring=chrstring,file.cpgbin=cpgfile,file.mrecpgbin=mrecpgfile,writefile=writefile,reportfile = reportfile)

datafile<-paste(dirwrite,"pval_H2034_H2080.bed",sep="")
writefile<-paste(dirwrite,"q_H2034_H2080.bed",sep="")
reportfile<-paste(dirwrite,"report_q_H2034_H2080.bed",sep="")
MnM.qvalue(datafile,writefile,reportfile)

file<-paste(dirwrite,"q_H2034_H2080.bed",sep="")
frames<-read.table(file, header=TRUE,sep="\t", as.is=TRUE)
DMR<-MnM.selectDMR(frames =frames , up = 1.45, down =1/1.45, p.value.MM = 0.01, p.value.SAGE = 0.01,q.value = 0.01,cutoff="q-value", quant= 0.6)
writefile<-paste(dirwrite,"DMR_e5_H2034_H2080.bed",sep="")
write.table(DMR, writefile,sep="\t", quote=FALSE,row.names=FALSE)


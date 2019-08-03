options(warn=-1)
library("coloc")

#cfile="VLDLR_coloc.txt"
#N=172925

args = commandArgs(trailingOnly=TRUE)
cfile=args[1]
N=as.numeric(args[2])

data = read.csv(cfile, sep=" ")
dataset1=list(beta=data$gtex.beta, varbeta=data$gtex.varbeta, type="quant", sdY=1)
dataset2=list(beta=data$gwas.beta, varbeta=data$gwas.varbeta, MAF=data$gwas.maf, N=N, type="quant")

myres = coloc.abf(dataset1, dataset2)
write.csv(myres[2], paste(cfile, ".coloc.txt", sep=""))
x=myres[2][[1]]
bestsnpdata=x[x$SNP.PP.H4==max(x$SNP.PP.H4),]
print(bestsnpdata)

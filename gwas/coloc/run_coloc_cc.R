options(warn=-1)

library("coloc")

args = commandArgs(trailingOnly=TRUE)
cfile=args[1]
S=as.numeric(args[2])

data = read.csv(cfile, sep=" ")
dataset1=list(beta=data$gtex.beta, varbeta=data$gtex.varbeta, type="quant", sdY=1)
dataset2=list(beta=data$gwas.beta, varbeta=data$gwas.varbeta, MAF=data$gwas.maf, s=S, type="cc")

myres = coloc.abf(dataset1, dataset2)
write.csv(myres[2], paste(cfile, ".coloc.txt", sep=""))
x=myres[2][[1]]
bestsnpdata=x[x$SNP.PP.H4==max(x$SNP.PP.H4),]
print(bestsnpdata)

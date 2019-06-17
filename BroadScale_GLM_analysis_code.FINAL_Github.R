#Code used to generate Supplementary Tables 2 and 3
#in Li et al. 2019 given source data
#Nick Altemose
#4 April 2019

library(MASS)
library(Hmisc)
library(scales)
library(corrplot)

#load the matrix of genome-wide predictor/response variable values into a data frame called df1
#binary r file below contains the same data as the source data file
load("GLM_input_data.rdata")

#note: there are
#2500 total crossovers
#1515 crossovers determined to be controlled by the Humanized PRDM9 allele
#768  crossovers determined to be controlled by the Cast PRDM9 allele
#
#1575 total non-crossovers
#878 non-crossovers determined to be controlled by the Humanized PRDM9 allele
#562 non-crossovers determined to be controlled by the Cast PRDM9 allele


########Neg Binom family GLMs for combined, CO, NCO
#first examine dispersion
#poisson and quasipoisson show overdispersion, neg binom seems like the best fit (see plot)
#use default log link for GLMs

#plot histogram of combined CO+NCO event counts
hist(breaks=seq(0,30,1),df1$COplusNCO,freq=F,col=rgb(1,0,0,0.5),ylim=c(0,0.15),main="Dist of CO+NCO counts in 5-Mb bins")
#add neg binom in blue
rnb = rnbinom(1000,size=3,mu=mean(df1$COplusNCO)) 
hist(rnb,breaks=seq(0,round(max(rnb))+1,1),add=T,col=rgb(0,0,1,0.5),freq=F)
#add poisson in green
rp = rpois(1000,lambda=mean(df1$COplusNCO)) 
hist(rp,breaks=seq(0,round(max(rp))+1,1),add=T,col=rgb(0,1,0,0.5),freq=F)
legend("topright",col=c('red','blue','green'),legend=c('Empirical','Neg Binom','Poisson'),pch=15)


###### plot correlation matrix (Supp Table 3 upper panel)

res2<-rcorr(as.matrix(df1))
corrplot(res2$r, type="lower", order="original", p.mat = res2$P, sig.level = 0.0001, insig = "pch",tl.col = "black", tl.srt = 45,addCoef.col='black',number.cex=0.7,addCoefasPercent=T,method="color",pch.col="gray")


############ run full models for each response variable, then plot results (Supp Table 3 lower table)

#GLM for combined CO+NCO both alleles
g = glm.nb(COplusNCO ~ DMC1 + H3K4me3 + H3K9me3 + ChIP_Input + Satellite + SNPs + Genes + GC + Compartment + CenDist + TelDist,data = df1)
summary(g)
#drop1(g, test="LRT")
pvals = coef(summary(g))[,4]
coeffs = g$coefficients

#GLM for COs both alleles
gCO = glm.nb(CO ~ DMC1 + H3K4me3 + H3K9me3 + ChIP_Input + Satellite + SNPs + Genes + GC + Compartment + CenDist + TelDist,data = df1)
summary(gCO)
#drop1(gCO, test="LRT")
pvals = cbind(pvals,coef(summary(gCO))[,4])
coeffs = cbind(coeffs,gCO$coefficients)

#GLM for NCOs both alleles
gNCO = glm.nb(NCO ~ DMC1 + H3K4me3 + H3K9me3 + ChIP_Input + Satellite + SNPs + Genes + GC + Compartment + CenDist + TelDist,data = df1)
summary(gNCO)
#drop1(gNCO, test="LRT")
pvals = cbind(pvals,coef(summary(gNCO))[,4])
coeffs = cbind(coeffs,gNCO$coefficients)


#GLM for combined CO+NCO Humanized allele
g = glm.nb(Hum_COplusNCO ~ Hum_DMC1 + H3K4me3 + H3K9me3 + ChIP_Input + Satellite + SNPs + Genes  + Genes + GC + Compartment + CenDist + TelDist,data = df1)
summary(g)
#drop1(g, test="LRT")
pvals = cbind(pvals,coef(summary(g))[,4])
coeffs = cbind(coeffs,g$coefficients)

#GLM for COs Humanized allele
gCO = glm.nb(Hum_CO ~ Hum_DMC1 + H3K4me3 + H3K9me3 + ChIP_Input + Satellite + SNPs + Genes + GC + Compartment + CenDist + TelDist,data = df1)
summary(gCO)
#drop1(gCO, test="LRT")
pvals = cbind(pvals,coef(summary(gCO))[,4])
coeffs = cbind(coeffs,gCO$coefficients)

#GLM for NCOs Humanized allele
gNCO = glm.nb(Hum_NCO ~ Hum_DMC1 + H3K4me3 + H3K9me3 + ChIP_Input + Satellite + SNPs + Genes + GC + Compartment+ CenDist + TelDist,data = df1)
summary(gNCO)
#drop1(gNCO, test="LRT")
pvals = cbind(pvals,coef(summary(gNCO))[,4])
coeffs = cbind(coeffs,gNCO$coefficients)


#GLM for combined CO+NCO Cast allele
g = glm.nb(Cast_COplusNCO ~ Cast_DMC1 + H3K4me3 + H3K9me3 + ChIP_Input + Satellite + SNPs + Genes   + GC + Compartment+  CenDist + TelDist,data = df1)
summary(g)
#drop1(g, test="LRT")
pvals = cbind(pvals,coef(summary(g))[,4])
coeffs = cbind(coeffs,g$coefficients)


#GLM for COs Cast allele
gCO = glm.nb(Cast_CO ~ Cast_DMC1 + H3K4me3 + H3K9me3 + ChIP_Input + Satellite + Genes + SNPs  + GC + Compartment+ CenDist + TelDist,data = df1)
summary(gCO)
#drop1(gCO, test="LRT")
pvals = cbind(pvals,coef(summary(gCO))[,4])
coeffs = cbind(coeffs,gCO$coefficients)


#GLM for NCOs Cast allele
gNCO = glm.nb(Cast_NCO ~ Cast_DMC1 + H3K4me3 + H3K9me3 + ChIP_Input + Satellite + Genes + SNPs + Genes + GC + Compartment+  CenDist + TelDist,data = df1)
summary(gNCO)
#drop1(gNCO, test="LRT")
pvals = cbind(pvals,coef(summary(gNCO))[,4])
coeffs = cbind(coeffs,gNCO$coefficients)


pvals2 = t(apply(pvals,1,function(x) as.numeric(scientific(x))))
pvals2[pvals2>0.1]=NA
coeff2=coeffs
coeff2[coeff2<0]=-1
coeff2[coeff2>0]=1
combo=round(-log(pvals2,base=10)*coeff2,2)
colnames(combo)=c("COs+NCOs ALL (n=4075)","COs ALL (n=2500)","NCOs ALL (n=1575)","COs+NCOs Hum (n=2393)","COs Hum (n=1515)","NCOs Hum (n=878)","COs+NCOs Cast (n=1330)","COs Cast (n=768)","NCOs Cast (n=562)")
corrplot(t(combo), type="full", order="original", is.corr=F,tl.col = "black",addCoef.col='black',number.cex=1,method="color",pch.col="gray",number.digits=5,na.label="square",na.label.col="white",tl.srt = 45,addgrid.col='gray',cl.pos='n')




############ trimming variables illustrates that including GC content as an explanatory variable reverses the direction of association of chromatin compartment A with recombination rate

g = glm.nb(COplusNCO ~ Compartment,data = df1)
summary(g)

g = glm.nb(COplusNCO ~ GC + Compartment,data = df1)
summary(g)

g = glm.nb(COplusNCO ~ Genes + GC + Compartment,data = df1)
summary(g)

g = glm.nb(COplusNCO ~ DMC1 + Genes + GC + Compartment,data = df1)
summary(g)

g = glm.nb(COplusNCO ~ DMC1 + Genes + GC + Compartment + TelDist,data = df1)
summary(g)


g = glm.nb(Hum_COplusNCO ~ Compartment,data = df1)
summary(g)

g = glm.nb(Hum_COplusNCO ~ GC + Compartment,data = df1)
summary(g)

g = glm.nb(Hum_COplusNCO ~ Genes + GC + Compartment,data = df1)
summary(g)

g = glm.nb(Hum_COplusNCO ~ Hum_DMC1 + Genes + GC + Compartment,data = df1)
summary(g)

g = glm.nb(Hum_COplusNCO ~ Hum_DMC1 + Genes + GC + Compartment + TelDist,data = df1)
summary(g)


g = glm.nb(Cast_COplusNCO ~ Compartment,data = df1)
summary(g)

g = glm.nb(Cast_COplusNCO ~ GC + Compartment,data = df1)
summary(g)

g = glm.nb(Cast_COplusNCO ~ Genes + GC + Compartment,data = df1)
summary(g)

g = glm.nb(Cast_COplusNCO ~ Cast_DMC1 + Genes + GC + Compartment,data = df1)
summary(g)

g = glm.nb(Cast_COplusNCO ~ Cast_DMC1 + Genes + GC + Compartment + TelDist,data = df1)
summary(g)





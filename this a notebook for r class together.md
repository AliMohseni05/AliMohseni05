
"""this a notebook for r class together"""
# jalese aval 
```
#varible and vector 
x<- 6
print(x)
x=c(2,5,6,9,10)
y<-2x
x*2
x+5
y=c(1.2,12,13,14,3)
y^2
mean(x)
sd(x)
summary(x,y)
var(x)
cor(x,y)
median(x)
median(y)
plot(1:10, type="1", col="blue")
plot(1:10, type="5", col="blue")


# Create a data frame with the predictor and response variables
data <- data.frame(x = c(151, 174, 138, 186, 128, 136, 179, 163, 152, 131),
                   y = c(63, 81, 56, 91, 47, 57, 76, 72, 62, 48))

# Fit a linear regression model
model <- lm(y ~ x, data = data)

# Print the summary of the model
summary(model)

# Predict the response variable for a new value of the predictor variable
new_x <- 20
new_y <- predict(model, newdata = data.frame(x = new_x))
print(new_y)

# Visualize the regression graphically
plot(data$x, data$y, main = "Linear Regression Example", xlab = "Predictor Variable", ylab = "Response Variable")
abline(model, col = "red")

#-----------------------------

#make data set
my_dataset<- data.frame(hx<-c(1,2,3,4,5,6,7,8,9,10),vy<-c(1,2,3,4,5,6,7,8,9,10))
#fit model from a liner reg and print the summary
reglin <- lm(hx ~ vy, data = my_dataset)
summary(reglin)

#show plate
plot(my_dataset$vy, my_dataset$hx, main = "titel", xlab = "v", ylab = "h")
#show line
abline(reglin, col = "blue", type="l", lwd=2, lty=1)
```
# jalse 2 
```
#-----------------------------
#jalse 2                     #
#-----------------------------
#make vector and remove NA data print them 
rawdata<- c(NA,21,23,23,45,NA,23,NA,98)
rna<-na.omit(rawdata)
cat(rawdata)
cat(rna)
v1<-c(1,2,3,5,NA)
v2<-c(12,42,42,54,44,44,3,3,3,4)
v3<-v1 + v2 # add to vector together in like numpuy in python

v4<-v2-v1
cat(v3)
cat(v4)
v1[v2]  #show v1 according list v2 
v1[-v2] #show v1 according some thing that excise in v2 
v2[v1]  # this one give error
?c()
x<-10.5
class(x)
x2<-12L
class(x2)
x3<-"this sring pyhone"
class(x3)
a<- 25+3i
class(a)
x4<-charToRaw("salam")
class(x4)
x5<-rawToChar(61 6c 69)
x5
my.ege<-32
my.ege
# we can use True and False as T and F 
#you can make coment by cont+shif+c  but it dosent work for me
#vector 1d the can be numice and other at same time 
vector_e<-c(12,19,19,19)
vector_e2<-c("hi","ali","m",2,3L,"other","حتی فارسی")
vector_e3<-c(TRUE,FALSE,T,F) #TRUE better 
vector_e3
class(vector_e3)
vector_e[1:5] 
v1<-c(1,2,4)
v2<-c("dna","rna","pro","RNAseq")
v2[v1]#this function form get number v1 
v1[v2]#this wouldent fine any thikng becusea it dont have any number 
v2[-v1]#show us somethink we dont have in v1 
newvec<-seq(1:10)
newvec2<- 1:20
newvec3<-rep(seq(1:10),11)
newvec3
vc<- rep(1:12,t=12) #this one can be by t or time or none it just need the 
#number of rep

s<- seq(1,12,by=0.1)
assign("A",c(100,12,14,15))#for make valiue by this way
A #A;vc for run two varialbe at onetime 
ad<-c(1,1,2,2)
ad2<-c(2,2,4,NA)
ad+ad2;ad*ad2;ad^ad2 #thay are must be same derae ot do opertion 
adsum<-sum(ad2,na.rm=T)
ad2[!is.na(ad2)]
my.weight<-c(12,NA,14)
mv<-na.omit(my.weight) #this one show us which one is missing or NA
print(mv) #and return us a funtion 
max(ad)

#loop
for(i in 1:10){print("helo world")}
for(x in 1:10){print(x)}
for(x in 1:10){print(2*x-1)}

install.packages("rmp") 
library(rmp) #class packet 
#na
###################################
#make matrix
vs1<-c(1,2,4)
vs2<-c(2,4,5)
vs3<-c(5,5,6)
rbind(vs1,vs2,vs3)#they should be equal 
cbind(vs1,vs2,vs3)#r= Row and c=colem
#vector bind by r of c 
mat<-rbind(vs1,vs2,vs3)
(mat)
t(mat) #transpose mat 
det(mat) #derteminal
mat2<-rbind(c(1,3),c(5,6)) #make matrix by another 
mat2
eigen(mat2) #eigen value thay use is pca analise
####
#practice mode 
mat3<-rbind(c(1,4),c(3,4))
mat3
det(mat3)
av<-c(1,0)
av2<-c(0,1)
mata<-rbind(av,av2)
mata
det(mata)
eigen(mata)

####
matc<- matrix(1:5,nrow=5,ncol=5, byrow=TRUE)
matc
det(matc)
matr<- matrix(runif(n=5,2,10),nrow=5,ncol=5, byrow=F,dimnames =list(c("a","b","c","e","f"),c("col1","clo2","col3","col4","col5")))
matr
det(matr)
#####
#data farame 
df<-data.frame(genom=c("rosha","srosh","ahamdi")
,tol= c(12,10,12),h= c(151, 174, 138),g=c(1,0,1))
# here use = importent
#df$mar2=c(1,1,3)

lost_row<-data.frame(genom= c("newmar"),tol= c(34),h= c(76),g=c(1))
pdf<-rbind(df,lost_row)

new_row= data.frame(genom=c("d"),tol= c(1),h= c(31),g=c(1)
pdf<-rbind(df,new_row) #row.name=FLASL 
pdf
```
# jalse 3.1 
```
#-----------------------------
#jalse 3.1                   #
#-----------------------------
organ<-data.frame(genid=c(1,3,6),organ=c("liver","kideny","hreat")
               ,experss= c(7.3,11.9,9.1),g=c(1,2,1),f=c(1,12,123))
(organ)
#add one col 
organ$mar2=c(1,1,3,4)

#for instal go this site 
#https://www.bioconductor.org/
#https://cran.r-project.org/
#for instal packagees use this one  
install.packages("ggplot2") #for install form repos =www.ulr.com
#a specic site 
#packegess form https://www.bioconductor.org/
install.packages(c("tidyr","readr"))#install more than one packages at ones
#if funtion be impty that would worke 
#and it will one new windon
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")


library(ggplot2)

#------
#3.2  #
#------
getwd() #where am i ?
setwd("G:/R CLASS practice/box") #go to here 
df<-read.csv("G://R CLASS practice/box/Boxplot_data.csv") # it one two // two ec
#c<- reasd.txt("G://R CLASS practice/box/data-WS015.txt")
df$id<- 1:nrow(df)
names(df)
head(df,10)
install.packages("reshape2")
library(ggplot2)
library(reshape2)
df.t<-melt(df,id.vars ="id",variable.name = "varlable",value.name = "valuus")

ggplot(df.t,aes(x=varlable,y=valuus,fill=varlable))+geom_boxplot(color= "black",notch = T)+
  scale_fill_manual(values = c("green","blue","orange","red"))+
  labs(title="my box",x="x data",y="y data")





#----------
# tamrin  #
#----------
#MAKE RAND NUMBER 
nomre<-runif(n=1,1,20)
if (nomre >=10) {print("you are passt")}
if (nomre <10) {print ("you are fall")}
nomre

#make range nomre
nomre<-runif(n=1,1,20)
if(nomre<=5){
  print("f")
} else if(nomre>5 & nomre<10 ){
  print("e")
} else if (nomre>=10 & nomre<=15){
  print("C")
} 
if(nomre>15){
  print("A")
}
nomre
#---------
melted.mind<-melt(mind,id=c("ID","BTW"), 
                  variable.name = "WTN", 
                  value.name = "Results")

head(melted.mind)   

##   ID BTW  WTN  Results
## 1  1   1 WTN1 6.939619
## 2  2   2 WTN1 5.978156
## 3  3   1 WTN1 4.515321
## 4  4   2 WTN1 3.502277
## 5  5   1 WTN1 4.972581
## 6  6   2 WTN1 3.275813
```
# aov 
```
#.....................................#
#                aov                  #
#.....................................#
https://www.scribbr.com/statistics/anova-in-r/

crop.data <- read.csv("G:/R CLASS practice/f/crop.data.csv", header = TRUE, 
                      colClasses = c("factor", "factor", "factor", "numeric"))
one.way <- aov(yield ~ fertilizer, data = crop.data)

summary(one.way)

two.way <- aov(yield ~ fertilizer + density, data = crop.data)

summary(two.way)

interaction <- aov(yield ~ fertilizer*density, data = crop.data)

summary(interaction)

blocking <- aov(yield ~ fertilizer + density + block, data = crop.data)

summary(blocking)
```
# LEA 
```
#.....................................#
#               LEA                   #
#.....................................#
getwd()
setwd("G:/R CLASS practice/lea2")
install.packages("LEA")
library(LEA)
#data("tutorial")
c<-read.geno( "CART-p-o-t-g.geno")
write.lfmm(c, "genotypes.lfmm")
write.geno(c, "genotypes.geno")
#write.env(tutorial.C, "gradients.env")
pc = pca("genotypes.lfmm", scale = TRUE)
tw = tracy.widom(pc)
plot(tw$percentage, pch = 19, col = "darkblue", cex = .8)
# main options
# K = number of ancestral populations
# entropy = TRUE computes the cross-entropy criterion,
# CPU = 4 is the number of CPU used (hidden input)
project = NULL
project = snmf("genotypes.geno",
               K = 1:10,
               entropy = TRUE,
               repetitions = 10,
               project = "new")
plot(project, col = "blue", pch = 19, cex = 1.2)
best = which.min(cross.entropy(project, K = 4))
my.colors <- c("tomato", "lightblue",
               "olivedrab", "gold")
barchart(project, K = 4, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)
#...................
p = snmf.pvalues(project,
                 entropy = TRUE,
                 ploidy = 2,
                 K = 4)
pvalues = p$pvalues
par(mfrow = c(2,1))
hist(pvalues, col = "orange")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .5)
#
#
# Missing genotype imputation using snmf 
#
# creation of a genotype matrix with missing genotypes
dat = as.numeric(tutorial.R)
dat[sample(1:length(dat), 100)] <- 9
dat <- matrix(dat, nrow = 50, ncol = 400)
write.lfmm(dat, "genoM.lfmm")
## [1] "genoM.lfmm"
project.missing = snmf("genoM.lfmm", K = 3,
                       entropy = TRUE, repetitions = 10,
                       project = "new")

####
# select the run with the lowest cross-entropy value
best = which.min(cross.entropy(project.missing, K = 3))
# Impute the missing genotypes
impute(project.missing, "genoM.lfmm",
       method = 'mode', K = 3, run = best)
## Missing genotype imputation for K = 4
## Missing genotype imputation for run = 1
## Results are written in the file: genoM.lfmm_imputed.lfmm
# Proportion of correct imputation results
dat.imp = read.lfmm("genoM.lfmm_imputed.lfmm")
mean( tutorial.R[dat == 9] == dat.imp[dat == 9] )
## [1] 0.83

project = NULL
project = lfmm("genotypes.lfmm",
               "gradients.env",
               K = 6,
               repetitions = 5,
               project = "new")

p = lfmm.pvalues(project, K = 6)
pvalues = p$pvalues

par(mfrow = c(2,1))
hist(pvalues, col = "lightblue")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .7)



####
for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  12
  print(paste("Expected FDR:", alpha))
  L = length(pvalues)
  # return a list of candidates with expected FDR alpha.
  # Benjamini-Hochberg's algorithm:
  w = which(sort(pvalues) < alpha * (1:L) / L)
  candidates = order(pvalues)[w]
  # estimated FDR and True Positive Rate
  Lc = length(candidates)
  estimated.FDR = sum(candidates <= 350)/Lc
  print(paste("Observed FDR:",
              round(estimated.FDR, digits = 2)))
  estimated.TPR = sum(candidates > 350)/50
  print(paste("Estimated TPR:",
              round(estimated.TPR, digits = 2)))
}
## [1] "Expected FDR: 0.05"
## [1] "Observed FDR: 0.05"
## [1] "Estimated TPR: 0.84"
## [1] "Expected FDR: 0.1"
## [1] "Observed FDR: 0.04"
## [1] "Estimated TPR: 0.86"
## [1] "Expected FDR: 0.15"
## [1] "Observed FDR: 0.06"
## [1] "Estimated TPR: 0.88"
## [1] "Expected FDR: 0.2"
## [1] "Observed FDR: 0.06"
## [1] "Estimated TPR: 0.92"


# load simulated data
data("offset_example")
# 200 diploid individuals genotyped at 510 SNP
Y <- offset_example$geno
# 4 environmental variables
X <- offset_example$env
mod.lfmm2 <- lfmm2(input = Y, env = X, K = 2)

# Simulate non-null effect sizes for 10 target loci
#individuals
n = 100
#loci
L = 1000
# Environmental variable
X = as.matrix(rnorm(n))
# effect sizes
B = rep(0, L)
target = sample(1:L, 10)
# GEA significance test
# showing the K = 2 estimated factors
plot(mod.lfmm2@U, col = "grey", pch = 20,
     xlab = "Factor 1",
     ylab = "Factor 2")

B[target] = runif(10, -10, 10)
# Create 3 hidden factors and their loadings
U = t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X)) +
  matrix(rnorm(3*n), ncol = 3)
V <- matrix(rnorm(3*L), ncol = 3)

pv <- lfmm2.test(object = mod.lfmm2,
                 input = Y,
                 env = X,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")

# Fitting an LFMM with K = 3 factors
mod <- lfmm2(input = Y, env = X, K = 3)

# Computing P-values and plotting their minus log10 values
pv <- lfmm2.test(object = mod,
                 input = Y,
                 env = X,
                 linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .6, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")

data("offset_example")
Y <- offset_example$geno
X <- offset_example$env
X.pred <- offset_example$env.pred

g.gap.scaled <- genetic.gap(input = Y,
                            env = X,
                            pred.env = X.pred,
                            scale = TRUE,
                            K = 2)
g.gap.scaled$vectors[,1:2]^2
## [,1] [,2]
## [1,] 0.513933202 0.2127500577
## [2,] 0.179886467 0.6453681204
## [3,] 0.302617773 0.0005067115
## [4,] 0.003562557 0.1413751104

par(mfrow = c(1,2))
barplot(g.gap.scaled$eigenvalues,
        col = "orange",
        xlab = "Axes",
        ylab = "Eigenvalues")
Delta = X[,1:2] - X.pred[,1:2]
squared.env.dist = rowSums(Delta^2)
plot(squared.env.dist, g.gap.scaled$offset, cex = .6)
```
# agricolae 
```
#.....................................#
#             agricolae               #
#.....................................#
install.packages("agricolae")
library(agricolae)
# 4 treatments and k=3 size block
trt<-c("A","B","C","D")
k<-3
outdesign<-design.bib(trt,k,serie=2,seed =41,kinds ="Super-Duper") # seed = 41
print(outdesign$parameters)
book<-outdesign$book
plots <-as.numeric(book[,1])
matrix(plots,byrow=TRUE,ncol=k)
print(outdesign$sketch)
# write in hard disk
write.csv(book,"book.csv", row.names=FALSE)
file.show("book.csv")
```
# ggcorrplo   
```
#.....................................#
#             ggcorrplo               #
#.....................................#

install.packages("ggcorrplot")
library(ggplot2)
library(ggcorrplot)
library(corrplot)

setwd("G://R CLASS practice/box/data-WS015.txt")
mydata <- read.table("G://R CLASS practice/box/data-WS015.txt" , header = TRUE)
M = cor(mydata)
p.mat <- cor_pmat(mydata)
head(p.mat[, 1:4])
# Compute a correlation matrix
corr <- round(cor(mydata), 1)
head(corr[, 1:6])

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(mydata)
head(p.mat[, 1:4])

# Visualize the correlation matrix
# --------------------------------
# method = "square" (default)
ggcorrplot(corr)

# Change colors and theme
# --------------------------------
# Argument colors
ggcorrplot(corr,method = "circle", hc.order = TRUE, type = "lower",
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"))
ggcorrplot(corr, hc.order = TRUE, type = "lower",
           lab = TRUE)


# Add correlation significance level
# --------------------------------
# Argument p.mat
# Barring the no significant coefficient
ggcorrplot(corr, hc.order = TRUE,
           type = "lower", p.mat = p.mat)



```

# blup
```
library(lme4)
library(dplyr)
#----------------------#
#       call data      #
#----------------------#
getwd()
setwd("G:/R CLASS practice/LS")
datak<-read.csv("blup4.csv",header =T)

str(datak)
head(datak)
tail(datak)

#-------------------------------------------------
rep=as.factor(datak$rep)
Year=as.factor(datak$Year)
Location=as.factor(datak$Location)
block=as.factor(datak$block)
plot=as.factor(datak$plot)
env=as.factor(datak$env)
G=as.factor(datak$G)
Y2 = as.numeric(datak$Y)
ph = as.numeric(datak$ph)
ssT = as.numeric(datak$ssT)
LOC=as.factor(datak$Location)
YEAR=as.factor(datak$Year)
REP=as.factor(datak$rep)
LINE=as.factor(datak$G)
#-----------------------------------------------
hist(Y2, col="gold")
boxplot(Y2~LOC, xlab="Location", ylab="Degrees Brix", main="Degrees Brix by Location", col="pink")

## BLUPS
# fit the model
brixmodel = lmer(Y2~ (1|LINE) + (1|LOC) + (1|YEAR) + (1|LINE:LOC) + (1|LINE:YEAR))


# estimate BLUPS
brixblup = ranef(brixmodel)
# look at output structure
str(brixblup)
# extract blup for line
brixlineblup = brixblup$LINE
# see the structure of the blup for each line
str(brixlineblup)
# save the brixlineblup output to a separate .csv file
write.csv(brixlineblup, file="BrixLineBLUPS.csv")

## Creating plots with the BLUPs
# Create a numeric vector with the BLUP for each line
LINEBLUP = brixlineblup[,1]
# Create a histogram with the BLUP for each line
hist(LINEBLUP, col="green")

## Compare BLUP to line averages on a scatterplot
lmean = tapply(Y2, LINE, na.rm=T, mean)
plot(LINEBLUP, lmean, col="black")

```


جمعت ریسل f7 
نوع جامعه  ریل 
مارک دارت :: غالب  است فاید 0 و 2 دارد 
فایل های فنوتیپ های  ژنتیک داریم 

در کویتل می توانیم فاصله ژنتیکی می توانمی به دست اوریم  در gwas  نمی شود
 
دو تا فایل 
فایل اول به اساس دستور عمل برنامه 
فایل دوم  بر اساس فایل مارک های شرک 
اسم مارک های و عدد های و ........ 
برگه دوم فقط شماره کروموز های می گذاریم برای na  0 می گذاریم  اسم مارک های به صورت عدد باشد. 
opee file 
*bin
choose data 
#yon cand see it 

#you can change format 
by missing rade 
read in artcel 
run it or bind it 

it will give new file two file 
a b  h for hybrid 

*map 
run iit in 3 steps 
and outputing 


pib* 

by QTL maping 

contle +shift and <- 
chosen all
!!!! 

-1 lost data 

icm
icm ip


start forme up top of  appitlication 



توالی هر مارک به دست به اوریم به .... در روی ژونم گندم بلست کردیم به جای مارک رو کرومزم به دست می آوریم. 
در سایت ncbi و توالی های خیلی شبیه هستند را انتخاب می کنیم. 
graingenes 
https://wheat.pw.usda.gov/GG3/
https://wheat.pw.usda.gov/blast/




https://knetminer.com/cereals/







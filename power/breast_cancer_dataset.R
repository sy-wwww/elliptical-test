file = getwd()
source(paste0(file,"/tests.R"))

# BiocManager::install("cancerdata")
library(cancerdata)
data("VIJVER")
a3 = assayData(VIJVER)
a3=a3$exprs

# Remove genes with missing values
a3_c = na.omit(a3)
dim(a3_c)


# Transform it into an n*p matrix
a3_c = t(a3_c)

# Centering 
data = scale(a3_c, center = TRUE, scale = FALSE)
data = as.matrix(data)


p = 200 #200, 300, 500
N = 100


re_elliptical = matrix(0,N,3)
re_normality = matrix(0,N,1)
for (i in 1:N) {
  set.seed(i)
  l = sample(22223,p)
  data1 = data[,l] #n*p
  
  re_normality[i,] = normality_test(data1)
  re_elliptical[i,] = elliptical_test(data1)
}


median(re_normality[,1])
median(re_elliptical[,3])

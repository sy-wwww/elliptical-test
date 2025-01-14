file = getwd()
source(paste0(file,"/tests.R"))


breastCancer = read.csv(paste0(file,"/elliptical_breast_cancer.csv"), header = TRUE) 
data = as.matrix(breastCancer)


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


mean(re_normality[,1])
mean(re_elliptical[,3])
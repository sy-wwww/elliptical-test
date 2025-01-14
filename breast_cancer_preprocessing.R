BiocManager::install("cancerdata")
library(cancerdata)
data("VIJVER")
a3 = assayData(VIJVER)
a3=a3$exprs
# dim(a3)
# rowSums(is.na(a3))
# colSums(is.na(a3))

# Remove genes with missing values
a3_c = na.omit(a3)
dim(a3_c)


# Transform it into an n*p matrix
a3_c = t(a3_c)

# Centering 
a3_f = scale(a3_c, center = TRUE, scale = FALSE)


write.csv(a3_f,'elliptical_breast_cancer.csv', row.names = TRUE)

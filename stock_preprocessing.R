file = getwd()
stock = read.csv(paste0(file,"/all_stocks_5yr.csv"), header = TRUE)

stock$time <-  as.POSIXct(stock$Date,format = "%Y-%m-%d")
stock$time <- format(stock$time, format = "%Y/%m/%d")
stock$time <- as.Date(stock$time)
startdate <- as.Date(format("2012-6-28", format = "%Y/%m/%d"))
stock <- stock[which(stock$time>startdate),-1]
stock$year<-format(stock$time, format="%Y")
stock$month<-format(stock$time, format="%m")
stock$day<-format(stock$time, format="%d")
enddate <- as.Date(format("2022-7-1", format = "%Y/%m/%d"))


## ranks of prices at the beginning of July 2012
o <- as.Date(format("2012-7-2", format = "%Y/%m/%d"))
o1 <- stock[which(stock$time==o),]
o2 <- o1[order(-o1$Close),]
o3 <- o2$Name
stock <- stock[which(stock$Date<enddate),]
temp1<-stock[order(stock$Name,stock$time),]

## ignore stocks having missing values
temp2<-aggregate(temp1$time,by=list(temp1$Name,temp1$year,
                                    temp1$month),FUN=max)
I <- unique(temp2$x)
which(I == "2019-02-13")
which(I == "2018-04-03")
which(I == "2018-10-10")
which(I == "2018-12-21")
which(I == "2022-06-07")
I = I[-c(18,38,102,123,64)]
stock <- stock[(stock$time %in% I),]

## labeled the stocks according to the ranks of their prices at the beginning of July 2012
data2={}
for (i in 1:485) {
  data2 <- rbind(data2,stock[stock$Name==o3[i],])
}

## ignore stocks that do not have complete histories
y = as.data.frame(table(stock$Name))
y2 = y[order(y$Freq),]
l = y[which(y$Freq==121),]
l1 = l$Var1
stock_d = data2[which(data2$Name %in% l1),]
stock_d1 = matrix(stock_d$Close,121,480)

o4 = o3[o3 %in% l1]
rowname = stock_d[which(stock_d$Name == 'NVR'),'Date']

## calculate stock market return 
data1 = {}
for (i in 1:120) {
  data1 = rbind(data1,log(stock_d1[i+1,]/stock_d1[i,]))
}

data1 = matrix(data1, ncol = 480)
rownames(data1) = rowname[2:121]
colnames(data1) = o4
## save the dataset
write.table(data1,file=paste0(file,"/stock_data_elliptical.csv"),row.names = TRUE,col.names = TRUE, sep=",")

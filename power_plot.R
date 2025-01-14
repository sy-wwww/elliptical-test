  file = getwd()
  
  
  source(paste0(file,"/funcs.r"))
  source(paste0(file,"/funcs2.r"))
  
  
  library(ggplot2)
  library(scales)
  
  m = 400
  rp = 0.5
  d = m*rp
  
  
  
  setting_sigma = "Sig4" 
  setting_distribution = 'dist1'
  
  
  
  
  k_l=c(0,0.05,0.1,0.15,0.2)
  k_l=c(0,0.2)
  
  l1 = length(k_l)
  
  
  
  
  
  r_norm = matrix(0,l1,2)
  r_ellp = rep(0,l1)
  
  f = 5
  k = 100
  
  for (i in 1:l1) {
    normality = matrix(0,f*k,2)
    elliptical = matrix(0,f*k,3)
    
    k_p = k_l[i]
    for (num in 1:f) {
      normality[(k*(num-1)+1):(k*num),] = as.matrix(read.csv(paste0(file,paste('/norm',k_p,m,d,setting_distribution,setting_sigma,num,".csv", sep="_")),header = TRUE))
      elliptical[(k*(num-1)+1):(k*num),] = as.matrix(read.csv(paste0(file,paste("/elliptical",k_p,m,d,setting_distribution,setting_sigma,num,".csv", sep="_")),header = TRUE))
  }
    aY = getpow(normality)
    r_norm[i,] = aY
    
    
    r_ellp[i] = mean(elliptical[,3]< 0.05)
  }
  
  
  
  r_norm = r_norm[,1]
  
  
  p1 = ggplot()+
    geom_line(data = NULL,aes(x = k_l,y = r_norm),color = "red",size=5)+
    geom_point(data = NULL,aes(x =  k_l,y = r_norm),color = "red",size=6)+
    geom_line(data = NULL,aes(x =  k_l,y = r_ellp),color = "blue",size=5)+
    geom_point(data = NULL,aes(x =  k_l,y = r_ellp),color = "blue",size=6)+
    geom_line(data = NULL,aes(x =  k_l,y = rep(0.05,length(k_l))),color = "black",size=3,linetype=2)+
    xlab("")+ylab("")
  
  
  p1 + theme_bw() + theme(
    strip.text.x = element_blank(),
    strip.background = element_rect(colour="white", fill="white"))+
    theme(legend.title= element_blank(),legend.text = element_text(size=28))+
    theme(axis.text=element_text(size=33),
          axis.title=element_text(size=30))+
    theme(legend.position=c(.1,.80),legend.text = element_text(size=0))+
    scale_y_continuous(breaks = c(0.15, 0.4, 0.6, 0.8,1)) +
    scale_x_continuous(breaks = k_l) + theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())+
    geom_text(show.legend = FALSE)+
    theme(legend.key.size = unit(4, 'lines'))+
    theme(legend.title.align = 0.5)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(legend.position = "none")
  

file = getwd()
source(paste0(file,"/tests.R"))

data = read.csv(paste0(file,"/stock_data_elliptical.csv"),header = FALSE)
data = as.matrix(data)


p_l = 100
p_u = 480
p_values = seq(p_l,p_u,10)

re_elliptical = matrix(0,length(p_values),3)
re_normality = matrix(0,length(p_values),1)


for (i in 1:length(p_values)) {
  b1 = data[1:120,1:p_values[i]]
  re_normality[i,] = normality_test(b1)
  re_elliptical[i,] = elliptical_test(b1)
}


library(ggplot2)
p1 = ggplot()+
  geom_line(data = NULL,aes(x = p_values,y = re_elliptical[,3],colour = "elliptical test"),size=3)+
  geom_point(data = NULL,aes(x = p_values,y = re_elliptical[,3],colour = "elliptical test"),size=4)+
  geom_line(data = NULL,aes(x = p_values,y = re_normality,colour = "normality test"),size=3)+
  geom_point(data = NULL,aes(x = p_values,y = re_normality,colour = "normality test"),size=4) +
  scale_color_manual(values = c("blue", "red"))


p1 + theme_bw() + theme(
  strip.text.x = element_blank(),
  strip.background = element_rect(colour="white", fill="white"))+
  # theme(legend.box.background=element_rect(),legend.box.margin=margin(9,9,9,9)) +
  theme(legend.title= element_blank(),legend.text = element_text(size=28))+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24))+
  theme(legend.position=c(.55,.50),legend.text = element_text(size=0))+
  scale_x_continuous(breaks = seq(100,480,50)) + theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())+
  geom_text(show.legend = FALSE)+
  theme(legend.key.size = unit(4, 'lines'))+
  theme(legend.title.align = 0.5)

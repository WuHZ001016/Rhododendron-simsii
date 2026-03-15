library(ggplot2)
da <- read.csv("C:\\Users\\25770\\Desktop\\A12-1 MPEA.csv")
da <- data.frame(TIME = da$Prompt.T.*1000,FL = da$Prompt.Fl.)
O <- da[da$TIME == 0.02,]
J <- da[da$TIME == 2,]
I <- da[da$TIME == 30,]
P <- da[max(da$FL),]

p <- ggplot(da,mapping = aes(x = TIME,y = FL)) +
  geom_point(colour = "darkgrey") +
  scale_x_continuous(limits = c(0.01,max(da$TIME)+1000),breaks = c(O$TIME,J$TIME,I$TIME,P$TIME)) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0,max(da$FL)+1000),breaks = c(O$FL,J$FL,I$FL,P$FL)) +
  geom_vline(xintercept = c(O$TIME,J$TIME,I$TIME,P$TIME),linetype = "dashed",colour = "grey") +
  stat_smooth() +
  labs(x = "Time(s)",y = "fluorescence",title = "A12-1 OJIP") +
  annotate(geom = "text",label = c("O","J","I","P"),x = c(O$TIME,J$TIME,I$TIME,P$TIME),y = c(O$FL,J$FL,I$FL,P$FL),size = 5,hjust = 1.2,vjust = c(-0.5,-1,-0.8,2)) +
  theme_classic()
print(p)
View(da)


max(da$FL)

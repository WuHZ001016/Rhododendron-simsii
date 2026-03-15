setwd("D:/ZheJiang A&F University/杜鹃文章/第一次收获/licor-6800数据/第4次")

library(readxl)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(dplyr)
library(multcompView)### 与下面的 agricolae 是两种字母显著性表示的方式
library(agricolae)### 见上文
library(tibble)
library(showtext)
library(plantecowrap)
library(plantecophys)
library(photosynthesis)
library(purrr)
library(tibble)
library(showtext)
library(car) ## 方差齐性检验 leveneTest
library(dunn.test)
library(patchwork)

#################################### 数据 #####
da.physiology <- read_excel("all.xlsx","physiology")
factor(da.physiology$co2)
factor(da.physiology$light)

da.physiology <- subset.data.frame(da.physiology, co2=='400ppm')

################################ 调用字体库 #####
font_add("Times", regular = "times.ttf", bold = "timesbd.ttf", italic = "timesi.ttf", bolditalic = "timesbi.ttf")
font_add("wryh", regular="msyh.ttc", bold="msyhbd.ttc")
font_add(family = "arial", regular = "arial.ttf", bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")
showtext_auto()


####################################### Ci ####################################
aov.Ci <- aov(Ci ~ light, data = da.physiology)
aov.Ci
summary(aov.Ci)

tukey <- aov(Ci ~ treatment, data = da.physiology) |> TukeyHSD()
tukey.Ci <- as.data.frame(tukey$treatment)
View(tukey.Ci)
write.csv(tukey.Ci,"Ci.csv")

# 生成字母标记 #
zimu <- HSD.test(aov.Ci,'light', alpha = 0.05)$groups
## 生成画图数据 ##
plot_Ci <- da.physiology |> 
  group_by(light) |> 
  summarise(
    mean = mean(Ci, na.rm = TRUE),
    n    = sum(Ci),
    sd   = sd(Ci, na.rm = TRUE),
    se   = sd / sqrt(n)
  ) |> 
left_join(zimu|> rownames_to_column("light"), by = "light")## 将字母加入数据画图

ggplot(plot_Ci, aes(x = light, y = mean)) +
  geom_col(fill = "grey", width = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 15, label = groups), size = 5) +
  labs(x = "treatment", y = "Ci") +
  theme_bw() +
  theme(panel.grid = element_blank())


################################ gs ###########################################
aov.gsw <- aov(gsw ~ treatment, data = da.physiology)
aov.gsw
summary(aov.gsw)
tukey <- aov(gsw ~ treatment, data = da.physiology) |> TukeyHSD()
tukey.gsw <- as.data.frame(tukey$treatment)

View(tukey.gsw)
write.csv(tukey.gsw,"gsw.csv")
################ 生成字母标记
letters <- multcompLetters4(aov.gsw, tukey)
letters_df <- as.data.frame.list(letters$treatment) # 关键：取出$treatment
letters_df$treatment <- rownames(letters_df)
colnames(letters_df)[1] <- "letters"
################### 生成画图数据
plot_gsw <- da.physiology  |> 
  group_by(treatment)  |> 
  summarise(
    mean = mean(gsw, na.rm = TRUE),
    n    = sum(!is.na(gsw)),
    sd   = sd(gsw, na.rm = TRUE),
    se   = sd / sqrt(n)
  )

plot_gsw <- left_join(plot_gsw, letters_df, by = "treatment")## 将字母加入数据画图

ggplot(plot_gsw, aes(x = treatment, y = mean)) +
  geom_col(fill = "grey40", width = 0.7)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.2, linewidth = 1) +
  geom_text(aes(y = mean + se + 0.003, label = letters), size = 5) +
  labs(x = "treatment", y = "gsw") +
  theme_bw() +
  theme(panel.grid = element_blank())

############################## PhiPS2 ##################################
aov.PhiPS2 <- aov(PhiPS2 ~ treatment, data = da.physiology)
aov.PhiPS2
summary(aov.PhiPS2)
tukey <- aov(PhiPS2 ~ treatment, data = da.physiology) |> TukeyHSD()
tukey.PhiPS2 <- as.data.frame(tukey$treatment)

View(tukey.PhiPS2)
write.csv(tukey.PhiPS2,"PhiPS2.csv")

# 生成字母标记 #
zimu <- HSD.test(aov.Ci,'light', alpha = 0.05)$groups
## 生成画图数据 ##
plot_Ci <- da.physiology |> 
  group_by(light) |> 
  summarise(
    mean = mean(Ci, na.rm = TRUE),
    n    = sum(Ci),
    sd   = sd(Ci, na.rm = TRUE),
    se   = sd / sqrt(n)
  ) |> 
  left_join(zimu|> rownames_to_column("light"), by = "light")## 将字母加入数据画图

ggplot(plot_Ci, aes(x = light, y = mean)) +
  geom_col(fill = "grey", width = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 15, label = groups), size = 5) +
  labs(x = "treatment", y = "Ci") +
  theme_bw() +
  theme(panel.grid = element_blank())

################################### ETR ##########################
capture.output(aov(ETR ~co2+light, da.physiology) |> summary(), file = "etr~co2+light.txt")

aov.ETR <- aov(ETR ~ treatment, data = da.physiology)
aov.ETR
summary(aov.ETR)
tukey <- aov(ETR ~ treatment, data = da.physiology) |> TukeyHSD()
tukey.ETR <- as.data.frame(tukey$treatment)
View(tukey.ETR)
write.csv(tukey.ETR,"ETR.csv")
# 生成字母标记 #
letters <- multcompLetters4(aov.ETR, tukey)
letters_df <- as.data.frame.list(letters$treatment) # 关键：取出$treatment
letters_df$treatment <- rownames(letters_df)
colnames(letters_df)[1] <- "letters"
## 生成画图数据 ##
plot_ETR <- da.physiology %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(ETR, na.rm = TRUE),
    n    = sum(!is.na(ETR)),
    sd   = sd(ETR, na.rm = TRUE),
    se   = sd / sqrt(n)
  )

plot_ETR <- left_join(plot_ETR, letters_df, by = "treatment")## 将字母加入数据画图

ggplot(plot_ETR, aes(x = treatment, y = mean)) +
  geom_col(fill = "grey40", width = 0.7)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.2, linewidth = 1) +
  geom_text(aes(y = mean + se + 3, label = letters), size = 5) +
  labs(x = "treatment", y = "ETR") +
  theme_bw() +
  theme(panel.grid = element_blank())


########################### 下面只考虑光周期 ##############################
##### WUE 瞬时水分利用效率 #####
wue.aov <- aov(WUE~light, da.physiology)
wue.aov |> summary()
wue.letters <- HSD.test(wue.aov, "light", alpha=0.05)$groups

wue <- da.physiology |> 
  group_by(light) |>
  summarise(
    mean=mean(WUE, na.rm=TRUE),
    n=n(),
    sd=sd(WUE, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(wue.letters |> rownames_to_column("light"), by="light")

ggplot(wue, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 0.3, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)", y = "WUE (μmol CO₂·mmol⁻¹ H₂O)", title="水分利用效率") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        plot.title = element_text(family="wryh", face="bold", size=35, hjust=0.5)
        )
ggsave(filename="水分利用率.png", width=8, height=8, units="cm", dpi=300, device="png")


############################ ETR #####

etr.aov <- aov(ETR~light, da.physiology)
etr.aov |> summary()
etr.letters <- HSD.test(etr.aov, "light", alpha=0.05)$groups

etr <- da.physiology |> 
  group_by(light) |>
  summarise(
    mean=mean(ETR, na.rm=TRUE),
    n=n(),
    sd=sd(ETR, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(etr.letters |> rownames_to_column("light"), by="light")

ggplot(etr, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 5, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)",
       y = expression(bold(paste('ETR (', mu, 'mol ', m^{-2}, ' ', s^{-1}, ')'))),
       title="电子传递效率") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        plot.title = element_text(family="wryh", face="bold", size=35, hjust=0.5)
  )
ggsave(filename="电子传递效率.png", width=8, height=8, units="cm", dpi=300, device="png")


############################ NPQ #####

npq.aov <- aov(NPQ~light, da.physiology)
npq.aov |> summary()
npq.letters <- HSD.test(npq.aov, "light", alpha=0.05)$groups

npq <- da.physiology |> 
  group_by(light) |>
  summarise(
    mean=mean(NPQ, na.rm=TRUE),
    n=n(),
    sd=sd(NPQ, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(npq.letters |> rownames_to_column("light"), by="light")

ggplot(npq, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 0.1, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)",
       y = "NPQ",
       title="非光化学猝灭") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        plot.title = element_text(family="wryh", face="bold", size=35, hjust=0.5))
ggsave(filename="非光化学淬灭.png", width=8, height=8, units="cm", dpi=300, device="png")


################################ PhiPS2 ###################################

phi.aov <- aov(PhiPS2~light, da.physiology)
phi.aov |> summary()
phi.letters <- HSD.test(phi.aov, "light", alpha=0.05)$groups

phi <- da.physiology |> 
  group_by(light) |>
  summarise(
    mean=mean(PhiPS2, na.rm=TRUE),
    n=n(),
    sd=sd(PhiPS2, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(phi.letters |> rownames_to_column("light"), by="light")

ggplot(phi, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 0.01, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)",
       y = "PhiPS2",
       title="光系统2实际量子效率") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        plot.title = element_text(family="wryh", face="bold", size=35, hjust=0.5))
ggsave(filename="光系统2实际量子效率PhiPS2.png", width=8, height=8, units="cm", dpi=300, device="png")

######################### qP 光化学淬灭系数（数值范围0-1）####################

qP.aov <- aov(qP~light, da.physiology)
qP.aov |> summary()
qP.letters <- HSD.test(qP.aov, "light", alpha=0.05)$groups

qP <- da.physiology |> 
  group_by(light) |>
  summarise(
    mean=mean(qP, na.rm=TRUE),
    n=n(),
    sd=sd(qP, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(qP.letters |> rownames_to_column("light"), by="light")

ggplot(qP, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 0.02, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)",
       y = "qP",
       title="光化学淬灭系数（0~1）") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        plot.title = element_text(family="wryh", face="bold", size=35, hjust=0.5))
ggsave(filename="光化学淬灭系数qP.png", width=8, height=8, units="cm", dpi=300, device="png")


######################## PhiCO2 CO2同化的量子效率 ############################

phi.aov <- aov(PhiCO2~light, da.physiology)
phi.aov |> summary()
phi.letters <- HSD.test(phi.aov, "light", alpha=0.05)$groups

phi <- da.physiology |> 
  group_by(light) |>
  summarise(
    mean=mean(PhiCO2, na.rm=TRUE),
    n=n(),
    sd=sd(PhiCO2, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(phi.letters |> rownames_to_column("light"), by="light")

ggplot(phi, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 0.005, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)",
       y = "PhiCO2",
       title="CO2同化的量子效率_光呼吸") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        plot.title = element_text(family="wryh", face="bold", size=35, hjust=0.5))
ggsave(filename="光化学淬灭系数qP.png", width=8, height=8, units="cm", dpi=300, device="png")

############################ Fo'光下最小荧光 ################################

fo.aov <- aov(`Fo'`~light, da.physiology)
fo.aov |> summary()
fo.letters <- HSD.test(fo.aov, "light", alpha=0.05)$groups

fo <- da.physiology |> 
  group_by(light) |>
  summarise(
    mean=mean(`Fo'`, na.rm=TRUE),
    n=n(),
    sd=sd(`Fo'`, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(fo.letters |> rownames_to_column("light"), by="light")

ggplot(fo, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 10, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)",
       y = "Fo'",
       title="光下最小荧光") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        plot.title = element_text(family="wryh", face="bold", size=35, hjust=0.5))
ggsave(filename="光下最小荧光Fo'.png", width=8, height=8, units="cm", dpi=300, device="png")

########### Fv'/Fm' 光适应下有热耗散等存在时开放的PSII反应中心激发能捕获效率（1-Fo'/Fm'） ################################

fvfm.aov <- aov(`Fv'/Fm'`~light, da.physiology)
fvfm.aov |> summary()
fvfm.letters <- HSD.test(fvfm.aov, "light", alpha=0.05)$groups

fvfm <- da.physiology |> 
  group_by(light) |>
  summarise(
    mean=mean(`Fv'/Fm'`, na.rm=TRUE),
    n=n(),
    sd=sd(`Fv'/Fm'`, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(fvfm.letters |> rownames_to_column("light"), by="light")

ggplot(fvfm, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 0.05, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)",
       y = "Fv'/Fm'",
       title="光适应下有热耗散等存在时开放的PSII反应中心激发能捕获效率") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        plot.title = element_text(family="wryh", face="bold", size=25, hjust=0.5))
ggsave(filename="光适应下有热耗散等存在时开放的PSII反应中心激发能捕获效率Fv'_Fm'.png", width=8, height=8, units="cm", dpi=300, device="png")

############################ 光下最大荧光 Fm' ################################

fm.aov <- aov(`Fm'`~light, da.physiology)
fm.aov |> summary()
fm.letters <- HSD.test(fm.aov, "light", alpha=0.05)$groups

fm <- da.physiology |> 
  group_by(light) |>
  summarise(
    mean=mean(`Fm'`, na.rm=TRUE),
    n=n(),
    sd=sd(`Fm'`, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(fm.letters |> rownames_to_column("light"), by="light")

ggplot(fm, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 22, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)",
       y = "Fm'",
       title="光下最大荧光") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        plot.title = element_text(family="wryh", face="bold", size=35, hjust=0.5))
ggsave(filename="光下最大光F'.png", width=8, height=8, units="cm", dpi=300, device="png")


################################ Ls 气孔限制值 ################################

ls.aov <- aov(Ls~light, da.physiology)
summary(ls.aov)
ls.letters <- HSD.test(ls.aov, "light", alpha=0.05)$groups

ls <- da.physiology |> 
  group_by(light) |>
  summarise(
    mean=mean(Ls, na.rm=TRUE),
    n=n(),
    sd=sd(Ls, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(ls.letters |> rownames_to_column("light"), by="light")

ggplot(ls, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 0.03, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)",
       y = "Limiting Value of Stomata",
       title="气孔限制值(Ca-Ci)/Ca") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        plot.title = element_text(family="wryh", face="bold", size=35, hjust=0.5))
ggsave(filename="气孔限制值Ls.png", width=8, height=8, units="cm", dpi=300, device="png")


################################## 光合速率A #####
da.physiology <- read_excel("all.xlsx","Sheet1")

### 正态检验 14h p-value = 0.0003829 ###
tapply(X = da.physiology$`A（qin1000）`, INDEX = da.physiology$light, FUN = "shapiro.test")
shapiro.test(da.physiology$`A（qin1000）`)
### 方差齐性检验 1.808e-07 *** ###
leveneTest(y = da.physiology$`A（qin1000）` ~ da.physiology$light, data = da.physiology)

ggboxplot(data = da.physiology, x = "light", y = "`A（qin1000）`",
          colour = 'black', fill = 'grey',
          add = "point", bxp.errorbar = TRUE, bxp.errorbar.width = 0.2,
          xlab = "光周期", ylab = "A_net") +
  stat_compare_means(paired = " wilcox.test", family = "Times", face = "bold", fize = 25) +
  theme_classic()+
  theme(axis.text = element_text(family = "Times", size = 20),
        axis.title = element_text(family = "Times", face = "bold", size = 25))



### 不能用ANOVA 使用Kruskal-Wallis检验 ###
An_Kruskal <- kruskal(y = da.physiology$A, trt = da.physiology$light, alpha = 0.05, p.adj = "none")

### 事后进行Dunn检验 ###
Dunn_A <- dun


An.aov <- aov(`A（qin1000）` ~ light, da.physiology)
summary(An.aov)
An.letters <- LSD.test(An.aov, "light", alpha=0.05)$groups

An <- da.physiology |> 
  group_by(light) |>
  summarise(
    mean=mean(`A（qin1000）`, na.rm=TRUE),
    n=n(),
    sd=sd(`A（qin1000）`, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(An.letters |> rownames_to_column("light"), by="light")

ggplot(An, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 0.4, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)",
       y = expression(bold(paste(A[net], " (", mu, "mol ", m^{-2}, s^{-1}, ")"))),
       title="净光合速率") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 30, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        plot.title = element_text(family="wryh", face="bold", size=35))

ggsave(filename="净光合速率Anet.png", width=8, height=8, units="cm", dpi=300, device="png")



















####################################


##### 气孔导度模型 #####

da.physiology <- read_excel("all.xlsx","physiology（Q=600）")
factor(da.physiology$co2)
factor(da.physiology$light)



da_light$RHcham <- da_light$RHcham/100
### A12 ###
gs12h <- fit_gs_model(data = subset(da_light, light == "12h"),
      varnames = list(A_net = "A", g_sw = "gsw", VPD = "VPDleaf", C_air = "Ca", RH = "RHcham"), model = c("BallBerry", "Leuning", "Medlyn_partial", "Medlyn_full"), D0 = 3)
gs12h$BallBerry$Parameters

models <- list(
  BallBerry       = gs12h$BallBerry$Parameters,
  Leuning         = gs12h$Leuning$Parameters,
  Medlyn_partial  = gs12h$Medlyn_partial$Parameters,
  Medlyn_full     = gs12h$Medlyn_full$Parameters
) |> bind_rows(.id = "model")
models
write.csv(models,"A12_gs_model(qin600).csv")
### A14 ###
gs14h <- fit_gs_model(data = subset(da_light, light == "14h"),
                      varnames = list(A_net = "A", g_sw = "gsw", VPD = "VPDleaf", C_air = "Ca", RH = "RHcham"), model = c("BallBerry", "Leuning", "Medlyn_partial", "Medlyn_full"), D0 = 3)
gs14h$BallBerry$Parameters
gs14h$Leuning$Parameters
gs14h$Medlyn_partial$Parameters
gs14h$Medlyn_full$Parameters
models <- list(
  BallBerry       = gs14h$BallBerry$Parameters,
  Leuning         = gs14h$Leuning$Parameters,
  Medlyn_partial  = gs14h$Medlyn_partial$Parameters,
  Medlyn_full     = gs14h$Medlyn_full$Parameters
) |> bind_rows(.id = "model")
models
write.csv(models,"A14_gs_model(qin600).csv")
### A16 ###
gs16h <- fit_gs_model(data = subset(da_light, light == "16h"),
                      varnames = list(A_net = "A", g_sw = "gsw", VPD = "VPDleaf", C_air = "Ca", RH = "RHcham"), model = c("BallBerry", "Leuning", "Medlyn_partial", "Medlyn_full"), D0 = 3)
gs16h$BallBerry$Parameters
gs16h$Leuning$Parameters
gs16h$Medlyn_partial$Parameters
gs16h$Medlyn_full$Parameters
models <- list(
  BallBerry       = gs16h$BallBerry$Parameters,
  Leuning         = gs16h$Leuning$Parameters,
  Medlyn_partial  = gs16h$Medlyn_partial$Parameters,
  Medlyn_full     = gs16h$Medlyn_full$Parameters
) |> bind_rows(.id = "model")
models
write.csv(models,"A16_gs_model(qin600).csv")




# gs12h$BallBerry$Graph +
#   labs(title="A12气孔导度模型（BallBerry）", y="gsw") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         plot.title = element_text(size=40, family="wryh", face="bold"),
#         axis.title = element_text(size=35, family="Times"),
#         axis.text = element_text(size=25, family="Times"))
# ggsave("A12 ballberry model.png", width=5, height=4, dpi=200,device="png")
# 
# gs12h$Leuning$Graph +
#   labs(title="A12气孔导度模型（Leuning）", y="gsw") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         plot.title = element_text(size=30, family="wryh", face="bold"),
#         axis.title = element_text(size=25, family="Times"),
#         axis.text = element_text(size=20, family="Times"))
# ggsave("A12 Leuning model.png", width=5, height=4, dpi=200,device="png")
# 
# gs12h$Medlyn_partial$Graph +
#   labs(title="A12气孔导度模型（Medlyn_partial）", y="gsw") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         plot.title = element_text(size=30, family="wryh", face="bold"),
#         axis.title = element_text(size=25, family="Times"),
#         axis.text = element_text(size=20, family="Times"))
# ggsave("A12 Medlyn_partial model.png", width=5, height=4, dpi=200,device="png")
# 
# gs12h$Medlyn_full$Graph +
#   labs(title="A12气孔导度模型（Medlyn_full）", y="gsw") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         plot.title = element_text(size=30, family="wryh", face="bold"),
#         axis.title = element_text(size=25, family="Times"),
#         axis.text = element_text(size=20, family="Times"))
# ggsave("A12 Medlyn_full model.png", width=5, height=4, dpi=200,device="png")
# 

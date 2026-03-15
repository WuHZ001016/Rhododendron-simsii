library(tidyverse)
library(ggplot2)
library(readxl)
library(purrr)
library(ggpubr)
library(ggsignif)
library(dplyr)
library(multcompView)### 与下面的agricolae是两种字母显著性表示的方式
library(agricolae)### 见上文
library(tibble)
library(showtext)
library(ggsci)
library(patchwork)
library(devEMF)

setwd("D:/ZheJiang A&F University/杜鹃文章/第一次收获/吴华政（20250904）全部数据")

##### 数据 #####
my_data <- read_excel("数据.xlsx","all")
da.light <- subset(my_data, co2=="400ppm")

##################################### 调用字体库 #####
font_add("Times", regular = "times.ttf", bold = "timesbd.ttf", italic = "timesi.ttf", bolditalic = "timesbi.ttf")
font_add("wryh", regular="msyh.ttc", bold="msyhbd.ttc")
font_add(family = "arial", regular = "arial.ttf", bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")
showtext_auto()

############################### 蛋白质含量 ##############################
aov1 <- aov(BCA~light, data=da.light)
summary(aov1) ### p=0.001

TukeyHSD(aov1)$light

tuk.letters <- HSD.test(aov1, "light", alpha=0.05)$groups

bca <- da.light |> 
  group_by(light) |> 
  summarise(
    mean=mean(BCA, na.rm=TRUE),
    n=n(),
    sd=sd(BCA, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(tuk.letters |> rownames_to_column("light"), by="light")

p_protein <- ggplot(bca, aes(x = light, y = mean)) +
  geom_col(mapping = aes(fill = light), width = 0.7)+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean + se + 3, label = groups), size = 8, family = "arial", fontface = "bold") +
  scale_y_continuous(limits=c(0,60), breaks=c(0,15,30,45,60), expand=0) +
  scale_x_discrete(labels = c("12h" = expression(12~h~d^{-1}),
                              "14h" = expression(14~h~d^{-1}),
                              "16h" = expression(16~h~d^{-1}))) +
  labs(x = "Photoperiod",
       y = expression(bold(`Protein Content (mg`~g^{-1}~`FW)`))
       ) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 20)
        ) +
  scale_fill_aaas(alpha = 0.8)
p_protein
ggsave(filename="蛋白质含量.tiff", device="tiff", dpi=300, width=5, height=5, units="cm")
ggsave(filename="蛋白质含量.png", device="png", dpi=300, width=5, height=5, units="cm")

################################ 葡萄糖 glucose ################################
aov2 <- aov(glucose~light, data=da.light)
summary(aov2) ## p=0.001

TukeyHSD(aov2)$light
glu.letters <- HSD.test(aov2, "light", alpha=0.05)$groups

glu <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(glucose, na.rm=TRUE),
    n=n(),
    sd=sd(glucose, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(glu.letters |> rownames_to_column("light"), by="light")

ggplot(glu, aes(x = light, y = mean, fill = light))+
  geom_col(width = 0.7)+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean + se + 0.25, label = groups), size = 8, family = "arial", fontface="bold") +
  scale_x_discrete(labels = c("12h" = expression(12~h~d^{-1}),
                              "14h" = expression(14~h~d^{-1}),
                              "16h" = expression(16~h~d^{-1}))) +
  scale_y_continuous(limits = c(0,6), breaks = c(0,1.5,3,4.5,6), expand = 0) +
  labs(x = "Photoperiod",
       y = expression(bold(Glucose~"(mg"~g^{-1}~"FW)")),
       title="(A)") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        title = element_text(face="bold", size=25)) +
  scale_fill_aaas(alpha = 0.8)

ggsave(filename="葡萄糖.tiff", device = "tiff", width=5, height=5, units="cm", dpi=300)
ggsave(filename="葡萄糖.png", device = "png", width=5, height=5, units="cm", dpi=300)

################################# 果糖 fructose ###################################
aov3 <- aov(fructose~light, data=da.light)
summary(aov3) ## p=0.002

TukeyHSD(aov3)$light
fru.letters <- HSD.test(aov3, "light", alpha=0.05)$groups

fru <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(fructose, na.rm=TRUE),
    n=n(),
    sd=sd(fructose, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(fru.letters |> rownames_to_column("light"), by="light")

ggplot(fru, aes(x=light, y=mean, fill=light))+
  geom_col(width = 0.7)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean + se + 0.1, label = groups), size = 8, family = "arial", fontface = "bold") +
  scale_x_discrete(labels = c("12h" = expression(12~h~d^{-1}),
                              "14h" = expression(14~h~d^{-1}),
                              "16h" = expression(16~h~d^{-1}))) +
  scale_y_continuous(limits = c(0,2.8), breaks = c(0,0.7,1.4,2.1,2.8), expand = 0) +
  labs(x = "Photoperiod",
       y = expression(bold(Fructose~"(mg"~g^{-1}~"FW)")),
       title = "(B)") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial", colour = "black"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 25)) +
  scale_fill_aaas(alpha = 0.8)

ggsave(filename="果糖.tiff", width=5, height=5, units="cm", dpi=300,device = 'tiff')
ggsave(filename="果糖.png", width=5, height=5, units="cm", dpi=300,device = 'png')


################################# 蔗糖 sucrose ##################################
aov4 <- aov(sucrose~light, data=da.light)
summary(aov4) ## p=0.002

TukeyHSD(aov4)$light
suc.letters <- HSD.test(aov4, "light", alpha=0.05)$groups

suc <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(sucrose, na.rm=TRUE),
    n=n(),
    sd=sd(sucrose, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(suc.letters |> rownames_to_column("light"), by="light")

ggplot(suc, aes(x=light, y=mean, fill = light))+
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean + se + 0.3, label = groups), size = 8, family = "arial", fontface = "bold") +
  scale_x_discrete(labels = c("12h" = expression(12~h~d^{-1}),
                              "14h" = expression(14~h~d^{-1}),
                              "16h" = expression(16~h~d^{-1}))) +
  scale_y_continuous(limits = c(0,8), breaks = c(0,2,4,6,8), expand = 0) +
  labs(x = "Photoperiod",
       y = expression(bold(Sucrose~"(mg"~g^{-1}~"FW)")),
       title="(C)") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        title = element_text(face = "bold", size = 25)) +
  scale_fill_aaas(alpha = 0.8)

ggsave(filename="蔗糖.tiff", width=5, height=5, units="cm", dpi=300, device="tiff")
ggsave(filename="蔗糖.png", width=5, height=5, units="cm", dpi=300, device="png")

################################## 麦芽糖 ##################################
aov5 <- aov(maltose~light, data=da.light)
summary(aov5) ## p=0.002

TukeyHSD(aov5)$light
mal.letters <- HSD.test(aov5, "light", alpha=0.05)$groups

mal <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(maltose, na.rm=TRUE),
    n=n(),
    sd=sd(maltose, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(mal.letters |> rownames_to_column("light"), by="light")

ggplot(mal, aes(x=light, y=mean))+
  geom_col(mapping = aes(fill = light), width = 0.7)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean + se + 0.05, label = groups), size = 8, family = "arial", fontface="bold") +
 scale_x_discrete(labels = c("12h" = expression(12~h~d^{-1}),
                             "14h" = expression(14~h~d^{-1}),
                             "16h" = expression(16~h~d^{-1}))) +
 scale_y_continuous(limits = c(0,1.2), breaks = c(0,0.3,0.6,0.9,1.2), expand = 0) +
  labs(x = "Photoperiod",
       y = expression(bold(paste("Maltose (mg ", g^{-1}, " FW)"))),
       title="mai ya tang") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        title = element_text(face = "bold", size = 25)) +
 scale_fill_aaas(alpha = 0.8)

ggsave(filename="麦芽糖.tiff", width=5, height=5, units="cm", dpi=300, device = 'tiff')
ggsave(filename="麦芽糖.png", width=5, height=5, units="cm", dpi=300, device = 'png')

################################## 淀粉 ######################################
aov6 <- aov(starch~light, data=da.light)
summary(aov6) ## p < 0.001

TukeyHSD(aov6)$light
sta.letters <- HSD.test(aov6, "light", alpha=0.05)$groups

sta <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(starch, na.rm=TRUE),
    n=n(),
    sd=sd(starch, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(sta.letters |> rownames_to_column("light"), by="light")

p_starch <- ggplot(sta, aes(x = light, y = mean)) +
  geom_col(mapping = aes(fill = light), width = 0.7, colour = "black", linewidth = 0.1) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean + se + 5, label = groups), size = 8, family = "arial", fontface="bold") +
  scale_x_discrete(labels = c("12h" = "12", "14h" = "14", "16h" = "16")) +
  scale_y_continuous(limits = c(0,60), breaks = c(0,15,30,45,60), expand = 0) +
  labs(x = expression(bold(Photoperiod~(h~d^{-1}))),
       y = expression(bold("Starch" ~ "(mg" ~ g^{-1} ~ "FW)"))
       ) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20)
        ) +
 scale_fill_brewer(palette = "YlOrRd")
p_starch
ggsave(filename = "淀粉.tiff", device = "tiff", width=5, height=5, units="cm", dpi=300)
ggsave(filename = "淀粉.png", device = "png", width=5, height=5, units="cm", dpi=300)

############################# 可溶性糖 soluble sugar ###############################
aov_kerongtang <- aov(`soluble sugar` ~ light, data=da.light)
summary(aov_kerongtang) ## p < 0.001

TukeyHSD(aov_kerongtang)$light
ss.letters <- HSD.test(aov_kerongtang, "light", alpha=0.05)$groups

ss <- da.light |> 
 group_by(light) |>
 summarise(
  mean=mean(`soluble sugar`, na.rm=TRUE),
  n=n(),
  sd=sd(`soluble sugar`, na.rm=TRUE),
  se=sd/sqrt(n)
 ) |> left_join(ss.letters |> rownames_to_column("light"), by="light")

p_ss <- ggplot(ss, aes(x=light, y=mean)) +
 geom_col(mapping = aes(fill = light), width = 0.7, colour = "black", linewidth = 0.1) +
 geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
               position = position_dodge(), width = 0.2, linewidth = 0.4) +
 geom_text(aes(y = mean + se + 1.2, label = groups), size = 8, family = "arial", fontface="bold") +
 scale_x_discrete(labels = c("12h" = "12", "14h" = "14", "16h" = "16")) +
 scale_y_continuous(limits = c(0,20), breaks = c(0,5,10,15,20), expand = 0) +
 labs(x = expression(bold(Photoperiod~(h~d^{-1}))),
      y = expression(bold(Soluble~Sugar~(mg~g^{-1}~FW)))
      ) +
 theme_classic() +
 theme(legend.position = "none",
       text = element_text(family = "arial"),
       axis.title = element_text(size = 25, face = "bold"),
       axis.text = element_text(size = 20)
       ) +
 scale_fill_brewer(palette = "YlOrRd")
p_ss
ggsave("可溶性糖.tiff", device = "tiff", width=5, height=5, units="cm", dpi=300)
ggsave("可溶性糖.png", device = "png", width=5, height=5, units="cm", dpi=300)


############################ 非结构性碳水化合物 NSC #########################
aov_nsc <- aov(NSC ~ light, data=da.light)
summary(aov_nsc) ## p < 0.001

TukeyHSD(aov_kerongtang)$light
nsc.letters <- HSD.test(aov_nsc, "light", alpha=0.05)$groups

nsc <- da.light |> 
 group_by(light) |>
 summarise(
  mean=mean(NSC, na.rm=TRUE),
  n=n(),
  sd=sd(NSC, na.rm=TRUE),
  se=sd/sqrt(n)
 ) |> left_join(nsc.letters |> rownames_to_column("light"), by="light")

p_nsc <- ggplot(nsc, aes(x=light, y=mean)) +
 geom_col(mapping = aes(fill=light), width = 0.7, colour = "black", linewidth = 0.1) +
 geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
 geom_text(aes(y=mean+se+4, label=groups), size=8, family="arial", fontface="bold") +
 scale_x_discrete(labels = c("12h" = "12", "14h" = "14", "16h" = "16")) +
 scale_y_continuous(limits = c(0,60), breaks = c(0,15,30,45,60), expand = 0) +
 labs(x = expression(bold(Photoperiod~(h~d^{-1}))),
      y = expression(bold(NSC~(mg~g^{-1}~FW)))
      ) +
 theme_classic() +
 theme(legend.position = "none",
       text = element_text(family = "arial"),
       axis.title = element_text(size = 25, face = "bold"),
       axis.text = element_text(size = 20)
       ) +
 scale_fill_brewer(palette = "YlOrRd")
p_nsc
ggsave("非结构性碳水化合物.tiff", device = "tiff", width=5, height=5, units="cm", dpi=300)
ggsave("非结构性碳水化合物.png", device = "png", width=5, height=5, units="cm", dpi=300)


######################## 腺苷二磷酸葡萄糖焦磷酸化酶 AGPase ############################
aov7 <- aov(AGPase~light, data=da.light)
summary(aov7) ## p < 0.001
TukeyHSD(aov7)$light
agp.letters <- HSD.test(aov7, "light", alpha=0.05)$groups

agp <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(AGPase, na.rm=TRUE),
    n=n(),
    sd=sd(AGPase, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(agp.letters |> rownames_to_column("light"), by="light")

ggplot(agp, aes(x=light, y=mean))+
  geom_col(mapping = aes(fill = light), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean + se + 45, label = groups), size = 8, family = "arial", fontface="bold") +
  scale_x_discrete(labels = c("12h" = expression(12~h~d^{-1}),
                              "14h" = expression(14~h~d^{-1}),
                              "16h" = expression(16~h~d^{-1}))) +
  scale_y_continuous(limits = c(0,1200), breaks = c(0,300,600,900,1200), expand = 0) +
  labs(x = "Photoperiod",
       y = expression(bold("AGPase (nmol" ~ min^{-1} ~ g^{-1} ~ "FW)")),
       title = "(A)") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        title = element_text(face = "bold", size = 25)) +
 scale_fill_aaas(alpha = 0.8)

ggsave(filename="AGPase.tiff", device='tiff',width=5, height=5, units="cm", dpi=300)
ggsave(filename="AGPase.png", device='png', width=5, height=5, units="cm", dpi=300)

################################ 可溶性淀粉合成酶 SSS #########################
aov8 <- aov(SSS~light, data=da.light)
summary(aov8)

TukeyHSD(aov8)$light
SSS.letters <- HSD.test(aov8, "light", alpha=0.05)$groups

SSS <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(SSS, na.rm=TRUE),
    n=n(),
    sd=sd(SSS, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(SSS.letters |> rownames_to_column("light"), by="light")

ggplot(SSS, aes(x=light, y=mean)) +
  geom_col(mapping = aes(fill = light), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean+se+2.5, label=groups), size=8, family="arial", fontface="bold") +
  scale_x_discrete(labels = c("12h" = expression(12~h~d^{-1}),
                             "14h" = expression(14~h~d^{-1}),
                             "16h" = expression(16~h~d^{-1}))) +
  scale_y_continuous(limits = c(0,60), breaks = c(0,15,30,45,60), expand = 0) +
  labs(x = "Photoperiod",
       y = expression(bold(paste("SSS (nmol ", min^{-1}, g^{-1}, " FW)"))),
       title = "(C)") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        title = element_text(face = "bold", size = 25)) +
  scale_fill_aaas(alpha = 0.8)

ggsave(filename="可溶性淀粉合成酶.tiff",device='tiff',width=5,height=5,units="cm",dpi=300)
ggsave(filename="可溶性淀粉合成酶.png",device='png',width=5,height=5,units="cm",dpi=300)

##################### 海藻糖-6-磷酸合成酶 TPS #########################
aov9 <- aov(TPS~light, data=da.light)
summary(aov9) ## p < 0.001

TukeyHSD(aov9)$light
tps.letters <- HSD.test(aov9, "light", alpha=0.05)$groups

tps <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(TPS, na.rm=TRUE),
    n=n(),
    sd=sd(TPS, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(tps.letters |> rownames_to_column("light"), by="light")

ggplot(tps, aes(x=light, y=mean)) +
  geom_col(mapping = aes(fill = light), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean + se + 2.5, label = groups), size=8, family="arial", fontface="bold") +
  scale_x_discrete(labels = c("12h" = expression(12~h~d^{-1}),
                              "14h" = expression(14~h~d^{-1}),
                              "16h" = expression(16~h~d^{-1}))) +
  scale_y_continuous(limits = c(0,60), breaks = c(0,15,30,45,60), expand = 0) +
  labs(x = "Photoperiod",
       y = expression(bold("TPS"~"(noml"~min^{-1}~g^{-1}~"FW)")),
       title="(A)") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        title = element_text(face = "bold", size = 25)) +
  scale_fill_aaas(alpha = 0.8)

ggsave(filename="海藻糖-6-磷酸合成酶.tiff",device='tiff',width=5,height=5,units="cm",dpi=300)
ggsave(filename="海藻糖-6-磷酸合成酶.png",device='png',width=5,height=5,units="cm",dpi=300)

##################### 核酮糖-1,5-二磷酸羧化酶 Rubisco #####################
### 正态检验 ###
tapply(X = da.light$`Rubisco_area`, INDEX = da.light$light, FUN = "shapiro.test")
aov10 <- aov(Rubisco_area~light, data=da.light)
summary(aov10) ## p>0.05

TukeyHSD(aov10)$light
rub.letters <- HSD.test(aov10, "light", alpha=0.05)$groups

rub <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(Rubisco_area, na.rm=TRUE),
    n=n(),
    sd=sd(Rubisco_area, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(rub.letters |> rownames_to_column("light"), by="light")

p_rubisco <- ggplot(rub, aes(x=light, y=mean))+
  geom_col(mapping = aes(fill = light), width = 0.7, colour = "black", linewidth = 0.1) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean + se + 5000, label = groups), size = 10, family = "arial", fontface="bold") +
  scale_x_discrete(labels = c("12h" = "12", "14h" = "14", "16h" = "16")) +
  scale_y_continuous(limits = c(0,10^5),
                     breaks = c(0,2.5*10^4,5*10^4,7.5*10^4,10^5),
                     labels = c(0,2.5,5,7.5,10),
                     expand = 0) +
  labs(x = expression(bold(Photoperiod~(h~d^{-1}))),
       y = expression(bold(paste(`Rubisco (nmol `,min^{-1},` `,m^{-2},`, ×`,10^5,`)`)))
       ) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        ) +
  scale_fill_brewer(palette = "YlOrRd")
p_rubisco
ggsave(filename="核酮糖-1,5-二磷酸羧化酶Rubisco area.tiff",device='tiff',width=5,height=5,units="cm",dpi=300)
ggsave(filename="核酮糖-1,5-二磷酸羧化酶Rubisco area.png",device='png',width=5,height=5,units="cm",dpi=300)

########################## beta粉酶 β-amylase #################################
aov11 <- aov(beta_amylase ~ light, data=da.light)
summary(aov11) ## p < 0.001
TukeyHSD(aov11)$light
beta.letters <- HSD.test(aov11, "light", alpha=0.05)$groups
beta <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(beta_amylase, na.rm=TRUE),
    n=n(),
    sd=sd(beta_amylase, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(beta.letters |> rownames_to_column("light"), by="light")

ggplot(beta, aes(x=light, y=mean)) +
  geom_col(mapping = aes(fill = light), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width=0.2, linewidth=0.4) +
  geom_text(aes(y = mean+se+0.7, label = groups), size=8, family="arial", fontface="bold") +
  scale_x_discrete(labels = c("12h" = expression(12~h~d^{-1}),
                             "14h" = expression(14~h~d^{-1}),
                             "16h" = expression(16~h~d^{-1}))) +
  scale_y_continuous(limits = c(0,16), breaks = c(0,4,8,12,16), expand = 0) +
  labs(x = "Photoperiod",
       y = expression(bold("β-amylase (mg"~min^{-1}~g^{-1}~"FW)")),
       title = "(B)") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        title = element_text(face="bold", size=25)) +
 scale_fill_aaas(alpha = 0.8)
ggsave(filename="β淀粉酶.tiff", width=5, height=5, units="cm", dpi=300, device="tiff")
ggsave(filename="β淀粉酶.png", width=5, height=5, units="cm", dpi=300, device="png")

###################### 叶绿素a chlorophyll a #####
aov12 <- aov(chlorophylla ~ light, data=da.light)
summary(aov12) ## p < 0.001
TukeyHSD(aov12)$light
chla.letters <- HSD.test(aov12, "light", alpha=0.05)$groups
chla <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(chlorophylla, na.rm=TRUE),
    n=n(),
    sd=sd(chlorophylla, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(chla.letters |> rownames_to_column("light"), by="light")

ggplot(chla, aes(x=light, y=mean)) +
  geom_col(mapping = aes(fill = light), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean + se + 0.03, label = groups), size = 8, family = "arial", fontface="bold") +
  scale_x_discrete(labels = c("12h" = expression(12~h~d^{-1}),
                              "14h" = expression(14~h~d^{-1}),
                              "16h" = expression(16~h~d^{-1}))) +
  scale_y_continuous(limits = c(0,0.8), breaks = c(0,0.2,0.4,0.6,0.8), expand = 0) +
  labs(x = "Photoperiod",
       y = expression(bold("Chlorophyll a"~"(mg"~g^{-1}~"FW)")),
       title = "Chl a") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        title = element_text(face="bold", size=25)) +
 scale_fill_aaas(alpha = 0.8)
ggsave(filename="叶绿素a.tiff",device="tiff",width=5,height=5,units="cm",dpi=300)
ggsave(filename="叶绿素a.png",device="png",width=5,height=5,units="cm",dpi=300)

########################### 叶绿素b chlorophyll b #####
aov13 <- aov(chlorophyllb ~ light, data=da.light)
summary(aov13) ## p < 0.001

TukeyHSD(aov13)$light
chlb.letters <- HSD.test(aov13, "light", alpha=0.05)$groups

chlb <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(chlorophyllb, na.rm=TRUE),
    n=n(),
    sd=sd(chlorophyllb, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(chlb.letters |> rownames_to_column("light"), by="light")

ggplot(chlb, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.1, linewidth = 0.8) +
  geom_text(aes(y = mean + se + 0.02, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)", y = "Chlorophyll b (mg/g fresh weight)", title="叶绿素b") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        title = element_text(family="wryh", face="bold", size=30))
ggsave(filename="叶绿素b.png", width=10, height=10, units="cm", dpi=300, device="png")

########################## 总叶绿素 chlorophyll #####
aov14 <- aov(chlorophyll ~ light, data=da.light)
summary(aov14) ## p < 0.001
TukeyHSD(aov14)$light
chl.letters <- HSD.test(aov14, "light", alpha=0.05)$groups
chl <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(chlorophyll, na.rm=TRUE),
    n=n(),
    sd=sd(chlorophyll, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(chl.letters |> rownames_to_column("light"), by="light")

p_chl <- ggplot(chl, aes(x=light, y=mean))+
  geom_col(aes(fill = light), width = 0.7, colour = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width=0.2, linewidth = 0.4) +
  geom_text(aes(y=mean+se+0.07, label=groups), size=10, family="arial", fontface="bold") +
  scale_x_discrete(labels = c("12h" = "12", "14h" = "14", "16h" = "16")) +
  scale_y_continuous(limits = c(0,1.2), breaks = c(0,0.3,0.6,0.9,1.2), expand = 0) +
  labs(x = expression(bold(Photoperiod~(h~d^{-1}))),
       y = expression(bold(`Chlorophyll (mg`~g^{-1}~"FW)"))
       ) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        ) +
  scale_fill_brewer(palette = "YlOrRd")
p_chl
ggsave(filename="叶绿素.tiff", width=5, height=5, units="cm", dpi=300, device="tiff")
ggsave(filename="叶绿素.png", width=5, height=5, units="cm", dpi=300, device="png")

########################### 蔗糖磷酸合成酶 SPS ################################
aov15 <- aov(SPS ~ light, data=da.light)
summary(aov15) ## p < 0.001
TukeyHSD(aov15)$light
sps.letters <- HSD.test(aov15, "light", alpha=0.05)$groups
sps <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(SPS, na.rm=TRUE),
    n=n(),
    sd=sd(SPS, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(sps.letters |> rownames_to_column("light"), by="light")

ggplot(sps, aes(x=light, y=mean))+
  geom_col(aes(fill = light), width = 0.7)+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean + se + 25, label = groups), size = 8, family = "arial", fontface="bold") +
 scale_x_discrete(labels = c("12h" = expression(12~h~d^{-1}),
                             "14h" = expression(14~h~d^{-1}),
                             "16h" = expression(16~h~d^{-1}))) +
  scale_y_continuous(limits = c(0,600), breaks = c(0,150,300,450,600), expand = 0) +
  labs(x = "Photoperiod",
       y = expression(bold(paste("SPS (μg ", min^{-1}, g^{-1}, "FW)"))),
       title= "(E)") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        title = element_text(face="bold", size=25)) +
  scale_fill_aaas(alpha = 0.8)
ggsave(filename="蔗糖磷酸合成酶.tiff", width=5, height=5, units="cm", dpi=300, device="tiff")
ggsave(filename="蔗糖磷酸合成酶.png", width=5, height=5, units="cm", dpi=300, device="png")

############################ 果糖-1,6-二磷酸酶 FBP ##################################
aov16 <- aov(FBP ~ light, data=da.light)
summary(aov16) ## p < 0.001

TukeyHSD(aov16)$light
fbp.letters <- HSD.test(aov16, "light", alpha=0.05)$groups

fbp <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(FBP, na.rm=TRUE),
    n=n(),
    sd=sd(FBP, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(fbp.letters |> rownames_to_column("light"), by="light")

ggplot(fbp, aes(x=light, y=mean))+
  geom_col(mapping = aes(fill = light), width = 0.7)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, linewidth=0.4) +
  geom_text(aes(y=mean+se+2, label=groups), size=8, family = "arial", fontface="bold") +
  scale_x_discrete(labels = c("12h" = expression(12~h~d^{-1}),
                              "14h" = expression(14~h~d^{-1}),
                              "16h" = expression(16~h~d^{-1}))) +
  scale_y_continuous(limits = c(0,40), breaks = c(0,10,20,30,40), expand = 0) +
  labs(x = "Photoperiod",
       y = expression(bold("FBP"~"(nmol"~min^{-1}~g^{-1}~"FW)")),
       title = "(D)") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        title = element_text(face="bold", size=25)) +
  scale_fill_aaas(alpha = 0.8)
ggsave(filename="果糖-1,6-二磷酸酶FBP.tiff", width=5, height=5, units="cm", dpi=300, device="tiff")
ggsave(filename="果糖-1,6-二磷酸酶FBP.png", width=5, height=5, units="cm", dpi=300, device="png")

################################# 碳酸酐酶 CA ####################################
aov17 <- aov(CA ~ light, data=da.light)
summary(aov17) ## p < 0.001

TukeyHSD(aov17)$light
ca.letters <- HSD.test(aov17, "light", alpha=0.05)$groups

ca <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(CA, na.rm=TRUE),
    n=n(),
    sd=sd(CA, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(ca.letters |> rownames_to_column("light"), by="light")

p_ca <- ggplot(ca, aes(x=light, y=mean))+
  geom_col(mapping = aes(fill = light), width = 0.7)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean + se + 0.5, label = groups), size = 8, family = "arial", fontface="bold") +
  scale_x_discrete(labels = c("12h" = expression(12~h~d^{-1}),
                              "14h" = expression(14~h~d^{-1}),
                              "16h" = expression(16~h~d^{-1}))) +
  scale_y_continuous(limits = c(0,12), breaks = c(0,3,6,9,12), expand = 0) +
  labs(x = "Photoperiod",
       y = expression(bold(paste(`CA (nmol `,min^{-1},` `,g^{-1},`)`))),
       ) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 20)
        ) +
  scale_fill_aaas(alpha = 0.8)
p_ca
ggsave(filename="碳酸酐酶.tiff", width=5, height=5, units="cm", dpi=300, device="tiff")
ggsave(filename="碳酸酐酶.png", width=5, height=5, units="cm", dpi=300, device="png")

########################### 肌醇-1-磷酸合成酶 MIPS1 #############################
aov18 <- aov(MIPS1 ~ light, data=da.light)
summary(aov18) ## p < 0.001
TukeyHSD(aov18)$light
mips1.letters <- HSD.test(aov18, "light", alpha=0.05)$groups
mips1 <- da.light |> 
  group_by(light) |>
  summarise(
    mean=mean(MIPS1, na.rm=TRUE),
    n=n(),
    sd=sd(MIPS1, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(mips1.letters |> rownames_to_column("light"), by="light")

p_mips1 <- ggplot(mips1, aes(x=light, y=mean))+
  geom_col(mapping = aes(fill = light), width = 0.7, colour = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean + se + 13, label = groups), size = 8, family = "arial", fontface="bold") +
  scale_x_discrete(labels = c("12h" = "12",
                              "14h" = "14",
                              "16h" = "16")) +
  scale_y_continuous(limits = c(0,200), breaks = c(0,50,100,150,200), expand = 0) +
  labs(x = expression(bold(Photoperiod~(h~d^{-1}))),
       y = expression(bold(paste("MIPS1 (U ",L^{-1},")")))
       ) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20)
        ) +
  scale_fill_brewer(palette = "YlOrRd")
p_mips1
ggsave(filename="肌醇-1-磷酸合成酶.tiff", width=5, height=5, units="cm", dpi=300, device="tiff")
ggsave(filename="肌醇-1-磷酸合成酶.png", width=5, height=5, units="cm", dpi=300, device="png")

############################ C、N #########################################
### 碳氮比 ###
tapply(X = da.light$`C_N`, INDEX = da.light$light, FUN = "shapiro.test") ###正态检验
aov_cn <- aov(`C_N` ~ light, data = da.light) ###ANOVA
summary(aov_cn)
cn_letters <- HSD.test(aov_cn, "light", alpha = 0.05)$groups ###多重比较

cn <- da.light |> 
 group_by(light) |>
 summarise(
  mean=mean(`C_N`, na.rm=TRUE),
  n=n(),
  sd=sd(`C_N`, na.rm=TRUE),
  se=sd/sqrt(n)
 ) |> left_join(cn_letters |> rownames_to_column("light"), by="light")

p_cn <- ggplot(cn, aes(x = light, y = mean))+
 geom_col(mapping = aes(fill = light), width = 0.7) +
 geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
 geom_text(aes(y = mean + se + 1.5, label = groups), size = 8, family = "arial", fontface = "bold") +
 scale_x_discrete(labels = c("12h" = expression(paste('12 h ',d^{-1})),
                             "14h" = expression(paste('14 h ',d^{-1})),
                             "16h" = expression(paste('16 h ',d^{-1})))) +
 scale_y_continuous(limits = c(0,40), breaks = c(0,10,20,30,40), expand = 0) +
 labs(x = "Photoperiod", y = "Carbon to Nitrogen Ratio") +
 theme_classic() +
 theme(legend.position = "none",
       text = element_text(family = "arial"),
       axis.title = element_text(size = 20, face = "bold"),
       axis.text = element_text(size = 20)) +
 scale_fill_aaas(alpha = 0.8)
p_cn
ggsave(filename="碳氮比.tiff", width=5, height=5, units="cm", dpi=300, device="tiff")
ggsave(filename="碳氮比.png", width=5, height=5, units="cm", dpi=300, device="png")

### 氮 ###
tapply(X = da.light$N, INDEX = da.light$light, FUN = "shapiro.test") 
aov_n <- aov(N ~ light, data = da.light)
summary(aov_n)
n_letters <- HSD.test(aov_n, "light", alpha = 0.05)$groups

n <- da.light |> 
 group_by(light) |>
 summarise(
  mean=mean(N, na.rm=TRUE),
  n=n(),
  sd=sd(N, na.rm=TRUE),
  se=sd/sqrt(n)
 ) |> left_join(n_letters |> rownames_to_column("light"), by="light")
### FIGURE ###
p_N <- ggplot(n, aes(x = light, y = mean))+
 geom_col(mapping = aes(fill = light), width = 0.7, colour = "black", linewidth = 0.1) +
 geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
 geom_text(aes(y=mean+se+0.2, label=groups), size=10, family = "arial", fontface="bold") +
 scale_x_discrete(labels = c("12h" = "12", "14h" = "14", "16h" = "16")) +
 scale_y_continuous(limits = c(0,4), breaks = c(0,1,2,3,4), expand = 0) +
 labs(x = expression(bold(Photoperiod~(h~d^{-1}))),
      y = "Nitrogen Content (%)"
      ) +
 theme_classic() +
 theme(legend.position = "none",
       text = element_text(family = "arial"),
       axis.title = element_text(size = 25, face = "bold"),
       axis.text = element_text(size = 20)
       ) +
 scale_fill_brewer(palette = "YlOrRd")
p_N
ggsave(filename="氮.tiff", width=5, height=5, units="cm", dpi=300, device="tiff")
ggsave(filename="氮.png", width=5, height=5, units="cm", dpi=300, device="png")

### 碳 ###
tapply(X = da.light$C, INDEX = da.light$light, FUN = "shapiro.test") 
aov_C <- aov(C ~ light, data = da.light)
summary(aov_C)
C_letters <- HSD.test(aov_C, "light", alpha = 0.05)$groups
C <- da.light |> 
 group_by(light) |>
 summarise(
  mean=mean(C, na.rm=TRUE),
  n=n(),
  sd=sd(C, na.rm=TRUE),
  se=sd/sqrt(n)
 ) |> left_join(C_letters |> rownames_to_column("light"), by="light")
### FIGURE ###
p_C <- ggplot(C, aes(x = light, y = mean))+
 geom_col(mapping = aes(fill = light), width = 0.7) +
 geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
 geom_text(aes(y = mean + se + 3, label = groups), size = 8, family = "arial", fontface = "bold") +
 scale_x_discrete(labels = c("12h" = expression(paste('12 h ',d^{-1})),
                             "14h" = expression(paste('14 h ',d^{-1})),
                             "16h" = expression(paste('16 h ',d^{-1})))) +
 scale_y_continuous(limits = c(0,80), breaks = c(0,20,40,60,80), expand = 0) +
 labs(x = "Photoperiod",
      y = "Total Carbon (%)"
      ) +
 theme_classic() +
 theme(legend.position = "none",
       text = element_text(family = "arial"),
       axis.title = element_text(size = 20, face = "bold"),
       axis.text = element_text(size = 20)
       ) +
 scale_fill_aaas(alpha = 0.8)
p_C
ggsave(filename="碳.tiff", width=5, height=5, units="cm", dpi=300, device="tiff")
ggsave(filename="碳.png", width=5, height=5, units="cm", dpi=300, device="png")


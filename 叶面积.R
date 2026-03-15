setwd("D:/ZheJiang A&F University/杜鹃文章/第一次收获/叶面积")

### 加载依赖包 ###
library(tidyverse)
library(ggplot2)
library(readxl)
library(plantecophys)
library(plantecowrap)
library(photosynthesis)
library(FitAQ)
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

#################### 调用字体库 ###################
font_add("Times", regular = "times.ttf", bold = "timesbd.ttf", italic = "timesi.ttf", bolditalic = "timesbi.ttf")
font_add("wryh", regular="msyh.ttc", bold="msyhbd.ttc")
font_add(family = "arial", regular = "arial.ttf", bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")
showtext_auto()


############################# data ####################################

da <- read_excel("叶面积.xlsx", "all")
factor(da$co2)
factor(da$light)
da_400ppm <- subset.data.frame(da, co2=="400ppm")

# da.400.12 <- subset(da, treatment == "400_12")
# da.400.14 <- subset(da, treatment == "400_14")
# da.400.16 <- subset(da, treatment == "400_16")
# da.800.12 <- subset(da, treatment == "800_12")
# da.800.14 <- subset(da, treatment == "800_12")



################################ 叶面积 ####################################
aov.leafarea <- aov(leafarea ~ light, data = da_400ppm)
aov.leafarea
summary(aov.leafarea)

# tukey <- aov(leafarea ~ treatment, data = da) |> TukeyHSD()
# tukey.leafarea <- as.data.frame(tukey$treatment)
# 
# view(tukey.leafarea)
# write.csv(tukey.leafarea,"leaf area.csv")
# # 生成字母标记 #
# letters <- multcompLetters4(aov.leafarea, tukey)
# letters_df <- as.data.frame.list(letters$treatment) # 关键：取出$treatment
# letters_df$treatment <- rownames(letters_df)
# colnames(letters_df)[1] <- "letters"

zimu <- HSD.test(aov.leafarea, "light", alpha=0.05)
zimu <- zimu$groups

area <- da_400ppm |> 
  group_by(light)  |> 
  summarise(
    mean = mean(leafarea, na.rm = TRUE),
    n    = n(),
    sd   = sd(leafarea, na.rm = TRUE),
    se   = sd / sqrt(n)
  ) |> left_join(zimu |> rownames_to_column("light"), by = "light")

p_leafarea <- ggplot(area, aes(x = light, y = mean)) +
  geom_col(mapping = aes(fill = light), width = 0.7)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean + se + 4, label = groups), size = 8,, fontface = "bold") +
  scale_x_discrete(labels = c("12h/d" = expression(paste('12 h ',d^{-1})),
                              "14h/d" = expression(paste('14 h ',d^{-1})),
                              "16h/d" = expression(paste('16 h ',d^{-1}))
                              )
                   ) +
  scale_y_continuous(limits = c(0,100), breaks = c(0,25,50,75,100), expand = 0) +
  labs(x = "Photoperiod",
       y = expression(bold(paste("Leaf Area", " (", mm^2, ")"))),
       title = "(D)") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size=25,face="bold"),
        axis.text = element_text(size=20),
        plot.title = element_text(face="bold", size=25)) +
 scale_fill_aaas(alpha = 0.8)
print(p_leafarea)
ggsave("叶面积.tiff", device = "tiff", width = 5, height = 5, unit = "cm", dpi = 300)
ggsave("叶面积.png", device = "png", width = 5, height = 5, unit = "cm", dpi = 300)

################################ 周长 perimeter ############################
aov.perimeter <- aov(leafperimeter ~ light, data=da_400ppm)
aov.perimeter
summary(aov.perimeter)
# tukey <- aov(leafperimeter ~ treatment, data = da) |> TukeyHSD()
# tukey.leafperimeter <- as.data.frame(tukey$treatment)
# 
# view(tukey.leafperimeter)
# write.csv(tukey.leafperimeter,"leaf perimeter.csv")
# # 生成字母标记 #
# letters <- multcompLetters4(aov.leafperimeter, tukey)
# letters_df <- as.data.frame.list(letters$treatment) # 关键：取出$treatment
# letters_df$treatment <- rownames(letters_df)
# colnames(letters_df)[1] <- "letters"

letter.per <- HSD.test(aov.perimeter, trt = "light", alpha = 0.05)$groups

perimeter <- da_400ppm  |> 
  group_by(light) |> 
  summarise(
    mean = mean(leafperimeter, na.rm = TRUE),
    n    = n(),
    sd   = sd(leafperimeter, na.rm = TRUE),
    se   = sd / sqrt(n),
    .groups = 'drop'
  ) |> 
  left_join(letter.per |> rownames_to_column("light"), by = "light")

ggplot(perimeter, aes(x = light, y = mean)) +
  geom_col(fill = "grey", colour="black", width = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.2, linewidth = 1) +
  geom_text(aes(y = mean + se + 5, label = letter.per), size = 5) +
  labs(x = "treatment", y = "平均叶周长 mm") +
  theme_bw() +
  theme(panel.grid = element_blank())


################################ 比叶面积 sla #################################
myda <- read_excel("叶面积.xlsx", "summary")
myda <- subset.data.frame(myda, co2=="400ppm")

aov.sla <- aov(sla ~ light, data = myda)
aov.sla
summary(aov.sla)
# tukey <- aov(sla ~ light, data = myda) |> TukeyHSD()
# tukey.sla <- as.data.frame(tukey$treatment)
# 
# view(tukey.sla)
# write.csv(tukey.sla,"比叶面积.csv")
# # 生成字母标记 #
# letters <- multcompLetters4(aov.sla, tukey)
# letters_df <- as.data.frame.list(letters$treatment) # 关键：取出$treatment
# letters_df$treatment <- rownames(letters_df)
# colnames(letters_df)[1] <- "letters"

sla_letters <- HSD.test(aov.sla, "light", alpha=0.05)$groups

sla <- myda  |> 
  group_by(light)  |> 
  summarise(
    mean = mean(sla, na.rm = TRUE),
    n    = n(),
    sd   = sd(sla, na.rm = TRUE),
    se   = sd / sqrt(n)
  ) |> 
left_join(sla_letters |> rownames_to_column("light"), by = "light")## 将字母加入数据画图

p_sla <- ggplot(sla, aes(x = light, y = mean, fill = light)) +
  geom_col(width = 0.7, linewidth = 0.1, colour = "black") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean+se+1000, label = groups), size = 10, family = "arial", fontface="bold") +
  scale_x_discrete(labels = c("12h/d" = "12", "14h/d" = "14", "16h/d" = "16")) +
  scale_y_continuous(limits = c(0,20000), breaks = c(0,5000,10000,15000,20000),
                     labels = c(0,5,10,15,20), expand = 0) +
  labs(x = expression(bold(Photoperiod~(h~d^{-1}))),
       y = expression(bold(paste(`SLA (`,mm^2,` `,g^{-1},`, ×`,10^3,")"))),
       ) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family="arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20)) +
 scale_fill_brewer(palette = "YlOrRd")
p_sla
ggsave(filename="比叶面积.tiff", width=5, height=5, units="cm", dpi=300, device="tiff")
ggsave(filename="比叶面积.png", width=5, height=5, units="cm", dpi=300, device="png")

setwd("D:/ZheJiang A&F University/杜鹃文章/第一次收获/质量、根冠比")

### 加载依赖包 ###
library(readxl)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(dplyr)
library(multcompView)### 与下面的agricolae是两种字母显著性表示的方式
library(agricolae)### 见上文
library(tibble)
library(showtext)
library(ggsci)
library(patchwork)


#### 你的数据 ####
dry <- read_excel("鲜重干重、根冠比等.xlsx", "Sheet1")
dry <- subset.data.frame(dry, co2 == "400ppm")

##### 调用字体库 #####
font_add("Times", regular = "times.ttf", bold = "timesbd.ttf", italic = "timesi.ttf", bolditalic = "timesbi.ttf")
font_add("wryh", regular="msyh.ttc", bold="msyhbd.ttc")
font_add(family = "arial", regular = "arial.ttf", bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")
showtext_auto()

#### sd、se、mean举例 #####
library(dplyr)
library(ggplot2)
# 假设你的数据框叫 df，包含两列：group（分组），y（数值）
# head(df)
# 1) 计算均值、SE、样本量 n
df_sum <- df %>%
  group_by(group) %>%
  summarise(
    mean = mean(y, na.rm = TRUE),
    n    = sum(!is.na(y)),
    sd   = sd(y, na.rm = TRUE),
    se   = sd / sqrt(n)
  )
ggplot(df_sum, aes(x = group, y = mean)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.15, linewidth = 0.7) +
  labs(x = NULL, y = "Mean ± SE") +
  theme_classic(base_size = 12)
# 假设 df 有三列：treatment, light, y
# head(df)
pos <- position_dodge(width = 0.7)  # 分组并排的关键
df_sum2 <- df %>%
  group_by(treatment, light) %>%
  summarise(
    mean = mean(y, na.rm = TRUE),
    n    = sum(!is.na(y)),
    sd   = sd(y, na.rm = TRUE),
    se   = sd / sqrt(n),
    .groups = "drop"
  )
ggplot(df_sum2, aes(x = treatment, y = mean, fill = light)) +
  geom_col(position = pos, width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = pos, width = 0.2, linewidth = 0.7) +
  labs(x = NULL, y = "Mean ± SE", fill = "Light") +
  theme_classic(base_size = 12)
font_add("Times", regular = "times.ttf",
         bold = "timesbd.ttf",
         italic = "timesi.ttf",
         bolditalic = "timesbi.ttf")

#######################################

#### R/T 方差分析、多重比较、字母标记 ####
# 单因素方差分析 #

fit1 <- aov(R_T~light,dry)
summary(fit1) # p valve合适，做多重比较 #
fit1.1 <- TukeyHSD(fit1)
result1.1 <- as.data.frame(fit1.1$light)
result1.1$comparison <- rownames(fit1.1)
result1.1
write.csv(result1.1,"R_T.csv")
# 添加字母标记 groups是画图的字母列 #
letters.RT <- HSD.test(fit1,"light",alpha=0.05)$groups

#### 带入你的数据 计算re、rd、mean，合并字母标记 ####
da.RT <- dry |> 
  group_by(light) |> 
  summarise(
    mean = mean(R_T, na.rm = TRUE),
    n = sum(!is.na(R_T)),
    sd = sd(R_T, na.rm = TRUE),
    se = sd / sqrt(n)
  ) |> 
  left_join(letters.RT |> rownames_to_column("light"), by = "light")

#### R_T画图 ####
p_rt <- ggplot(da.RT, aes(x = light, y = mean, fill = light)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.1) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean+se+0.12, label = groups), size = 10, family = "arial", fontface = "bold") +
  scale_x_discrete(labels = c("12h/d" = "12", "14h/d" = "14", "16h/d" = "16")) +
  scale_y_continuous(limits = c(0,2.4), breaks = c(0,0.6,1.2,1.8,2.4), expand = 0) +
  labs(x = expression(bold(Photoperiod~(h~d^{-1}))),
       y = "R/S ratio") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size=25,face="bold"),
        axis.text = element_text(size=20)
        ) +
  scale_fill_brewer(palette = "YlOrRd")
p_rt
ggsave("根冠比.tiff", device="tiff", dpi=300, width=5, height=5, units="cm")
ggsave("根冠比.png", device="png", dpi=300, width=5, height=5, units="cm")




#### root ratio ####
# 单因素方差分析 #
fit2 <- aov(root_ratio ~ treatment, dry)
summary(fit2)
# p valve合适，做多重比较 #
fit2.1 <- TukeyHSD(fit2)
result2.1 <- as.data.frame(fit2.1$treatment)
write.csv(result2.1, "root ratio.csv")
# 赋予字母标记 #
letters <- multcompLetters4(fit2, fit2.1)
letters_df <- as.data.frame.list(letters$treatment) # 关键：取出$treatment
letters_df$treatment <- rownames(letters_df)
colnames(letters_df)[1] <- "letters"
# 带入数据，合并re、rd、mean、字母标记 #
da_sum <- dry |> 
  group_by(treatment) |> 
  summarise(
    mean = mean(root_ratio, na.rm = TRUE),
    n = sum(!is.na(root_ratio)),
    sd = sd(root_ratio, na.rm = TRUE),
    se = sd / sqrt(n)
  ) |> 
  left_join(letters_df, by = "treatment")
# 画图 #
ggplot(da_sum, aes(x = treatment, y = mean)) +
  geom_col(fill = "grey40", width = 0.7)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.2, linewidth = 1) +
  geom_text(aes(y = mean + se + 0.03, label = letters), size = 5) +
  labs(x = "treatment", y = "root ratio") +
  theme_bw() +
  theme(panel.grid = element_blank())

#### stem ratio ####
# 单因素方差分析 #
fit3 <- aov(stem_ratio ~ treatment, dry)
summary(fit3)
# p valve合适，做多重比较 #
fit3.1 <- TukeyHSD(fit3)
result3.1 <- as.data.frame(fit3.1$treatment)
write.csv(result3.1, "stem ratio.csv")
# 赋予字母标记 #
letters <- multcompLetters4(fit3, fit3.1)
letters_df <- as.data.frame.list(letters$treatment) # 关键：取出$treatment
letters_df$treatment <- rownames(letters_df)
colnames(letters_df)[1] <- "letters"
# 带入数据，合并re、rd、mean、字母标记 #
da_sum <- dry |> 
  group_by(treatment) |> 
  summarise(
    mean = mean(stem_ratio, na.rm = TRUE),
    n = sum(!is.na(stem_ratio)),
    sd = sd(stem_ratio, na.rm = TRUE),
    se = sd / sqrt(n)
  ) |> 
  left_join(letters_df, by = "treatment")
# 画图 #
ggplot(da_sum, aes(x = treatment, y = mean)) +
  geom_col(fill = "grey40", width = 0.7)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.2, linewidth = 1) +
  geom_text(aes(y = mean + se + 0.02, label = letters), size = 5) +
  labs(x = "treatment", y = "stem ratio") +
  theme_bw() +
  theme(panel.grid = element_blank())


#### leaf ratio ####
# 单因素方差分析 #
fit4 <- aov(leaf_ratio ~ treatment, dry)
summary(fit4)
# p valve合适，做多重比较 #
fit4.1 <- TukeyHSD(fit4)
result4.1 <- as.data.frame(fit4.1$treatment)
write.csv(result4.1, "leaf ratio.csv")
# 赋予字母标记 #
letters <- multcompLetters4(fit4, fit4.1)
letters_df <- as.data.frame.list(letters$treatment) # 关键：取出$treatment
letters_df$treatment <- rownames(letters_df)
colnames(letters_df)[1] <- "letters"
# 带入数据，合并re、rd、mean、字母标记 #
da_sum <- dry |> 
  group_by(treatment) |> 
  summarise(
    mean = mean(leaf_ratio, na.rm = TRUE),
    n = sum(!is.na(leaf_ratio)),
    sd = sd(leaf_ratio, na.rm = TRUE),
    se = sd / sqrt(n)
  ) |> 
  left_join(letters_df, by = "treatment")
# 画图 #
ggplot(da_sum, aes(x = treatment, y = mean)) +
  geom_col(fill = "grey40", width = 0.7)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.2, linewidth = 1) +
  geom_text(aes(y = mean + se + 0.02, label = letters), size = 5) +
  labs(x = "treatment", y = "leaf ratio") +
  theme_bw() +
  theme(panel.grid = element_blank())



#### 双因素 ####
dry_sum2 <- dry %>%
  group_by(co2, light, treatment) %>%
  summarise(
    mean = mean(R_T, na.rm = TRUE),
    n    = sum(!is.na(R_T)),
    sd   = sd(R_T, na.rm = TRUE),
    se   = sd / sqrt(n),
    .groups = "drop"
  )
fit5 <- aov(R_T~ co2 + light, dry)
summary(fit5)

# p valve合适，做多重比较 #
fit5.1 <- TukeyHSD(fit5)
result5.1 <- as.data.frame(fit5.1$co2)
write.csv(result1.1,"R_T co2.csv")
result5.1 <- as.data.frame(fit5.1$light)
write.csv(result1.1,"R_T light.csv")

# 赋予字母标记 #
letters <- multcompLetters4(fit5, fit5.1)
letters_df <- as.data.frame.list(letters$co2, letters$light) # 关键：取出 $co2 和 $light
letters_df$treatment <- rownames(letters_df)
colnames(letters_df)[1] <- "letters"

ggplot(dry_sum2, aes(x = treatment, y = mean)) +
  geom_col(fill = "grey40", width = 0.7)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.2, linewidth = 1) +
  geom_text(aes(y = mean + se + 0.1, label = letters), size = 5) +
  labs(x = "treatment", y = "root top ratio (R/T)", title = "根冠比") +
  theme_bw() +
  theme(panel.grid = element_blank())






#### 干物质含量 ####
da <- subset(dry, co2=="400ppm")
fit.aov <- aov(LDMC~light, data = da)
summary(fit.aov)## p = 0.00229
LDMC.letter <- HSD.test(fit.aov, "light", alpha=0.05)$groups
da.ldmc <- da |> 
  group_by(light) |> 
  summarise(
    mean = mean(LDMC, na.rm = TRUE),
    N = n(),
    sd = sd(LDMC, na.rm = TRUE),
    se = sd/sqrt(N)
  ) |> 
  left_join(LDMC.letter |> rownames_to_column("light"), by = "light")

ggplot(da.ldmc, aes(x = light, y = mean)) +
  geom_col(fill = "grey", colour="black", width = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean+se+0.02, label=groups),size =10,family="Times",face="bold") +
  labs(x="Photoperiod (h/d)", y="LDMC", title="干物质含量") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(family = "Times"),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 20)
        )
ggsave("干物质含量.png",device="png",width=8,height=8,dpi=300,units="cm")


####################################### 生物量 ##############################
da <- subset.data.frame(dry, co2=="400ppm")
tapply(X = da$biomass, INDEX = da$light, FUN = "shapiro.test")

fit.aov <- aov(biomass~light, data = da)
summary(fit.aov)## p = 0.00229
TukeyHSD(fit.aov)
tk.letter <- HSD.test(fit.aov, "light", alpha=0.05)$groups
da.bm <- da |> 
  group_by(light) |> 
  summarise(
    mean = mean(biomass, na.rm = TRUE),
    N = n(),
    sd = sd(biomass, na.rm = TRUE),
    se = sd/sqrt(N)
  ) |> 
  left_join(tk.letter |> rownames_to_column("light"), by = "light")

p_biomass <- ggplot(da.bm, aes(x = light, y = mean, fill = light)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean+se+10, label = groups), size = 10, family = "arial", fontface = "bold") +
  scale_x_discrete(labels = c("12h/d" = "12",
                              "14h/d" = "14",
                              "16h/d" = "16")) +
  scale_y_continuous(limits = c(0,200), breaks = c(0,50,100,150,200), expand = 0) +
  labs(x = expression(bold(Photoperiod~(h~d^{-1}))),
       y = "Total Biomass (g)"
       ) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20)
        ) +
 scale_fill_brewer(palette = "YlOrRd")
p_biomass
ggsave(filename="生物量.tiff", width=5, height=5, device="tiff", dpi=300, units="cm")
ggsave(filename="生物量.png", width=5, height=5, device="png", dpi=300, units="cm")



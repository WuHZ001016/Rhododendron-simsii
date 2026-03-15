# 加载所需的库
library(readxl)  # 读取Excel文件
library(ggplot2) # 绘制箱线图
library(dplyr)   # 数据处理
library(ggsignif) # 显著性标记

library(readxl)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(dplyr)
library(multcompView)### 与下面的 agricolae 是两种字母显著性表示的方式
library(agricolae)### 见上文
library(tibble)
library(showtext)
library(RColorBrewer)
library(ggsci)
library(patchwork)

setwd("D:/ZheJiang A&F University/杜鹃文章/第一次收获/株高地径")

##### 调用字体库 #####
font_add("Times", regular = "times.ttf", bold = "timesbd.ttf", italic = "timesi.ttf", bolditalic = "timesbi.ttf")
font_add("wryh", regular="msyh.ttc", bold="msyhbd.ttc")
font_add(family = "arial", regular = "arial.ttf", bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")
showtext_auto()

###### 株高 #########

# 读取数据
data <- read_excel("编号  时间.棵树.房间(2).xlsx", sheet = "Sheet3")

# 查看数据结构
str(data)

# 筛选A12、A14、A16、E12、E14的数据
# 这里假设“处理”列存储了不同的处理（A12, A14, A16, E12, E14）
data_subset <- data %>%
  filter(处理 %in% c("A12", "A14", "A16", "E12", "E14"))

# A12、A14、A16 之间的株高差箱线图和方差分析
data_A <- data_subset %>% filter(处理 %in% c("A12", "A14", "A16"))
ggplot(data_A, aes(x = 处理, y = 株高差, fill = 处理)) +
  geom_boxplot() +
  stat_compare_means(method = "anova") +
  theme_classic() +
  ggtitle("A12, A14, A16 株高差") +
  theme(legend.position = "none")

# E12、E14 之间的株高差箱线图和方差分析
data_E <- data_subset %>% filter(处理 %in% c("E12", "E14"))
ggplot(data_E, aes(x = 处理, y = 株高差, fill = 处理)) +
  geom_boxplot() +
  stat_compare_means(method = "anova") +
  theme_classic() +
  ggtitle("E12, E14 株高差") +
  theme(legend.position = "none")

# A12 和 E12 之间的株高差箱线图和方差分析
data_AE <- data_subset %>% filter(处理 %in% c("A12", "E12"))
ggplot(data_AE, aes(x = 处理, y = 株高差, fill = 处理)) +
  geom_boxplot() +
  stat_compare_means(method = "anova") +
  theme_classic() +
  ggtitle("A12 和 E12 株高差") +
  theme(legend.position = "none")

# A14 和 E14 之间的株高差箱线图和方差分析
data_AE2 <- data_subset %>% filter(处理 %in% c("A14", "E14"))
ggplot(data_AE2, aes(x = 处理, y = 株高差, fill = 处理)) +
  geom_boxplot() +
  stat_compare_means(method = "anova") +
  theme_classic() +
  ggtitle("A14 和 E14 株高差") +
  theme(legend.position = "none")



###### 地径 ########

# 加载所需的库
library(readxl)   # 读取Excel文件
library(ggplot2)  # 绘制箱线图
library(dplyr)    # 数据处理
library(ggsignif) # 显著性标记
library(rstatix)  # 统计分析和多重比较

# 读取数据
data <- read_excel("编号  时间.棵树.房间(2).xlsx", sheet = "Sheet3")

# 查看数据结构
str(data)

# 筛选A12、A14、A16、E12、E14的数据
# 假设“处理”列存储了不同的处理（A12, A14, A16, E12, E14）
data_subset <- data %>%
  filter(处理 %in% c("A12", "A14", "A16", "E12", "E14"))

# A12、A14、A16之间的地径差箱线图和方差分析
data_A <- data_subset %>% filter(处理 %in% c("A12", "A14", "A16"))
# 单因素方差分析
anova_A <- aov(地径差 ~ 处理, data = data_A)
summary(anova_A)
# 多重比较
pairwise_A <- data_A %>% 
  anova_test(dv = 地径差, between = 处理) %>%
  adjust_pvalue(method = "bonferroni")

ggplot(data_A, aes(x = 处理, y = 地径差, fill = 处理)) +
  geom_boxplot() +
  stat_compare_means(method = "anova") +
  stat_compare_means(comparisons = list(c("A12", "A14"), c("A12", "A16"), c("A14", "A16")),
                     label = "p.signif") +
  theme_classic() +
  ggtitle("A12, A14, A16 地径差") +
  theme(legend.position = "none")

# E12、E14之间的地径差箱线图和方差分析
data_E <- data_subset %>% filter(处理 %in% c("E12", "E14"))
# 单因素方差分析
anova_E <- aov(地径差 ~ 处理, data = data_E)
summary(anova_E)
# 多重比较
pairwise_E <- data_E %>% 
  anova_test(dv = 地径差, between = 处理) %>%
  adjust_pvalue(method = "bonferroni")

ggplot(data_E, aes(x = 处理, y = 地径差, fill = 处理)) +
  geom_boxplot() +
  stat_compare_means(method = "anova") +
  stat_compare_means(comparisons = list(c("E12", "E14")),
                     label = "p.signif") +
  theme_classic() +
  ggtitle("E12, E14 地径差") +
  theme(legend.position = "none")

# A12 和 E12 之间的地径差箱线图和方差分析
data_AE <- data_subset %>% filter(处理 %in% c("A12", "E12"))
# 单因素方差分析
anova_AE <- aov(地径差 ~ 处理, data = data_AE)
summary(anova_AE)
# 多重比较
pairwise_AE <- data_AE %>% 
  anova_test(dv = 地径差, between = 处理) %>%
  adjust_pvalue(method = "bonferroni")

ggplot(data_AE, aes(x = 处理, y = 地径差, fill = 处理)) +
  geom_boxplot() +
  stat_compare_means(method = "anova") +
  stat_compare_means(comparisons = list(c("A12", "E12")),
                     label = "p.signif") +
  theme_classic() +
  ggtitle("A12, E12 地径差") +
  theme(legend.position = "none")

# A14 和 E14 之间的地径差箱线图和方差分析
data_AE2 <- data_subset %>% filter(处理 %in% c("A14", "E14"))
# 单因素方差分析
anova_AE2 <- aov(地径差 ~ 处理, data = data_AE2)
summary(anova_AE2)
# 多重比较
pairwise_AE2 <- data_AE2 %>% 
  anova_test(dv = 地径差, between = 处理) %>%
  adjust_pvalue(method = "bonferroni")

ggplot(data_AE2, aes(x = 处理, y = 地径差, fill = 处理)) +
  geom_boxplot() +
  stat_compare_means(method = "anova") +
  stat_compare_means(comparisons = list(c("A14", "E14")),
                     label = "p.signif") +
  theme_classic() +
  ggtitle("A14，E14 地径差") +
  theme(legend.position = "none")


###### 肌醇1磷酸合成酶 ######
# 加载包
library(ggplot2)
library(multcompView)  # 用于显著性字母标注
library(agricolae)     # 也可以做 LSD.test
library(dplyr)

# 读取数据
data <- read.csv("肌醇1磷酸合成酶.csv")

# 确认处理列为因子
data$处理 <- as.factor(data$处理)

# -------- 分析函数：给定分组名，画箱线图+ANOVA+多重比较 --------
compare_groups <- function(df, groups, title){
  sub <- df %>% filter(处理 %in% groups)
  
  # 单因素方差分析
  model <- aov(MIPS1 ~ 处理, data = sub)
  print(summary(model))
  
  # 多重比较（这里用 TukeyHSD）
  tuk <- TukeyHSD(model)
  print(tuk)
  
  # 提取分组显著性字母
  tuk.cld <- multcompLetters4(model, tuk)
  letters <- as.data.frame.list(tuk.cld$`处理`)
  letters$处理 <- rownames(letters)
  
  # 合并数据用于绘图
  plot_data <- sub %>%
    left_join(letters, by = "处理")
  
  # 画图
  p <- ggplot(plot_data, aes(x=处理, y=MIPS1, fill=处理)) +
    geom_boxplot() +
    geom_jitter(width=0.2, alpha=0.6) +
    geom_text(aes(label=Letters, y=max(MIPS1)+5), vjust=0) +
    theme_bw() +
    labs(title = title, y="MIPS1", x="处理")
  
  print(p)
}

# -------- 分别比较 --------
compare_groups(data, c("A12","A14","A16"), "A12 vs A14 vs A16")
compare_groups(data, c("E12","E14"), "E12 vs E14")
compare_groups(data, c("A12","E12"), "A12 vs E12")
compare_groups(data, c("A14","E14"), "A14 vs E14")


########################### 株高、地径 new ###############################
### new shoot growth ###
da <- read_excel("data.xlsx","Sheet3")
A_da <- subset.data.frame(da, co2 == "400ppm")
tapply(X = da$株高差, INDEX = da$light, FUN = "shapiro.test")
height.aov <- aov(株高差~light, A_da)
height.aov |> summary()
height.letters <- HSD.test(height.aov, "light", alpha=0.05)$groups
height <- A_da |> 
  group_by(light) |>
  summarise(
    mean=mean(株高差, na.rm=TRUE),
    n=n(),
    sd=sd(株高差, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(height.letters |> rownames_to_column("light"), by="light")
p_height <- ggplot(height, aes(x=light, y=mean, fill=light)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.1) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean+se+1, label = groups),
            size = 10, family = "arial", fontface="bold") +
  labs(x = expression(bold(Photoperiod~(h~d^{-1}))),
       y = "New Shoot Growth (cm)") +
  scale_y_continuous(limits=c(0,20), breaks=c(0,5,10,15,20), expand=0) +
  scale_x_discrete(labels = c("12h" = "12", "14h" = "14", "16h" = "16")) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20)
        ) +
 scale_fill_brewer(palette = "YlOrRd")
p_height
ggsave(filename="新梢生长.tiff", width=5, height=5, units="cm", dpi=300, device="tiff")
ggsave(filename="新梢生长.png", width=5, height=5, units="cm", dpi=300, device="png")

### 地径 ###
ground_diameter.aov <- aov(地径差~light, A_da)
ground_diameter.aov |> summary()
ground_diameter.letters <- HSD.test(ground_diameter.aov, "light", alpha=0.05)$groups

ground_diameter <- A_da |> 
  group_by(light) |>
  summarise(
    mean=mean(地径差, na.rm=TRUE),
    n=n(),
    sd=sd(地径差, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(ground_diameter.letters |> rownames_to_column("light"), by="light")

ggplot(ground_diameter, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 0.3, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)", y = "Plant ground diameter Difference (mm)", title="地径差") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        scale_y_continuous(limits=c(0,22), breaks=c(0,5,10,15,20)),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        title = element_text(family="wryh", face="bold", size=30))
ggsave(filename="地径差.png", width=8, height=8, units="cm", dpi=300, device="png")

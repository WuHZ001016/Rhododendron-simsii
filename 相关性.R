setwd("D:/ZheJiang A&F University/杜鹃文章/第一次收获/相关性分析")

library(readxl)
library(ggplot2)
library(ggcorrplot)
library(ggpubr)
library(corrplot)
library(pheatmap)
library(GGally)
library(showtext)
library(dplyr)
library(purrr)
library(tidyr)
library(ggsci)
library(patchwork)
library(ggtext)
library(openxlsx)

##### 调用字体库 #####
font_add("times", regular = "times.ttf", bold = "timesbd.ttf", italic = "timesi.ttf", bolditalic = "timesbi.ttf")
font_add("wryh", regular="msyh.ttc", bold="msyhbd.ttc")
font_add(family = "arial", regular = "arial.ttf", bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")
showtext_auto()

############################### biomass #####################################

##### 数据 #####
da <- read_excel("生物量相关.xlsx", "Sheet1")
da$photoperiod <- factor(da$photoperiod, levels = c("12h/d", "14h/d", "16h/d"))
da$id <- factor(da$id)

tapply(X = da$MIPS1, INDEX = da$photoperiod, FUN = "shapiro.test")
shapiro.test(da$MIPS1)

cand_vars <- c(
  "glucose", "fructose", "sucrose", "maltose", "starch",  # 糖 + 淀粉
  "AGPase", "SSS", "TPS", "RUBP", "beta_amylase",         # 淀粉代谢/光合相关酶
  "SPS", "FBP", "CA", "MIPS1"                             # 糖代谢/光呼吸相关
  # 如果你想把 BCA、R_T 也算进去，可以在这里加： "R_T", "BCA"
)

##### 相关分析 #####

get_cor <- function(df, yvar, xvars, method = "pearson") {
  map_dfr(xvars, \(x) {
    ct <- cor.test(df[[yvar]], df[[x]], method = method)
    tibble(
      y        = yvar,
      x        = x,
      r        = unname(ct$estimate),   # 相关系数
      p        = ct$p.value,            # p 值
      n        = sum(complete.cases(df[[yvar]], df[[x]]))  # 样本数
    )
  })
}

### pearson相关系数 ###
cor_biomass <- get_cor(da, "biomass", cand_vars) %>%
  arrange(p)

### spearman相关系数 ###
cor_biomass <- get_cor(da, "biomass", cand_vars, method = "spearman") %>%
  arrange(p)

cor_biomass

##### 画图 #####
### 选出画图的参数 ###
vars_plot <- c("glucose", "sucrose", "AGPase", "MIPS1", "starch", "SSS")

da_long <- da %>%
  select(photoperiod, biomass, all_of(vars_plot)) %>%
  pivot_longer(cols = all_of(vars_plot),
               names_to = "variable",
               values_to = "value")

ggplot(da_long, aes(x = value, y = biomass, color = photoperiod)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ variable, scales = "free_x") +
  theme_classic() +
  labs(x = NULL, y = "Biomass", color = "Photoperiod") +
  theme(text = element_text(family = "times"),
        axis.title.y = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        legend.title = element_text(family="times", face = "bold", size = 20),
        legend.text = element_text(family = "times", size = 20),
        strip.text.x = element_text(size = 25, face = "bold"))
ggsave("biomass_pearson.png", width = 18, height = 12, device = "png", units = 'cm', dpi = 300)


############################# 净光合速率 Pn、Anet #####
##### 数据 #####
da <- read_excel("光合速率相关.xlsx", "Sheet1")
da$photoperiod <- factor(da$photoperiod, levels = c("12h/d", "14h/d", "16h/d"))
da$id <- factor(da$id)
tapply(da$A, INDEX = da$photoperiod, FUN = "shapiro.test")

##### 确定相关性数据 #####
cand_vars_Pn <- c(
  "glucose", "fructose", "sucrose", "maltose", "starch",  # 糖 + 淀粉
  "AGPase", "SSS", "TPS", "RUBP", "beta_amylase",         # 代谢/淀粉酶/光合相关
  "chlorophylla", "chlorophyllb", "chlorophyll",          # 叶绿素
  "SPS", "FBP", "CA", "MIPS1",                            # 糖代谢/光周期相关酶
  "C", "N", "C_N"                                         # C、N 及其比值
)

##### 相关性分析 #####

get_cor <- function(df, yvar, xvars, method = "pearson") {
  map_dfr(xvars, \(x) {
    ct <- cor.test(df[[yvar]], df[[x]], method = method)
    tibble(
      y = yvar,
      x = x,
      r = unname(ct$estimate),                          # 相关系数
      p = ct$p.value,                                   # p 值
      n = sum(complete.cases(df[[yvar]], df[[x]]))      # 有效样本数
    )
  })
}

cor_Pn <- get_cor(da, "A", cand_vars_Pn, method = "pearson") %>%
  arrange(p)
cor_Pn
write.csv(cor_Pn, "Anet_pearson.csv")
cor_Pn <- get_cor(da, "A", cand_vars_Pn, method = "spearman") %>%
  arrange(p)
cor_Pn
write.csv(cor_Pn, "Anet_spearman.csv")

##### 画图 #####
# 选几个最有代表性的指标 #
vars_Pn_plot <- c("sucrose", "maltose", "glucose", "MIPS1", "chlorophylla", "chlorophyll")

da_long_Pn <- da %>%
  select(photoperiod, A, all_of(vars_Pn_plot)) %>%
  pivot_longer(cols = all_of(vars_Pn_plot),
               names_to = "variable",
               values_to = "value")

# 多面板图（推荐，论文里好放）#
ggplot(da_long_Pn, aes(x = value, y = A, color = photoperiod)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ variable, scales = "free_x") +
  theme_classic() +
  labs(x = NULL,
       y = expression(bold(paste("A (μmol ", CO[2], m^{-2}, s^{-1}, ")"))),
       color = "Photoperiod") +
  theme(text = element_text(family = "times"),
        axis.title.y = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        legend.title = element_text(family="times", face = "bold", size = 25),
        legend.text = element_text(family = "times", size = 20),
        strip.text.x = element_text(size = 25, face = "bold"))
ggsave("Anet.png", width = 18, height = 12, device = "png", units = 'cm', dpi = 300)



##### 用矩阵热图分析 #######
## ================================
## 1. 加载包 & 导入数据
## ================================
# 如尚未安装，请先运行：
# install.packages(c("tidyverse", "Hmisc"))

library(tidyverse)
library(Hmisc)

# 读入数据
dat_raw <- read.csv("biomass.csv", header = TRUE)

# 只保留数值型变量，并去掉含 NA 的行（避免 cor.test 报错）
dat_num <- dat_raw %>%
  dplyr::select(where(is.numeric)) %>%
  na.omit()

# 变量名顺序
vars <- colnames(dat_num)

## ================================
## 2. 分别计算 Pearson & Spearman 的 r 与 p
## ================================
# Hmisc::rcorr 会同时返回相关系数矩阵 r 和 p 值矩阵 P
pearson_res  <- rcorr(as.matrix(dat_num), type = "pearson")
spearman_res <- rcorr(as.matrix(dat_num), type = "spearman")

r_pearson <- pearson_res$r
p_pearson <- pearson_res$P

r_spear   <- spearman_res$r
p_spear   <- spearman_res$P

# 有时候 rcorr 的 p 矩阵上三角是 NA，这里把它补成对称矩阵
p_pearson[upper.tri(p_pearson)] <- t(p_pearson)[upper.tri(p_pearson)]
p_spear[upper.tri(p_spear)]    <- t(p_spear)[upper.tri(p_spear)]

## ================================
## 3. 整理成长表（长格式），构造上三角 Pearson、下三角 Spearman
## ================================
# 构造所有 (var1, var2) 的组合
grid <- expand.grid(
  var1 = vars,
  var2 = vars,
  stringsAsFactors = FALSE
)

# 给每个组合一个索引 i, j（对应在矩阵中的行列）
grid$i <- match(grid$var1, vars)
grid$j <- match(grid$var2, vars)

grid2 <- grid %>%
  mutate(
    # 上三角（i < j）使用 Pearson，下三角（i > j）使用 Spearman
    method = dplyr::case_when(
      i < j ~ "Pearson",
      i > j ~ "Spearman",
      TRUE ~ "Diag"            # 对角线
    ),
    # 按 method 选择相关系数 r
    r = dplyr::case_when(
      method == "Pearson"  ~ r_pearson[cbind(i, j)],
      method == "Spearman" ~ r_spear[cbind(i, j)],
      TRUE ~ 1            # 对角线自己和自己相关为 1
    ),
    # 按 method 选择 p 值
    p = dplyr::case_when(
      method == "Pearson"  ~ p_pearson[cbind(i, j)],
      method == "Spearman" ~ p_spear[cbind(i, j)],
      TRUE ~ NA_real_
    ),
    # 显著性标记
    signif = dplyr::case_when(
      is.na(p)      ~ "",
      p < 0.01      ~ "***",
      p < 0.05      ~ "**",
      p < 0.1       ~ "*",
      TRUE          ~ ""
    ),
    # 圆圈大小与 p 值相关：这里用 -log10(p) 表示显著性（p 越小，值越大，圆越大）
    size_val = ifelse(is.na(p), NA_real_, -log10(p + 1e-10))
  )

# 如果不想画对角线，可以去掉 method == "Diag"
plot_df <- grid2 %>%
  filter(method != "Diag", !is.na(r))

# 控制坐标轴顺序，使矩阵图对称、上三角在右上
plot_df <- plot_df %>%
  mutate(
    var1 = factor(var1, levels = vars),         # x 轴
    var2 = factor(var2, levels = rev(vars))     # y 轴反转，视觉上更像相关矩阵
  )

## ================================
## 4. 绘图：上三角 Pearson，下三角 Spearman
## ================================
library(ggplot2)
library(grid)   # 用于 unit()

p_corr <- ggplot(plot_df, aes(x = var1, y = var2)) +
  # 先画圆：大小 = -log10(p)，颜色 = r
  geom_point(
    aes(size = size_val, fill = r),
    shape  = 21,             # 空心圆，可单独控制描边和填充
    color  = "grey60",       # 浅灰色描边，让圆在白底上也能看见
    stroke = 0.3             # 描边宽度
  ) +
  # 再叠加显著性星号
  geom_text(
    aes(label = signif),
    size = 8
  ) +
  # 颜色：负相关蓝，接近 0 为浅灰，正相关红
  scale_fill_gradient2(
    low  = "#3B4992FF",      # 强负相关
    mid  = "#F0F0F0",        # 用浅灰代替纯白，避免“隐身”
    high = "#BB0021FF",      # 强正相关
    limits = c(-1, 1),
    name  = "Pearson'r\nSpearman'ρ",
    guide = guide_colorbar(
      # 控制色条高度/宽度，使其和图本身更协调
      barheight = unit(3, "cm"),
      barwidth  = unit(0.5, "cm"),
      # 单独给颜色图例的标题设置样式，关键是 lineheight
      title.theme = element_text(
        family    = "Times",
        size      = 30,
        face      = "bold",
        lineheight = 0.3      # <--- 减小两行之间间距的关键
      )
    )
  ) +
  # 圆的大小范围 + 大小图例的样式
  scale_size_continuous(
    name  = "-log10 (p)",
    range = c(2, 7),
    breaks = 0:5,
    guide = guide_legend(
      title.theme = element_text(
        family = "Times",
        size   = 30,
        face   = "bold"
      ),
      label.theme = element_text(
        family = "Times",
        size   = 25
      ),
      # 图例示例圆的外观
      override.aes = list(
        fill  = "grey80",
        color = "grey60"
      )
    )
  ) +
  coord_equal() +
  labs(title = "Pearson correlation",
       y = "Spearman correlation") +
  theme_bw() +
  theme(
    text = element_text(family = "Times"),
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white"),
    
    plot.title = element_text(size = 35, face = "bold", hjust = 1),
    
    axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    axis.title.y = element_text(size = 35, face = "bold", hjust = 0),
    axis.text.x  = element_text(angle = 40, hjust = 1, size = 25),
    axis.text.y  = element_text(angle = -30, vjust = 1, size = 25),
    
    legend.position   = "right",
    # 这两个控制整体图例块大小
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(0.5, "cm"),
    legend.spacing.y  = unit(0.3, "cm"),
    
    # 全局的图例标题再统一一下（对单行的 -log10(p) 也生效）
    legend.title = element_text(
      size      = 30,
      face      = "bold",
      lineheight = 1.2
    ),
    legend.text = element_text(size = 25)
  )

print(p_corr)

# 建议图整体稍微放大一点，这样图和图例视觉更平衡
ggsave("biomass相关性.png",
       plot   = p_corr,
       device = "png",
       width  = 15,  # 原来是 10 cm，可以适当放大
       height = 15,
       units  = "cm",
       dpi    = 300)









##################################### 生物量 相关指数 #####################################
# 1. 设置工作目录 (根据您的实际文件位置修改，如果文件在当前目录则不需要)
# setwd("C:/Your/Path/Here")

# 2. 读取数据文件
# 请确保文件名与您的实际文件名一致
data <- read_excel("生物量相关.xlsx", "Sheet1") # check.names=FALSE 防止列名中的特殊字符被修改

# 3. 数据预处理
# 排除非数值列 (根据您的描述，前两列 'photoperiod' 和 'id' 不需要参与计算)
# 如果您的数据列顺序不同，请调整这里的索引
analysis_data <- data[, -c(1, 2)] 

# 确保所有列都是数值型 (以防万一)
analysis_data <- as.data.frame(lapply(analysis_data, as.numeric))

# 4. 初始化结果容器
results <- data.frame(
 Indicators = character(),
 Correlation = numeric(),
 P_value = numeric(),
 stringsAsFactors = FALSE
)

# 获取参与分析的列名
col_names <- colnames(analysis_data)
n <- length(col_names)

# 5. 循环计算每一对指标的相关性
# 使用双重循环遍历所有唯一的组合
for (i in 1:(n-1)) {
 for (j in (i+1):n) {
  var1 <- col_names[i]
  var2 <- col_names[j]
  
  # 执行 Pearson 相关性检验
  # use = "complete.obs" 处理缺失值（如果有）
  test_result <- cor.test(analysis_data[[var1]], analysis_data[[var2]], 
                          method = "spearman", use = "complete.obs")
  
  # 提取相关系数和 p 值
  cor_value <- as.numeric(test_result$estimate)
  p_val <- as.numeric(test_result$p.value)
  
  # 组合指标名，例如 "biomass-R_T"
  pair_name <- paste(var1, var2, sep = "-")
  
  # 将结果添加到数据框
  results <- rbind(results, data.frame(
   Indicators = pair_name,
   Correlation = cor_value,
   P_value = p_val
  ))
 }
}

# 6. 重命名列以符合您的要求
colnames(results) <- c("指标对", "相关性指数", "p值")

# 7. 输出结果到 CSV 文件
output_filename <- "生物量有关的spearman.csv"
write.csv(results, output_filename, row.names = FALSE)

print(paste("计算完成！结果已保存为:", output_filename))


####################################### 光合速率 相关指数 ####################################
# install.packages("openxlsx")  # 若未安装请先安装
library(openxlsx)

infile  <- "光合速率相关原始数据 对应light res.xlsx"
outfile <- "光合速率相关性结果 对应light res.xlsx"

# dat <- read.csv(infile, header = TRUE, stringsAsFactors = FALSE)
dat <- read_excel(infile, "光合速率相关原始数据 对应light res")
# 仅保留数值型列（如 photoperiod 为字符/因子会自动剔除）
num_dat <- dat[, sapply(dat, is.numeric), drop = FALSE]

calc_pairwise_cor <- function(df, method = c("pearson", "spearman")) {
 method <- match.arg(method)
 vars <- colnames(df)
 comb <- combn(vars, 2, simplify = FALSE)
 
 res <- lapply(comb, function(v) {
  x <- df[[v[1]]]
  y <- df[[v[2]]]
  ok <- complete.cases(x, y)
  n  <- sum(ok)
  
  if (n < 3) {
   r <- NA_real_
   p <- NA_real_
  } else {
   ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method))
   r  <- unname(ct$estimate)
   p  <- ct$p.value
  }
  
  data.frame(
   指标1 = v[1],
   指标2 = v[2],
   n = n,
   相关系数 = r,
   p_value = p,
   stringsAsFactors = FALSE
  )
 })
 
 do.call(rbind, res)
}

pearson_res  <- calc_pairwise_cor(num_dat, method = "pearson")
spearman_res <- calc_pairwise_cor(num_dat, method = "spearman")

# 为了区分列名，可重命名“相关系数”列
names(pearson_res)[names(pearson_res) == "相关系数"]  <- "Pearson_r"
names(spearman_res)[names(spearman_res) == "相关系数"] <- "Spearman_rho"

# 写入一个工作簿的两个工作表
wb <- createWorkbook()

addWorksheet(wb, "Pearson相关性")
writeData(wb, "Pearson相关性", pearson_res)

addWorksheet(wb, "Spearman相关性")
writeData(wb, "Spearman相关性", spearman_res)

saveWorkbook(wb, outfile, overwrite = TRUE)




############ biomass 点图+回归 ###########

### 生长-光合限制格局 biomass vs vcmax/jmax ###
p1 <- ggplot(data=da, mapping=aes(x=biomass, y=`Vcmax/Jmax`)) +
 geom_point(size=1, colour="grey30") +
 geom_smooth(colour="blue", method="glm", linewidth=1) +
 labs(x = "Biomass",
      y = expression(bold(italic(V)[bold(cmax)] / italic(J)[bold(max)]))) +
 theme_bw() +
 stat_cor(
  aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~~")),
  method = "spearman", label.x = 80, label.y = 0.37,
  fontface = "bold", size = 5, family = "arial"
 ) +
 theme(panel.grid = element_blank(), 
       text = element_text(family="arial"),
       axis.title = element_text(size=20, face="bold"),
       axis.text = element_text(size=15))
p1

### "糖积累—生长"的耦合 biomass vs SS ###
p2 <- ggplot(data=da, mapping=aes(x=biomass, y=`soluble sugar`)) +
 geom_point(size=1, colour="grey30") +
 geom_smooth(colour="blue", method="lm", linewidth=1) +
 labs(x = "Biomass", y = "Soluble sugar") + theme_bw() +
 stat_cor(
  aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~~")),
  method = "spearman", label.x = 80, label.y = 8.5,
  fontface = "bold", size = 5, family = "arial"
 ) +
 theme(panel.grid=element_blank(), 
       text=element_text(family="arial"),
       axis.text=element_text(size=15),
       axis.title = element_text(size = 20, face = "bold")
       )
p2

### 电子传递能力-叶绿素（把 Jmax的变化落到光系统/色素层面） vcmax vs chl ###
p3 <- ggplot(data=da, mapping=aes(x=Jmax, y=chlorophyll)) +
 geom_point(size=1, colour="grey30") +
 geom_smooth(colour="blue", method="lm", linewidth=1) +
 labs(x=NULL, y=NULL, title="Jmax ~ Chlorophyll") + theme_bw() +
 stat_cor(
  aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~~")),
  method = "spearman", label.x = 70, label.y = 0.55,
  fontface = "bold", size = 8, family = "arial"
 ) +
 theme(panel.grid=element_blank(), 
       text=element_text(family="arial"),
       plot.title=element_text(size=20, face="bold", hjust=0.5),
       axis.text=element_text(size=15))
p3

### 碳库-叶绿素（体现“储碳—光合结构”的权衡） NSC vs chl ###
p4 <- ggplot(data=da, mapping=aes(x=NSC, y=chlorophyll)) +
 geom_point(size=1, colour="grey30") +
 geom_smooth(colour="blue", method="lm", linewidth=1) +
 labs(x=NULL, y=NULL, title="NSC ~ Chlorophyll") + theme_bw() +
 stat_cor(
  aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~~")),
  method = "spearman", label.x = 35.5, label.y = 0.5,
  fontface = "bold", size = 8, family = "arial"
 ) +
 theme(panel.grid=element_blank(), 
       text=element_text(family="arial"),
       plot.title=element_text(size=20, face="bold", hjust=0.5),
       axis.text=element_text(size=15))
p4

### 含氮组分-储碳（蛋白与淀粉/NSC的“此消彼长”） protein vs starch ###
p5 <- ggplot(data=da, mapping=aes(x=protein, y=starch)) +
 geom_point(size=1, colour="grey30") +
 geom_smooth(colour="blue", method="lm", linewidth=1) +
 labs(x=NULL, y=NULL, title="Protein ~ Starch") + theme_bw() +
 stat_cor(
  aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~~")),
  method = "spearman", label.x = 25, label.y = 24,
  fontface = "bold", size = 8, family = "arial"
 ) +
 theme(panel.grid=element_blank(), 
       text=element_text(family="arial"),
       plot.title=element_text(size=20, face="bold", hjust=0.5),
       axis.text=element_text(size=15))
p5

### 生物量 vs mips1 ###
p6 <- ggplot(data=da, mapping=aes(x=biomass, y=MIPS1)) +
 geom_point(size=1, colour="grey30") +
 geom_smooth(colour="blue", method="lm", linewidth=1) +
 stat_cor(
  aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~~")),
  method = "spearman", label.x = 90, label.y = 80,
  fontface = "bold", size = 5, family = "arial"
  ) +
 labs(x = "Biomass", y = "MIPS1") +
 theme_bw() + 
 theme(panel.grid=element_blank(),
       axis.title=element_text(size=20, face="bold", hjust=0.5),
       text=element_text(family="arial"),
       axis.text=element_text(size=15))
p6

### 糖 vs mips1 ###
p7 <- ggplot(data=da, mapping=aes(x=`soluble sugar`, y=MIPS1)) +
 geom_point(size=1, colour="grey30") +
 geom_smooth(colour="blue", method="lm", linewidth=1) +
 stat_cor(
  aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~~")),
  method = "pearson", label.x = 9, label.y = 65,
  fontface = "bold", size = 5, family = "arial") +
 labs(x = "Soluble sugar", y = "MIPS1") +
 theme_bw() + 
 theme(panel.grid=element_blank(),
       axis.title=element_text(size=20, face="bold", hjust=0.5),
       text=element_text(family="arial"),
       axis.text=element_text(size=15))
p7






#### MIPS1 vs Vc/J ####
p3 <- ggplot(data=da, mapping = aes(x= MIPS1, y = `Vcmax/Jmax`)) +
 geom_point(size=1, colour="grey30") +
 geom_smooth(colour="blue", method="glm", linewidth=1) +
 labs(x=NULL, y=NULL, title="MIPS1 ~ Vcmax/Jmax") + theme_bw() +
 stat_cor(
  aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~~")),
  method = "spearman", label.x = 80, label.y = 0.4,
  fontface = "bold", size = 8, family = "arial"
 ) +
 # annotate("text", label=("r = -0.96\np = 0.02"), x=100, y=0.375, family="arial", fontface="bold.italic", size=8, lineheight=0.3) +
 theme(panel.grid=element_blank(), 
       text=element_text(family="arial"),
       plot.title=element_text(size=20, face="bold", hjust=0.5),
       axis.text=element_text(size=15))
p3

### 生物量 and MIPS1 vs 可溶糖 ###
p4 <- ggplot(data=da, mapping=aes(x=biomass, y=`soluble sugar`)) +
 geom_point(size=1, colour="grey30") +
 geom_smooth(colour="blue", method="glm", linewidth=1) +
 labs(title="Biomass ~ Soluble Sugar", x=NULL, y=NULL) + theme_bw() +
 stat_cor(
  aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~~")),
  method = "spearman", label.x = 80, label.y = 8.5,
  fontface = "bold", size = 8, family = "arial"
 ) +
 # annotate("text", label="r = -0.83\np = 0.01", x=100, y=8.5, family="arial", size=8, fontface="bold.italic", lineheight=0.3) +
 theme(panel.grid=element_blank(), 
       text=element_text(family="arial"),
       plot.title=element_text(size=20, face="bold", hjust=0.5),
       axis.text=element_text(size=15))
p4

p5 <- ggplot(data=da, mapping=aes(x=MIPS1, y=`soluble sugar`)) +
 geom_point(size=1, colour="grey30") +
 geom_smooth(colour="blue", method="glm", linewidth=1) +
 labs(title="MIPS1 ~ Soluble Sugar", x=NULL, y=NULL) + theme_bw() +
 stat_cor(
  aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~~")),
  method = "spearman", label.x = 80, label.y = 9.2,
  fontface = "bold", size = 8, family = "arial"
 ) +
 # annotate("text", label="r = -0.83\np = 0.01", x=100, y=8.5, family="arial", size=8, fontface="bold.italic", lineheight=0.3) +
 theme(panel.grid=element_blank(), 
       text=element_text(family="arial"),
       plot.title=element_text(size=20, face="bold", hjust=0.5),
       axis.text=element_text(size=15))
p5





#################### mips1入手 点图+回归 ################

p6 <- ggplot(data=da, mapping=aes(x=MIPS1, y=SLA)) +
 geom_point(size=2, colour="grey30") +
 geom_smooth(colour="blue", method="glm", linewidth=1.5) +
 scale_y_continuous(limits=c(3000,20000)) +
 labs(x="MIPS1", y="SLA") + theme_bw() + 
 annotate("text", label="r = -0.92, p < 0.001", x=95, y=5000, family="arial", size=10, fontface="bold.italic") +
 theme(panel.grid=element_blank(), 
       text=element_text(family="arial"),
       axis.title=element_text(size=25,face="bold"),
       axis.text=element_text(size=15))
p6

p7 <- ggplot(data=da, mapping=aes(x=MIPS1, y=`soluble sugar`)) +
 geom_point(size=2, colour="grey30") +
 geom_smooth(colour="blue", method="glm", linewidth=1.5) +
 scale_y_continuous(limits=c(8.5,12)) +
 labs(x="MIPS1", y="Soluble Sugar") + theme_bw() + 
 annotate("text", label="r = -0.9, p < 0.001", x=95, y=9.1, family="arial", size=10, fontface="bold.italic") +
 theme(panel.grid=element_blank(), 
       text=element_text(family="arial"),
       axis.title=element_text(size=25,face="bold"),
       axis.text=element_text(size=15))
p7

p8 <- ggplot(data=da, mapping=aes(x=MIPS1, y=Vcmax)) +
 geom_point(size=2, colour="grey30") +
 geom_smooth(colour="blue", method="glm", linewidth=1.5) +
 scale_y_continuous(limits=c(29,36)) +
 labs(x="MIPS1", y="Vcmax") + theme_bw() + 
 annotate("text", label="r = -0.92, p = 0.001", x=95, y=30, family="arial", size=10, fontface="bold.italic") +
 theme(panel.grid=element_blank(), 
       text=element_text(family="arial"),
       axis.title=element_text(size=25,face="bold"),
       axis.text=element_text(size=15))
p8

p9 <- ggplot(data=da, mapping=aes(x=MIPS1, y=Vcmax)) +
 geom_point(size=2, colour="grey30") +
 geom_smooth(colour="blue", method="glm", linewidth=1.5) +
 scale_y_continuous(limits=c(29,36)) +
 labs(x="MIPS1", y="Vcmax") + theme_bw() + 
 annotate("text", label="r = -0.92, p = 0.001", x=95, y=30, family="arial", size=10, fontface="bold.italic") +
 theme(panel.grid=element_blank(), 
       text=element_text(family="arial"),
       axis.title=element_text(size=25,face="bold"),
       axis.text=element_text(size=15))
p9


p <- p1+p2+p3+p4+p5 +
 plot_layout(design=c("BBAAADD
                      CCAAAEE"))
print(p)
ggsave("biomass.tiff", p, width=12, height=8, device="tiff",units="cm", dpi=300)
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
library(scales)
library(zoo)

setwd("D:/ZheJiang A&F University/杜鹃文章/第一次收获/mpea数据/收获")
da.fl <- read_excel("data.xlsx", "all_")
da.fl.A <- subset.data.frame(da.fl, co2=="400ppm")

#################### 调用字体库 ###################
font_add("Times", regular = "times.ttf", bold = "timesbd.ttf", italic = "timesi.ttf", bolditalic = "timesbi.ttf")
font_add("wryh", regular="msyh.ttc", bold="msyhbd.ttc")
font_add(family = "arial", regular = "arial.ttf", bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")
showtext_auto()

################################### OJIP #####################################

df_cleaned <- da.fl.A %>%
 group_by(light) %>%
 arrange(`Prompt T.(ms)`) %>%
 # 1. 首先剔除那些已知的、受干扰严重的特定时间点（扩大范围）
 filter(!(`Prompt T.(ms)` %in% c(2.0, 3.0, 4.0, 30.0))) %>%
 # 2. 对荧光值进行小窗口滚动平均，消除单点毛刺
 mutate(`Prompt Fl.` = rollapply(`Prompt Fl.`, width = 3, FUN = mean, fill = "extend")) %>%
 ungroup()

# 注意：filter 掉之后，ggplot 会自动根据 scale_x_log10 连接相邻点，
# 这比单纯用均值替换一个点效果更好，因为 log 坐标系下点非常密集。

p_ojip <- ggplot(data = df_cleaned, mapping = aes(x = `Prompt T.(ms)`, y = `Prompt Fl.`)) +
  geom_point(mapping = aes(colour = light), alpha = 0.8, size = 0.5) +
  # geom_smooth(mapping = aes(colour = light), method = "gam", show.legend = FALSE) +
  annotation_logticks(sides = "b") +
  scale_x_log10(
    breaks = c(0.1, 1, 10, 100, 1000), 
    labels = trans_format("log10", math_format(10^.x)), # 自动生成 10的次幂格式
    limits = c(NA, 3000) # 根据数据适当调整上限
  ) +
  scale_y_continuous(limits = c(4000,25000),
                     breaks = c(5000,10000,15000,20000,25000),
                     labels = c(5,10,15,20,25)) +
  scale_colour_manual(
    name = "Photoperiod",
    values = c("12h/d" = "#FFEDA0",
               "14h/d" = "#FD8D3C",
               "16h/d" = "#BD0026"),
    labels = c("12h/d" = expression(12~h~d^{-1}),
               "14h/d" = expression(14~h~d^{-1}),
               "16h/d" = expression(16~h~d^{-1})),
    breaks = c("12h/d", "14h/d", "16h/d")) +
  labs(x = "Time (ms)",
       y = expression(bold(paste("Fluorescence (×", bold(10^3), ")"))),
       ) +
  theme_classic() +
  theme(text = element_text(family = "arial"),
        legend.position = c(0.25, 0.7), 
        legend.background = element_blank(), # 移除图例背景框
        legend.key = element_blank(),
        legend.key.height = unit(0.05, "cm"),
        legend.key.width = unit(0.1, "cm"),
        legend.key.spacing.x = unit(0.1, "cm"),
        legend.title = element_text(face = "bold", size = 20),
        axis.title = element_text(size = 25, face="bold"),
        axis.text = element_text(size = 20),
        axis.line = element_line(linewidth = 0.3, colour = "black"),
        axis.ticks = element_line(linewidth = 0.3, colour = "black"),
        axis.ticks.length = unit(0.1, "cm"), # 刻度线稍长一点更专业
        legend.text = element_text(size=20))
p_ojip
ggsave("ojip_400ppm.tiff", device="tiff", dpi=300, width=7, height=5, units = "cm")
ggsave("ojip_400ppm.png", device="png", dpi=300, width=7, height=5, units = "cm")



############################## 暗叶绿素荧光参数数据 ##############################

setwd("D:/ZheJiang A&F University/杜鹃文章/第一次收获/mpea数据")

fl <- read_excel("D:/ZheJiang A&F University/杜鹃文章/第一次收获/mpea数据/荧光.xlsx", "flr")
fl <- subset.data.frame(fl, co2=="400ppm")

############################### 初始荧光Fo ####################################

fo.aov <- aov(`Fo`~light, fl)
fo.aov |> summary()
fo.letters <- HSD.test(fo.aov, "light", alpha=0.05)$groups

fo <- fl |> 
  group_by(light) |>
  summarise(
    mean=mean(`Fo`, na.rm=TRUE),
    n=n(),
    sd=sd(`Fo`, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(fo.letters |> rownames_to_column("light"), by="light")

ggplot(fo, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 300, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)",
       y = "Fo",
       title="初始荧光") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        plot.title = element_text(family="wryh", face="bold", size=35, hjust=0.5))
ggsave(filename="初始荧光Fo.png", width=8, height=8, units="cm", dpi=300, device="png")

############################## 最大荧光 Fm ####################################

fm.aov <- aov(`Fm`~light, fl)
fm.aov |> summary()
fm.letters <- HSD.test(fm.aov, "light", alpha=0.05)$groups

fm <- fl |> 
  group_by(light) |>
  summarise(
    mean=mean(`Fm`, na.rm=TRUE),
    n=n(),
    sd=sd(`Fm`, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(fm.letters |> rownames_to_column("light"), by="light")

ggplot(fm, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 1300, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)",
       y = "Fm",
       title="最大荧光") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        plot.title = element_text(family="wryh", face="bold", size=35, hjust=0.5))
ggsave(filename="最大荧光Fm.png", width=8, height=8, units="cm", dpi=300, device="png")

################################ 可变荧光 Fv ##################################

fv.aov <- aov(`Fv`~light, fl)
fv.aov |> summary()
fv.letters <- HSD.test(fv.aov, "light", alpha=0.05)$groups

fv <- fl |> 
  group_by(light) |>
  summarise(
    mean=mean(`Fv`, na.rm=TRUE),
    n=n(),
    sd=sd(`Fv`, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(fv.letters |> rownames_to_column("light"), by="light")

ggplot(fv, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 1300, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)",
       y = "Fv",
       title="可变荧光") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        plot.title = element_text(family="wryh", face="bold", size=35, hjust=0.5))
ggsave(filename="可变荧光Fv.png", width=8, height=8, units="cm", dpi=300, device="png")

################################ Fv/Fm ##################################

fvfm.aov <- aov(`Fv/Fm`~light, fl)
fvfm.aov |> summary()
fvfm.letters <- HSD.test(fvfm.aov, "light", alpha=0.05)$groups

fvfm <- fl |> 
  group_by(light) |>
  summarise(
    mean=mean(`Fv/Fm`, na.rm=TRUE),
    n=n(),
    sd=sd(`Fv/Fm`, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(fvfm.letters |> rownames_to_column("light"), by="light")

p_fv_fm <- ggplot(fvfm, aes(x=light, y=mean))+
  geom_col(aes(fill = light), width = 0.7, colour = "black", linewidth = 0.1)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = mean + se + 0.06, label = groups), size = 10, family = "arial", fontface="bold") +
  scale_x_discrete(labels = c("12h/d" = "12", "14h/d" = "14", "16h/d" = "16")) +
  scale_y_continuous(limits = c(0,1.2), breaks = c(0,0.3,0.6,0.9,1.2), expand = 0) +
  labs(x = expression(bold(Photoperiod~(h~d^{-1}))),
       y = expression(bold(paste(F[v],"/",F[m]))),
       ) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "arial"),
        axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20)
        ) +
  scale_fill_brewer(palette = "YlOrRd")
p_fv_fm
ggsave(filename="最大光化学效率Fv_Fm.tiff", width=5, height=5, units="cm", dpi=300, device="tiff")
ggsave(filename="最大光化学效率Fv_Fm.png", width=5, height=5, units="cm", dpi=300, device="png")


########################### 初始荧光/最大荧光 Fv/Fo ###########################

fvfo.aov <- aov(`Fv/Fo`~light, fl)
fvfo.aov |> summary()
fvfo.letters <- HSD.test(fvfo.aov, "light", alpha=0.05)$groups

fvfo <- fl |> 
  group_by(light) |>
  summarise(
    mean=mean(`Fv/Fo`, na.rm=TRUE),
    n=n(),
    sd=sd(`Fv/Fo`, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(fvfo.letters |> rownames_to_column("light"), by="light")

ggplot(fvfo, aes(x=light, y=mean))+
  geom_col(fill = "grey", width = 0.8, colour = "black")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.3, linewidth = 1) +
  geom_text(aes(y = mean + se + 0.3, label = groups), size = 15, family = "Times", fontface="bold") +
  labs(x = "Photoperiod (h/d)",
       y = "Fv/Fo",
       title="初始荧光/最大荧光") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(family="Times", size = 25, face = "bold"),
        axis.text = element_text(family="Times", size = 25),
        plot.title = element_text(family="wryh", face="bold", size=35, hjust=0.5))
ggsave(filename="初始荧光_最大荧光Fv_Fo.png", width=8, height=8, units="cm", dpi=300, device="png")

























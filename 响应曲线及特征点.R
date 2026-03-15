library(ggplot2)
library(readxl)
library(plantecophys)
library(plantecowrap)
library(photosynthesis)
library(FitAQ)
library(purrr)
library(ggpubr)
library(ggsignif)
library(ggsci)
library(dplyr)
library(multcompView)### 与下面的agricolae是两种字母显著性表示的方式
library(agricolae)### 见上文
library(tibble)
library(showtext)
library(patchwork)
library(car)
library(scales)

setwd("D:/ZheJiang A&F University/杜鹃文章/第一次收获/licor-6800数据/第4次")

##### 数据 #####
mydaco2 <- read_excel("all.xlsx","co2")
mydalight <- read_excel("all.xlsx","light")

##### 调用字体库 #####
font_add("Times", regular = "times.ttf", bold = "timesbd.ttf", italic = "timesi.ttf", bolditalic = "timesbi.ttf")
font_add("wryh", regular="msyh.ttc", bold="msyhbd.ttc")
font_add(family = "arial", regular = "arial.ttf", bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")
showtext_auto()

########################### 光响应 ###############################
# 1) 读入与筛选
dalight <- read_excel("all.xlsx", sheet = "light") %>%
  mutate(replication = as.factor(replication)) %>%
  filter(treatment == "400_12")

# 2) 按重复分别拟合（NRH 模型；A_net = A，PPFD = Qin）
#    若 Qin 为入射光，可将 usealpha_Q = TRUE（alpha_Q=0.84 默认）
fits <- dalight %>%
  split(.$replication) %>%
  map(~ fit_photosynthesis(
    .data     = .x,
    .photo_fun= "aq_response",
    .vars     = list(.A = A, .Q = Qin), # A_net = "A", PPFD = "Qin"
    .method   = "ls",
    usealpha_Q = FALSE,  # 入射光→吸收光: TRUE；这里按你的字段直接拟合
    quiet     = TRUE
  ))

# 3) 提取参数表
params <- fits %>%
  imap_dfr(~ {
    b <- coef(.x)
    tibble(replication = .y,
           k_sat   = b["k_sat"],
           phi_J   = b["phi_J"],
           theta_J = b["theta_J"],
           Rd      = b["Rd"])
  })

# 4) 生成平滑拟合曲线用于绘图
Qmax <- max(dalight$Qin, na.rm = TRUE)
preds <- fits %>%
  imap_dfr(~ {
    b <- coef(.x)
    grid <- tibble(Qin = seq(0, Qmax, length.out = 300))
    # 若 usealpha_Q = TRUE，改为 Q_abs = 0.84 * Qin
    grid %>%
      mutate(Ahat = marshall_biscoe_1980(
        Q_abs = Qin,                # 若用吸收光则改为 0.84*Qin
        k_sat = b["k_sat"],
        phi_J = b["phi_J"],
        theta_J = b["theta_J"]
      ) - b["Rd"],
      replication = .y)
  })

# 5) 叠加绘图（6条曲线一张图，按重复分色）
ggplot() +
  geom_point(data = dalight, aes(Qin, A, color = replication), alpha = 1) +
  geom_line(data = preds, aes(Qin, Ahat, color = replication), linewidth = 1.5, alpha = 0.5) +
  labs(title = "Light-response curves (A12)",
       x = expression(paste("PPFD (", mu, "mol", "·", m^{-2}, "·", s^{-1}, ")")),
       y = expression(paste(A[net]," (", mu, "mol·", m^{-2}, ")"))) +
  annotate("text", x = Inf, y = -Inf,
           label = "Lcp[mean] == 15.457", parse = TRUE,
           hjust = 1, vjust = -2.5, size = 4) +
  annotate("text", x = Inf, y = -Inf,
           label = "Lsp[mean] == 334.114", parse = TRUE,
           hjust = 1, vjust = -1, size = 4) +
  theme_classic()


###光补偿点 nls ###
FitLCP(FitAQ(data = subset.data.frame(mydalight, treatment == "400_12" & replication == "2ND"), A = A, Q = Qin, provide.model = T))
###光饱和点 nls ###
FitSat(FitAQ(data = subset.data.frame(mydalight, treatment == "400_12" & replication == "2ND"), A = A, Q = Qin, provide.model = T))



##### E14_6 #####
# 依赖
library(minpack.lm)

# 你的数据子集
d <- subset(mydalight, treatment == "800_14" & replication == "6TH")
d <- subset(d, is.finite(A) & is.finite(Qin))
d <- d[order(d$Qin), ]

# —— 启动值（关键）——
Rd0    <- max(0, -min(d$A[d$Qin <= 50], na.rm = TRUE))  # 用近0光的观测估 Rd
Asat0  <- max(d$A, na.rm = TRUE)                  
Amax0  <- Asat0 + Rd0          
alpha0 <- try(coef(lm(A ~ 0 + Qin, data=d, subset = Qin <= 150))[["Qin"]], silent=TRUE)
if (inherits(alpha0, "try-error") || !is.finite(alpha0) || alpha0 <= 0)
  alpha0 <- (Amax0 + Rd0) / max(d$Qin, na.rm=TRUE)
theta0 <- 0.7

# —— 非矩形双曲线：优先模型 —— 
fit <- try(nlsLM(
  A ~ ((alpha*Qin + Amax - sqrt((alpha*Qin + Amax)^2 - 4*theta*alpha*Qin*Amax)) / (2*theta)) - Rd,
  data = d,
  start = list(alpha = alpha0, Amax = Amax0, theta = theta0, Rd = Rd0),
  lower = c(0,   0,    0.3, 0),
  upper = c(1.0, Inf,  1.0, Amax0 + Rd0),
  control = nls.lm.control(maxiter = 200, ftol = 1e-10, ptol = 1e-10)
), silent = TRUE)

# —— 若仍失败，退回矩形双曲线 —— 
if (inherits(fit, "try-error")) {
  fit <- nlsLM(
    A ~ (alpha*Qin)/(1 + (alpha*Qin/Amax)) - Rd,
    data = d,
    start = list(alpha = alpha0, Amax = Amax0, Rd = Rd0),
    lower = c(0, 0, 0),
    control = nls.lm.control(maxiter = 200, ftol = 1e-10, ptol = 1e-10)
  )
}

summary(fit)
co <- coef(fit)
Asat <- unname(co["Amax"] - co["Rd"])

# —— 求“达到 90% Asat 的光强 Q90”，带自适应括区，避免 f.lower 报错 —— 
Afun <- function(Q) {
  if ("theta" %in% names(co)) {
    ((co["alpha"]*Q + co["Amax"] - sqrt((co["alpha"]*Q + co["Amax"])^2 -
                                          4*co["theta"]*co["alpha"]*Q*co["Amax"]))/(2*co["theta"])) - co["Rd"]
  } else {
    (co["alpha"]*Q)/(1 + (co["alpha"]*Q/co["Amax"])) - co["Rd"]
  }
}
target <- 0.9 * Asat
upper <- max(d$Qin, na.rm=TRUE) * 3
g <- function(Q) Afun(Q) - target
while (is.finite(upper) && upper < 1e4 && g(upper) < 0) upper <- upper * 1.5
Q90 <- uniroot(g, c(0, upper))$root

list(params = co, Asat = Asat, Q90 = Q90)






#### vcmax ####
da.光合指数 <- read_excel("all.xlsx", "光合指数")

capture.output(aov(Vcmax ~co2+light, da.光合指数) |> summary(), file = "Vcmax~co2+light.csv")

aov.Vcmax <- aov(Vcmax ~ treatment, data = da.光合指数)
aov.Vcmax
summary(aov.Vcmax)
tukey <- aov(Vcmax ~ treatment, data = da.光合指数) |> TukeyHSD()
tukey.Vcmax <- as.data.frame(tukey$treatment)

View(tukey.Vcmax)
write.csv(tukey.Vcmax,"Vcmax.csv")
# 生成字母标记 #
letters <- multcompLetters4(aov.Vcmax, tukey)
letters_df <- as.data.frame.list(letters$treatment) # 关键：取出$treatment
letters_df$treatment <- rownames(letters_df)
colnames(letters_df)[1] <- "letters"
## 生成画图数据 ##
plot_Vcmax <- da.光合指数 %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(Vcmax, na.rm = TRUE),
    n    = sum(!is.na(Vcmax)),
    sd   = sd(Vcmax, na.rm = TRUE),
    se   = sd / sqrt(n)
  )

plot_Vcmax <- left_join(plot_Vcmax, letters_df, by = "treatment")## 将字母加入数据画图

ggplot(plot_Vcmax, aes(x = treatment, y = mean)) +
  geom_col(fill = "grey40", width = 0.7)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.2, linewidth = 1) +
  geom_text(aes(y = mean + se + 1.5, label = letters), size = 5) +
  labs(x = "treatment", y = "Vcmax") +
  theme_bw() +
  theme(panel.grid = element_blank())

#### Jmax ####
capture.output(aov(Jmax ~co2+light, da.光合指数) |> summary(), file = "Jmax~co2+light.csv")

aov.Jmax <- aov(Jmax ~ treatment, data = da.光合指数)
aov.Jmax
summary(aov.Jmax)
tukey <- aov(Jmax ~ treatment, data = da.光合指数) |> TukeyHSD()
tukey.Jmax <- as.data.frame(tukey$treatment)

View(tukey.Jmax)
write.csv(tukey.Jmax,"Jmax.csv")
# 生成字母标记 #
letters <- multcompLetters4(aov.Jmax, tukey)
letters_df <- as.data.frame.list(letters$treatment) # 关键：取出$treatment
letters_df$treatment <- rownames(letters_df)
colnames(letters_df)[1] <- "letters"
## 生成画图数据 ##
plot_Jmax <- da.光合指数 %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(Jmax, na.rm = TRUE),
    n    = sum(!is.na(Jmax)),
    sd   = sd(Jmax, na.rm = TRUE),
    se   = sd / sqrt(n)
  )

plot_Jmax <- left_join(plot_Jmax, letters_df, by = "treatment")## 将字母加入数据画图

ggplot(plot_Jmax, aes(x = treatment, y = mean)) +
  geom_col(fill = "grey40", width = 0.7)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.2, linewidth = 1) +
  geom_text(aes(y = mean + se + 3.5, label = letters), size = 5) +
  labs(x = "treatment", y = "Jmax") +
  theme_bw() +
  theme(panel.grid = element_blank())

#### Lcp ####

capture.output(aov(Lcp ~co2+light, da.光合指数) |> summary(), file = "Lcp~co2+light.csv")

aov.Lcp <- aov(Lcp ~ treatment, data = da.光合指数)
aov.Lcp
summary(aov.Lcp)
tukey <- aov(Lcp ~ treatment, data = da.光合指数) |> TukeyHSD()
tukey.Lcp <- as.data.frame(tukey$treatment)

View(tukey.Lcp)
write.csv(tukey.Lcp,"Lcp.csv")
# 生成字母标记 #
letters <- multcompLetters4(aov.Lcp, tukey)
letters_df <- as.data.frame.list(letters$treatment) # 关键：取出$treatment
letters_df$treatment <- rownames(letters_df)
colnames(letters_df)[1] <- "letters"
## 生成画图数据 ##
plot_Lcp <- da.光合指数 %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(Lcp, na.rm = TRUE),
    n    = sum(!is.na(Lcp)),
    sd   = sd(Lcp, na.rm = TRUE),
    se   = sd / sqrt(n)
  )

plot_Lcp <- left_join(plot_Lcp, letters_df, by = "treatment")## 将字母加入数据画图

ggplot(plot_Lcp, aes(x = treatment, y = mean)) +
  geom_col(fill = "grey40", width = 0.7)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.2, linewidth = 1) +
  geom_text(aes(y = mean + se + 3.5, label = letters), size = 5) +
  labs(x = "treatment", y = "Lcp") +
  theme_bw() +
  theme(panel.grid = element_blank())
#### Lsp ####
capture.output(aov(Lsp ~co2+light, da.光合指数) |> summary(), file = "Lsp~co2+light.csv")

aov.Lsp <- aov(Lsp ~ treatment, data = da.光合指数)
aov.Lsp
summary(aov.Lsp)
tukey <- aov(Lsp ~ treatment, data = da.光合指数) |> TukeyHSD()
tukey.Lsp <- as.data.frame(tukey$treatment)

View(tukey.Lsp)
write.csv(tukey.Lsp,"Lsp.csv")
# 生成字母标记 #
letters <- multcompLetters4(aov.Lsp, tukey)
letters_df <- as.data.frame.list(letters$treatment) # 关键：取出$treatment
letters_df$treatment <- rownames(letters_df)
colnames(letters_df)[1] <- "letters"
## 生成画图数据 ##
plot_Lsp <- da.光合指数 %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(Lsp, na.rm = TRUE),
    n    = sum(!is.na(Lsp)),
    sd   = sd(Lsp, na.rm = TRUE),
    se   = sd / sqrt(n)
  )

plot_Lsp <- left_join(plot_Lsp, letters_df, by = "treatment")## 将字母加入数据画图

ggplot(plot_Lsp, aes(x = treatment, y = mean)) +
  geom_col(fill = "grey40", width = 0.7)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(), width = 0.2, linewidth = 1) +
  geom_text(aes(y = mean + se + 25, label = letters), size = 5) +
  labs(x = "treatment", y = "Lsp") +
  theme_bw() +
  theme(panel.grid = element_blank())









############################### 新的开始  ###############################

setwd("D:\\ZheJiang A&F University\\杜鹃文章\\第一次收获\\licor-6800数据\\第4次")
lightres.da <- read.csv("light response.csv", fileEncoding = "GBK")
co2res.da <- read.csv("co2 response.csv", fileEncoding = "GBK")

##### co2 response #####
################################# 12 h d^-1 ########################
co2_a12 <- subset(co2res.da, treatment == '400_12') |> 
fitaci(varnames=list(ALEAF='A', Tleaf='TleafCnd', Ci='Ci', PPFD='Qin'),
       fitmethod='default') 
aci_a12 <- ggplot(co2_a12$df, aes(x = Ci)) +
  # 绘制实测数据点
  geom_point(aes(y = Ameas,color = "Measured"), alpha = 1, size = 1.5) +
  # 绘制模型拟合曲线
  geom_line(aes(y = Amodel, color = "Model Fit"), linewidth = 1) +
  # 绘制Ac和Aj曲线 - 使用颜色映射
  geom_line(aes(y = Ac, color = "Ac (Rubisco-limited)"), linetype = "dashed", alpha = 1) +
  geom_line(aes(y = Aj, color = "Aj (Light-limited)"), linetype = "dashed", alpha = 1) +
 annotate("text", label = expression(bold("12"~h~d^{-1})), x=300, y=25, size=10, family="arial", parse=TRUE) +
  # 手动设置颜色
  scale_color_manual(
    name = "Curves",
    labels = c("Ac (Rubisco-limited)" = "Ac",
               "Aj (Light-limited)" = "Aj"
               ),
    values = c(
      "Measured" = "darkblue",
      "Model Fit" = "red",
      "Ac (Rubisco-limited)" = "darkgreen",
      "Aj (Light-limited)" = "orange")
  ) +
  # 添加主题和标签
  theme_classic() +
  labs(x = expression(bold(paste("Ci ","(μmol ",mol^{-1},")"))),
       y = expression(bold(paste("A (μmol ",m^{-2}," ",s^{-1},")")))) +
  # 自定义主题
  theme(
    text = element_text(family = "arial"),
    axis.title = element_text(size = 25,face = "bold"),
    axis.text = element_text(size = 20),
    legend.position = "none",
    legend.title = element_text(size = 20,face = "bold"),
    legend.text = element_text(size = 20)
    )
aci_a12
ggsave("a12 aci.tiff",device="tiff",width=7.5,height=5,units="cm",dpi=300)
ggsave("a12 aci.png",device="png",width=7.5,height=5,units="cm",dpi=300)

################################# 14 h d^-1 ########################

co2_a14 <- subset(co2res.da, treatment == "400_14") |>
  fitaci(varnames=list(ALEAF='A', Tleaf='TleafCnd', Ci='Ci', PPFD='Qin'),
         fitmethod='default')

aci_a14 <- ggplot(co2_a14$df, aes(x = Ci)) +
  # 绘制实测数据点
  geom_point(aes(y = Ameas,color = "Measured"), alpha = 1, size = 1.5) +
  # 绘制模型拟合曲线
  geom_line(aes(y = Amodel, color = "Model Fit"), linewidth = 1) +
  # 绘制Ac和Aj曲线 - 使用颜色映射
  geom_line(aes(y = Ac, color = "Ac (Rubisco-limited)"), linetype = "dashed", alpha = 1) +
  geom_line(aes(y = Aj, color = "Aj (Light-limited)"), linetype = "dashed", alpha = 1) +
  annotate("text", label = expression(bold("14"~h~d^{-1})), x=300, y=25, size=10, family="arial", parse=TRUE) +
  # 手动设置颜色
  scale_color_manual(
    name = "Curves",
    labels = c("Ac (Rubisco-limited)" = "Ac",
               "Aj (Light-limited)" = "Aj"
    ),
    values = c(
      "Measured" = "darkblue",
      "Model Fit" = "red",
      "Ac (Rubisco-limited)" = "darkgreen",
      "Aj (Light-limited)" = "orange")
  ) +
  scale_y_continuous(breaks = c(0,10,20,30),
                     limits = c(-3,30),
                     expand = 0) +
  # 添加主题和标签
  theme_classic() +
  labs(x = expression(bold(paste("Ci ","(μmol ",mol^{-1},")"))),
       y = expression(bold(paste("A (μmol ",m^{-2}," ",s^{-1},")")))) +
  # 自定义主题
  theme(
    text = element_text(family = "arial"),
    axis.title = element_text(size = 25,face = "bold"),
    axis.text = element_text(size = 20),
    legend.position = "none",
    legend.title = element_text(size = 20,face = "bold"),
    legend.text = element_text(size = 20))
aci_a14
ggsave("a14 aci.tiff",device="tiff",width=7.5,height=5,units="cm",dpi=300)
ggsave("a14 aci.png",device="png",width=7.5,height=5,units="cm",dpi=300)

#################################### 16 h d^-1 ################################

co2_a16 <- subset(co2res.da, treatment == "400_16") |>
  fitaci(varnames=list(ALEAF='A', Tleaf='TleafCnd', Ci='Ci', PPFD='Qin'),
         fitmethod='default')

aci_a16 <- ggplot(co2_a16$df, aes(x = Ci)) +
  # 绘制实测数据点
  geom_point(aes(y = Ameas,color = "Measured"), alpha = 1, size = 1.5) +
  # 绘制模型拟合曲线
  geom_line(aes(y = Amodel, color = "Model Fit"), linewidth = 1) +
  # 绘制Ac和Aj曲线 - 使用颜色映射
  geom_line(aes(y = Ac, color = "Ac (Rubisco-limited)"), linetype = "dashed", alpha = 1) +
  geom_line(aes(y = Aj, color = "Aj (Light-limited)"), linetype = "dashed", alpha = 1) +
 annotate("text", label = expression(bold("16"~h~d^{-1})), x=300, y=25, size=10, family="arial", parse=TRUE) +
  # 手动设置颜色
  scale_color_manual(
    name = "Curves",
    labels = c("Ac (Rubisco-limited)" = "Ac",
               "Aj (Light-limited)" = "Aj"),
    values = c(
      "Measured" = "darkblue",
      "Model Fit" = "red",
      "Ac (Rubisco-limited)" = "darkgreen",
      "Aj (Light-limited)" = "orange")
  ) +
  scale_y_continuous(limits = c(-2,30), breaks = c(0,10,20,30)) +
  # 添加主题和标签
  theme_classic() +
  labs(x = expression(bold(paste("Ci ","(μmol ",mol^{-1},")"))),
       y = expression(bold(paste("A (μmol ",m^{-2}," ",s^{-1},")")))) +
  theme(
    text = element_text(family = "arial"),
    axis.title = element_text(size = 25,face = "bold"),
    axis.text = element_text(size = 20),
    legend.position = "none",
    legend.title = element_text(size = 20,face = "bold"),
    legend.text = element_text(size = 20)
    )
aci_a16
ggsave("a16 aci.tiff",device="tiff",width=7.5,height=5,units="cm",dpi=300)
ggsave("a16 aci.png",device="png",width=7.5,height=5,units="cm",dpi=300)


### E12 ###

subset(co2res.da, treatment == "800_12") |>
  fitaci(varnames=list(ALEAF='A', Tleaf='TleafCnd', Ci='Ci', PPFD='Qin'),
         fitmethod='default') |> plot()

### E14 ###

subset(co2res.da, treatment == "800_14") |>
  fitaci(varnames=list(ALEAF='A', Tleaf='TleafCnd', Ci='Ci', PPFD='Qin'),
         fitmethod='default') |> plot()


################## light response ##############
ggplot(subset(lightres.da, co2 == 400), aes(x=Qin,y=A,colour=treatment))+
  geom_point()+
  stat_smooth(method="glm")

lightres_A <- subset(lightres.da, co2 == 400)

fits <- lightres_A %>%
  split(.$treatment) %>%
  map(~ fit_photosynthesis(
    .data     = .x,
    .photo_fun= "aq_response",
    .vars     = list(.A = A, .Q = Qin), # A_net = "A", PPFD = "Qin"
    .method   = "ls",
    usealpha_Q = FALSE,  # 入射光→吸收光: TRUE；这里按你的字段直接拟合
    quiet     = TRUE
  ))

# 3) 提取参数表
params <- fits %>%
  imap_dfr(~ {
    b <- coef(.x)
    tibble(replication = .y,
           k_sat   = b["k_sat"],
           phi_J   = b["phi_J"],
           theta_J = b["theta_J"],
           Rd      = b["Rd"])
  })

# 4) 生成平滑拟合曲线用于绘图
Qmax <- max(lightres_A$Qin, na.rm = TRUE)
preds <- fits %>%
  imap_dfr(~ {
    b <- coef(.x)
    grid <- tibble(Qin = seq(0, Qmax, length.out = 300))
    # 若 usealpha_Q = TRUE，改为 Q_abs = 0.84 * Qin
    grid %>%
      mutate(Ahat = marshall_biscoe_1980(
        Q_abs = Qin,                # 若用吸收光则改为 0.84*Qin
        k_sat = b["k_sat"],
        phi_J = b["phi_J"],
        theta_J = b["theta_J"]
      ) - b["Rd"],
      treatment = .y)
  })

# 5) 叠加绘图（3条曲线一张图，按重复分色）
p_lightres <- ggplot(data = lightres_A, aes(x = Qin, y = A, color = treatment)) +
  geom_point(alpha = 1, size = 1.2, show.legend = TRUE) +
  geom_line(data = preds, aes(Qin, Ahat, color = treatment), linewidth = 1, alpha = 0.7, show.legend = FALSE) +
  labs(x = expression(bold(paste("PPFD (μmol ", m^{-2}, s^{-1},")"))),
       y = expression(bold(paste("A (μmol ", m^{-2}, ` `, s^{-1}, ")"))),
       colour = "Photoperiod") +
  scale_y_continuous(limits = c(-2.5,8),
                     breaks = c(-2,0,2,4,6,8),
                     expand = 0
                     ) +
  scale_x_continuous(limits = c(0,1601),
                     breaks = c(0,400,800,1200,1600)
                     ) +
  scale_color_manual(labels = c("400_12" = expression(12~h~d^{-1}),
                            "400_14" = expression(14~h~d^{-1}),
                            "400_16" = expression(16~h~d^{-1})
                            ),
                     values = c("400_12" = "#FFEDA0",
                              "400_14" = "#FD8D3C",
                              "400_16" = "#BD0026")
                   ) +
  theme_classic() +
  theme(
    text = element_text(family = "arial", colour = "black"),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    legend.position = "inside",
    legend.background = element_blank(),
    legend.position.inside = c(0.6,0.4),
    legend.key.height = unit(3, "pt"),
    legend.key.width  = unit(5, "pt"),
    legend.key.justification = c(0.5,0.5)
        )
p_lightres
ggsave(filename="light res A.tiff", device="tiff", dpi=300, width=6, height=5, units="cm")
ggsave(filename="light res A.png", device="png", dpi=300, width=6, height=5, units="cm")

###########################

################################# 修改后Vcmax ######################

shuju <- read_excel("all.xlsx", "co2指数（改）")

### 正态检验 ###
tapply(X = shuju$Vcmax, INDEX = shuju$light, FUN = "shapiro.test")
### 方差齐性检验 ###
leveneTest(y = Vcmax~light, data = shuju)
### 方差分析 无差异 ###
vc.aov <- aov(Vcmax~light, shuju)
summary(vc.aov)
### 事后多重比较 无差异 ###
vc.letters <- HSD.test(vc.aov, "light", alpha=0.05)$groups
vc.letters <- LSD.test(vc.aov, "light", alpha=0.05)$groups

vc <- shuju |> 
  group_by(light) |>
  summarise(
    mean=mean(Vcmax, na.rm=TRUE),
    n=n(),
    sd=sd(Vcmax, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(vc.letters |> rownames_to_column("light"), by="light")
write.csv(vc,"修改后Vcmax.csv")

# ggplot(vc, aes(x=light, y=mean))+
#   geom_col(fill = "grey", width = 0.8, colour = "black")+
#   geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
#                 position = position_dodge(), width = 0.3, linewidth = 1) +
#   geom_text(aes(y = mean + se + 0.4, label = groups), size = 15, family = "Times", fontface="bold") +
#   labs(x = "Photoperiod (h/d)",
#        y = expression(bold(paste(A[net], " (", mu, "mol ", m^{-2}, s^{-1}, ")"))),
#        title="净光合速率") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.title = element_text(family="Times", size = 25, face = "bold"),
#         axis.text = element_text(family="Times", size = 25),
#         plot.title = element_text(family="wryh", face="bold", size=35, hjust=0.5))
# ggsave(filename="净光合速率Anet.png", width=8, height=8, units="cm", dpi=300, device="png")


################################# 修改后Jmax #################################
### 正态检验 ###
tapply(X = shuju$Jmax, INDEX = shuju$light, FUN = "shapiro.test")
### 方差齐性检验 ###
leveneTest(y = Jmax~light, data = shuju)
### anova ###
jm.aov <- aov(Jmax~light, shuju)
summary(jm.aov)

jm.letters <- HSD.test(jm.aov, "light", alpha=0.05)$groups

jm <- shuju |> 
  group_by(light) |>
  summarise(
    mean=mean(Jmax, na.rm=TRUE),
    n=n(),
    sd=sd(Jmax, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(jm.letters |> rownames_to_column("light"), by="light")
jm
write.csv(jm,"修改后jmax.csv")

###################### Vcmax/Jmax ##############

### 正态检验 ###
tapply(X = shuju$`Vc/J`, INDEX = shuju$light, FUN = "shapiro.test")
### 方差齐性检验 ###
leveneTest(y = `Vc/J`~light, data = shuju)
### anova ###
VcJ.aov <- aov(`Vc/J`~light, shuju)
summary(VcJ.aov)

VcJ.letters <- HSD.test(VcJ.aov, "light", alpha=0.05)$groups

VcJ <- shuju |> 
  group_by(light) |>
  summarise(
    mean=mean(`Vc/J`, na.rm=TRUE),
    n=n(),
    sd=sd(`Vc/J`, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(VcJ.letters |> rownames_to_column("light"), by="light")
VcJ

write.csv(VcJ,"修改后Vc_J.csv")

###################### 修改后Lcp #########################
shuju <- read_excel("all.xlsx", "光指数（改）")
### 正态检验 ###
tapply(X = shuju$Lcp, INDEX = shuju$light, FUN = "shapiro.test")
### 方差齐性检验 ###
leveneTest(y = Lcp ~ light, data = shuju)
### ANOVA ###
lcp.aov <- aov(Lcp~light, shuju)
summary(lcp.aov)
lcp.letters <- HSD.test(lcp.aov, "light", alpha=0.05)$groups

lcp <- shuju |> 
  group_by(light) |>
  summarise(
    mean=mean(Lcp, na.rm=TRUE),
    n=n(),
    sd=sd(Lcp, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(lcp.letters |> rownames_to_column("light"), by="light")

lcp
write.csv(lcp,"修改后lcp.csv")

################################ 修改后Lsp ###################################
shuju <- read_excel("all.xlsx", "光指数（改）")
### 正态检验 ###
tapply(X = shuju$Lsp, INDEX = shuju$light, FUN = "shapiro.test")
### 方差齐性检验 ###
leveneTest(y = Lsp ~ light, data = shuju)
### ANOVA ###
lsp.aov <- aov(Lsp~light, shuju)
summary(lsp.aov)
lsp.letters <- HSD.test(lsp.aov, "light", alpha=0.05)$groups

lsp <- shuju |> 
  group_by(light) |>
  summarise(
    mean=mean(Lsp, na.rm=TRUE),
    n=n(),
    sd=sd(Lsp, na.rm=TRUE),
    se=sd/sqrt(n)
  ) |> left_join(lsp.letters |> rownames_to_column("light"), by="light")

lsp
write.csv(lsp,"修改后lsp.csv")






#############################################################################

############ 12h/d Amax from curve #######
### A12-1 id=4 ###
a12_1 <- subset.data.frame(mydalight, treatment == "400_12" & replication == "4TH")
fit <- fit_photosynthesis(
 .data = a12_1,
 .photo_fun = "aq_response",
 .vars = list(.A = A, .Q = "Qin")
)
coef(fit)[["k_sat"]]

### A12-2 id=5 ###
a12_2 <- subset.data.frame(mydalight, treatment == "400_12" & replication == "5TH")
fit <- fit_photosynthesis(
 .data = a12_2,
 .photo_fun = "aq_response",
 .vars = list(.A = A, .Q = "Qin")
)
coef(fit)[["k_sat"]]

### A12-3 id=6 ###
a12_3 <- subset.data.frame(mydalight, treatment == "400_12" & replication == "6TH")
fit <- fit_photosynthesis(
 .data = a12_3,
 .photo_fun = "aq_response",
 .vars = list(.A = A, .Q = "Qin")
)
coef(fit)[["k_sat"]]

############ 14h/d Amax from curve #######
### A14-1 id=1 ###
a14_1 <- subset.data.frame(mydalight, treatment == "400_14" & replication == "1ST")
fit <- fit_photosynthesis(
 .data = a14_1,
 .photo_fun = "aq_response",
 .vars = list(.A = A, .Q = "Qin")
)
coef(fit)[["k_sat"]]

### A14-2 id=4 ###
a14_2 <- subset.data.frame(mydalight, treatment == "400_14" & replication == "4TH")
fit <- fit_photosynthesis(
 .data = a14_2,
 .photo_fun = "aq_response",
 .vars = list(.A = A, .Q = "Qin")
)
coef(fit)[["k_sat"]]

### A14-3 id=6 ###
a14_3 <- subset.data.frame(mydalight, treatment == "400_14" & replication == "6TH")
fit <- fit_photosynthesis(
 .data = a14_3,
 .photo_fun = "aq_response",
 .vars = list(.A = A, .Q = "Qin")
)
coef(fit)[["k_sat"]]

############ 16h/d Amax from curve #######
### A16-1 id=4 ###
a16_1 <- subset.data.frame(mydalight, treatment == "400_16" & replication == "4TH")
fit <- fit_photosynthesis(
 .data = a16_1,
 .photo_fun = "aq_response",
 .vars = list(.A = A, .Q = "Qin")
)
coef(fit)[["k_sat"]]

### A16-2 id=5 ###
a16_2 <- subset.data.frame(mydalight, treatment == "400_16" & replication == "5TH")
fit <- fit_photosynthesis(
 .data = a16_2,
 .photo_fun = "aq_response",
 .vars = list(.A = A, .Q = "Qin")
)
coef(fit)[["k_sat"]]

### A16-3 id=6 ###
a16_3 <- subset.data.frame(mydalight, treatment == "400_16" & replication == "6TH")
fit <- fit_photosynthesis(
 .data = a16_3,
 .photo_fun = "aq_response",
 .vars = list(.A = A, .Q = "Qin")
)
coef(fit)[["k_sat"]]



############

################################### Rd ###################################

## 几种Rd的计算方法 详情看help文件 ##
  # fit_r_light_kok  Kok(1956)的方法 #
  # fit_r_light_WalkerOrt  Walker和Ort方法(Laisk) #
  # fit_r_light_yin  Yin方法(2011) 该方法只适用于改进的Rd计算 #

## 必要的参数：A、PPFD、Ci、PhiPS2(Yin方法) ##

# fit_photosynthesis(
#   .data,
#   .photo_fun,   .photo_fun = "r_light"
#   .model = "default",
#   .vars = NULL,
#   .method = "ls",
#   ...,
#   quiet = FALSE,
#   brm_options = NULL
# )

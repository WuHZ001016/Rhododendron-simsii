# 安装一次即可
# install.packages(c("readr","stringr","dplyr","ggplot2","scales"))

library(readr); library(stringr); library(dplyr); library(ggplot2); library(scales)

file <- "E:\\mpea\\收获\\A12-1.scs"  # <- 修改为你的 .scs 文件路径

# 1) 读取原始行
lines <- readr::read_lines(file)

# 2) 找到包含列名的那一行（m-PEA 的总表头，后半段是时间点）
hdr_idx <- which(str_detect(lines, "^Acq\\. No,Time,Date,Trace Colour"))
stopifnot(length(hdr_idx) >= 1)
header_tokens <- str_split(lines[hdr_idx[1]], ",", simplify = FALSE)[[1]]

# 3) 时间列从“第一个看起来像数字”的 token 开始
is_num <- function(x) grepl("^\\s*[-+]?\\d*\\.?\\d+(e[-+]?\\d+)?\\s*$", x, ignore.case = TRUE)
time_start <- which(sapply(header_tokens, is_num))[1]
stopifnot(!is.na(time_start))

# 4) 抽取 PF（Prompt Fluorescence）这一行的数据
## pf_idx <- which(str_detect(lines, ",PF\\s*,")) ##
library(stringr)

# 帮你看一眼文件里有哪些候选的“ID”标签（第6列），便于确认：
peek_ids <- function(lines, n = 10) {
  hits <- sapply(seq_along(lines), function(i){
    toks <- str_split(lines[i], ",", simplify = TRUE)
    if (ncol(toks) >= 6) trimws(toks[1,6]) else NA_character_
  })
  unique(na.omit(hits))[1:n]
}

ids_preview <- peek_ids(lines, n = 20)
print(ids_preview)  # 这里你会看到诸如 "01-01--PF", "01-01--Ch2", "01-01--DF000" 等

# 更粗鲁地找到 PF 行：检查第6列（ID）是否包含 "PF"（忽略大小写与空格），且不包含 "DF"
is_pf_line <- function(line) {
  toks <- str_split(line, ",", simplify = TRUE)
  if (ncol(toks) < 6) return(FALSE)
  id <- trimws(toks[1,6])
  grepl("PF", id, ignore.case = TRUE) && !grepl("\\bDF\\b", id, ignore.case = TRUE)
}

pf_hits <- which(vapply(lines, is_pf_line, logical(1)))
stopifnot(length(pf_hits) >= 1)   # 如果这里还报错，ids_preview 的输出发我一下，我们对照实际标签再匹配

pf_idx <- pf_hits[1]  # 取第一个 PF 行

stopifnot(length(pf_idx) >= 1)
pf_tokens <- str_split(lines[pf_idx[1]], ",", simplify = FALSE)[[1]]

# 5) 取出时间与荧光值
time_vec <- suppressWarnings(as.numeric(header_tokens[time_start:length(header_tokens)]))
F_vec    <- suppressWarnings(as.numeric(pf_tokens[time_start:length(header_tokens)]))

# 如果存在 DF（delayed fluorescence）行，官方注记要求 ÷100 才是正确幅度；
# 本段是 PF，不需要做 ÷100。若你日后处理 DF，可在取值后 F_vec <- F_vec/100。

# 6) 计算 Fo, Fm, Fv/Fm，并定位 J、I 典型时间点
Fo <- F_vec[1]
Fm <- max(F_vec, na.rm = TRUE)
FvFm <- (Fm - Fo) / Fm

# 典型时间点（秒）：O≈20 µs, J≈2 ms, I≈30 ms；P 取 Fm 所在时刻
closest_idx <- function(target, x) which.min(abs(x - target))
idx_O <- 1
idx_J <- closest_idx(0.002, time_vec)
idx_I <- closest_idx(0.03,  time_vec)
idx_P <- which.max(F_vec)

marks <- tibble(
  label = c("O","J (~2 ms)","I (~30 ms)","P (=Fm)"),
  t     = c(time_vec[idx_O], time_vec[idx_J], time_vec[idx_I], time_vec[idx_P]),
  F     = c(F_vec[idx_O],    F_vec[idx_J],   F_vec[idx_I],   F_vec[idx_P])
)

# 7) 绘图（对数时间轴）
df <- tibble(time = time_vec, F = F_vec)

p <- ggplot(df, aes(time, F)) +
  geom_line() +
  scale_x_log10(
    breaks = c(2e-5, 2e-4, 2e-3, 2e-2, 2e-1, 1),
    labels = label_number()
  ) +
  labs(
    x = "Time (s, log scale)",
    y = "Fluorescence (a.u.)",
    title = "OJIP Transient (PF) from m-PEA .scs"
  ) +
  theme_minimal(base_size = 13)

# 标注 O/J/I/P 点与水平/竖直参考线
p <- p +
  geom_point(data = marks, aes(t, F)) +
  geom_text(data = marks, aes(t, F, label = label), vjust = -0.7, size = 4) +
  geom_vline(xintercept = c(time_vec[idx_O], time_vec[idx_J], time_vec[idx_I], time_vec[idx_P]),
             linetype = "dashed") +
  annotate("label", x = max(time_vec, na.rm = TRUE), y = Inf, hjust = 1, vjust = 1,
           label = sprintf("Fo = %.0f\nFm = %.0f\nFv/Fm = %.3f", Fo, Fm, FvFm))

print(p)

# 8) 如需导出图片：
# ggsave("OJIP_A12-1.png", p, width = 7, height = 5, dpi = 300)

# ===== 可选：把结果保存为 CSV，便于后续统计 =====
# write_csv(df, "OJIP_A12-1.csv")

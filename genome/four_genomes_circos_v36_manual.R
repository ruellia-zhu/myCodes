suppressPackageStartupMessages({
  library(scales)
})

root <- "~/work/circos"
sizedir <- file.path(root, "01_sizes")
densdir <- file.path(root, "04_density")
outdir  <- file.path(root, "07_plot")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# V3.6 重新设计版
# 说明：
# 1. 不再使用共享 sector 的 circlize 布局来承载所有物种。
# 2. 每个物种单独按自身染色体长度分配弧长，因此同一编号染色体在不同物种中长度可以不同。
# 3. chr14 与 chr01 之间保留大 gap，系统树插入该 gap。
# =========================================================

# =========================
# 0. 可调参数
# =========================
FONT_FAMILY <- "Times"   # 常用值: "Times", "serif"

# 角度系统：0 度在正上方，顺时针增加
GAP_CENTER_DEG <- 315     # 大 gap 中心位置。常用 300~330，左上一般取 310~320
GAP_BIG_DEG    <- 48      # 大 gap 角度。常用 40/45/50/55
GAP_SMALL_DEG  <- 1.2     # 染色体间小 gap。常用 0.8/1.0/1.2/1.5

# 轨道平滑。1 = 不平滑。
GENE_SMOOTH_K   <- 1      # 常用 1/2/3
REPEAT_SMOOTH_K <- 1      # 常用 1/2/3

# 极端值截断，防止极端窗口把峰顶拉坏
GENE_CAP_QUANTILE   <- 0.97  # 常用 0.95/0.97/0.99
REPEAT_CAP_QUANTILE <- 0.98  # 常用 0.98/0.99

# 半径布局
R_OUTER       <- 0.92    # 最外圈半径
SPECIES_GAP_R <- 0.012   # 物种层间距
GROUP_THICK   <- 0.150   # 每个物种大层总厚度

# 每个物种大层内部结构
CHR_BAND_THICK   <- 0.028   # 染色体色带厚度
CHR_AXIS_OFFSET  <- 0.010   # 长度标尺相对染色体带外缘偏移
GENE_BASE_OFFSET <- 0.044   # gene 轨道下边界相对 species_outer 的偏移
GENE_TRACK_THICK <- 0.026   # gene 峰最大厚度
REPEAT_BASE_OFFSET <- 0.079 # repeat 轨道下边界相对 species_outer 的偏移
REPEAT_TRACK_THICK <- 0.026 # repeat 峰最大厚度
BG_INNER_OFFSET <- 0.110    # 物种背景带内缘相对 species_outer 的偏移

# 标尺设置
AXIS_MAJOR_MB      <- 25     # 主刻度 Mb 间隔。常用 10/20/25/50
AXIS_MINOR_TICKS   <- 2      # 每主刻度间的小刻度数。常用 0/1/2/4
AXIS_LABEL_CEX     <- 0.40   # 外圈标尺文字大小。这里是 base 图 cex。
AXIS_SHOW_END_LABEL <- TRUE  # 是否补终点刻度
CHR_LABEL_CEX      <- 0.55   # 染色体名称字号

# 系统树位置，0~1 画布坐标
TREE_BOX <- c(0.11, 0.30, 0.67, 0.84)  # left, right, bottom, top

# 图例位置，0~1 画布坐标
TRACK_LEGEND_X <- 0.03
TRACK_LEGEND_Y <- 0.95
SPECIES_LEGEND_X <- 0.80
SPECIES_LEGEND_Y <- 0.95
LEGEND_TITLE_CEX <- 0.95
LEGEND_ITEM_CEX  <- 0.78

# 输出清晰度
PNG_WIDTH  <- 7000
PNG_HEIGHT <- 7000
PNG_RES    <- 600
PNG_TYPE   <- if (capabilities("cairo")) "cairo" else getOption("bitmapType")

# 标题
TITLE <- "Four genomes circular comparison"
TITLE_CEX <- 1.25
TITLE_Y <- 0.985

# 物种顺序：外 -> 内
species_info <- data.frame(
  key = c("PP", "CW", "SJO", "SJ"),
  latin = c("P. platycarpum", "C. wilsonii", "S. japonicum f. oligophyllum", "S. japonicum"),
  stringsAsFactors = FALSE
)

# 配色
cols <- list(
  gene_fill   = "#C47AB8",
  gene_line   = "#9A4E93",
  repeat_fill = "#8DC17B",
  repeat_line = "#5E9E4A",
  species_bg = c(
    "P. platycarpum" = "#F1DFC2",
    "C. wilsonii" = "#D8ECD7",
    "S. japonicum f. oligophyllum" = "#D9E4F5",
    "S. japonicum" = "#EFD7E8"
  ),
  chr_cols = c(
    "#8FB7AA", "#D6B27C", "#AEB9DA", "#C8A7C9",
    "#9DB8C8", "#D1C38E", "#AFC8B0", "#D8AFA6",
    "#88B0BC", "#C7B59E", "#B8C9E2", "#CDB8D8",
    "#B8D3A2", "#E0B8A6"
  ),
  tree_tip = c(
    "P. platycarpum" = "#9B6A2A",
    "C. wilsonii" = "#4D8B4D",
    "S. japonicum f. oligophyllum" = "#4B79BD",
    "S. japonicum" = "#9A4E93"
  )
)

# =========================
# 1. 读取数据
# =========================
read_sizes <- function(sp) {
  x <- read.table(file.path(sizedir, paste0(sp, ".chrom.sizes")), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(x) <- c("chr_raw", "len")
  x$chr <- sub("^[A-Z]+chr", "chr", x$chr_raw)
  x$chr <- factor(x$chr, levels = sprintf("chr%02d", 1:14))
  x <- x[order(x$chr), c("chr", "len")]
  rownames(x) <- NULL
  x
}

read_gene_density <- function(sp) {
  x <- read.table(file.path(densdir, paste0(sp, ".gene_density.500k.tsv")), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(x) <- c("chr_raw", "start", "end", "value")
  x$chr <- sub("^[A-Z]+chr", "chr", x$chr_raw)
  x[, c("chr", "start", "end", "value")]
}

read_repeat_density <- function(sp) {
  x <- read.table(file.path(densdir, paste0(sp, ".repeat_density_frac.500k.tsv")), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(x) <- c("chr_raw", "start", "end", "value")
  x$chr <- sub("^[A-Z]+chr", "chr", x$chr_raw)
  x[, c("chr", "start", "end", "value")]
}

size_list <- setNames(lapply(species_info$key, read_sizes), species_info$key)
gene_list <- setNames(lapply(species_info$key, read_gene_density), species_info$key)
rep_list  <- setNames(lapply(species_info$key, read_repeat_density), species_info$key)

all_gene_values <- unlist(lapply(gene_list, function(x) x$value), use.names = FALSE)
all_rep_values  <- unlist(lapply(rep_list,  function(x) x$value), use.names = FALSE)
gene_cap <- as.numeric(quantile(all_gene_values, GENE_CAP_QUANTILE, na.rm = TRUE))
rep_cap  <- as.numeric(quantile(all_rep_values,  REPEAT_CAP_QUANTILE, na.rm = TRUE))

# =========================
# 2. 基础工具函数
# =========================
polar_xy <- function(theta_deg, r) {
  th <- theta_deg * pi / 180
  cbind(x = r * sin(th), y = r * cos(th))
}

annulus_polygon <- function(theta1, theta2, r_in, r_out, n = 120) {
  th_outer <- seq(theta1, theta2, length.out = n)
  th_inner <- seq(theta2, theta1, length.out = n)
  outer <- polar_xy(th_outer, r_out)
  inner <- polar_xy(th_inner, r_in)
  rbind(outer, inner)
}

arc_line <- function(theta1, theta2, r, n = 150) {
  polar_xy(seq(theta1, theta2, length.out = n), r)
}

radial_line <- function(theta, r1, r2) {
  polar_xy(c(theta, theta), c(r1, r2))
}

smooth_numeric <- function(v, k = 1) {
  if (length(v) < k || k <= 1) return(v)
  out <- stats::filter(v, rep(1 / k, k), sides = 2)
  out <- as.numeric(out)
  out[is.na(out)] <- v[is.na(out)]
  out
}

fmt_mb <- function(bp_vec) {
  mb <- bp_vec / 1e6
  ifelse(abs(mb - round(mb)) < 1e-8,
         as.character(as.integer(round(mb))),
         formatC(mb, format = "f", digits = 1))
}

# =========================
# 3. 每个物种独立角度布局
# =========================
make_species_layout <- function(size_df) {
  chr_order <- as.character(size_df$chr)
  lens <- size_df$len
  names(lens) <- chr_order

  total_data_angle <- 360 - GAP_BIG_DEG - GAP_SMALL_DEG * (length(chr_order) - 1)
  angle_per_bp <- total_data_angle / sum(lens)

  gap_start <- GAP_CENTER_DEG - GAP_BIG_DEG / 2
  gap_end   <- GAP_CENTER_DEG + GAP_BIG_DEG / 2

  cur <- gap_end
  out <- vector("list", length(chr_order))

  for (i in seq_along(chr_order)) {
    chr <- chr_order[i]
    w <- lens[i] * angle_per_bp
    th1 <- cur
    th2 <- cur + w
    out[[i]] <- data.frame(chr = chr, len = lens[i], theta1 = th1, theta2 = th2, theta_mid = (th1 + th2) / 2)
    cur <- th2
    if (i < length(chr_order)) cur <- cur + GAP_SMALL_DEG
  }

  do.call(rbind, out)
}

layout_list <- setNames(lapply(size_list, make_species_layout), species_info$key)

# =========================
# 4. 物种半径布局
# =========================
make_radii <- function(n_species) {
  res <- vector("list", n_species)
  outer <- R_OUTER
  for (i in seq_len(n_species)) {
    species_outer <- outer
    species_inner <- outer - GROUP_THICK
    res[[i]] <- list(
      species_outer = species_outer,
      species_inner = species_inner,
      bg_outer = species_outer,
      bg_inner = species_outer - BG_INNER_OFFSET,
      chr_outer = species_outer,
      chr_inner = species_outer - CHR_BAND_THICK,
      gene_base = species_outer - GENE_BASE_OFFSET,
      gene_top  = species_outer - GENE_BASE_OFFSET + GENE_TRACK_THICK,
      rep_base  = species_outer - REPEAT_BASE_OFFSET,
      rep_top   = species_outer - REPEAT_BASE_OFFSET + REPEAT_TRACK_THICK
    )
    outer <- species_inner - SPECIES_GAP_R
  }
  names(res) <- species_info$key
  res
}

radii_list <- make_radii(nrow(species_info))

# =========================
# 5. 绘图部件
# =========================
open_canvas <- function(bg = "white") {
  par(mar = c(0, 0, 0, 0), family = FONT_FAMILY, xpd = NA, bg = bg)
  plot(c(-1.18, 1.18), c(-1.18, 1.18), type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1)
}

draw_annulus <- function(theta1, theta2, r_in, r_out, col, border = NA, lwd = 1, n = 120) {
  xy <- annulus_polygon(theta1, theta2, r_in, r_out, n = n)
  polygon(xy[,1], xy[,2], col = col, border = border, lwd = lwd)
}

draw_arc <- function(theta1, theta2, r, col = "black", lwd = 1, n = 150) {
  xy <- arc_line(theta1, theta2, r, n = n)
  lines(xy[,1], xy[,2], col = col, lwd = lwd)
}

draw_radial <- function(theta, r1, r2, col = "black", lwd = 1) {
  xy <- radial_line(theta, r1, r2)
  lines(xy[,1], xy[,2], col = col, lwd = lwd)
}

draw_species_background <- function(sp_key, sp_label) {
  lay <- layout_list[[sp_key]]
  rr  <- radii_list[[sp_key]]
  bg_col <- alpha(cols$species_bg[sp_label], 0.78)
  for (i in seq_len(nrow(lay))) {
    draw_annulus(lay$theta1[i], lay$theta2[i], rr$bg_inner, rr$bg_outer,
                 col = bg_col, border = alpha(cols$species_bg[sp_label], 1), lwd = 0.5)
  }
}

draw_species_chrband <- function(sp_key) {
  lay <- layout_list[[sp_key]]
  rr  <- radii_list[[sp_key]]
  for (i in seq_len(nrow(lay))) {
    draw_annulus(lay$theta1[i], lay$theta2[i], rr$chr_inner, rr$chr_outer,
                 col = cols$chr_cols[i], border = "white", lwd = 0.7)

    # 染色体标签
    xy_lab <- polar_xy(lay$theta_mid[i], rr$chr_outer + 0.055)
    text(xy_lab[1,1], xy_lab[1,2], labels = lay$chr[i], cex = CHR_LABEL_CEX,
         srt = -(lay$theta_mid[i]), family = FONT_FAMILY)

    # 标尺
    len_bp <- lay$len[i]
    major_step <- AXIS_MAJOR_MB * 1e6
    major_rel <- seq(0, floor(len_bp / major_step) * major_step, by = major_step)
    if (AXIS_SHOW_END_LABEL && tail(major_rel, 1) < len_bp) {
      major_rel <- c(major_rel, len_bp)
    }
    major_rel <- unique(round(major_rel))

    # 主刻度
    for (bp in major_rel) {
      frac <- bp / len_bp
      th <- lay$theta1[i] + frac * (lay$theta2[i] - lay$theta1[i])
      draw_radial(th, rr$chr_outer + CHR_AXIS_OFFSET, rr$chr_outer + CHR_AXIS_OFFSET + 0.010,
                  col = "grey25", lwd = 0.5)

      lab <- fmt_mb(bp)
      xy_ticklab <- polar_xy(th, rr$chr_outer + CHR_AXIS_OFFSET + 0.023)
      text(xy_ticklab[1,1], xy_ticklab[1,2], labels = lab,
           cex = AXIS_LABEL_CEX, family = FONT_FAMILY,
           srt = -(th), col = "grey20")
    }

    # 小刻度
    if (AXIS_MINOR_TICKS > 0 && length(major_rel) >= 2) {
      for (j in seq_len(length(major_rel) - 1)) {
        s <- major_rel[j]
        e <- major_rel[j + 1]
        mids <- seq(s, e, length.out = AXIS_MINOR_TICKS + 2)[-c(1, AXIS_MINOR_TICKS + 2)]
        for (bp in mids) {
          frac <- bp / len_bp
          th <- lay$theta1[i] + frac * (lay$theta2[i] - lay$theta1[i])
          draw_radial(th, rr$chr_outer + CHR_AXIS_OFFSET, rr$chr_outer + CHR_AXIS_OFFSET + 0.006,
                      col = "grey35", lwd = 0.4)
        }
      }
    }
  }
}

draw_density_track <- function(sp_key, df, cap_value, smooth_k, base_r, top_r, fill_col, line_col) {
  lay <- layout_list[[sp_key]]
  size_df <- size_list[[sp_key]]

  for (i in seq_len(nrow(lay))) {
    chr <- as.character(lay$chr[i])
    chr_len <- lay$len[i]
    d <- df[df$chr == chr, ]
    if (nrow(d) == 0) next
    d <- d[order(d$start, d$end), ]
    vals <- pmin(d$value, cap_value) / cap_value
    vals <- smooth_numeric(vals, k = smooth_k)

    # 先画基线
    draw_arc(lay$theta1[i], lay$theta2[i], base_r, col = alpha("grey45", 0.6), lwd = 0.5)

    # 每个窗口一个小 polygon，保留真实锯齿
    for (j in seq_len(nrow(d))) {
      frac1 <- d$start[j] / chr_len
      frac2 <- d$end[j] / chr_len
      th1 <- lay$theta1[i] + frac1 * (lay$theta2[i] - lay$theta1[i])
      th2 <- lay$theta1[i] + frac2 * (lay$theta2[i] - lay$theta1[i])
      r_top <- base_r + vals[j] * (top_r - base_r)
      poly_xy <- annulus_polygon(th1, th2, base_r, r_top, n = 10)
      polygon(poly_xy[,1], poly_xy[,2], col = alpha(fill_col, 0.74), border = NA)
    }

    # 峰顶折线
    mids <- (d$start + d$end) / 2
    thm  <- lay$theta1[i] + (mids / chr_len) * (lay$theta2[i] - lay$theta1[i])
    rtop <- base_r + vals * (top_r - base_r)
    xy_top <- polar_xy(thm, rtop)
    lines(xy_top[,1], xy_top[,2], col = line_col, lwd = 0.45)
  }
}

draw_all_rings <- function() {
  for (idx in seq_len(nrow(species_info))) {
    sp_key <- species_info$key[idx]
    sp_label <- species_info$latin[idx]
    draw_species_background(sp_key, sp_label)
    draw_species_chrband(sp_key)
    rr <- radii_list[[sp_key]]
    draw_density_track(sp_key, gene_list[[sp_key]], gene_cap, GENE_SMOOTH_K,
                       base_r = rr$gene_base, top_r = rr$gene_top,
                       fill_col = cols$gene_fill, line_col = cols$gene_line)
    draw_density_track(sp_key, rep_list[[sp_key]], rep_cap, REPEAT_SMOOTH_K,
                       base_r = rr$rep_base, top_r = rr$rep_top,
                       fill_col = cols$repeat_fill, line_col = cols$repeat_line)
  }
}

# =========================
# 6. 树与图例
# =========================
draw_tree_unit <- function() {
  left <- TREE_BOX[1]; right <- TREE_BOX[2]; bottom <- TREE_BOX[3]; top <- TREE_BOX[4]

  # 用户树的分支长度
  x_root <- 0
  x_pp   <- 0.0244259
  x_n1   <- 0.0244259
  x_cw   <- 0.0244259 + 0.0246292
  x_n2   <- 0.0244259 + 0.0300592
  x_sjo  <- x_n2 + 0.00399534
  x_sj   <- x_n2 + 0.00400791
  x_max  <- max(x_pp, x_cw, x_sjo, x_sj)
  sx <- function(x) left + (x / x_max) * (right - left) * 0.78

  # 纵向位置与物种层次一致，外->内
  y_pp  <- top - 0.03
  y_cw  <- top - 0.09
  y_sjo <- top - 0.15
  y_sj  <- top - 0.21
  y_n2 <- (y_sjo + y_sj) / 2
  y_n1 <- (y_cw + y_n2) / 2

  seg <- function(x1,y1,x2,y2, col="grey30", lwd=1.3) segments(x1,y1,x2,y2,col=col,lwd=lwd,xpd=NA)
  txt <- function(x,y,lab,col,cex=0.62) text(x,y,lab,col=col,cex=cex,font=3,adj=c(0,0.5),xpd=NA,family=FONT_FAMILY)

  seg(sx(x_root), y_pp, sx(x_root), y_n1)
  seg(sx(x_root), y_pp, sx(x_pp), y_pp)
  seg(sx(x_root), y_n1, sx(x_n1), y_n1)
  seg(sx(x_n1), y_cw, sx(x_n1), y_n2)
  seg(sx(x_n1), y_cw, sx(x_cw), y_cw)
  seg(sx(x_n1), y_n2, sx(x_n2), y_n2)
  seg(sx(x_n2), y_sjo, sx(x_n2), y_sj)
  seg(sx(x_n2), y_sjo, sx(x_sjo), y_sjo)
  seg(sx(x_n2), y_sj, sx(x_sj), y_sj)

  dx <- 0.010
  txt(sx(x_pp)+dx, y_pp,  "P. platycarpum", cols$tree_tip["P. platycarpum"], 0.60)
  txt(sx(x_cw)+dx, y_cw,  "C. wilsonii", cols$tree_tip["C. wilsonii"], 0.60)
  txt(sx(x_sjo)+dx, y_sjo, "S. japonicum f. oligophyllum", cols$tree_tip["S. japonicum f. oligophyllum"], 0.56)
  txt(sx(x_sj)+dx, y_sj,  "S. japonicum", cols$tree_tip["S. japonicum"], 0.60)
}

draw_legends_unit <- function() {
  # 左上 Tracks
  x0 <- TRACK_LEGEND_X; y0 <- TRACK_LEGEND_Y
  text(x0, y0, "Tracks", cex = LEGEND_TITLE_CEX, font = 2, adj = c(0,0.5), xpd=NA, family=FONT_FAMILY)
  rect(x0, y0-0.060, x0+0.018, y0-0.043, col = cols$gene_fill, border = cols$gene_line, xpd=NA)
  text(x0+0.026, y0-0.0515, "Gene density", cex=LEGEND_ITEM_CEX, adj=c(0,0.5), xpd=NA, family=FONT_FAMILY)
  rect(x0, y0-0.104, x0+0.018, y0-0.087, col = cols$repeat_fill, border = cols$repeat_line, xpd=NA)
  text(x0+0.026, y0-0.0955, "Repeat density", cex=LEGEND_ITEM_CEX, adj=c(0,0.5), xpd=NA, family=FONT_FAMILY)
  rect(x0, y0-0.148, x0+0.018, y0-0.131, col = cols$chr_cols[1], border = NA, xpd=NA)
  text(x0+0.026, y0-0.1395, "Chromosome band", cex=LEGEND_ITEM_CEX, adj=c(0,0.5), xpd=NA, family=FONT_FAMILY)

  # 右上 species，使用拉丁名缩写
  x1 <- SPECIES_LEGEND_X; y1 <- SPECIES_LEGEND_Y
  text(x1, y1, "Species groups", cex = LEGEND_TITLE_CEX, font = 2, adj = c(0,0.5), xpd=NA, family=FONT_FAMILY)
  ys <- c(y1-0.050, y1-0.086, y1-0.122, y1-0.158)
  labs <- species_info$latin
  for (i in seq_along(labs)) {
    rect(x1, ys[i]-0.009, x1+0.018, ys[i]+0.009,
         col = alpha(cols$species_bg[labs[i]], 0.92), border = "grey70", xpd=NA)
    text(x1+0.026, ys[i], labs[i], cex = 0.58, adj = c(0,0.5), xpd=NA, family=FONT_FAMILY)
  }
}

draw_title_unit <- function() {
  text(0.5, TITLE_Y, TITLE, cex = TITLE_CEX, font = 2, adj = c(0.5,0.5), xpd=NA, family=FONT_FAMILY)
}

# =========================
# 7. 输出设备
# =========================
open_pdf <- function(path, bg = "white") {
  pdf(path, width = 12, height = 12, bg = bg, family = FONT_FAMILY, useDingbats = FALSE, onefile = TRUE)
}

open_png <- function(path, bg = "white") {
  png(path, width = PNG_WIDTH, height = PNG_HEIGHT, res = PNG_RES, bg = bg,
      type = PNG_TYPE, antialias = "subpixel")
}

# =========================
# 8. 绘图包装函数
# =========================
draw_full <- function() {
  open_canvas(bg = "white")
  draw_all_rings()
  par(new = TRUE)
  plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i", family = FONT_FAMILY)
  draw_title_unit()
  draw_tree_unit()
  draw_legends_unit()
}

draw_circos_only <- function(bg = "white") {
  open_canvas(bg = bg)
  draw_all_rings()
}

draw_tree_only <- function(bg = "transparent") {
  par(mar = c(0,0,0,0), family = FONT_FAMILY, xpd = NA, bg = bg)
  plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  draw_tree_unit()
}

draw_legends_only <- function(bg = "transparent") {
  par(mar = c(0,0,0,0), family = FONT_FAMILY, xpd = NA, bg = bg)
  plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  draw_legends_unit()
}

draw_title_only <- function(bg = "transparent") {
  par(mar = c(0,0,0,0), family = FONT_FAMILY, xpd = NA, bg = bg)
  plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  draw_title_unit()
}

# =========================
# 9. 输出文件
# =========================
open_pdf(file.path(outdir, "four_genomes_circos_v36_full_manual.pdf"), bg = "white")
draw_full()
dev.off()

open_png(file.path(outdir, "four_genomes_circos_v36_full_manual.png"), bg = "white")
draw_full()
dev.off()

open_pdf(file.path(outdir, "four_genomes_circos_v36_circos_only_manual.pdf"), bg = "white")
draw_circos_only(bg = "white")
dev.off()

open_png(file.path(outdir, "four_genomes_circos_v36_circos_only_transparent_manual.png"), bg = "transparent")
draw_circos_only(bg = "transparent")
dev.off()

open_png(file.path(outdir, "four_genomes_circos_v36_tree_only_transparent_manual.png"), bg = "transparent")
draw_tree_only(bg = "transparent")
dev.off()

open_png(file.path(outdir, "four_genomes_circos_v36_legends_only_transparent_manual.png"), bg = "transparent")
draw_legends_only(bg = "transparent")
dev.off()

open_png(file.path(outdir, "four_genomes_circos_v36_title_only_transparent_manual.png"), bg = "transparent")
draw_title_only(bg = "transparent")
dev.off()

cat("Done.\n")
cat(file.path(outdir, "four_genomes_circos_v36_full_manual.pdf"), "\n")
cat(file.path(outdir, "four_genomes_circos_v36_full_manual.png"), "\n")
cat(file.path(outdir, "four_genomes_circos_v36_circos_only_manual.pdf"), "\n")
cat(file.path(outdir, "four_genomes_circos_v36_circos_only_transparent_manual.png"), "\n")
cat(file.path(outdir, "four_genomes_circos_v36_tree_only_transparent_manual.png"), "\n")
cat(file.path(outdir, "four_genomes_circos_v36_legends_only_transparent_manual.png"), "\n")
cat(file.path(outdir, "four_genomes_circos_v36_title_only_transparent_manual.png"), "\n")
cat("Gene cap (", GENE_CAP_QUANTILE * 100, "th percentile) = ", gene_cap, "\n", sep = "")
cat("Repeat cap (", REPEAT_CAP_QUANTILE * 100, "th percentile) = ", rep_cap, "\n", sep = "")
cat("PNG type = ", PNG_TYPE, ", resolution = ", PNG_RES, " dpi\n", sep = "")

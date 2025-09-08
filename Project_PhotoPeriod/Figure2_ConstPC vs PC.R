pacman::p_load(
  Ipaper, data.table, dplyr, lubridate, 
  ggplot2, gg.layers, ggrepel
)

# f = "./OUTPUT/PMLV2China_flux37_LAI_glass,ConstPC_gof.csv"
# d = fread(f)
indir = "./Project_PhotoPeriod/OUTPUT/"

sites = c("儋州", "哀牢山", "西双版纳-Rubber", "西双版纳-EBF", "鼎湖山")
# d[site %in% sites, ]

# dat[GLASS <= -2]
set_font()

plot_Figure1 <- function(VAR="GPP", type_lai = "whit", 
  x = "ConstPC", y = "WithPC", 
  lims = c(0.2, 1), thr = 0.02) {
  
  dat_x = fread(glue("{indir}/PMLV2China_flux37_LAI_{type_lai},{x}_gof.csv"))
  dat_y = fread(glue("{indir}/PMLV2China_flux37_LAI_{type_lai},{y}_gof.csv"))

  df = listk(x = dat_x, y = dat_y) %>% melt_list("type")
  dat = df %>%
    select(-n_valid) %>%
    melt(c("type", "site", "var")) %>%
    dcast(site + var + variable ~ type, value.var = "value") %>%
    .[variable %in% c("NSE", "KGE", "R2"), ]

  brks <- c(thr, Inf) %>% c(-rev(.), .)
  pdat = dat[var == VAR & x >= lims[1] & y >= lims[1]] %>%
    mutate(diff = y - x, diff_lev = cut(diff, brks))

  p <- ggplot(pdat, aes(x, y)) +
    geom_point() +
    geom_point(data = pdat[abs(diff) >= thr], aes(color = diff_lev), alpha = 0.8) +
    # geom_point(data = pdat[(WithPC - NonPC) <= -0.02,], color = "red", alpha = 0.8) +
    # geom_text_repel(data = pdat[abs(NonPC - WithPC) >= 0.02,],
    #   aes(label = site), hjust = 0, vjust = 0, color = "blue") +
    geom_text_repel(
      data = pdat[abs(diff) >= thr], aes(color = diff_lev, label = site), min.segment.length = 0.2, 
      hjust = 0, vjust = 0, alpha = 0.8) +
    # geom_text_repel(data = pdat[(WithPC - NonPC) <= -0.02,],
    #   aes(label = site), hjust = 0, vjust = 0, color = "red", alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "red", alpha = 0.8) +
    coord_cartesian(xlim = lims, ylim = lims, expand = FALSE) +
    scale_color_manual(values = c("red", "darkgreen")) +
    theme(legend.position = "none") +
    labs(x = x, y = y) +
    # lims(x = lims, y = lims) +
    facet_wrap(~variable, scales = "free") 
  prefix = toupper(type_lai)
  write_fig(p, glue("Figure3_{x}&{y}_{prefix}_{VAR}.pdf"), 10, 4.5)
}

## LAI
# type_lai = "whit"
type_lai = "glass"
# x = "ConstPC"
x = "NonPC"
y = "ConstPC"
# y = "WithPC"

plot_Figure1("GPP", type_lai, thr = 0.02, x = x, y = y)
plot_Figure1("ET", type_lai, lims = c(0.2, 0.9), thr = 0.02, x = x, y = y)

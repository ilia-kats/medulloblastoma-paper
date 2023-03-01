import matplotlib as mpl
import matplotlib.pyplot as plt
import plotnine
from plotnine import *

mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "Helvetica"
mpl.rcParams["font.size"] = mpl.rcParams["axes.labelsize"] = mpl.rcParams[
    "axes.titlesize"
] = mpl.rcParams["xtick.labelsize"] = mpl.rcParams["ytick.labelsize"] = mpl.rcParams[
    "legend.fontsize"
] = mpl.rcParams[
    "figure.titlesize"
] = 7

plotnine.options.base_family = "Helvetica"
th = theme_bw(base_size=7, base_family="Helvetica") + theme(
    line=element_line(size=0.5),
    rect=element_rect(size=0.5),
    panel_grid_minor=element_blank(),
    panel_border=element_line(color="black"),
    axis_ticks=element_line(color="black"),
    axis_ticks_minor = element_blank(),
    axis_text=element_text(color="black", size=7),
    strip_background=element_blank(),
    strip_text=element_text(color="black", size=7),
    legend_text=element_text(size=7),
    legend_key=element_blank(),
    plot_title=element_text(ha="center"),
    aspect_ratio=1,
)
theme_set(th)

lfs_sporadic_color_scale = scale_color_manual(
    values=("#1b9e77", "#e6ab02"), breaks=("LFS", "sporadic")
)

lfs_sporadic_fill_scale = scale_fill_manual(
    values=("#1b9e77", "#e6ab02"), breaks=("LFS", "sporadic")
)

human_sample_map = {
    "1.1": "LFS1",
    "1.2": "LFS2",
    "1.3": "LFS3",
    "3.4": "LFS8",
    "2.1": "LFS4",
    "2.2": "LFS5",
    "2.4": "LFS6",
    "3.1": "LFS7",
    "4.1": "S1",
    "4.3": "S2",
    "5.2": "S3",
    "5.3": "S4",
    "5.4": "S5"
}

pdx_sample_map = {
    "1.2": "PDX B",
    "1.3": "PDX D1",
    "3.4": "PDX E1",
    "1.4": "PDX E2",
    "2.4": "PDX D2",
    "3.3": "PDX E3",
    "3.2": "PDX E4",
    "2.1": "PDX C1",
    "2.2": "PDX C2",
    "2.3": "PDX C3",
    "1.1": "PDX A"
}

import pandas as pd
import scanpy as sc

from plotnine import *

from plot_settings import human_sample_map, lfs_sporadic_fill_scale

def make_qc_df(samples):
    df = []
    for sname, sample in samples.items():
        cdf = {"sample": sname, "spots": sample.n_obs}
        for stat, name in (("total_counts", "total counts"), ("n_genes_by_counts", "detected genes")):
            cdf[f"median {name}"] = sample.obs[stat].median()
            cdf[f"mean {name}"] = sample.obs[stat].mean()
            cdf[f"stddev {name}"] = sample.obs[stat].std()
        df.append(pd.DataFrame(cdf, index=[0]))
    return pd.concat(df, axis=0).sort_values("sample").reset_index(drop=True)

def plot_gene_violin(samples, gene, plot_samples=False):
    xvar = "sample" if plot_samples else "type"
    xlab = "patient #" if plot_samples else None
    p = (
        ggplot(
            pd.concat(
                [
                    sc.get.obs_df(adata, [gene]).assign(
                        sample=human_sample_map[sample],
                        type=("LFS" if (int(sample[0]) <= 3) else "sporadic"),
                    )
                    for sample, adata in samples.items() if gene in adata.var_names
                ],
                axis=0,
            )
            .reset_index(drop=True)
            .melt(
                id_vars=["sample", "type"],
                var_name="gene",
                value_name="norm_expression",
            ),
            aes(xvar, "norm_expression", fill="type"),
        )
        + geom_violin(
            scale="width",
            color="black",
            n=4096,
            width=0.7,
            draw_quantiles=[0.25, 0.5, 0.75],
        )
        + scale_y_continuous(limits=(0, None), expand=(0, 0, 0.04, 0))
        + lfs_sporadic_fill_scale
        + labs(x=xlab, y="norm. expression", fill="")
        + theme(axis_ticks_major_x=element_blank(), aspect_ratio="auto")
    )
    if plot_samples:
        p += scale_x_discrete(
            labels=lambda x: [k[-1] for k in x.keys()]
        )
    else:
        p += guides(fill=None)
    return p

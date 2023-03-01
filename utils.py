import pandas as pd

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

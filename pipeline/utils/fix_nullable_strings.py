import pandas as pd
from anndata import AnnData

def fix_nullable_strings(adata: AnnData) -> AnnData:
    def fix_index(idx):
        return pd.Index([str(x) for x in idx], dtype=object, name=idx.name)

    def fix_df(df: pd.DataFrame):
        df.index = fix_index(df.index)

        for col in df.columns:
            s = df[col]

            if isinstance(s.dtype, pd.StringDtype):
                df[col] = s.astype(object)

            elif isinstance(s.dtype, pd.CategoricalDtype):
                cats = pd.Index(list(s.cat.categories), dtype=object)
                df[col] = pd.Categorical.from_codes(
                    s.cat.codes,
                    categories=cats,
                    ordered=s.cat.ordered,
                )

    fix_df(adata.obs) #type: ignore
    fix_df(adata.var) #type: ignore

    return adata
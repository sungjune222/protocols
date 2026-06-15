import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import Optional, Dict, List, Union
from anndata import AnnData

def spatial_visualization(
    adata: AnnData, 
    mask_dict: Dict[str, np.ndarray], 
    color_dict: Dict[str, str],
    filter_dict: Optional[Dict[str, Union[str, List[str]]]] = None
):
    if set(mask_dict.keys()) != set(color_dict.keys()):
        raise ValueError("Keys of mask_dict and color_dict must match")
    
    filter_desc = "" 
    processed_mask_dict = mask_dict.copy()

    if filter_dict:
        filter_mask = pd.Series(True, index=adata.obs.index)
        
        for col, values in filter_dict.items():
            if col not in adata.obs.columns:
                raise ValueError(f"Column '{col}' not found in adata.obs")
            
            if isinstance(values, str):
                values = [values]
            filter_mask &= adata.obs[col].isin(values)
            
        adata = adata[filter_mask]
        filter_mask_np = filter_mask.to_numpy()
        for k, v in processed_mask_dict.items():
            processed_mask_dict[k] = np.asarray(v)[filter_mask_np]
        
        filter_desc = " | ".join([f"{k}: {v}" for k, v in filter_dict.items()])

    mask_matrix = np.vstack(list(processed_mask_dict.values())).astype(bool)
    overlap_counts = mask_matrix.sum(axis=0)
    
    if np.any(overlap_counts > 1):
        raise ValueError("Overlapping regions found in mask_dict")
    
    if np.sum(overlap_counts == 0) > 0:
        print(f"Warning: {np.sum(overlap_counts == 0)} cells are not covered by any mask in mask_dict")
    
    if len(overlap_counts) != adata.n_obs:
        raise ValueError("Length of masks in mask_dict must match number of observations in adata")

    _, ax = plt.subplots(figsize=(10, 10)) 
    
    assert isinstance(adata.obs, pd.DataFrame)
    spatial_coords = adata.obs[['x_axis', 'y_axis']].to_numpy()

    for label, mask in processed_mask_dict.items():
        mask = mask.astype(bool)
        if mask.sum() > 0:
            ax.scatter(
                spatial_coords[mask, 0], 
                spatial_coords[mask, 1], 
                color=color_dict[label], alpha=1, s=3, edgecolors="none", 
                label=f"{label} (n={mask.sum()})" 
            )
    ax.set_aspect("equal")

    mask_names = ", ".join(list(mask_dict.keys()))
    main_title = f"Spatial Expression of {mask_names}"
    if filter_desc:
        ax.set_title(f"{main_title}\n(Filtered by: {filter_desc})", fontsize=14)
    else:
        ax.set_title(main_title, fontsize=16)
    
    ax.legend(
        loc="center left", bbox_to_anchor=(1.02, 0.5), 
        markerscale=4, fontsize=10, title="Mask Groups"
    )
    
    plt.tight_layout()
    plt.show()
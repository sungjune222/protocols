import os
import pickle
import pyranges as pr
from pipeline.utils.env import find_env_dir

def get_ensg_to_symbol(force_update=False):
    root_dir = find_env_dir("ROOT_DIR")
    cache_dir = os.path.join(root_dir, "references", "processed")
    cache_path = os.path.join(cache_dir, "ensg_to_symbol.pkl")
    
    gtf_path = os.path.join(root_dir, "references", "raw", "Homo_sapiens.GRCh38.115.gtf")

    if os.path.exists(cache_path) and not force_update:
        with open(cache_path, "rb") as f:
            return pickle.load(f)

    print("Parsing GTF to create gene mapping (this may take a while)...")
    
    gtf = pr.read_gtf(gtf_path, as_df=True)
    genes = gtf[gtf["Feature"] == "gene"].copy()
    genes["gene_id_clean"] = genes["gene_id"].astype(str).str.split(".", n=1).str[0]
    mapping = dict(zip(genes["gene_id_clean"], genes["gene_name"]))

    os.makedirs(cache_dir, exist_ok=True) 
    with open(cache_path, "wb") as f:
        pickle.dump(mapping, f)
    
    return mapping

def get_ensmusg_to_symbol(force_update=False):
    root_dir = find_env_dir("ROOT_DIR")
    cache_dir = os.path.join(root_dir, "references", "processed")
    cache_path = os.path.join(cache_dir, "ensmusg_to_symbol.pkl")
    
    gtf_path = os.path.join(root_dir, "references", "raw", "Mus_musculus.GRCm39.115.gtf")

    if os.path.exists(cache_path) and not force_update:
        with open(cache_path, "rb") as f:
            return pickle.load(f)

    print("Parsing GTF to create mouse gene mapping (this may take a while)...")
    
    gtf = pr.read_gtf(gtf_path, as_df=True)
    genes = gtf[gtf["Feature"] == "gene"].copy()
    genes["gene_id_clean"] = genes["gene_id"].astype(str).str.split(".", n=1).str[0]
    mapping = dict(zip(genes["gene_id_clean"], genes["gene_name"]))

    os.makedirs(cache_dir, exist_ok=True) 
    with open(cache_path, "wb") as f:
        pickle.dump(mapping, f)
    
    return mapping
import os
import re
from pipeline.utils.env import find_env_dir

def extract_mt_genes(species, mt_seqnames=["NC_005089.1", "MT", "chrM"]):
    root_dir = find_env_dir("ROOT_DIR")
    if species.lower() == "mouse":
        gtf_path = os.path.join(root_dir, "references", "raw", "GCF_000001635.27_GRCm39_genomic.gtf")
    elif species.lower() == "human":
        gtf_path = os.path.join(root_dir, "references", "raw", "GCF_000001405.40_GRCh38.p14_genomic.gtf")
    else:
        raise ValueError(f"Unsupported species: {species}. Supported species are 'mouse' and 'human'.")
    
    mt_genes = set()
    nuclear_genes = set()
    
    with open(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            parts = line.strip().split('\t')

            if len(parts) < 9 or parts[2] != 'gene':
                continue
                
            seqname = parts[0]
            attributes = parts[8]
            
            match = re.search(r'(?:gene|gene_name)\s+"([^"]+)"', attributes)
            if match:
                gene_name = match.group(1)
                if seqname in mt_seqnames:
                    mt_genes.add(gene_name)
                else:
                    nuclear_genes.add(gene_name)

    ambiguous_genes = mt_genes.intersection(nuclear_genes)
    
    if ambiguous_genes:
        error_msg = (
            f"Ambiguous genes detected! {len(ambiguous_genes)} gene(s) found in "
            f"both MT and nuclear genomes: {', '.join(ambiguous_genes)}"
        )
        raise ValueError(error_msg)

    return list(mt_genes)
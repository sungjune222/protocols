# Batch sizes
# This settings are based on 64GB RAM and RTX 5070ti GPU with 16GB VRAM
DOUBLET_REMOVAL_VAE_BATCH_SIZE = 1000
DOUBLET_REMOVAL_SOLO_BATCH_SIZE = 4000
SINGLE_CELL_VAE_BATCH_SIZE = 5000
# You could increase this value if you have more RAM (not VRAM), default value is 1024
DE_ANALYSIS_BATCH_SIZE = 100

# Number of CPU cores you want to use for parallel processing
CPU_CORE_COUNT = 16

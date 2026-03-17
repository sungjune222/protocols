import scvi
import torch
from pipeline.utils.env import find_env

CPU_CORE_COUNT = int(find_env("N_THREADS"))
torch.set_float32_matmul_precision("medium")
scvi.settings.num_threads = 8

# If you have enough RAM, you can set pin_memory to True for faster data transfer to GPU
DataLoader = dict(
    num_workers=8,
    pin_memory=True,
)

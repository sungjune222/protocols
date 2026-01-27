import scvi
import torch

from pipeline.config.constants import CPU_CORE_COUNT

# **If you are using NVIDIA GPU** this will adjust `float32_matmul_precision` to 'high' (TF32, default 'highest which is FP32),
# which accelerates matrix multiplication.
torch.set_float32_matmul_precision("high")
scvi.settings.num_threads = CPU_CORE_COUNT

# If you have enough RAM, you can set pin_memory to True for faster data transfer to GPU
DataLoader = dict(
    pin_memory=False,
)

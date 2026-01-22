import scvi
import torch

# **If you are using NVIDIA GPU** this will adjust `float32_matmul_precision` to 'high' (TF32, default 'highest which is FP32),
# which accelerates matrix multiplication.
torch.set_float32_matmul_precision("high")
scvi.settings.num_threads = 12

DataLoader = dict(
    pin_memory=True,
)

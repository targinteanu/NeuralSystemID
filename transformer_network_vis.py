import math
import numpy as np
from scipy.signal import firwin, filtfilt
import torch
import matplotlib.pyplot as plt
from myPytorchModels import TimeSeriesConv as myNetMdl
from csv2numpy import prepTimeSeqData
from torchinfo import summary
from torchviz import make_dot

# set params -------------------------------------------------------------------------------------
hzn = .2 # EVALUATION sample time, s
groupsize=15
numgroups=5
numgroupsunpaired=2
#fc = np.array([4,10,27,60,90]) # freq band center freqs
fc = np.array([4,10,27]) # freq band center freqs
netfile = "neural_network_pytorch_093523772db8074d4eac32acee3d4d2943fb7193_3.pth"
dt_target = 0.01 # model sample time, s
seq_len = 128 # model transformer samples
hzn_len = math.ceil(hzn / dt_target)  # horizon as multiple of MODEL Ts, NOT data Ts 
filtorder = 201

# define mdl struct ====================================================================
model = myNetMdl(dim_in=165, dim_out=165, time_len=seq_len, group_size=groupsize, num_groups=numgroups, tuple_size=3, numGrpUnpaired=numgroupsunpaired)
if netfile:
    model.load_state_dict(torch.load(netfile, map_location=torch.device('cpu'))) # load model weights
else:
    # bias the model frequency prediction to the center of each band 
    fcenter = torch.tensor(fc, dtype=torch.float32)
    fbias = fcenter.repeat_interleave(groupsize)
    fbias = fbias*dt_target*2*math.pi # scale by model sample time and 2pi to convert to radians
    fbias = torch.atanh(fbias / (math.pi + 1e-5)) # apply inverse tanh to get bias in pre-tanh space
    with torch.no_grad():
        model.fcoFreq.bias.copy_(fbias)

# Print all parameter information
total_params = 0

for name, param in model.named_parameters():
    num_params = param.numel()
    total_params += num_params

    print(f"{name}")
    print(f"  shape: {tuple(param.shape)}")
    print(f"  requires_grad: {param.requires_grad}")
    print(f"  num params: {num_params}")
    print(f"  mean: {param.data.mean():.6f}")
    print(f"  std:  {param.data.std():.6f}")
    print()

print(f"Total parameters: {total_params:,}")

x = torch.randn(1, seq_len, 165) # dummy input to visualize model structure
u = torch.randn(1, seq_len, 1)
y = model(x, u)

summary(model, input_data=[x, u])

#dot = make_dot(y[0,0,0], params=dict(model.named_parameters()))
#dot.render("transformer_network_vis", format="png")
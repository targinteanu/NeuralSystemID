import torch
import onnx
import onnxscript
from myPytorchModels import TimeSeriesConv as myNetMdl
from csv2numpy import prepTimeSeqData

netfile = "neural_network_pytorch_f246691b075e327a9a373ae3f43229bf0b6ad320(1)" # .pth -> .onnx
numtap = 128
groupsize=15
numgroups=5
numgroupsunpaired=2
dt = 0.01 # model sample time, s
hzn = 1 # samples

fs, feature_names, _, Xs, Ys, _, _, _, YsRaw, _, _, Us, _ = prepTimeSeqData(
    seq_len=numtap, hzn_len=hzn, dt_target=dt, drawFromStart=False, maxNumel=1e8)

Xs = torch.tensor(Xs, dtype=torch.float32)
Ys = torch.tensor(Ys, dtype=torch.float32)
Us = torch.tensor(Us, dtype=torch.float32)

model = myNetMdl(dim_in=Xs.shape[-1], dim_out=Ys.shape[-1], dim_u=Us.shape[-1], time_len=numtap, group_size=groupsize, num_groups=numgroups, tuple_size=3, numGrpUnpaired=numgroupsunpaired)
model.load_state_dict(torch.load(netfile+".pth"))
model.eval()

torch.onnx.export(
    model, 
    (Xs,Us), 
    netfile+".onnx", 
    input_names=["X", "U"], 
    output_names=["Y"], 
    opset_version=18, 
    do_constant_folding=True
)
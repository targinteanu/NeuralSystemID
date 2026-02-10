import math
import torch
import torch.nn as nn
import torch.nn.functional as F


class PairedDualNet(nn.Module):
    """
    One timepoint input, one timepoint output predictor without time sequence, 
    convolution, or transformer.
    Input features are assumed to be in pairs (i.e. real and imaginary), 
    groups (i.e. frequency band), and threads (i.e. channels).
    """

    def __init__(self, input_dim=101, output_dim=100, group_size=7, num_groups=7):
        super().__init__()

        C1 = 32
        C2 = 32

        self.num_groups = num_groups
        self.group_size = group_size
        self.num_pairs = num_groups * group_size  
        self.used_for_pairing = self.num_pairs * 2
        self.leftover_dim = input_dim - self.used_for_pairing
        self.pair_output = 8

        # Stage 1: pairwise linear 
        self.pair_fc1 = nn.Linear(2, C1)  
        self.pair_fc2 = nn.Linear(C1, self.pair_output)  

        # Stage 2: linear over flattened group
        self.group_fcA = nn.Linear(group_size * self.pair_output, C2)
        self.group_fcB = nn.Linear(num_groups * self.pair_output, C2)
        # MLP input = 56 + 56 + leftover(3) = 115
        mlp_in = (num_groups * C2) + (group_size * C2) + self.leftover_dim

        self.fc1 = nn.Linear(mlp_in, 256)
        self.fc2 = nn.Linear(256, 256)
        self.fc3 = nn.Linear(256, 256)
        self.fc4 = nn.Linear(256, 256)
        self.fc5 = nn.Linear(256, 256)
        self.fc6 = nn.Linear(256, 256)
        self.fc7 = nn.Linear(256, output_dim)

    def forward(self, x):
        B = x.size(0)

        x_left = x[:, self.used_for_pairing:]  # (B,n)
        x_used = x[:, :self.used_for_pairing]  # (B,2N)
        x_pairs = torch.stack(
            (
                x_used[:, :self.num_pairs],            # first half
                x_used[:, self.num_pairs:],            # second half
            ),
            dim=2
        )                                              # (B, N, 2)

        # Stage 1
        p = F.gelu(self.pair_fc1(x_pairs))     # (B,N,C1)
        p = F.gelu(self.pair_fc2(p))           # (B,N,pair_output)

        # Stage 2A: groups
        p_groups = p.view(B, self.num_groups, self.group_size * self.pair_output)  # (B,num_groups,...)
        a = F.gelu(self.group_fcA(p_groups))                      # (B,7,8)
        a_flat = a.view(B, -1)                                     # (B,56)

        # Stage 2B: threads
        p_threads = p.view(B, self.num_groups, self.group_size, self.pair_output)   # (B,7,7,2)
        p_threads = p_threads.permute(0,2,1,3).contiguous()           # (B,7,7,2)
        p_threads = p_threads.view(B, self.group_size, -1)            # (B,7,14)
        b = F.gelu(self.group_fcB(p_threads))                         # (B,7,8)
        b_flat = b.view(B, -1)                                        # (B,56)

        # MLP
        h = torch.cat([a_flat, b_flat, x_left], dim=1)        # (B,115)
        h = F.gelu(self.fc1(h))
        h = F.gelu(self.fc2(h))
        h = F.gelu(self.fc3(h))
        h = F.gelu(self.fc4(h))
        h = F.gelu(self.fc5(h))
        h = F.gelu(self.fc6(h))
        return self.fc7(h)


def elu_feature_map(x):
    # FAVOR+ feature map
    return F.elu(x) + 1


class LinearAttention(nn.Module):
    def __init__(self, dim, num_heads):
        super().__init__()
        assert dim % num_heads == 0
        self.num_heads = num_heads
        self.head_dim = dim // num_heads

        self.q_proj = nn.Linear(dim, dim)
        self.k_proj = nn.Linear(dim, dim)
        self.v_proj = nn.Linear(dim, dim)
        self.out_proj = nn.Linear(dim, dim)

    def forward(self, x):
        """
        x: (B, T, D)
        """
        B, T, D = x.shape

        # project to Q, K, V
        q = self.q_proj(x)
        k = self.k_proj(x)
        v = self.v_proj(x)

        # reshape to multi-head
        q = q.view(B, T, self.num_heads, self.head_dim)
        k = k.view(B, T, self.num_heads, self.head_dim)
        v = v.view(B, T, self.num_heads, self.head_dim)

        # apply kernel feature map
        q = elu_feature_map(q)
        k = elu_feature_map(k)

        # normalize to avoid numerical issues
        eps = 1e-6

        # Compute KV = sum_t (phi(K_t) * V_t)
        # (B, num_heads, head_dim, head_dim)
        kv = torch.einsum("bthd,bthm->bhmd", k, v)

        # Compute normalizer: z = 1 / (sum_t phi(K_t) * 1)
        # (B, num_heads, head_dim)
        z = 1 / (torch.einsum("bthd,bhd->bth", q, k.sum(dim=1)) + eps)

        # Compute output: y_t = (phi(Q_t) * KV) * z_t
        # (B, T, num_heads, head_dim)
        out = torch.einsum("bthd,bhmd->bthm", q, kv)
        out = out * z.unsqueeze(-1)

        # merge heads
        out = out.reshape(B, T, D)

        return self.out_proj(out)


class LinearTransformerEncoderLayer(nn.Module):
    """
    Meant to be a faster, less memory-intensive replacement for nn.TransformerEncoderLayer
    using linear attention. Not yet tested. 
    """

    def __init__(self, dim_model, num_heads, dim_ff=128, dropout=0.1):
        super().__init__()

        self.attn = LinearAttention(dim_model, num_heads)
        self.norm1 = nn.LayerNorm(dim_model)
        self.dropout1 = nn.Dropout(dropout)

        self.ff = nn.Sequential(
            nn.Linear(dim_model, dim_ff),
            nn.ReLU(),
            nn.Linear(dim_ff, dim_model),
        )
        self.norm2 = nn.LayerNorm(dim_model)
        self.dropout2 = nn.Dropout(dropout)

    def forward(self, x):
        # Attention block
        attn_out = self.attn(x)
        x = self.norm1(x + self.dropout1(attn_out))

        # Feedforward block
        ff_out = self.ff(x)
        x = self.norm2(x + self.dropout2(ff_out))

        return x


class SinusoidalPositionalEncoding(nn.Module):
    def __init__(self, dim_model, max_len=5000):
        super().__init__()

        # Create matrix of shape (max_len, dim_model)
        pe = torch.zeros(max_len, dim_model)
        position = torch.arange(0, max_len).unsqueeze(1)  # (max_len, 1)

        # Divide by log-based frequencies
        div_term = torch.exp(torch.arange(0, dim_model, 2) * (-math.log(10000.0) / dim_model))

        # Apply sin to even indices, cos to odd
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)

        # Register as buffer so it's saved with model but not trainable
        self.register_buffer('pe', pe.unsqueeze(0))  # (1, max_len, dim_model)

    def forward(self, x):
        """
        x: (B, T, dim_model)
        returns: x + positional_encoding[:, :T, :]
        """
        T = x.size(1)
        return x + self.pe[:, :T, :]


class TimeSeriesTransformer(nn.Module):
    """
    This is an encoder-only transformer that first acts on some pre-processed 
    input features which are assumed to be in pairs (i.e. real and imaginary), 
    groups (i.e. frequency band), and threads (i.e. channels). These features 
    are first processed together in a time-independent manner and then fed into 
    the transformer. A single transformer output is processed to predict the 
    final output. 
    """

    def __init__(self, dim_in, dim_out, dim_u=1, time_len=64, group_size=7, num_groups=7, tuple_size=2):
        super().__init__()

        # transformer properties
        dim_model=32
        len_model=time_len # i.e. "T"
        num_heads=4
        num_layers=8
        dim_ff=128

        # preprocessing features -------------------------------------------
        
        C1 = 32 # stage 1 hidden dim
        C2 = 32 # stage 2 hidden dim

        self.num_groups = num_groups
        self.group_size = group_size
        self.tuple_size = tuple_size
        self.num_pairs = num_groups * group_size  
        self.used_for_pairing = self.num_pairs * tuple_size
        self.leftover_dim = dim_in - self.used_for_pairing
        self.pair_output = 8

        # Stage 1: pairwise linear 
        self.pair_fc1 = nn.Linear(tuple_size, C1)  
        self.pair_fc2 = nn.Linear(C1, self.pair_output)  

        # Stage 2: linear over flattened group
        self.group_fcA = nn.Linear(group_size * self.pair_output, C2)
        self.group_fcB = nn.Linear(num_groups * self.pair_output, C2)
        mlp_in = (num_groups * C2) + (group_size * C2) + self.leftover_dim

        # stage 3B: feature-only processing to get to dim_model
        """
        self.fc1 = nn.Linear(mlp_in, 256)
        #self.fc2 = nn.Linear(256, 256)
        self.fc3 = nn.Linear(256, 128)
        self.fc4 = nn.Linear(128, 64)
        self.fc5 = nn.Linear(64, dim_model)
        """
        self.fc1 = nn.Linear(mlp_in, dim_model)

        # ---------------------------------------------------------------------------------------

        # Positional embedding
        self.pos_emb = SinusoidalPositionalEncoding(dim_model, max_len=len_model)

        # Transformer Encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=dim_model,
            nhead=num_heads,
            dim_feedforward=dim_ff,
            activation='gelu',
            batch_first=True  # lets inputs be (B, T, D)
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)

        # MLP latent dynamics 
        self.fc2 = nn.Linear(dim_model + dim_u, 512)
        self.fc3 = nn.Linear(512, 512)
        self.fc4 = nn.Linear(512, 512)
        self.fc5 = nn.Linear(512, dim_model)

        # Output head for next-step prediction ------------------------------------------------
        self.fco1 = nn.Linear(dim_model, 128)
        #self.fco2 = nn.Linear(64, 64)
        #self.fco3 = nn.Linear(64, 64)
        self.fco4 = nn.Linear(128, dim_out)

    def forward(self, x, u_seq):
        """
        x: (B, T, dim_in)
        """
        if x.ndim != 3:
            raise ValueError(f"Expected input ndim=3, got {x.ndim}")
        T = x.size(1)
        B = x.size(0)
        rollout = u_seq.size(1)

        x_left = x[:, :, self.used_for_pairing:]  # (B,T,n)
        x_used = x[:, :, :self.used_for_pairing]  # (B,T,kN)
        xy_skip = x_used[:,-1,:]  # for skip connection at output
        """
        x_pairs = torch.stack( 
            (
                x_used[:, :, :self.num_pairs],             # first half
                x_used[:, :, self.num_pairs:],             # second half
            ),
            dim=3
        )  # (B,T,N,2)
        """
        x_pairs = x_used.view(x_used.shape[0], x_used.shape[1], self.tuple_size, -1).permute(0,1,3,2).contiguous() # (B,T,N,k)

        # Stage 1
        p = F.gelu(self.pair_fc1(x_pairs))     # (B,T,N,C1)
        p = F.gelu(self.pair_fc2(p))           # (B,T,N,pair_output)
        # Stage 2A: groups
        p_groups = p.view(B, T, self.num_groups, self.group_size * self.pair_output)  # (B,T,num_groups,...)
        a = F.gelu(self.group_fcA(p_groups))                      
        a_flat = a.view(B, T, -1)                                 

        # Stage 2B: threads
        p_threads = p.view(B, T, self.num_groups, self.group_size, self.pair_output)   
        p_threads = p_threads.permute(0,1,3,2,4).contiguous()           
        p_threads = p_threads.view(B, T, self.group_size, -1)            # (B,T,group_size,...)
        b = F.gelu(self.group_fcB(p_threads))                         
        b_flat = b.view(B, T, -1)                                        
        h = torch.cat([a_flat, b_flat, x_left], dim=2)        # (B,T,mlp_in)

        # Stage 3 MLP
        h = F.gelu(self.fc1(h))
        #h = F.gelu(self.fc2(h))
        """
        h = F.gelu(self.fc3(h))
        h = F.gelu(self.fc4(h))
        h = F.gelu(self.fc5(h))
        """

        # Transformer --------------------------------------------------------------------

        h = self.pos_emb(h)
        h = self.encoder(h)

        # latent dynamics 
        z = h[:, -1, :]  # (B, dim_model)
        zskip = z
        for r in range(rollout):
            u = u_seq[:, r, :] # (B, dim_u)
            z = torch.cat([z, u], dim=1) # (B, dim_model+dim_u)
            z = F.gelu(self.fc2(z))
            z = F.gelu(self.fc3(z))
            z = F.gelu(self.fc4(z))
            z = F.gelu(self.fc5(z))
            z = z + zskip # skip connection
            zskip = z

        # Output head --------------------------------------------------------------------
        y = F.gelu(self.fco1(z))
        #y = F.gelu(self.fco2(y))
        #y = F.gelu(self.fco3(y))
        y = self.fco4(y)  # (B, dim_out)
        out = y + xy_skip # skip connection
        return out


class TimeSeriesConvTransformer(nn.Module):
    """
    Should work like TimeSeriesTransformer but with added time-convolutional
    layers before the transformer to capture local time patterns; therefore, a 
    longer time sequence (time_len) can be used, but the memory requirements of 
    the transformer (len_model) remain smaller. Could not get this to work well yet.
    TO DO: try combining stage 3A and 3B into a single conv processing stream without groups.
    """

    def __init__(self, dim_in, dim_out, time_len=256, group_size=7, num_groups=7, tuple_size=2):
        super().__init__()

        # transformer properties
        dim_model=32
        len_model=64 # i.e. "T"
        num_heads=4
        num_layers=8
        dim_ff=128

        # time conv properties
        # These should all be kept smaller than time_len
        K1 = 8
        K2 = 64
        #K3 = 64

        # preprocessing features -------------------------------------------
        
        C1 = 32
        C2 = 32

        self.num_groups = num_groups
        self.group_size = group_size
        self.tuple_size = tuple_size
        self.num_pairs = num_groups * group_size  
        self.used_for_pairing = self.num_pairs * tuple_size
        self.leftover_dim = dim_in - self.used_for_pairing
        self.pair_output = 8

        # Stage 1: pairwise linear 
        self.pair_fc1 = nn.Linear(tuple_size, C1) 
        self.pair_fc2 = nn.Linear(C1, self.pair_output) 

        # Stage 2: linear over flattened group
        self.group_fcA = nn.Linear(group_size * self.pair_output, C2)
        self.group_fcB = nn.Linear(num_groups * self.pair_output, C2)
        # MLP input = 56 + 56 + leftover(3) = 115
        mlp_in = (num_groups * C2) + (group_size * C2) + self.leftover_dim

        # stage 3A: time-only processing to get to len_model
        self.time_conv1 = nn.Conv1d(mlp_in, mlp_in, groups=mlp_in, kernel_size=K1)
        self.time_conv2 = nn.Conv1d(mlp_in, mlp_in, groups=mlp_in, kernel_size=K2)
        #self.time_conv3 = nn.Conv1d(mlp_in, mlp_in, groups=mlp_in, kernel_size=K3)
        #self.time_fc = nn.Linear(time_len - K1 - K2 - K3 + 3, len_model)
        self.time_fc = nn.Linear(time_len - K1 - K2 + 2, len_model)

        # stage 3B: feature-only processing to get to dim_model
        self.fc1 = nn.Linear(mlp_in, 256)
        #self.fc2 = nn.Linear(256, 256)
        self.fc3 = nn.Linear(256, 128)
        self.fc4 = nn.Linear(128, 64)
        self.fc5 = nn.Linear(64, dim_model)

        # ---------------------------------------------------------------------------------------

        # Positional embedding
        self.pos_emb = SinusoidalPositionalEncoding(dim_model, max_len=len_model)

        # Transformer Encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=dim_model,
            nhead=num_heads,
            dim_feedforward=dim_ff,
            activation='gelu',
            batch_first=True  # lets inputs be (B, T, D)
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)

        # Output head for next-step prediction ------------------------------------------------
        self.fco1 = nn.Linear(dim_model, 64)
        #self.fco2 = nn.Linear(64, 64)
        #self.fco3 = nn.Linear(64, 64)
        self.fco4 = nn.Linear(64, dim_out)

    def forward(self, x):
        """
        x: (B, T, dim_in)
        """
        if x.ndim != 3:
            raise ValueError(f"Expected input ndim=3, got {x.ndim}")
        T = x.size(1)
        #if T > self.pos_emb.size(1):
        #    raise RuntimeError(f"Sequence length T={T} exceeds time_len={self.pos_emb.size(1)}. "
        #                       "Either increase time_len or ensure inputs have smaller T.")
        #x = self.input_proj(x) + self.pos_emb[:, :T, :]
        B = x.size(0)

        x_left = x[:, :, self.used_for_pairing:]  # (B,T,n)
        x_used = x[:, :, :self.used_for_pairing]  # (B,T,kN)
        """
        x_pairs = torch.stack( 
            (
                x_used[:, :, :self.num_pairs],            # first half
                x_used[:, :, self.num_pairs:],            # second half
            ),
            dim=3
        )  # (B,T,N,2)
        """
        x_pairs = x_used.view(x_used.shape[0], x_used.shape[1], self.tuple_size, -1).permute(0,1,3,2).contiguous() # (B,T,N,k)

        # Stage 1
        p = F.gelu(self.pair_fc1(x_pairs))     # (B,T,N,C1)
        p = F.gelu(self.pair_fc2(p))           # (B,T,N,pair_output)

        # Stage 2A: groups
        p_groups = p.view(B, T, self.num_groups, self.group_size * self.pair_output)  # (B,T,num_groups,...)
        a = F.gelu(self.group_fcA(p_groups))                      
        a_flat = a.view(B, T, -1)                                 

        # Stage 2B: threads
        p_threads = p.view(B, T, self.num_groups, self.group_size, self.pair_output)   
        p_threads = p_threads.permute(0,1,3,2,4).contiguous()           
        p_threads = p_threads.view(B, T, self.group_size, -1)            # (B,T,group_size,...)
        b = F.gelu(self.group_fcB(p_threads))                        
        b_flat = b.view(B, T, -1)                                        
        h = torch.cat([a_flat, b_flat, x_left], dim=2)       # (B,T,mlp_in)

        # Stage 3A time processing
        h = h.permute(0,2,1) # (B, dim_model, T)
        h = F.gelu(self.time_conv1(h))
        h = F.gelu(self.time_conv2(h))
        #h = F.gelu(self.time_conv3(h))
        h = F.gelu(self.time_fc(h))
        h = h.permute(0,2,1) # (B, T, dim_model)

        # Stage 3B MLP
        h = F.gelu(self.fc1(h))
        #h = F.gelu(self.fc2(h))
        h = F.gelu(self.fc3(h))
        h = F.gelu(self.fc4(h))
        h = F.gelu(self.fc5(h))

        # Transformer --------------------------------------------------------------------

        h = self.pos_emb(h)
        z = self.encoder(h)

        # Output head --------------------------------------------------------------------
        y = z[:, -1, :]  # (B, dim_model)
        y = F.gelu(self.fco1(y))
        #y = F.gelu(self.fco2(y))
        #y = F.gelu(self.fco3(y))
        out = self.fco4(y)  # (B, dim_out)
        return out


# ----------------------------
# architecture for seq2seq: 
# encoder-decoder transformer for time series forecasting; not yet operational
# ----------------------------

# Helper: causal mask for decoder
def generate_square_subsequent_mask(sz: int, device: torch.device):
    """Upper-triangular mask with 0 on and below diagonal, -inf above to mask future tokens."""
    mask = torch.triu(torch.ones((sz, sz), device=device), diagonal=1).bool()
    # nn.Transformer modules expect float mask with -inf for masked positions if using add_mask
    # but when using boolean mask arguments (src_key_padding_mask / tgt_key_padding_mask) behavior differs.
    # We'll return a float mask suitable for use in `attn_mask` (additive).
    attn_mask = torch.full((sz, sz), float('-inf'), device=device)
    attn_mask[~mask] = 0.0
    return attn_mask  # shape (sz, sz)

# Seq2Seq Transformer Model
class Seq2SeqTimeSeriesTransformer(nn.Module):
    def __init__(
        self,
        dim_in,         # number of input features (all features present in encoder input)
        dim_out,        # number of features to predict (subset of dim_in)
        seq_len=64,
        horizon=32,
        dim_model=32,
        num_heads=2,
        num_encoder_layers=2,
        num_decoder_layers=2,
        dim_ff=64,
        dropout=0.1,
        #max_positional_len=5000,
        device=torch.device('cpu'),
    ):
        super().__init__()
        self.device = device
        self.seq_len = seq_len
        self.horizon = horizon
        self.dim_in = dim_in
        self.dim_out = dim_out
        self.dim_model = dim_model

        # 1) Input projection: map raw input features -> model dimension
        self.input_proj = nn.Linear(dim_in, dim_model)

        # 2) Positional embeddings (learned)
        #self.pos_emb_enc = nn.Parameter(torch.randn(1, max_positional_len, dim_model))
        #self.pos_emb_dec = nn.Parameter(torch.randn(1, max_positional_len, dim_model))
        # sinusoidal positional embeddings
        self.pos_emb_enc = SinusoidalPositionalEncoding(dim_model, max_len=seq_len)
        self.pos_emb_dec = SinusoidalPositionalEncoding(dim_model, max_len=horizon)

        # 3) Transformer encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=dim_model,
            nhead=num_heads,
            dim_feedforward=dim_ff,
            batch_first=True,
            dropout=dropout,
            activation='gelu'
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_encoder_layers)

        # 4) Transformer decoder
        decoder_layer = nn.TransformerDecoderLayer(
            d_model=dim_model,
            nhead=num_heads,
            dim_feedforward=dim_ff,
            batch_first=True,
            dropout=dropout,
            activation='gelu'
        )
        self.decoder = nn.TransformerDecoder(decoder_layer, num_layers=num_decoder_layers)

        # 5) Output projection: map model dimension -> predicted features
        self.output_proj = nn.Linear(dim_model, dim_out)

        # optional: small linear to initialize first decoder input from encoder summary, if needed
        self.dec_init_proj = nn.Linear(dim_model, dim_model)

        self._reset_parameters()
        self.to(device)

    def _reset_parameters(self):
        # Initialization similar to PyTorch Transformer
        for p in self.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p)

    def encode(self, src):
        """
        src: (batch, seq_len, dim_in)
        returns: memory (batch, seq_len, dim_model)
        """
        b, t, _ = src.shape
        x = self.input_proj(src)                       # (b, t, dim_model)
        #x = x + self.pos_emb_enc[:, :t, :]             # add positional embedding
        x = self.pos_emb_enc(x)                        # add positional embedding
        memory = self.encoder(x)                       # (b, t, dim_model)
        return memory

    def decode(self, tgt, memory, tgt_mask=None, tgt_key_padding_mask=None, memory_key_padding_mask=None):
        """
        tgt: (batch, tgt_len, dim_model)  -- already projected + pos emb
        memory: encoder outputs (batch, src_len, dim_model)
        """
        out = self.decoder(
            tgt=tgt,
            memory=memory,
            tgt_mask=tgt_mask,
            tgt_key_padding_mask=tgt_key_padding_mask,
            memory_key_padding_mask=memory_key_padding_mask,
        )  # (b, tgt_len, dim_model)
        return out

    def forward(self, src, tgt_in):
        """
        Forward pass used in training with teacher forcing.

        src: (batch, seq_len, dim_in)
        tgt_in: (batch, tgt_len, dim_out) -- ground truth future sequence shifted right
                Example: if horizon = H and dim_out=F,
                tgt_in[:, 0, :] should be a start token (e.g., zeros)
                tgt_in[:, 1:, :] are ground-truth up to H-1 steps (teacher forcing).
        Returns:
            logits: (batch, tgt_len, dim_out)
        """
        device = src.device
        b, src_len, _ = src.shape
        _, tgt_len, _ = tgt_in.shape
        #assert src_len <= self.pos_emb_enc.shape[1], "increase max_positional_len"
        #assert tgt_len <= self.pos_emb_dec.shape[1], "increase max_positional_len"

        # 1) Encode
        memory = self.encode(src)  # (b, src_len, dim_model)

        # 2) Prepare decoder input embeddings: we project the provided tgt_in (dim_out) into dim_model
        #    Common pattern: during training, feed the true previous outputs (teacher forcing).
        #    Because tgt_in contains actual feature values (not model embeddings), we map them
        #    into model space with a small linear layer (re-using input_proj for simplicity if dims match).
        #    Here, to keep model flexible, we'll use a small linear (input_proj_dec) implemented on the fly.
        # Simple approach: expand dim_out -> dim_model with a linear
        # For simplicity re-use input_proj if dim_in == dim_out; otherwise make a quick linear on the fly:
        if self.dim_in == self.dim_out:
            # reuse input projection if dims align
            tgt_emb = self.input_proj(tgt_in)            # (b, tgt_len, dim_model)
        else:
            # lightweight linear mapping for decoder inputs (not saved as parameter to keep API simple)
            # but better to declare a module if used heavily; here we allocate on the fly for clarity
            # (we'll do a proper parametrized projection to avoid re-creating params each forward)
            if not hasattr(self, 'tgt_input_proj'):
                self.tgt_input_proj = nn.Linear(self.dim_out, self.dim_model).to(self.device)
            tgt_emb = self.tgt_input_proj(tgt_in)

        # add positional embeddings
        #tgt_emb = tgt_emb + self.pos_emb_dec[:, :tgt_len, :]
        tgt_emb = self.pos_emb_dec(tgt_emb)

        # 3) Create causal mask for decoder self-attention
        tgt_mask = generate_square_subsequent_mask(tgt_len, device=device)  # (tgt_len, tgt_len)

        # 4) Run decoder: decoder expects target embeddings + encoder memory
        dec_out = self.decode(
            tgt=tgt_emb,
            memory=memory,
            tgt_mask=tgt_mask,
        )  # (b, tgt_len, dim_model)

        # 5) Project decoder outputs to predicted feature space
        logits = self.output_proj(dec_out)  # (b, tgt_len, dim_out)
        return logits

    @torch.no_grad()
    def generate(self, src, horizon=None, start_token=None):
        """
        Autoregressive generation / inference.

        src: (batch, seq_len, dim_in)
        horizon: number of steps to generate (defaults to self.horizon)
        start_token: (dim_out,) or (batch, dim_out) initial token fed to decoder as t=0.
                     If None, zeros are used.

        Returns:
            generated: (batch, horizon, dim_out)
        """
        self.eval()
        device = src.device
        b = src.shape[0]
        horizon = self.horizon if horizon is None else horizon

        # 1) encode
        memory = self.encode(src)  # (b, src_len, dim_model)

        # 2) prepare initial decoder input (t=0)
        if start_token is None:
            cur_input = torch.zeros((b, 1, self.dim_out), device=device)  # (b, 1, dim_out)
        else:
            # allow vector or batch vector
            st = torch.tensor(start_token, device=device, dtype=src.dtype)
            if st.dim() == 1:
                st = st.unsqueeze(0).expand(b, -1)
            cur_input = st.unsqueeze(1)  # (b, 1, dim_out)

        generated = []
        # autoregressive loop
        for t in range(horizon):
            # build tgt_in as the sequence of already generated tokens for this batch
            if len(generated) == 0:
                tgt_in = cur_input  # (b, 1, dim_out)
            else:
                # stack previously generated tokens
                prev = torch.cat(generated, dim=1)  # (b, t, dim_out)
                tgt_in = torch.cat([cur_input, prev], dim=1)  # (b, t+1, dim_out)

            # forward through model (teacher forcing not applied)
            logits = self.forward(src, tgt_in)  # (b, t+1, dim_out)
            # take the last step's predictions as next token
            next_token = logits[:, -1:, :]  # (b, 1, dim_out)
            generated.append(next_token)

            # note: optional sampling or temperature can be applied here for stochastic outputs

        gen = torch.cat(generated, dim=1)  # (b, horizon, dim_out)
        return gen
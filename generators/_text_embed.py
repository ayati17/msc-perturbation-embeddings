"""
Shared text -> vector step for text_struct and text_semantic. 
env: hf
"""
from __future__ import annotations
import os
import numpy as np
import torch
from tqdm import tqdm
from transformers import AutoModel, AutoTokenizer
from chemembed.registry import TextModel

BATCH_SIZE = 8


def embed_texts(model_spec: TextModel, texts: list[str]) -> np.ndarray:
    tokenizer, net, device = _load(model_spec)
    prepared = [model_spec.prefix + t for t in texts]

    out = []
    for i in tqdm(range(0, len(prepared), BATCH_SIZE), desc=model_spec.hf_id, leave=False):
        batch = prepared[i:i + BATCH_SIZE]
        encoded = tokenizer(batch, padding=True, truncation=True,
                            max_length=model_spec.max_length, return_tensors="pt")
        _warn_if_truncated(encoded, model_spec, i)
        encoded = {k: v.to(device) for k, v in encoded.items()}
        if model_spec.pooling == "last":
            mask = encoded["attention_mask"]
            encoded["position_ids"] = (mask.cumsum(-1) - 1).clamp(min=0)
        with torch.no_grad():
            hidden = net(**encoded).last_hidden_state
        pooled = _pool(hidden, encoded["attention_mask"], model_spec.pooling)
        out.append(pooled.float().cpu().numpy())

    return np.vstack(out)


def _load(model_spec: TextModel):
    device = "cuda" if torch.cuda.is_available() else "cpu"
    dtype = torch.float32
    if device == "cuda":
        dtype = torch.bfloat16 if torch.cuda.get_device_capability()[0] >= 8 else torch.float16

    token = os.environ.get("HF_TOKEN")          
    tokenizer = AutoTokenizer.from_pretrained(model_spec.hf_id, token=token)
    
    net = AutoModel.from_pretrained(model_spec.hf_id, dtype=dtype, token=token)
    net = net.to(device).eval()

    tokenizer.padding_side = model_spec.padding_side
    vocab_size = getattr(net.config, "vocab_size", None)

    pad = tokenizer.pad_token or tokenizer.eos_token
    if pad is None or (vocab_size and tokenizer.convert_tokens_to_ids(pad) >= vocab_size):
        vocab = tokenizer.get_vocab()
        pad = next((t for t in ("[PAD]", "<pad>", "<|pad|>") if t in vocab), None)
        pad = pad or tokenizer.convert_ids_to_tokens(0)
        print(f"[embed] no registered pad token; using {pad!r} "
              f"(id {tokenizer.convert_tokens_to_ids(pad)})")
    tokenizer.pad_token = pad
    if vocab_size and (tokenizer.pad_token_id or 0) >= vocab_size:
        fallback = tokenizer.unk_token or tokenizer.convert_ids_to_tokens(0)
        tokenizer.pad_token = fallback
        print(f"[embed] pad id outside vocab ({vocab_size}); using {fallback!r}")

    print(f"[embed] {model_spec.hf_id} on {device} ({dtype}), "
          f"pooling={model_spec.pooling}, padding={model_spec.padding_side}")
    return tokenizer, net, device


def _pool(hidden: torch.Tensor, mask: torch.Tensor, how: str) -> torch.Tensor:
    if how == "mean":
        m = mask.unsqueeze(-1).to(hidden.dtype)
        return (hidden * m).sum(1) / m.sum(1).clamp(min=1e-9)
    if how == "last":
        if not (mask[:, -1] == 1).all():
            raise ValueError("last-token pooling needs left padding; check padding_side")
        return hidden[:, -1]
    if how == "cls":
        return hidden[:, 0]
    raise ValueError(f"unknown pooling {how!r}")


def _warn_if_truncated(encoded, model_spec: TextModel, offset: int) -> None:
    lengths = encoded["attention_mask"].sum(1)
    if (lengths >= model_spec.max_length).any():
        n = int((lengths >= model_spec.max_length).sum())
        print(f"[embed] WARNING: {n} text(s) hit max_length={model_spec.max_length} "
              f"in batch starting at {offset} - content was truncated")
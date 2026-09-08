"""
Runs each spec as a subprocess in its own environment.
compounds/all.csv in, artifact dir out.
"""

from __future__ import annotations
import os
import subprocess
import sys
from chemembed.artifacts import artifact_dir, exists
from chemembed.config import ENVS, ROOT
from chemembed.registry import Spec


def env_python(env: str) -> str:
    candidate = ENVS / env / "bin" / "python"
    if candidate.exists():
        return str(candidate)
    print(f"[warn] no env. at {candidate}, falling back to current interpreter")
    return sys.executable


def generate(specs: list[Spec | str], *, force: bool = False) -> None:
    specs = [s if isinstance(s, Spec) else Spec.parse(s) for s in specs]
    for spec in specs:
        if exists(spec) and not force:
            print(f"[skip] {spec} -> {artifact_dir(spec)}")
            continue
        cmd = [env_python(spec.env), str(ROOT / spec.script), "--spec", str(spec)]
        print(f"[run ] {spec}  ({spec.env})")
        env = os.environ | {"PYTHONPATH": str(ROOT)}
        result = subprocess.run(cmd, env = env)
        if result.returncode != 0:
            print(f"[FAIL] {spec} exited {result.returncode}")
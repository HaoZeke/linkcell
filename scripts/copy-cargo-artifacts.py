#!/usr/bin/env python3
"""Build linkcell with Cargo and copy its C ABI artifacts for Meson."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cargo", required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--profile", choices=("debug", "release"), required=True)
    parser.add_argument("--shared", required=True)
    parser.add_argument("--static", required=True)
    parser.add_argument("--out-shared", type=Path, required=True)
    parser.add_argument("--out-static", type=Path, required=True)
    parser.add_argument("--implib")
    parser.add_argument("--out-implib", type=Path)
    args = parser.parse_args()
    if (args.implib is None) != (args.out_implib is None):
        parser.error("--implib and --out-implib must be supplied together")
    return args


def main() -> None:
    args = parse_args()
    command = [
        args.cargo,
        "rustc",
        "--lib",
        "--manifest-path",
        str(args.source / "Cargo.toml"),
        "--target-dir",
        str(args.target),
        "--features",
        "capi",
    ]
    if (args.source / "Cargo.lock").is_file():
        command.append("--locked")
    if args.profile == "release":
        command.append("--release")

    rustc = ["--crate-type=cdylib", "--crate-type=staticlib"]
    if sys.platform.startswith("linux"):
        rustc.append(f"-Clink-arg=-Wl,-soname,{args.shared}")
    elif sys.platform == "darwin":
        rustc.append(f"-Clink-arg=-Wl,-install_name,@rpath/{args.shared}")
    command.extend(["--", *rustc])
    subprocess.check_call(command)

    built = args.target / args.profile
    shutil.copy2(built / args.shared, args.out_shared)
    shutil.copy2(built / args.static, args.out_static)
    if args.implib is not None:
        shutil.copy2(built / args.implib, args.out_implib)


if __name__ == "__main__":
    main()

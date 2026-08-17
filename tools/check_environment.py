#!/usr/bin/env python3
"""在运行主流程前检查 Python 依赖和 Rscript 是否可用。"""
from __future__ import annotations

import importlib.util
import argparse
import shutil
import sys
from pathlib import Path

import yaml


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/project.yaml", help="YAML configuration file (default: config/project.yaml).")
    arguments = parser.parse_args()
    packages = ["numpy", "pandas", "scipy", "statsmodels", "matplotlib", "yaml"]
    missing = [name for name in packages if importlib.util.find_spec(name) is None]
    if missing:
        raise SystemExit("Missing Python packages: " + ", ".join(missing))
    config_path = Path(arguments.config).expanduser().resolve()
    if not config_path.is_file():
        raise SystemExit(f"Configuration file not found: {config_path}")
    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    rscript = str(config.get("runtime", {}).get("rscript", "Rscript"))
    if shutil.which(rscript) is None and not Path(rscript).is_file():
        raise SystemExit(f"Configured Rscript executable is unavailable: {rscript}")
    print(f"Environment check passed (Rscript: {rscript}).")


if __name__ == "__main__":
    main()

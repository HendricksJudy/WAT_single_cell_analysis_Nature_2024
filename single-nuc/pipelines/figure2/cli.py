"""Command-line interface for the Figure 2 single-nucleus pipeline."""

from __future__ import annotations

import argparse
from importlib import util as importlib_util
from pathlib import Path
from typing import Optional, Sequence

try:  # pragma: no cover - allows running as a module
    from .pipeline import Figure2Config, run_pipeline
except ImportError:  # pragma: no cover - executed when run as a script
    MODULE_PATH = Path(__file__).resolve().with_name("pipeline.py")
    spec = importlib_util.spec_from_file_location("figure2_pipeline", MODULE_PATH)
    if spec is None or spec.loader is None:  # pragma: no cover
        raise ImportError("Unable to load figure2 pipeline module")
    pipeline = importlib_util.module_from_spec(spec)
    spec.loader.exec_module(pipeline)
    Figure2Config = pipeline.Figure2Config  # type: ignore[attr-defined]
    run_pipeline = pipeline.run_pipeline  # type: ignore[attr-defined]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the Figure 2 single-nucleus analysis workflow")
    parser.add_argument("config", type=Path, help="Path to a YAML or JSON configuration file")
    parser.add_argument("--integrated-h5ad", type=Path, default=None, help="Override the integrated .h5ad input path")
    parser.add_argument("--output-dir", type=Path, default=None, help="Override the output directory")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    cfg = Figure2Config.from_file(args.config)
    if args.integrated_h5ad is not None:
        cfg.integrated_h5ad = args.integrated_h5ad
    if args.output_dir is not None:
        cfg.output_dir = args.output_dir
    run_pipeline(cfg)


if __name__ == "__main__":
    main()

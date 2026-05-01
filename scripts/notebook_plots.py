#!/usr/bin/env python3
"""Extract plots from executed notebooks and save as PNGs for LLM review.

Usage:
    python scripts/notebook_plots.py notebooks/smoke_test.ipynb
    python scripts/notebook_plots.py notebooks/          # all notebooks

Outputs PNG files to notebooks/_plots/ directory.
Each plot gets a descriptive filename based on notebook name and cell index.
"""

import json
import os
import sys
from pathlib import Path


def extract_plots(notebook_path: str, output_dir: str = None) -> list:
    """Extract all matplotlib figures from an executed notebook.

    Looks for PNG image data in cell outputs (standard Jupyter output format).
    """
    nb_path = Path(notebook_path)
    with open(nb_path) as f:
        nb = json.load(f)

    if output_dir is None:
        output_dir = nb_path.parent / "_plots"
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    plots = []
    nb_name = nb_path.stem

    for cell_idx, cell in enumerate(nb["cells"]):
        if cell["cell_type"] != "code":
            continue
        for output in cell.get("outputs", []):
            if output.get("output_type") != "display_data":
                continue
            # Jupyter stores images as base64 in data["image/png"]
            data = output.get("data", {})
            png_b64 = data.get("image/png")
            if png_b64:
                plot_path = output_dir / f"{nb_name}_cell{cell_idx}.png"
                import base64

                with open(plot_path, "wb") as f:
                    f.write(base64.b64decode(png_b64))
                plots.append(str(plot_path))

    return plots


def main():
    if len(sys.argv) < 2:
        print("Usage: python scripts/notebook_plots.py <notebook_or_dir> ...")
        print("  Pass a .ipynb file or a directory (processes all .ipynb files)")
        sys.exit(1)

    all_plots = []
    for arg in sys.argv[1:]:
        p = Path(arg)
        if p.suffix == ".ipynb":
            all_plots.extend(extract_plots(str(p)))
        elif p.is_dir():
            for nb in sorted(p.glob("*.ipynb")):
                if nb.name == "conftest.py":
                    continue
                all_plots.extend(extract_plots(str(nb)))

    if not all_plots:
        print("No plots found.")
        sys.exit(0)

    print(f"Extracted {len(all_plots)} plot(s):")
    for p in all_plots:
        print(f"  {p}")
    print()
    print("To review with LLM, send the PNG paths using @ notation:")
    print("  @ " + " @ ".join(all_plots))


if __name__ == "__main__":
    main()

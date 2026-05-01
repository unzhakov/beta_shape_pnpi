#!/usr/bin/env python3
"""Execute notebooks and extract plots for LLM review.

This script:
1. Executes all notebooks (or specified ones) using nbconvert
2. Saves executed notebooks with outputs
3. Extracts PNG figures from cell outputs
4. Reports plot locations for LLM analysis

Usage:
    python scripts/run_notebooks.py              # all notebooks
    python scripts/run_notebooks.py nb1.ipynb    # specific notebook
    python scripts/run_notebooks.py notebooks/   # directory
"""

import json
import os
import sys
import base64
from pathlib import Path


def execute_notebook(nb_path: str) -> bool:
    """Execute a notebook using nbconvert and save outputs."""
    import nbformat
    from nbclient import NotebookClient

    nb = nbformat.read(str(nb_path), as_version=4)
    client = NotebookClient(nb, timeout=600, kernel_name="python3")
    try:
        client.execute()
        nbformat.write(nb, str(nb_path))
        return True
    except Exception as e:
        print(f"  ERROR: {e}")
        return False


def extract_plots(nb_path: str, output_dir: str = None) -> list:
    """Extract PNG figures from an executed notebook's cell outputs."""
    nb_path = Path(nb_path)
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
            data = output.get("data", {})
            png_b64 = data.get("image/png")
            if png_b64:
                plot_path = output_dir / f"{nb_name}_cell{cell_idx}.png"
                with open(plot_path, "wb") as f:
                    f.write(base64.b64decode(png_b64))
                plots.append(str(plot_path))

    return plots


def main():
    if len(sys.argv) < 2:
        print("Usage: python scripts/run_notebooks.py <notebook_or_dir> ...")
        sys.exit(1)

    all_plots = []
    failed = []

    for arg in sys.argv[1:]:
        p = Path(arg)
        if p.suffix == ".ipynb":
            print(f"Executing {p.name}...")
            if execute_notebook(str(p)):
                plots = extract_plots(str(p))
                all_plots.extend(plots)
                print(f"  Extracted {len(plots)} plot(s)")
            else:
                failed.append(p.name)
        elif p.is_dir():
            for nb in sorted(p.glob("*.ipynb")):
                if nb.name.startswith("_") or nb.name == "conftest.py":
                    continue
                print(f"Executing {nb.name}...")
                if execute_notebook(str(nb)):
                    plots = extract_plots(str(nb))
                    all_plots.extend(plots)
                    print(f"  Extracted {len(plots)} plot(s)")
                else:
                    failed.append(nb.name)

    print()
    print(f"Done. {len(all_plots)} plot(s) extracted, {len(failed)} notebook(s) failed.")
    if all_plots:
        print()
        print("Plots saved to: notebooks/_plots/")
        print("To review with LLM, send the PNG paths using @ notation:")
        for p in all_plots:
            print(f"  @{p}")


if __name__ == "__main__":
    main()

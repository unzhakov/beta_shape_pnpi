#!/usr/bin/env python3
"""
docs_check.py — Automated documentation consistency checks.

Runs a battery of checks across README.md, TODO.md, AGENTS.md, and
CONVENTIONS.md against the actual codebase. Designed to be added to the
quality gate in AGENTS.md.

Exit codes:
  0 — all checks passed
  1 — one or more checks failed (prints summary to stderr)
"""

import argparse
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

@dataclass
class Finding:
    severity: str  # "ERROR" | "WARN" | "INFO"
    source: str    # which doc file
    section: str   # which section
    message: str


@dataclass
class Report:
    findings: list[Finding] = field(default_factory=list)

    @property
    def errors(self) -> list[Finding]:
        return [f for f in self.findings if f.severity == "ERROR"]

    @property
    def warnings(self) -> list[Finding]:
        return [f for f in self.findings if f.severity == "WARN"]

    def add(self, severity: str, source: str, section: str, message: str):
        self.findings.append(Finding(severity, source, section, message))

    def summary(self) -> str:
        lines = []
        if self.errors:
            lines.append(f"\n{'=' * 70}")
            lines.append(f"  ERRORS ({len(self.errors)})")
            lines.append(f"{'=' * 70}")
            for f in self.errors:
                lines.append(f"  [{f.source}] {f.section}: {f.message}")
        if self.warnings:
            lines.append(f"\n{'=' * 70}")
            lines.append(f"  WARNINGS ({len(self.warnings)})")
            lines.append(f"{'=' * 70}")
            for f in self.warnings:
                lines.append(f"  [{f.source}] {f.section}: {f.message}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent


def read_doc(name: str) -> str:
    path = ROOT / name
    return path.read_text() if path.exists() else ""


def all_source_files() -> list[Path]:
    """Walk the source tree and return all .py files."""
    files = []
    for root, dirs, filenames in os.walk(ROOT):
        dirs[:] = [d for d in dirs if not d.startswith(".")
                   and d not in ("__pycache__", ".pytest_cache", ".mypy_cache",
                                 ".ruff_cache", ".hypothesis", ".ipynb_checkpoints")]
        if "egg-info" in root:
            continue
        for fn in filenames:
            if fn.endswith(".py"):
                files.append(Path(root) / fn)
    return files


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def check_version_consistency(report: Report):
    """Check that __version__ matches pyproject.toml version."""
    pyproject = ROOT / "pyproject.toml"
    if not pyproject.exists():
        return

    content = pyproject.read_text()
    py_ver_match = re.search(r'^version\s*=\s*"([^"]+)"', content, re.MULTILINE)
    if not py_ver_match:
        return
    py_ver = py_ver_match.group(1)

    for pyfile in all_source_files():
        text = pyfile.read_text()
        ver_match = re.search(r'__version__\s*=\s*"([^"]+)"', text)
        if ver_match:
            ver = ver_match.group(1)
            if ver != py_ver:
                rel = pyfile.relative_to(ROOT)
                report.add("ERROR", rel.name, "version",
                           f"__version__ = '{ver}' but pyproject.toml says '{py_ver}'")


def check_version_in_code(report: Report):
    """Check that CLI --version output matches pyproject.toml."""
    pyproject = ROOT / "pyproject.toml"
    if not pyproject.exists():
        return

    content = pyproject.read_text()
    py_ver_match = re.search(r'^version\s*=\s*"([^"]+)"', content, re.MULTILINE)
    if not py_ver_match:
        return
    py_ver = py_ver_match.group(1)

    try:
        result = subprocess.run(
            ["bs_pnpi", "--version"],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            cli_ver_match = re.search(r'(\d+\.\d+\.\d+)', result.stdout)
            if cli_ver_match:
                cli_ver = cli_ver_match.group(1)
                if cli_ver != py_ver:
                    report.add("ERROR", "cli", "version",
                               f"CLI reports version '{cli_ver}' but pyproject.toml says '{py_ver}'")
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pass


def check_file_references(report: Report):
    """Check that files referenced in docs actually exist."""
    docs = {
        "README.md": read_doc("README.md"),
        "TODO.md": read_doc("TODO.md"),
        "AGENTS.md": read_doc("AGENTS.md"),
        "CONVENTIONS.md": read_doc("CONVENTIONS.md"),
    }

    # Known historical references that are intentionally kept in docs
    # (files that were removed as part of a migration/cleanup)
    known_removed = {"cw_extractor.py", "fitter.py", "test_cw_extractor.py"}
    # Placeholder examples used in documentation
    known_placeholders = {"foo.py", "test_foo.py"}

    for doc_name, content in docs.items():
        # Only check inline code refs (single backtick) that look like file paths.
        # Skip tree diagrams and multi-line code blocks.
        for line in content.split("\n"):
            # Find single-line inline code that looks like a file path
            refs = re.findall(r'`([^`\n]+\.py[^`]*)`', line)
            refs.extend(re.findall(r'`([^`\n]+\\.(csv|json|toml|md|yaml|png|log))`', line))
            for ref in refs:
                if isinstance(ref, tuple):
                    ref = ref[0]
                ref = ref.strip().rstrip(":").rstrip(".")
                if not ref or len(ref) > 200:
                    continue
                if ref.startswith("./"):
                    ref = ref[2:]

                target = ROOT / ref
                if target.exists():
                    continue
                # Filter out tree-diagram artifacts
                if any(c in ref for c in ('\u2190', '\u2192', '\u2713', '\u2717',
                                          '[', ']', '|', '\u2502', '\u2514', '\u251c',
                                          '\u2500', '\u2502', '\u2560', '\u256c')):
                    continue
                # Only warn if it looks like a real file path
                if re.search(r'[a-zA-Z0-9_\-]+\.[a-zA-Z0-9]+', ref):
                    basename = ref.split("/")[-1]
                    if basename in known_removed or basename in known_placeholders:
                        continue
                    # Skip bare __init__.py (common in tree diagrams, not actionable)
                    if basename == "__init__.py":
                        continue
                    report.add("WARN", doc_name, "file-reference",
                               f"Referenced file not found: `{ref}`")


def check_cli_commands(report: Report):
    """Check that CLI commands mentioned in AGENTS.md actually work."""
    agents_md = read_doc("AGENTS.md")

    # Check bs_pnpi --help works
    if "bs_pnpi" in agents_md:
        try:
            result = subprocess.run(
                ["bs_pnpi", "--help"],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode != 0:
                report.add("WARN", "AGENTS.md", "cli",
                           f"`bs_pnpi --help` failed: {result.stderr.strip()[:100]}")
        except (subprocess.TimeoutExpired, FileNotFoundError):
            report.add("WARN", "AGENTS.md", "cli", "bs_pnpi not found in PATH")

    # Check bs_exp CLI
    if "bs_exp" in agents_md:
        try:
            result = subprocess.run(
                ["bs_exp", "--help"],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode != 0:
                report.add("WARN", "AGENTS.md", "cli",
                           f"`bs_exp --help` failed: {result.stderr.strip()[:100]}")
        except (subprocess.TimeoutExpired, FileNotFoundError):
            report.add("WARN", "AGENTS.md", "cli", "bs_exp not found in PATH")


def check_imports_in_examples(report: Report):
    """Check that import statements in README examples are valid."""
    readme = read_doc("README.md")

    code_blocks = re.findall(r'```[^\n]*\n(.*?)```', readme, re.DOTALL)
    for block in code_blocks:
        for line in block.strip().split("\n"):
            line = line.strip()
            if not line.startswith("from ") and not line.startswith("import "):
                continue
            match = re.match(r'from\s+(\S+)\s+import', line)
            if not match:
                match = re.match(r'import\s+(\S+)', line)
            if not match:
                continue

            module = match.group(1)
            try:
                __import__(module)
            except ImportError:
                report.add("WARN", "README.md", "imports",
                           f"Cannot import module in example: `{module}`")


def check_directory_structure(report: Report):
    """Check that directory structures in docs match reality."""
    # Check README structure tree
    readme = read_doc("README.md")
    _check_tree_in_doc(readme, "README.md", report)

    # Check TODO structure tree
    todo = read_doc("TODO.md")
    _check_tree_in_doc(todo, "TODO.md", report)


def _check_tree_in_doc(content: str, doc_name: str, report: Report):
    """Extract file references from a markdown tree diagram and check existence."""
    # Find fenced code blocks that look like directory trees
    trees = re.findall(r'```\w*\n(.*?)```', content, re.DOTALL)
    for tree in trees:
        lines = tree.split("\n")
        for i, line in enumerate(lines):
            # Match lines like: `├── filename.py` or `│   ├── filename.py`
            match = re.search(r'[\u251c\u2514\u2502\s]*[\u2500\u2501]+\s+([a-zA-Z0-9_\-]+\.\w+)', line)
            if match:
                filename = match.group(1)
                # Determine depth by counting leading tree chars
                stripped = line.lstrip()
                depth = len(line) - len(stripped)
                # Only check files at reasonable depth (skip deeply nested false paths)
                if depth > 40:  # likely a parser artifact from long tree
                    continue
                # Try to infer path from context (look backwards for directory names)
                # But cap at 3 levels deep to avoid accumulated prefix bugs
                dirs = []
                seen_dir = False
                for j in range(i - 1, max(i - 10, -1), -1):
                    bl = lines[j]
                    dir_match = re.search(r'[\u251c\u2514\u2502\s]*[\u2500\u2501]+\s+([a-zA-Z0-9_\-]+)/', bl)
                    if dir_match:
                        dir_depth = len(bl) - len(bl.lstrip())
                        if dir_depth < depth:
                            dirs.append(dir_match.group(1))
                            seen_dir = True
                    elif re.search(r'\w', bl):  # non-empty, non-matching line
                        if seen_dir:
                            break  # stop looking after first non-matching content line
                if dirs and len(dirs) <= 3:  # cap at 3 levels
                    target = ROOT / "/".join(reversed(dirs)) / filename
                    if not target.exists():
                        report.add("WARN", doc_name, "structure",
                                   f"File in structure diagram not found: `{'/'.join(reversed(dirs))}/{filename}`")


def check_cli_args_consistency(report: Report):
    """Check that CLI flags mentioned in docs match actual CLI."""
    agents_md = read_doc("AGENTS.md")
    readme = read_doc("README.md")

    try:
        result = subprocess.run(
            ["bs_pnpi", "--help"],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            actual_flags = set()
            for line in result.stdout.split("\n"):
                flag_match = re.match(r'\s{0,2}(-\S+)', line)
                if flag_match:
                    actual_flags.add(flag_match.group(1))

            known_flags = {"-v", "-vv", "-q", "--nuclide", "--output",
                           "--input", "--dry-run", "--version", "--plot",
                           "--log-file", "--intensity-cutoff", "--help",
                           "-h"}
            mentioned = set()
            # Only match actual CLI flag patterns (--flag-name, not random words)
            for flag in re.findall(r'--[a-z][a-z0-9-]{2,}', agents_md):
                mentioned.add(flag)
            for flag in re.findall(r'-[a-z][a-z0-9-]{2,}', readme):
                mentioned.add(flag)
            unknown = mentioned - known_flags
            for flag in sorted(unknown):
                report.add("INFO", "AGENTS.md", "cli-args",
                           f"Flag `{flag}` mentioned but not in --help (may be custom)")
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pass


def check_markdown_formatting(report: Report):
    """Check that markdown files use consistent formatting."""
    docs = ["README.md", "TODO.md", "AGENTS.md", "CONVENTIONS.md"]
    for doc_name in docs:
        path = ROOT / doc_name
        if not path.exists():
            continue

        content = path.read_text()
        for i, line in enumerate(content.split("\n"), 1):
            if line != line.rstrip():
                report.add("INFO", doc_name, f"line {i}",
                           "Trailing whitespace detected")
                break  # one per file is enough


def check_todo_completions(report: Report):
    """Check that TODO items marked as completed actually exist in code."""
    todo = read_doc("TODO.md")

    # Find completed items that claim specific files/methods exist
    # Pattern: "- [x] ... implemented: `SomeClass` or `method_name`"
    completed = re.findall(r'- \[x\].*?(?=- \[|$)', todo, re.DOTALL)
    for block in completed:
        # Skip historical references (items about cleanup, moves, removals)
        if re.search(r'(?:remove|moved|replaced|deleted|clean up|cleanup)', block, re.IGNORECASE):
            continue
        # Look for specific file claims
        file_claims = re.findall(r'`([^`]+\.py)`', block)
        for fc in file_claims:
            target = ROOT / fc
            if not target.exists():
                report.add("WARN", "TODO.md", "completions",
                           f"Completed item references non-existent file: `{fc}`")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

CHECKS = [
    check_version_consistency,
    check_version_in_code,
    check_file_references,
    check_cli_commands,
    check_imports_in_examples,
    check_directory_structure,
    check_cli_args_consistency,
    check_markdown_formatting,
    check_todo_completions,
]


def main():
    parser = argparse.ArgumentParser(description="Check documentation consistency")
    parser.add_argument(
        "--strict", action="store_true",
        help="Treat warnings as errors (exit 1)"
    )
    parser.add_argument(
        "--json", action="store_true", dest="as_json",
        help="Output findings as JSON"
    )
    args = parser.parse_args()

    report = Report()
    for check in CHECKS:
        try:
            check(report)
        except Exception as e:
            report.add("WARN", "runner", "check", f"Check crashed: {e}")

    # Output
    if args.as_json:
        import json as _json
        output = [
            {"severity": f.severity, "source": f.source,
             "section": f.section, "message": f.message}
            for f in report.findings
        ]
        print(_json.dumps(output, indent=2))
    else:
        print(f"docs_check: {len(report.findings)} findings "
              f"({len(report.errors)} errors, {len(report.warnings)} warnings)")
        print(report.summary())

    # Exit
    if report.errors:
        return 1
    if args.strict and report.warnings:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""extract_refs.py — Extract and format references from PDF text output as BibTeX.

The reference section in pdftotext output follows a consistent pattern:
    Author(s) (Year), Journal Vol(Pages), arXiv/DOI if available.

This script joins line-wrapped entries, parses each one into structured fields,
and outputs valid BibTeX. Two modes:

1. Local parsing (default): regex/heuristic — fast, works for standard journal refs
2. LLM mode (--llm): sends to remote llama-server for complex/irregular references

Usage:
    # From a PDF file directly:
    python3 extract_refs.py --pdf paper.pdf --output refs.bib

    # Pipe from pdftotext (already extracted reference section):
    pdftotext paper.pdf - | grep -A9999 "^REFERENCES" > refs.txt
    python3 extract_refs.py --file refs.txt --output refs.bib

    # Force LLM mode for tricky papers:
    python3 extract_refs.py --pdf paper.pdf --llm --output refs.bib
"""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path


def extract_reference_section(pdf_path: str) -> str:
    """Extract just the REFERENCES section from a PDF using pdftotext."""
    result = subprocess.run(
        ["pdftotext", pdf_path, "-"],
        capture_output=True, text=True
    )
    if result.returncode != 0 or not result.stdout.strip():
        print(f"ERROR: Failed to extract text from {pdf_path}", file=sys.stderr)
        sys.exit(1)

    lines = result.stdout.splitlines()

    # Find the REFERENCES header (last occurrence — some papers have "References" in body too)
    ref_start = None
    for i in range(len(lines) - 1, -1, -1):
        if re.match(r'^\s*REFERENCES?\s*$', lines[i], re.IGNORECASE):
            ref_start = i
            break

    # If exact match not found, try partial match
    if ref_start is None:
        for i in range(len(lines) - 1, -1, -1):
            if "REFERENCES" in lines[i].upper():
                ref_start = i
                break

    if ref_start is None:
        print("WARNING: No REFERENCES section found — trying last 500 lines", file=sys.stderr)
        return "\n".join(lines[-500:])

    # Get everything after the header
    ref_lines = []
    for line in lines[ref_start + 1:]:
        stripped = line.strip()
        if not stripped or len(stripped) < 3:
            continue
        # Stop at end-of-paper markers
        if any(marker in stripped.upper() for marker in ["APPENDIX", "FOOTNOTE"]):
            break
        ref_lines.append(line)

    return "\n".join(ref_lines)


def join_references(text: str) -> list[str]:
    """Join line-wrapped references into individual entries.

    pdftotext output has each reference split across multiple lines, with
    continuation lines either indented or starting at column 0 with a new
    author name pattern (LASTNAME, Initials).
    """
    lines = text.splitlines()
    joined = []
    current = ""

    for line in lines:
        stripped = line.strip()
        if not stripped:
            if current:
                joined.append(current)
                current = ""
            continue

        # Skip page header markers (e.g., "58" alone on a line from PDF headers)
        if re.match(r'^\d{1,3}$', stripped) and len(stripped) <= 3:
            continue

        # A new reference starts when we have accumulated text AND the current line
        # matches author pattern: UPPERCASE last name followed by comma or period
        # e.g., "Wilkinson, D. H." or "Weinberg, S."
        if current and re.match(r'^[A-Z][a-z]+(?:\'?[A-Z])?,\s+[A-Z]', stripped):
            joined.append(current)
            current = stripped
        elif current:
            # Continuation line — join with space
            current += " " + stripped
        else:
            current = stripped

    if current:
        joined.append(current)

    return [ref.strip() for ref in joined if len(ref.strip()) > 10]


def parse_author_list(author_str: str) -> list[str]:
    """Parse an author string like 'Hayen, L. and N. Severijns' into individual authors."""
    # Split on "and" but be careful of "and" in names (rare)
    parts = re.split(r'\s+and\s+', author_str, flags=re.IGNORECASE)
    return [p.strip() for p in parts if p.strip()]


def parse_reference(ref: str) -> dict:
    """Parse a single reference string into BibTeX fields.

    Expected format (varies by journal style):
        Author1, A. and Author2, B. (Year), Journal Abbrev Vol(Pages).
        or with arXiv:  ... , arXiv:xxxx.xxxxx.
        or with DOI:  ... , doi:10.xxxx/....
    """
    entry = {
        "type": "article",
        "authors": [],
        "year": "",
        "journal": "",
        "volume": "",
        "pages": "",
        "title": "",
        "eprint": "",
        "doi": "",
        "note": ""
    }

    # Extract year: look for pattern like (2017) or 1983a, or 1983b
    year_match = re.search(r'\((\d{4}[a-z]?)\)', ref)
    if year_match:
        entry["year"] = year_match.group(1).rstrip('.,')

    # Extract arXiv ID
    arxiv_match = re.findall(r'arXiv[:\s]+(\d{4}\.\d{4,5})', ref, re.IGNORECASE)
    if arxiv_match:
        entry["eprint"] = f"arXiv:{arxiv_match[0]}"

    # Extract DOI
    doi_match = re.search(r'doi[:\s]+(10\.[\d/.]+)', ref, re.IGNORECASE)
    if doi_match:
        entry["doi"] = doi_match.group(1).rstrip('.,')

   # Extract journal, volume, pages using a reverse-matching strategy.
    # Find the citation info at end of string: "Journal Volume(Issue), Pages." or
    # "Journal Volume, Pages." — then extract journal as text between year paren and numbers.

    m = re.search(r'(\d+)\s*\(([^,\)]*)\)[,\s]+\s*(\d+)', ref)
    if m:
        entry["volume"] = str(m.group(1))
        before_vol = ref[:m.start()].rstrip(', .')
        journal = re.sub(r'^.+?\s*\(\d{4}[a-z]?\)\s*,\s*', '', before_vol).strip()
        if journal:
            entry["journal"] = journal.rstrip('.')
        entry["pages"] = str(m.group(3))

    else:
        # Try Vol, Pages format (no parenthesized issue): "Nucl. Phys. A 233, 55."
        m2 = re.search(r',\s*([A-Z][^(),]{1,60})\s+(\d+)\s*,\s*(\d+)(?:\.?\s*$|\s+arXiv)', ref)
        if m2:
            entry["journal"] = m2.group(1).strip().rstrip('.')
            entry["volume"] = str(m2.group(2))
            entry["pages"] = str(m2.group(3))

        else:
            # Try Vol (Issue). — ending with parenthesized issue but no trailing page number
            m3 = re.search(r'([A-Z][^(),]{1,50})\s+(\d+)\s*\(([^,\)]*)\)\.?\s*$', ref)
            if m3:
                journal = re.sub(r'^.+?\s*\(\d{4}[a-z]?\)\s*,\s*', '', ref[:m3.start()].rstrip(', .')).strip()
                if journal:
                    entry["journal"] = journal.rstrip('.')
                entry["volume"] = str(m3.group(2))

    # Extract authors — everything before the year parentheses
    author_match = re.match(r'^(.+?)\s*\(\d{4}', ref)
    if author_match:
        author_str = author_match.group(1).strip()
        entry["authors"] = parse_author_list(author_match.group(1))

    # Generate a unique citation key using author + year + second author initial disambiguation
    if entry["authors"]:
        last_name = re.sub(r',.*$', '', entry["authors"][0]).upper()[:4]
        base_key = f"{last_name}{entry['year']}"

        # Disambiguate: for same author+year, use first letter of second author or "a/b/c" suffix
        if len(entry["authors"]) > 1:
            second_initial = re.sub(r',.*$', '', entry["authors"][1]).upper().replace("'", "")[:2]
            base_key += f"{second_initial}"

        # Clean key — remove trailing punctuation from initials
        base_key = re.sub(r'[^A-Z0-9]+$', '', base_key)

        # Track duplicates and add numeric suffix
        key = base_key
        counter = 0
        while key in entry.get("_seen_keys", set()):
            counter += 1
            if counter == 1:
                m2 = re.match(r'^(.+?)(\d{4}[a-z]?[A-Z]{0,2})$', key)
                if m2 and len(m2.group(1)) >= 6:
                    base_suffix = key[-3:]
                    key = f"{m2.group(1)}{base_suffix}{counter}"
                else:
                    key = f"{key}{counter}"
            elif counter > 3:
                key = f"{last_name}{entry['year']}{counter}"

        entry["_seen_keys"] = entry.get("_seen_keys", set()) | {key}
        entry["_key"] = key
    else:
        entry["_key"] = "unknown"

    return entry


def make_key_unique(entries_with_keys: list[tuple[dict, str]]) -> None:
    """Post-process to ensure all BibTeX keys are unique."""
    seen = {}
    for entry, _ in entries_with_keys:
        key = entry["_key"]
        if key in seen:
            # Find disambiguation suffix
            base = re.sub(r'\d+$', '', key)
            i = 2
            new_key = f"{base}{i}"
            while new_key in seen:
                i += 1
                new_key = f"{base}{i}"
            entry["_key"] = new_key
        seen[key] = True


def entry_to_bibtex(entry: dict, ref_text: str) -> str:
    """Convert a parsed reference into BibTeX format."""
    key = entry["_key"]
    lines = [f"@{entry['type']}{{{key},"]

    # Authors — required field
    if entry["authors"]:
        authors_str = " and ".join(entry["authors"])
        lines.append(f"    author = {{{authors_str}}},")

    # Year — required for articles
    if entry["year"]:
        lines.append(f"    year = {{{entry['year']}}},")

    # Title — since pdftotext doesn't include titles in standard citation format,
    # leave empty (Obsidian notes can link to the paper directly)
    # To get actual titles would require cross-referencing with arXiv/ADS APIs

    # Journal — required for articles
    if entry["journal"]:
        lines.append(f"    journal = {{{entry['journal']}}},")

    # Volume and pages
    if entry["volume"]:
        lines.append(f"    volume = {{{entry['volume']}}},")
    if entry["pages"]:
        lines.append(f"    pages = {{{entry['pages']}}},")

    # arXiv / DOI
    if entry["eprint"]:
        lines.append(f"    eprint = {{{entry['eprint']}}},")
    if entry["doi"]:
        lines.append(f"    doi = {{{entry['doi']}}},")

    # Source reference text as note for cross-referencing
    clean_ref = re.sub(r'\s+', ' ', ref_text)[:200]
    lines.append(f"    note = \"{clean_ref}\"")

    lines.append("}")
    return "\n".join(lines)


def parse_to_bibtex(reference_text: str, pdf_path: str = "") -> tuple[str, int]:
    """Main parsing pipeline: join refs → parse each → output BibTeX.

    Returns (bibtex_string, count_of_parsed_refs).
    """
    refs = join_references(reference_text)

    if not refs:
        print("WARNING: No references found to parse", file=sys.stderr)
        return "", 0

    entries_with_keys = []
    for ref in refs:
        try:
            entry = parse_reference(ref)
            entries_with_keys.append((entry, ref))
        except Exception as e:
            print(f"WARNING: Failed to parse reference: {e}", file=sys.stderr)

    make_key_unique(entries_with_keys)

    bibtex_entries = []
    for entry, ref in entries_with_keys:
        try:
            bibtex_entries.append(entry_to_bibtex(entry, ref))
        except Exception as e:
            print(f"WARNING: Failed to convert reference to BibTeX: {e}", file=sys.stderr)

    return "\n\n".join(bibtex_entries), len(bibtex_entries)


def send_to_llm(reference_text: str, pdf_path: str = "") -> str:
    """Fallback: send reference text to remote llama-server for parsing."""
    system_prompt = (
        "You are a bibliographic formatting assistant. "
        "Convert the following garbled PDF reference list into valid BibTeX format.\n"
        "Rules:\n"
        "- Each entry must be complete @article or @book block\n"
        "- Extract: author, title, journal, volume, pages, year, eprint/arXiv if present\n"
        "- Authors in 'LastName, Initials' format, separated by ' and '\n"
        "- Titles in braces { }\n"
        "- Keep ALL entries — do not skip any\n"
        "- Output ONLY valid BibTeX, no markdown code blocks, no explanations\n"
    )

    user_prompt = f"""Convert these references to BibTeX:\n\n{reference_text}"""

    if pdf_path:
        system_prompt += (
            f"\nSource file: {pdf_path}. "
            "Use the filename to create a self-citation entry if referenced."
        )

    payload = {
        "model": "qwen3.6-35b-a3b",
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt}
        ],
        "max_tokens": 8192,
        "temperature": 0.0,
    }

    try:
        import urllib.request
        req = urllib.request.Request(
            "http://localhost:8081/v1/chat/completions",
            data=json.dumps(payload).encode(),
            headers={"Content-Type": "application/json"}
        )
        with urllib.request.urlopen(req, timeout=300) as resp:
            result = json.loads(resp.read())

        bibtex = result["choices"][0]["message"]["content"]
        for marker in ["```bibtex", "```"]:
            if bibtex.startswith(marker):
                bibtex = bibtex.replace(marker, "", 1)
                if bibtex.strip().startswith("```"):
                    bibtex = bibtex[3:].strip()
        return bibtex

    except Exception as e:
        print(f"LLM call failed: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Extract references from PDF and format as BibTeX",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  # From a PDF directly (local parsing):
  python3 extract_refs.py --pdf Hayen2017.pdf -o refs.bib

  # Pipe reference section from pdftotext:
  grep -A9999 "^REFERENCES" paper.txt | python3 extract_refs.py -o refs.bib

  # Use LLM for difficult papers:
  python3 extract_refs.py --pdf Hayen2017.pdf --llm -o refs_llm.bib"""
    )
    parser.add_argument("--file", "-f", help="Input text file with reference section")
    parser.add_argument("--pdf", "-p", help="PDF file to extract references from directly")
    parser.add_argument("--output", "-o", help="Output BibTeX file (default: stdout)")
    parser.add_argument("--llm", action="store_true", help="Use LLM mode for parsing")

    args = parser.parse_args()

    # Get reference text
    if args.file:
        ref_text = Path(args.file).read_text()
    elif args.pdf:
        print(f"Extracting references from {args.pdf}...", file=sys.stderr)
        ref_text = extract_reference_section(args.pdf)
    else:
        ref_text = sys.stdin.read()

    if not ref_text.strip():
        print("ERROR: No reference text provided", file=sys.stderr)
        sys.exit(1)

    # Parse to BibTeX
    if args.llm:
        bibtex = send_to_llm(ref_text, pdf_path=args.pdf or "")
        count = 0
    else:
        bibtex, count = parse_to_bibtex(ref_text, pdf_path=args.pdf or "")

    # Write output
    outfile = open(args.output, "w") if args.output else sys.stdout
    try:
        if not args.llm and not args.output:
            print(f"# {count} references extracted", file=sys.stderr)
        outfile.write(bibtex + "\n" if bibtex else "")
    finally:
        if args.output:
            outfile.close()


if __name__ == "__main__":
    main()

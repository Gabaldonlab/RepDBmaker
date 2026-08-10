#!/usr/bin/env python3
"""Remap `taxonkit create-taxdump`'s taxids to a small, dense, sequential
range, in place. This is because MMseqs doesn't handle sparse taxid ranges well because of
memory issues.

Usage: compact_taxdump_ids.py <taxdump_dir> [taxid_map_file]
"""
import sys
from pathlib import Path


def _fields(line):
    return line.rstrip("\n").split("\t")


def _remap_file(path, taxid_cols, mapping):
    """Rewrite `path` in place, replacing each `mapping`-known taxid at the
    given (0-indexed) columns. Columns not in `mapping` (shouldn't happen for
    a self-consistent taxdump, but not fatal) are left untouched."""
    if not path.exists():
        return
    lines = path.read_text().splitlines(keepends=True)
    if not lines:
        return
    with open(path, "w") as fh:
        for line in lines:
            fields = _fields(line)
            for col in taxid_cols:
                if fields[col] in mapping:
                    fields[col] = str(mapping[fields[col]])
            fh.write("\t".join(fields) + "\n")


def main():
    taxdump_dir = Path(sys.argv[1])
    taxid_map_file = Path(sys.argv[2]) if len(sys.argv) > 2 else None

    nodes_path = taxdump_dir / "nodes.dmp"

    # Collect every taxid that appears (own + parent column), preserving
    # first-seen order so the mapping is deterministic for identical input.
    # Root (parent of itself) lands on whatever number it's first seen at -
    # its actual value doesn't matter, only self-consistency across files.
    mapping = {}
    with open(nodes_path) as fh:
        for line in fh:
            fields = _fields(line)
            for taxid in (fields[0], fields[2]):
                if taxid not in mapping:
                    mapping[taxid] = len(mapping) + 1  # 1-based, dense

    if len(mapping) <= 1:
        return  # 0 or 1 taxon: already trivially dense, nothing to do

    _remap_file(nodes_path, (0, 2), mapping)
    _remap_file(taxdump_dir / "names.dmp", (0,), mapping)
    _remap_file(taxdump_dir / "merged.dmp", (0, 2), mapping)
    _remap_file(taxdump_dir / "delnodes.dmp", (0,), mapping)

    if taxid_map_file and taxid_map_file.exists():
        lines = taxid_map_file.read_text().splitlines()
        with open(taxid_map_file, "w") as fh:
            for line in lines:
                code, taxid = line.split("\t")
                fh.write(f"{code}\t{mapping.get(taxid, taxid)}\n")


if __name__ == "__main__":
    main()

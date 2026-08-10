"""Rename + clean proteomes into a single gzipped RepDB FASTA, in parallel.

Drop-in replacement for parse_gnm_v2.py. Same inputs/outputs and the *same*
sequence IDs (`{taxid}_{code}_AA{padded-counter}`) and accession map, but:

  * work is split across `snakemake.threads` worker processes (v2 was
    effectively single-threaded despite `threads: 112`);
  * each worker streams its share into ONE gzip chunk opened once, instead of
    re-opening the output gzip in append mode for every proteome (v2 did
    ~150k open/write/close cycles, each restarting the deflate stream -> the
    dominant cost);
  * chunks are concatenated (concatenated gzip members = a valid gzip) into the
    final FASTA and accession map.

The per-file IDs depend only on the file and the within-file order, so they are
identical no matter how the proteomes are distributed across workers.

Parallelism uses an explicit *fork* context with raw ``Process`` objects: fork
inherits the worker function from the parent, so nothing is pickled. That is
what lets it run under Snakemake's ``script:`` execution (where a Pool would try
to pickle the worker by qualified name and fail), independent of the
interpreter's default start method.
"""

import gzip
import json
import multiprocessing
import os
import shutil
import sys

# J, B, O, Z, U -> X ; * (stop), whitespace and '\r' deleted. The whitespace/'\r'
# deletion is NOT in v2: some source proteomes (block-formatted or
# Windows-line-ended FASTA) embed spaces/tabs/carriage returns *within*
# sequence lines, which `line.strip()` below does not touch (it only trims
# line edges) - those characters used to pass straight through into the final
# RepDB FASTA and only surface much later as diamond's cryptic "Invalid
# character in sequence: ' '", deep inside a multi-GB gzipped file.
_CLEANING_TABLE = str.maketrans("JBOZU", "XXXXX", "* \t\r")
WRAP = 60          # sequence line width, as in v2
COMPRESSLEVEL = 6  # v2 used gzip's default (9); 6 is much faster, ~same size


def open_stream(filename: str, mode: str = "rt"):
    """Transparently handle gzipped or plain-text files."""
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    return open(filename, mode)


def read_fasta(fafile: str) -> "dict[str, str]":
    """{primary_header: cleaned_sequence}. Primary header = token before the
    first space/tab (matches v2, which also deduplicates by that token)."""
    seqs: "dict[str, str]" = {}
    hdr = None
    parts: "list[str]" = []
    with open_stream(fafile) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line[0] == ">":
                if hdr is not None:
                    seqs[hdr] = "".join(parts).translate(_CLEANING_TABLE)
                    parts = []
                toks = line[1:].split(maxsplit=1)
                hdr = toks[0].split("\t", 1)[0] if toks else ""
            else:
                parts.append(line.strip())
        if hdr is not None:
            seqs[hdr] = "".join(parts).translate(_CLEANING_TABLE)
    return seqs


def chunk_paths(tmpdir: str, bin_id: int):
    base = os.path.join(tmpdir, f"chunk_{bin_id:05d}")
    return base + ".fa.gz", base + ".map", base + ".stat"


def process_bin(bin_id: int, jobs, tmpdir: str):
    """Write one worker's proteomes to a gzip chunk + a map chunk, plus a small
    JSON stat file (protein count + skipped codes) the parent reads back."""
    fa_path, map_path, stat_path = chunk_paths(tmpdir, bin_id)
    n_proteins = 0
    skipped: "list[str]" = []

    with gzip.open(fa_path, "wt", compresslevel=COMPRESSLEVEL) as fo, \
            open(map_path, "w") as mo:
        for path, code, taxid in jobs:
            if taxid is None:  # no taxid for this code -> cannot form the ID
                skipped.append(code)
                continue

            seqs = read_fasta(path)
            width = len(str(len(seqs))) if seqs else 1

            fa_lines: "list[str]" = []
            map_lines: "list[str]" = []
            for idx, (orig_id, seq) in enumerate(seqs.items()):
                protid = "AA" + str(idx).zfill(width)   # v2-compatible ID tail
                new_id = f"{taxid}_{code}_{protid}"
                fa_lines.append(f">{new_id}\n")
                fa_lines.extend(seq[l:l + WRAP] + "\n"
                                for l in range(0, len(seq), WRAP))
                map_lines.append(f"{orig_id}\t{new_id}\n")
                n_proteins += 1

            fo.write("".join(fa_lines))
            mo.write("".join(map_lines))

    with open(stat_path, "w") as s:                     # written last == success
        json.dump({"n": n_proteins, "skipped": skipped}, s)


def main():
    taxdump = snakemake.input.taxdump           # noqa: F821
    table = snakemake.input.table               # noqa: F821
    out_fa = snakemake.output.fa                 # noqa: F821
    out_map = snakemake.output.idmap            # noqa: F821
    threads = max(1, int(snakemake.threads))    # noqa: F821

    # code -> taxid (written by taxonkit create-taxdump -A1)
    taxid_map_path = os.path.join(taxdump, "taxid.map")
    taxids: "dict[str, str]" = {}
    with open(taxid_map_path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                taxids[p[0]] = p[1]

    # jobs: (path, code, taxid, size); size drives load balancing
    jobs = []
    with open(table) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 2:
                continue
            path, code = parts[0], parts[1]
            try:
                size = os.path.getsize(path)
            except OSError:
                size = 0
            jobs.append((path, code, taxids.get(code), size))

    if not jobs:
        open(out_fa, "w").close()
        open(out_map, "w").close()
        print("parse_gnm_v3: empty genome table, nothing to do.", file=sys.stderr)
        return

    # greedy longest-processing-time bin-packing by file size, so the workers
    # finish together despite very uneven proteome sizes.
    nbins = min(threads, len(jobs))
    bins = [[] for _ in range(nbins)]
    loads = [0] * nbins
    for path, code, taxid, size in sorted(jobs, key=lambda x: x[3], reverse=True):
        b = loads.index(min(loads))
        bins[b].append((path, code, taxid))
        loads[b] += size

    tmpdir = os.path.join(os.path.dirname(out_fa) or ".", ".parse_tmp")
    os.makedirs(tmpdir, exist_ok=True)
    for fn in os.listdir(tmpdir):           # clear stale chunks from a prior run
        os.remove(os.path.join(tmpdir, fn))

    # fork context: workers inherit process_bin (no pickling), so this works
    # under Snakemake's exec'd script regardless of the default start method.
    # (fork is safe here: the process is single-threaded at this point.)
    ctx = multiprocessing.get_context("fork")
    procs, ran = [], []
    for i, bin_jobs in enumerate(bins):
        if not bin_jobs:
            continue
        p = ctx.Process(target=process_bin, args=(i, bin_jobs, tmpdir))
        p.start()
        procs.append(p)
        ran.append(i)
    for i, p in zip(ran, procs):
        p.join()
        if p.exitcode != 0:
            raise RuntimeError(f"parse_gnm_v3: worker for bin {i} failed "
                               f"(exit code {p.exitcode})")

    # concatenate chunks in bin order (gzip members concatenate into a valid
    # gzip); read back each worker's stat file (missing == the worker crashed).
    total_proteins = 0
    all_skipped: "list[str]" = []
    with open(out_fa, "wb") as fo:
        for i in ran:
            fa_path, _, stat_path = chunk_paths(tmpdir, i)
            with open(fa_path, "rb") as ci:
                shutil.copyfileobj(ci, fo, 1024 * 1024)
            with open(stat_path) as s:
                st = json.load(s)
            total_proteins += st["n"]
            all_skipped.extend(st["skipped"])
    with open(out_map, "w") as mo:
        for i in ran:
            _, map_path, _ = chunk_paths(tmpdir, i)
            with open(map_path) as ci:
                shutil.copyfileobj(ci, mo, 1024 * 1024)

    for i in ran:
        for p in chunk_paths(tmpdir, i):
            try:
                os.remove(p)
            except OSError:
                pass
    try:
        os.rmdir(tmpdir)
    except OSError:
        pass

    if all_skipped:
        shown = ", ".join(all_skipped[:20]) + (" ..." if len(all_skipped) > 20 else "")
        print(f"WARNING: {len(all_skipped)} proteome(s) had no taxid in "
              f"{taxid_map_path} and were SKIPPED: {shown}", file=sys.stderr)
    print(f"parse_gnm_v3: wrote {total_proteins} proteins from "
          f"{len(jobs) - len(all_skipped)} proteomes across {len(ran)} "
          f"workers -> {out_fa}", file=sys.stderr)


if __name__ == "__main__":
    main()

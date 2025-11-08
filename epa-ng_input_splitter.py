#!/usr/bin/env python3
"""
SVL to per-study(/batch) FASTAs + dedup CSV with sequences + paired sv_long splits.

Inputs:
  - svl.csv with columns: specimen,study,sv,count  (optionally a batch column)
  - FASTA file with all SV sequences (headers match sv IDs)

Notes:
  - If the 'study' column is missing in svl.csv, the script will create it with value 'all' so you can still split into chunks.
  - If you provide a batch column name and the column is missing, it will be ignored with a warning.

Outputs:
  - out/dedup_svl_with_seq.csv                   (first occurrence per sv + sequence)
  - out/study_fastas/{group}.fasta               (non-chunked)
  - out/study_fastas/{group}.partXXX.fasta       (when chunked)
  - out/svl_groups/{group}.sv_long.csv           (non-chunked) OR out/svl_groups/{group}.partXXX.sv_long.csv (when chunked)
    where {group} = study or study__batch (sanitized)

Usage:
  python epa-ng_input_splitter.py --svl svl.csv --fasta all_svs.fasta --out out
  python epa-ng_input_splitter.py --svl svl.csv --fasta all_svs.fasta --out out --chunk-size 5000
  # With batch column named "batch":
  python epa-ng_input_splitter.py --svl svl.csv --fasta all_svs.fasta --out out --batch-col batch
"""

import argparse
import os
import re
import sys
from collections import OrderedDict, defaultdict

import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser(description="Create per-study FASTAs and a dedup CSV with sequences.")
    ap.add_argument("--svl", required=True, help="Path to svl.csv (columns: specimen,study,sv,count)")
    ap.add_argument("--fasta", required=True, help="FASTA with all SV sequences (headers match 'sv')")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--chunk-size", type=int, default=0,
                    help="If >0, split each study FASTA into parts with this many sequences per file.")
    ap.add_argument("--study-col", default="study", help="Column name for study (default: 'study')")
    ap.add_argument("--batch-col", default=None, help="Optional column name for batch/cohort (default: None)")
    ap.add_argument("--sv-col", default="sv", help="Column name for SV ID (default: 'sv')")
    ap.add_argument("--sep", default=",", help="CSV delimiter for svl.csv (default: ',')")
    return ap.parse_args()


def read_fasta_to_dict(fp):
    """
    Simple FASTA parser: returns OrderedDict{id: sequence}, preserves file order.
    Header parsing: takes the first whitespace-delimited token after '>'.
    """
    seqs = OrderedDict()
    current_id = None
    buf = []
    with open(fp, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                # save previous
                if current_id is not None:
                    seqs[current_id] = "".join(buf)
                # parse new id
                header = line[1:].strip()
                sv_id = header.split()[0]
                current_id = sv_id
                buf = []
            else:
                buf.append(line)
        # flush last
        if current_id is not None:
            seqs[current_id] = "".join(buf)
    return seqs


def sanitize_name(name):
    """Safe-ish filename from study name."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name))


def write_fasta(out_path, id_list, id_to_seq):
    missing = 0
    with open(out_path, "wt") as out_h:
        for sv in id_list:
            seq = id_to_seq.get(sv)
            if not seq:
                missing += 1
                continue
            out_h.write(f">{sv}\n")
            # Wrap lines to 80 chars for readability (optional)
            for i in range(0, len(seq), 80):
                out_h.write(seq[i:i+80] + "\n")
    return missing


def chunk_iterable(iterable, size):
    """Yield chunks (lists) of given size from iterable."""
    chunk = []
    for x in iterable:
        chunk.append(x)
        if len(chunk) >= size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    out_align_dir = os.path.join(args.out, "study_fastas")
    os.makedirs(out_align_dir, exist_ok=True)
    out_svl_dir = os.path.join(args.out, "svl_groups")
    os.makedirs(out_svl_dir, exist_ok=True)

    # Load data
    df = pd.read_csv(args.svl, sep=args.sep, dtype=str)
    # ensure expected cols; 'sv' is required, 'study' is optional
    if args.sv_col not in df.columns:
        print(f"[ERROR] Missing required column in svl: '{args.sv_col}'", file=sys.stderr)
        sys.exit(2)
    if args.study_col not in df.columns:
        # synthesize a single cohort to enable chunking without study
        df[args.study_col] = 'all'
    # handle optional batch column
    use_batch = False
    if args.batch_col:
        if args.batch_col in df.columns:
            df[args.batch_col] = df[args.batch_col].fillna('NA')
            use_batch = True
        else:
            print(f"[WARN] batch column '{args.batch_col}' not found; proceeding without batch grouping.", file=sys.stderr)
            args.batch_col = None

    # Keep only first occurrence of each SV (global)
    df_first = df.drop_duplicates(subset=[args.sv_col], keep="first").copy()

    # Load FASTA into dict
    id_to_seq = read_fasta_to_dict(args.fasta)

    # Add sequence column to dedupbed df
    df_first["sequence"] = df_first[args.sv_col].map(id_to_seq.get)

    # Report missing sequences (present in svl but not in FASTA)
    missing_mask = df_first["sequence"].isna()
    n_missing = int(missing_mask.sum())
    if n_missing > 0:
        print(f"[WARN] {n_missing} SV(s) in svl not found in FASTA. They will be omitted where needed.", file=sys.stderr)

    # Write deduped CSV with sequence
    dedup_csv = os.path.join(args.out, "dedup_svl_with_seq.csv")
    df_first.to_csv(dedup_csv, index=False)
    print(f"[OK] Wrote dedup CSV with sequences: {dedup_csv}")

    # Build per-group SV lists from the ORIGINAL df (not dedupbed)
    group_cols = [args.study_col] + ([args.batch_col] if args.batch_col else [])
    group_to_svs = defaultdict(list)
    for _, row in df.iterrows():
        key = tuple(row[c] for c in group_cols)
        sv = row[args.sv_col]
        group_to_svs[key].append(sv)

    def key_to_stem(key_tuple):
        parts = []
        for col, val in zip(group_cols, key_tuple):
            parts.append(f"{col}-{sanitize_name(val)}")
        return "__".join(parts) if parts else "all"

    # Create per-group FASTAs and paired sv_long splits (optionally chunked)
    summary = []
    for key, svs in group_to_svs.items():
        # Keep first occurrence order within group to reduce duplicates in output
        seen = set()
        ordered_unique_svs = []
        for s in svs:
            if s not in seen:
                seen.add(s)
                ordered_unique_svs.append(s)

        stem = key_to_stem(key)
        # Build boolean mask across all grouping columns once
        mask = True
        for col, val in zip(group_cols, key):
            mask = mask & (df[col] == val)
        df_sub = df.loc[mask].copy()

        if args.chunk_size and args.chunk_size > 0:
            # Write chunked FASTAs and matching per-chunk sv_long CSVs
            total_missing = 0
            part_idx = 0
            for part_idx, chunk in enumerate(chunk_iterable(ordered_unique_svs, args.chunk_size), start=1):
                # FASTA per chunk
                out_fp = os.path.join(out_align_dir, f"{stem}.part{part_idx:03d}.fasta")
                missing = write_fasta(out_fp, chunk, id_to_seq)
                total_missing += missing
                # sv_long per chunk (only rows for SVs in this chunk)
                chunk_svl = df_sub[df_sub[args.sv_col].isin(chunk)].copy()
                chunk_svl_path = os.path.join(out_svl_dir, f"{stem}.part{part_idx:03d}.sv_long.csv")
                chunk_svl.to_csv(chunk_svl_path, index=False)
            summary.append((stem, len(ordered_unique_svs), part_idx, total_missing))
            print(f"[OK] {stem}: wrote {part_idx} chunk(s), {len(ordered_unique_svs)} SVs, missing seqs: {total_missing}; sv_long: per-chunk in {out_svl_dir}")
        else:
            # Single FASTA and single sv_long per group
            out_fp = os.path.join(out_align_dir, f"{stem}.fasta")
            missing = write_fasta(out_fp, ordered_unique_svs, id_to_seq)
            svl_path = os.path.join(out_svl_dir, f"{stem}.sv_long.csv")
            df_sub.to_csv(svl_path, index=False)
            summary.append((stem, len(ordered_unique_svs), 1, missing))
            print(f"[OK] {stem}: wrote 1 file, {len(ordered_unique_svs)} SVs, missing seqs: {missing}; sv_long: {svl_path}")

    # Write a small summary TSV
    summary_tsv = os.path.join(args.out, "study_fasta_summary.tsv")
    with open(summary_tsv, "wt") as h:
        h.write("group\tn_svs\tfiles_written\tmissing_seqs\n")
        for group, n_svs, files_written, missing in summary:
            h.write(f"{group}\t{n_svs}\t{files_written}\t{missing}\n")
    print(f"[OK] Summary: {summary_tsv}")


if __name__ == "__main__":
    main()
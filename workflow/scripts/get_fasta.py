#!/usr/bin/env python3

import subprocess
import polars as pl
import os
import shutil
import argparse
import sys

def subset_fasta_mmseqs(id_file, repdb_mmseqs, output_fasta):
    # 1. Prepare workspace
    out_dir = os.path.dirname(os.path.abspath(output_fasta))
    tmp_dir = os.path.join(out_dir, "mmseqs_tmp")

    os.makedirs(tmp_dir, exist_ok=True)

    sub_db_path = os.path.join(tmp_dir, "subset_db")
    ids_to_extract = os.path.join(tmp_dir, "ids_to_extract.txt")

    lookup_file = f"{repdb_mmseqs}.lookup"
    
    if not os.path.exists(lookup_file):
        print(f"Error: Lookup file {lookup_file} not found. Is '{repdb_mmseqs}' a valid MMseqs DB prefix?")
        sys.exit(1)
    
    try:
        with open(id_file, 'r') as f:
            sequence_names = set(line.strip() for line in f)

        # 2. Extract specific IDs using Polars
        # Adjust logic to match your specific TSV structure
        data_df = pl.read_csv(lookup_file, separator='\t', has_header=False, 
                              new_columns=["id", "sequence_name", "other_column"])
        out_df = data_df.filter(pl.col("sequence_name").is_in(sequence_names))        
        # In this example, we assume we want the 'id' column for the MMseqs subdb
        out_df.select("id").write_csv(ids_to_extract, include_header=False)

        # 3. Create the subset database
        # --subdb-mode 1: Create a subdb with actual sequence data
        # --id-mode 0: Match by the entry names/accessions provided in your list
        subprocess.run([
            "mmseqs", "createsubdb", 
            "--subdb-mode", "1", 
            "--id-mode", "0",
            ids_to_extract, 
            repdb_mmseqs, 
            sub_db_path
        ], check=True)

        # 4. Convert back to FASTA
        subprocess.run(["mmseqs", "convert2fasta", sub_db_path, output_fasta], check=True)

        print(f"Successfully created subset FASTA: {output_fasta}")

    finally:
        # 5. Cleanup intermediate MMseqs binary files
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)

def main():
    parser = argparse.ArgumentParser(
        description="Subset an MMseqs2 database into a FASTA file using its lookup table."
    )
    parser.add_argument(
        "id_file", 
        help="File with the ids you want to extract"
    )
    
    parser.add_argument(
        "output_fasta", 
        help="Path where the resulting subset FASTA will be saved"
    )
    parser.add_argument(
        "repdb_mmseqs", 
        help="Path/prefix to the MMseqs database (e.g., 'results/db_mmseqs')"
    )

    args = parser.parse_args()

    # We don't use the file_exists helper on repdb_mmseqs directly 
    # because it's a prefix, not necessarily a single file. 
    # We validate the lookup file inside the function instead.

    subset_fasta_mmseqs(args.id_file, args.repdb_mmseqs, args.output_fasta)

if __name__ == "__main__":
    main()

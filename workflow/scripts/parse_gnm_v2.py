import string
import gzip
import sys # For flushing output

# Define the translation table once at the module level
# J, B, O, U, Z are replaced by X
# * (stop codon) is deleted (mapped to None)
_CLEANING_TABLE = str.maketrans('JBOZU', 'XXXXX', '*')

def clean_seq(seq: str) -> str:
    """Uses a single str.translate operation for maximum speed."""
    return seq.translate(_CLEANING_TABLE)

def open_file_stream(filename: str, mode: str = 'rt'):
    """Transparently handles GZIP and plain text files."""
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    return open(filename, mode)

def read_fasta(fafile: str) -> dict[str, str]:
    """Reads FASTA into a dictionary {header: clean_sequence}."""
    sequences = {}
    current_header = None
    seq_lines = []

    with open_file_stream(fafile) as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                # Process the previous sequence
                if current_header is not None:
                    sequences[current_header] = clean_seq("".join(seq_lines))
                    seq_lines = []

                # Store the new header (only the part before the first space/tab)
                header_parts = line[1:].split(maxsplit=1)[0].split('\t', maxsplit=1)
                current_header = header_parts[0]
            else:
                # Accumulate sequence segments
                seq_lines.append(line)

        # Handle the last sequence in the file
        if current_header is not None and seq_lines:
            sequences[current_header] = clean_seq("".join(seq_lines))

    return sequences


ALPH = string.ascii_uppercase # A-Z
ALPH_LEN = len(ALPH)

def generate_protid_generator(total_seqs: int):
    """
    A generator to produce the ProtID segment (e.g., A A 0001) efficiently.
    """
    k_padding_width = len(str(total_seqs))
    k = 0
    i = 0
    j = 0
    
    while True:
        # e.g., 'A', 'A', '0001'
        protid = ALPH[i] + ALPH[j] + str(k).zfill(k_padding_width)
        yield protid
        
        # Advance the counters
        k += 1
        # Simplified, robust iteration over the two letters
        if k >= total_seqs: # Stop condition if needed, otherwise just let the caller handle it
            pass 
        if j >= ALPH_LEN:
            j = 0
            i += 1
            if i >= ALPH_LEN:
                i = 0
        if k % total_seqs == 0: # Simple way to advance the letters after a full run (if the alphabet is used as part of the total count)
            j += 1


def write_fasta(seqs: dict[str, str], renaming_data: dict, seqlen: int = 60):
    """
    Writes sequences with new IDs to FASTA and MAP files.
    """
    ofilenm = renaming_data['output_fasta']
    mapfile = renaming_data['output_map']
    filename_code = renaming_data['file_code']
    taxid = renaming_data.get('taxid')
    taxid_dict = renaming_data.get('taxid_dict') # For virus/special cases
    
    fasta_buffer = []
    map_buffer = []

    protid_gen = generate_protid_generator(len(seqs))

    for original_id, seq_data in seqs.items():
        protid = next(protid_gen)
        primary_original_id = original_id.split(' ')[0]
        
        # --- ID Logic ---
        if taxid:
            new_id = f'{taxid}_{filename_code}_{protid}'
        elif taxid_dict:
            # Use a dictionary lookup for taxid based on the sequence ID
            tax_id_lookup = taxid_dict.get(primary_original_id, 'UNKNOWN')
            new_id = f'{tax_id_lookup}_{filename_code}_{protid}'
        else:
            # Fallback if no renaming metadata is provided
            new_id = original_id 

        # --- Buffering Output ---
        fasta_buffer.append(f'>{new_id}\n')
        
        # Sequence splitting with a list comprehension
        fasta_buffer.extend([
            seq_data[l:l+seqlen] + '\n' 
            for l in range(0, len(seq_data), seqlen)
        ])
        
        # Map file line
        map_buffer.append(f'{primary_original_id}\t{new_id}\n')
        
    # --- Final Writing (Append Mode) ---
    # Write FASTA (gzipped)
    if ofilenm:
        with gzip.open(ofilenm, "at") as ofile:
            ofile.write("".join(fasta_buffer))

    # Write Map File (plain text)
    if mapfile:
        with open(mapfile, 'a') as omapfile:
            omapfile.write("".join(map_buffer))


def main_process(table_file: str, taxid_map_file: str, output_fasta: str, output_map: str):
    """The main entry point, processing the table of files."""

    # 1. Load the Tax ID map ONCE
    taxids = {}
    with open(taxid_map_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                taxids[parts[0]] = parts[1]
    
    # 2. Clear output files ONCE
    open(output_fasta, 'w').close()
    open(output_map, 'w').close()

    # 3. Process the input table
    with open(table_file) as table:
        for line in table:
            parts = line.split()
            if len(parts) < 2:
                continue
                
            input_fasta_file = parts[0]
            file_code = parts[1]
            
            print(f'Processing: {input_fasta_file}', file=sys.stderr, flush=True) # Use stderr for logging

            # A. Read sequences
            sequences = read_fasta(input_fasta_file)
            
            # B. Prepare renaming metadata
            renaming_data = {
                'output_fasta': output_fasta,
                'output_map': output_map,
                'file_code': file_code
            }

            if file_code == "C-RVDB":
                # Special case: Tax ID lookup based on sequence header
                renaming_data['taxid_dict'] = taxids
            else:
                # Standard case: Fixed Tax ID for the entire file
                tax_id_for_file = taxids.get(file_code)
                if tax_id_for_file:
                    renaming_data['taxid'] = tax_id_for_file
                else:
                    print(f"Warning: Tax ID not found for file code '{file_code}'. Skipping renaming.", file=sys.stderr)
                    # Could skip the file or write without renaming here

            # C. Write results (handles all buffering and I/O)
            write_fasta(sequences, renaming_data)

# Example usage (would replace snakemake calls)
# if __name__ == '__main__':
#     # Hardcoded or command line arguments would replace the snakemake access
#     main_process(
#         table_file='data/meta/broaddb_genome_table.tsv',
#         taxid_map_file='data/taxdump/broaddb_taxdump/taxid.map',
#         output_fasta='test/output.fa.gz',
#         output_map='test/id.map'
#     )

#!/usr/bin/env python3

# import argparse
import string
import gzip
import concurrent.futures

alph = string.ascii_uppercase

def clean_seq(seq):
    """Replace ambiguous amino acids and remove stop codons."""
    return seq.replace('J', 'X').replace('B', 'X').replace('O', 'X') \
              .replace('U', 'X').replace('Z', 'X').replace('*', '')

def openfile(filename):
    """Open plain or gzipped FASTA file."""
    return gzip.open(filename, 'rt') if filename.endswith('.gz') else open(filename, 'r')


def read_fasta(fafile):
    """Read sequences from a FASTA file into a dict."""
    seqs = {}
    seqid = ''
    s = []
    with openfile(fafile) as infile:
        for line in infile:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if seqid:
                    seqs[seqid] = clean_seq(''.join(s))
                seqid = line.split(' ')[0].split('\t')[0][1:]
                s = []
            else:
                s.append(line)
        if seqid and s:
            seqs[seqid] = clean_seq(''.join(s))
    return seqs


def write_fasta(seqs, rename=False, taxid=None, taxiddic = None, virus=False, filename=None,
                seqlen=60, ofilenm=None, append=True, mapfile=None, ogfile=None):
    ostr = ''
    if rename:
        i = 0
        j = 0
        k = 0
    dict_ids = {}
    for seq in seqs:
        if rename:
            if taxid is not None:
                protid = alph[i] + alph[j] + str(k).zfill(len(str(len(seqs))))
                new_id = ('%s_%s_%s' % (taxid, filename, protid))
            elif virus and taxiddic is not None:
                protid = alph[i] + alph[j] + str(k).zfill(len(str(len(seqs))))
                virus = seq.split(' ')[0]
                new_id = ('%s_%s_%s' % (taxiddic[virus], filename, protid))
            ostr += ('>%s\n' % (new_id))
            dict_ids[seq] = new_id


            j += 1
            if i > len(alph) - 1:
                i = 0
                if j > len(alph) - 1:
                    j = 0
            elif j > len(alph) - 1:
                j = 0
                i += 1
                if i > len(alph) - 1:
                    i = 0
            k += 1

        else:
            ostr += ('>%s\n' % seq)
        for l in range(0, len(seqs[seq]), seqlen):
            ostr += seqs[seq][l:l+seqlen] + '\n'
    
    if ofilenm is not None:
        with gzip.open(ofilenm, "at") as ofile:
            ofile.write(ostr)
            ofile.close()
        #     if append:
        #     ofile = open(ofilenm, 'a')
        # else:
        #     ofile = open(ofilenm, 'w')

    if mapfile is not None:
        with open(mapfile, 'a') as omapfile:
            for seq in dict_ids:
                id_string = ('%s\t%s\n') % (seq.split(' ')[0], dict_ids[seq])
                omapfile.write(id_string)

    return(ostr)


# TaxID, filename, prot_id (internal)

def process_line(line, taxids, snakemake_output0, snakemake_output1):
    file, filenm = line.split()[0], line.split()[1]
    print(f'Processing: {file}', flush=True)
    seqs = read_fasta(file)
    write_fasta(seqs, True,
                taxid=taxids[filenm],
                filename=filenm, ofilenm=snakemake_output0, append=True,
                mapfile=snakemake_output1, ogfile=file)

if __name__ == '__main__':
    # args = parser.parse_args()
    mapfile = f"{snakemake.input[1]}/taxid.map"
    taxids = {x.split('\t')[0]: x.split('\t')[1].replace('\n', '') for x in open(mapfile)}
    num_threads = snakemake.threads
    open(snakemake.output[0], 'w').close()
    open(snakemake.output[1], 'w').close()
    lines = []
    with open(snakemake.input[0]) as table:
        for line in table:
            lines.append(line)
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(process_line, line, taxids, snakemake.output[0], snakemake.output[1]) for line in lines]
        concurrent.futures.wait(futures)

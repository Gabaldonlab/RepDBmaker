#!/usr/bin/env python3

# import argparse
import string
import gzip

alph = string.ascii_uppercase

def clean_seq(seq):
    oseq = seq.replace('J', 'X').replace('B', 'X').replace('O', 'X')
    oseq = oseq.replace('U', 'X').replace('Z', 'X').replace('*', '')

    return oseq

def openfile(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt') 
    else:
        return open(filename, 'r')

def read_fasta(fafile):
    seqs = {}
    seqid = ''

    with openfile(fafile) as infile:
        for line in infile:
            line = line.replace('\n', '')
            if '>' in line:
                if seqid == '':
                    seqid = line.split(' ')[0].split('\t')[0].replace('>', '')
                    s = []
                else:
                    seqs[seqid] = clean_seq(''.join(s))
                    seqid = line.split(' ')[0].split('\t')[0].replace('>', '')
                    s = []
            else:
                s.append(line)
        if len(s) > 0:
            seqs[seqid] = clean_seq(''.join(s))

    return(seqs)

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
        if append:
            ofile = open(ofilenm, 'a')
        else:
            ofile = open(ofilenm, 'w')
        ofile.write(ostr)
        ofile.close()

    if mapfile is not None:
        with open(mapfile, 'a') as omapfile:
            for seq in dict_ids:
                id_string = ('%s\t%s\n') % (seq.split(' ')[0], dict_ids[seq])
                omapfile.write(id_string)

    return(ostr)


# TaxID, filename, prot_id (internal)



# Parse data
# parser = argparse.ArgumentParser(
#     description="Create a fasta file for BroadDB with coded ids starting from the input table"
# )

# parser.add_argument("-i", "--input", default="None", 
#                     help="input table",required=True)
# parser.add_argument("-t", "--taxid", default="None", 
#                     help="taxid map",required=True)
# parser.add_argument("-o", "--out", default="None", 
#                     help="output fasta fil",required=True)
# parser.add_argument("-m", "--map", default="None", 
#                     help="output map file",required=True)

# input_table = 'data/meta/broaddb_genome_table.tsv'
# broaddb_map = 'data/taxdump/broaddb_taxdump/taxid.map'
# ofilenm = 'test/test.fa'
# ofilemap = 'test/id.map'

if __name__ == '__main__':
    # args = parser.parse_args()
    mapfile = snakemake.input[1]+"/taxid.map"
    taxids = {x.split('\t')[0]: x.split('\t')[1].replace('\n', '') for x in open(mapfile)}

    open(snakemake.output[0], 'w').close()
    open(snakemake.output[1], 'w').close()

    with open(snakemake.input[0]) as table:
        for line in table:
            file = line.split()[0]
            filenm = line.split()[1]
            print('Processing: %s' % file, flush=True)
            seqs = read_fasta(file)
            if filenm == "C-RVDB":
                write_fasta(seqs, True, 
                            taxiddic = taxids, virus=True,
                            filename=filenm, ofilenm=snakemake.output[0], append=True, 
                            mapfile=snakemake.output[1], ogfile=file)
            else:
                write_fasta(seqs, True, 
                            taxid=taxids[filenm],
                            filename=filenm, ofilenm=snakemake.output[0], append=True, 
                            mapfile=snakemake.output[1], ogfile=file)

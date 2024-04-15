#!/usr/bin/env python3

import pandas as pd
import argparse
import subprocess as sp

def parse_args():
    parser=argparse.ArgumentParser(description="Get a new taxdump including new close proteomes")
    parser.add_argument('-i', '--input', dest='input', required=True, type=str,
                    help='Input file: mnemo\ttaxonomy'),
    parser.add_argument('-o', '--outdir', dest='out', required=True, type=str,
                    help='output taxdump directory'),
    parser.add_argument('-m', '--map', dest='map', required=True, type=str,
                    help='output taxid map file'),
    args=parser.parse_args()
    return args
# taxdump_edit.pl -names names.dmp -nodes nodes.dmp -taxa NAME -parent XXX -rank NAME -division X

keydict = {'k': 'superkingdom',
           'p': 'phylum', 'c': 'class','o': 'order',
           'f': 'family', 'g': 'genus', 's': 'species'}

ranks = list(keydict.keys())


if __name__ == '__main__':

    inputs=parse_args()

    seen = {}
    taxid_map = {}
    cmdl = list()

    input_data = pd.read_csv(inputs.input, sep='\t', header=None, names=["mnemo"] + ranks)

    # Loop trough species
    for i, row in input_data.iterrows():
        # Loop through clades from bigger to smaller
        for i, col in enumerate(input_data.columns[1:]):
            name = row.loc[col]
            rank = keydict[ranks[i]]
            # this is because in some case two consecutive clades are the same
            seen_clade = name+"_"+rank
            # this is to reduce taxonkit calls
            if seen_clade not in seen:
                # look for clade in old taxdump
                command = ("echo %s | taxonkit name2taxid -r --data-dir %s | cut -f2,3 | grep -w %s" %
                          (name , inputs.out, keydict[col]))
                result = sp.run(command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
                clade = result.stdout.strip()
                # if it's already there save it
                if clade:
                    taxid = int(clade.split("\t")[0])
                    seen[seen_clade] = taxid
                    print("found %s in taxdump with %s" % (name, taxid))
                    # if it's a species you'll need it for the taxid map
                    if col == "s":
                        taxid_map[row.loc["mnemo"]] = taxid
                else:
                    # The parent should ideally always be there as there are no new superkingdoms
                    # and we go from bigger to smaller
                    rank = keydict[ranks[i-2]]
                    seen_clade_parent = row.iloc[i-1]+"_"+rank
                    parent = seen[seen_clade_parent]
                    cmd = ('perl scripts/taxdump_edit.pl -names %s/names.dmp -nodes %s/nodes.dmp '
                    '-taxa \"%s\" -parent %s '
                    '-rank %s -division 8' % (inputs.out, inputs.out, name, parent, rank))
                    cmdl.append(cmd)
                    result = sp.run(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
                    taxid = result.stdout.split("\n")[0].split("=")[1].split(".")[0].strip()
                    print("%s not found, added with taxid:%s" % (name, taxid))

                    seen[seen_clade] = taxid

                    if col == "s":
                        taxid_map[row.loc["mnemo"]] = int(taxid)
            elif col == "s":
                taxid_map[row.loc["mnemo"]] = seen[seen_clade]
    print(seen)
    outmap = pd.DataFrame(list(taxid_map.items()))
    outmap.to_csv(inputs.map, sep='\t', index=False, header=None)


    # cmdl = list(set(cmdl))
    # for item in cmdl:
    #     print(item)

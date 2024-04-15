#!/usr/bin/env python

from pathlib import Path

taxid_mapping_dict = {}

count = 1

with open(Path("test/broaddb").joinpath("names.dmp")) as fin:
    for line in fin:
        taxid = line.split("\t")[0]
        taxid_mapping_dict[taxid] = str(count)
        count += 1

with open(Path("test/broaddb").joinpath("names.dmp")) as fin, open(
    Path("test/broaddb").joinpath("names.dmp.new"), "w"
) as fout:
    for line in fin:
        line = line.strip().split("\t")
        line[0] = taxid_mapping_dict[line[0]]
        line = "\t".join(line)
        fout.write(f"{line}\n")

with open(Path("test/broaddb").joinpath("nodes.dmp")) as fin, open(
    Path("test/broaddb").joinpath("nodes.dmp.new"), "w"
) as fout:
    for line in fin:
        line = line.strip().split("\t")
        line[0] = taxid_mapping_dict[line[0]]
        line[2] = taxid_mapping_dict[line[2]]
        line = "\t".join(line)
        fout.write(f"{line}\n")

Path("test/broaddb").joinpath("names.dmp").unlink()
Path("test/broaddb").joinpath("nodes.dmp").unlink()
Path("test/broaddb").joinpath("nodes.dmp.new").rename("test/broaddb/nodes.dmp")
Path("test/broaddb").joinpath("names.dmp.new").rename("test/broaddb/names.dmp")
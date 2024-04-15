#!/usr/bin/env python3

'''
Proteome renaming
Copyright (C) 2023  Moisès Bernabeu

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

import logging
import sys
from utils import *

ranks = {'d': 'domain', 'p': 'phylum',
         'c': 'class', 'o': 'order', 'f': 'family',
         'g': 'genus', 's': 'species'}


def read_lineage(lineage, logging):
    # Genome map dictionary
    gnm_map = {}

    # Parents dictionary
    taxa_parent = {v: {} for k, v in ranks.items()}

    for line in open(lineage, 'r'):
        try:
            gnm, lineage = line.strip().split('\t')
        except ValueError:
            logging.error(f'Incorrect line: {line}')
        
        # The first node is the root
        parent = 'root'

        # Parse taxa in order
        for i, taxon in enumerate(lineage.split(';')):
            code, taxon = taxon.split('__')

            try:
                rank = ranks[code]
            except KeyError:
                logging.error(f'The key {code} for the rank does not exist: {line}.')
                sys.exit(1)
            
            if taxon in taxa_parent[rank]:
                if taxa_parent[rank][taxon] != parent:
                    logging.error(f'Inconsistent parent: {taxon} {rank} {line}.')
                    sys.exit(1)
            else:
                taxa_parent[rank][taxon] = parent
            
            # Update parent
            parent = taxon
        
        if rank != 'species':
            logging.error(f'Incomplete record: {line}.')
            sys.exit(1)
        # Here parent means species as it is the last one
        gnm_map[gnm] = parent

    return {'gnm_map': gnm_map, 'taxa_parent': taxa_parent}


def create_taxtree(taxa_parent, logging):
    taxtree = {'1': {'parent': '1', 'rank': 'no rank', 'name': 'root'}}

    # Current tax id
    cid = 1

    # Assign taxIDs to taxa
    parent_id = {'root': cid}
    species_id = {}
    for rank in ranks.values():
        taxon_id = {}
        for taxon, parent in taxa_parent[rank].items():
            cid += 1
            taxtree[str(cid)] = {'parent': parent_id[parent],
                                 'rank': rank,
                                 'name': taxon}
            taxon_id[taxon] = cid

            if rank == 'species':
                species_id[cid] = taxon

        parent_id = {k: v for k, v in taxon_id.items()}
    
    logging.debug('Taxa tree generated.')

    return {'taxtree': taxtree, 'parent_id': parent_id,
            'species_id': species_id}


def write_taxdump(gnm_map, taxtree, parent_id, outdir):
    # Generating the taxid.map file
    create_folder(f'{outdir}/taxdump/')
    with open(f'{outdir}/taxdump/taxid.map', 'w') as f:
        for g, species in sorted(gnm_map.items()):
            f.write(f'{g}\t{parent_id[species]}\n')
    logging.debug(f'{outdir}/taxdump/taxid.map created.')

    # Generating the nodes.dmp and names.dmp files
    nodes = open(f'{outdir}/taxdump/nodes.dmp', 'w')
    names = open(f'{outdir}/taxdump/names.dmp', 'w')

    for id_, rec in sorted(taxtree.items(), key=lambda x: int(x[0])):
        print(id_, rec)
        nodes.write(f"{id_}\t|\t{rec['parent']}\t|\t{rec['rank']}\t|\n")
        names.write(f"{id_}\t|\t{rec['name']}\t|\t\t|\tscientific name\t|\n")
    
    logging.debug(f'{outdir}/taxdump/nodes.dmp and names.dmp created.')

    return 0


def read_taxdump(nodes, names, taxmap, logging):
    names_id = {}
    for line in open(names, 'r'):
        line = line.replace('|', '').strip().replace('\t\t', '\t').split('\t')
        names_id[line[0]] = line[1]

    taxtree = {}
    taxa_parent = {v: {} for k, v in ranks.items()}
    for line in open(nodes, 'r'):
        taxid, parent, rank = line.replace('\t|', '').strip().split('\t')
        taxtree[taxid] = {'parent': parent,
                          'rank': rank,
                          'name': names_id[taxid]}
        if rank != 'no rank':
            taxa_parent[rank][names_id[taxid]] = names_id[parent]
    
    gnm_map = {}
    for line in open(taxmap):
        gnm, taxid = line.split()
        gnm_map[gnm] = taxtree[taxid]['name']
    
    logging.debug('Taxdump read.')

    return {'taxtree': taxtree, 'taxa_parent': taxa_parent, 'gnm_map': gnm_map}


def check_lineage(lineage, nodes, names, taxmap, logging):
    lng = read_lineage(lineage, logging)
    txd = read_taxdump(nodes, names, taxmap, logging)

    lng_taxa_parent = lng['taxa_parent']
    txd_taxa_parent = txd['taxa_parent']
    
    for rank in ranks.values():
        for taxa, parent in lng_taxa_parent[rank].items():
            try:
                txd_parent = txd_taxa_parent[rank][taxa]
                if txd_parent == parent:
                    print(f'{taxa} has the same parent.')
                else:
                    print(f'{taxa} does not have the same parent.')
            except KeyError:
                logging.error(f'{taxa} not in the taxdump')

    return 0


def update_taxdump():

    return 0


def main():
    odir = '../outputs/'
    prefix = 'bdb'
    update = True

    logfile = f'{odir}/{prefix}.log'


    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[logging.FileHandler(logfile), logging.StreamHandler()]
    )

    # lng = read_lineage('../data/lng.txt', logging)
    # txtr = create_taxtree(lng['taxa_parent'])
    # write_taxdump(lng['gnm_map'], txtr['taxtree'], txtr['parent_id'], '../outputs/')
    # print(txtr['species_id'])

    read_taxdump('../outputs/taxdump/nodes.dmp',
                 '../outputs/taxdump/names.dmp',
                 '../outputs/taxdump/taxid.map', logging)

    check_lineage('../data/lng.txt',
                  '../outputs/taxdump/nodes.dmp',
                  '../outputs/taxdump/names.dmp',
                  '../outputs/taxdump/taxid.map', logging)

    return 0


if __name__ == '__main__':
    main()

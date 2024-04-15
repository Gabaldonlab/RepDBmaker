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

from utils import *
from glob import glob
import argparse as ap
import pandas as pd
import logging


def main():
    args = ap.ArgumentParser()
    args.add_argument('-i', '--input', dest='input', required=True,
                      help=('Input genomes table.'),
                      metavar='<path/to/table>')
    args.add_argument('-o', '--output', dest='odir',
                      help='Output directory',
                      metavar='<path/to/output>',
                      default='outputs')
    args.add_argument('-p', '--prefix', dest='prefix',
                      help='Output prefix',
                      metavar='<STR>',
                      default='output')
    args = args.parse_args()

    table = args.input
    odir = args.odir
    prefix = args.prefix
    logfile = '%s/%s.log' % (odir, prefix)

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[logging.FileHandler(logfile), logging.StreamHandler()]
    )

    ofastadir = ('%s/%s_proteomes/' % (odir, prefix)).replace('//', '/')
    corresdir = ('%s/%s_correspondences/' % (odir, prefix)).replace('//', '/')

    genomes = [basename(x, '.fa') for x in glob(ofastadir + '/*/*')]

    if all([line.split('\t')[1] in genomes for line in open(table, 'r') if line != '\n']):
        logging.error('All the genomes are already parsed or you used extant mnemonic.')
        sys.exit(0)

    create_folder(odir)
    create_folder(ofastadir)
    create_folder(corresdir)

    maxfolder = glob(ofastadir + '/*')
    if len(maxfolder) > 0:
        maxfolder.sort()
        maxfolder = maxfolder[-1]
        if len(glob(maxfolder + '/*')) < 500:
            bnm = basename(maxfolder)
            prefix = bnm.split('_', 1)[0]
            i = int(bnm.split('_', 1)[1])
            ofastasubdir = '%s/%s_%s' % (ofastadir, prefix, str(i).zfill(4))
            corressubdir = '%s/%s_%s' % (corresdir, prefix, str(i).zfill(4))
        else:
            bnm = basename(maxfolder)
            prefix = bnm.split('_', 1)[0]
            i = int(bnm.split('_', 1)[1]) + 1
            ofastasubdir = '%s/%s_%s' % (ofastadir, prefix, str(i).zfill(4))
            corressubdir = '%s/%s_%s' % (corresdir, prefix, str(i).zfill(4))
            create_folder(ofastasubdir)
            create_folder(corressubdir)
    else:
        ofastasubdir = '%s/%s_%s' % (ofastadir, prefix, '0000')
        corressubdir = '%s/%s_%s' % (corresdir, prefix, '0000')
        create_folder(ofastasubdir)
        create_folder(corressubdir)
        i = 0


    for line in open(table, 'r'):
        line = line.strip()
        if line != '':
            if len(glob(ofastasubdir + '/*')) < 500:
                line = line.strip().split()
                if line[1] in genomes:
                    logging.debug('%s already parsed' % line[1])
                    continue
                seqs = read_fasta(line[0])
                oseqsfile = '%s/%s.fa' % (ofastasubdir, line[1])
                rnmd = write_fasta(seqs, ofilenm=oseqsfile, rename=True,
                                   mnemo=line[1], taxid=line[2])

                logging.debug('%s parsed' % line[1])

                df = pd.DataFrame(rnmd['correspondence']).transpose()
                df.to_csv('%s/%s.tsv' % (corressubdir, line[1]), sep='\t', header=False)
                logging.debug('%s correspondence written' % line[1])
            else:
                line = line.strip().split()
                if line[1] in genomes:
                    logging.debug('%s already parsed' % line[1])
                    continue
                seqs = read_fasta(line[0])
                oseqsfile = '%s/%s.fa' % (ofastasubdir, line[1])
                rnmd = write_fasta(seqs, ofilenm=oseqsfile, rename=True,
                                   mnemo=line[1], taxid=line[2])
                
                logging.debug('%s parsed' % line[1])

                df = pd.DataFrame(rnmd['correspondence']).transpose()
                df.to_csv('%s/%s.tsv' % (corressubdir, line[1]), sep='\t', header=False)
                logging.debug('%s correspondence written' % line[1])

                i += 1
                ofastasubdir = '%s/%s_%s' % (ofastadir, prefix, str(i).zfill(4))
                corressubdir = '%s/%s_%s' % (corresdir, prefix, str(i).zfill(4))
                create_folder(ofastasubdir)
                create_folder(corressubdir)

    return 0


if __name__ == '__main__':
    main()
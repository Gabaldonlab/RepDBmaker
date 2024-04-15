#!/usr/bin/env python3

'''
Useful function for processing the BroadDB files
Copyright (C) 2023  Moisès Bernabeu & Giacomo Mutti

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

import subprocess as sp
import sys
import os
import string

alph = string.ascii_uppercase

def create_folder(name):
    '''
    This function checks the existance of a folder and if it is not created
    it creates

    Args:
        name (str): string with the folder name

    Returns:
        int: 0
    '''

    if not os.path.exists(name):
        try:
            os.mkdir(name)
        except Exception:
            print('Unable to create the directory: %s' % name)

    return 0


def read_fasta(fasta_file, delimiter=''):
    '''
    It reads a fasta file

    The function reads a plain text line by line in fasta format, if the line
    contains the character ">" first, writes the previous sequence in the
    dictionary, then restarts the name and sequence variables, then write the
    last sequence.

    Args:
        fasta_file (str): fasta file path
        delimiter (str): the delimiter in the fasta header, if it exists

    Returns:
        dict: dictionary with the sequences where the keys are the sequence
        name and the value is the string with the sequence
    '''

    # Initialising variables
    seqs = dict()
    name = ''
    s = list()

    for line in open(fasta_file):
        # Removing spaces before and after the line
        line = line.strip()

        if '>' in line:
            # Writing previous sequence in the dictionary with the older name
            if name != '':
                seqs[name] = ''.join(s)

            # Restarting the name variable with the new sequence name checking
            # the existence of a delimiter, if it exists, we will use the
            # string before the delimiter
            if delimiter == '':
                name = line.replace('>', '')
            else:
                name = line.replace('>', '').split(delimiter)[0]

            # Restarting the sequence list
            s = list()
        else:
            # Appending sequences lines to the sequence list
            s.append(line.upper())

    # Writing the last sequence where there are not new lines with ">"
    # here the sequences are cleaned of unusual characters in amino acid code
    if len(s) != 0:
        seqs[name] = clean_seq(''.join(s))

    return seqs


def write_fasta(seqs, ofilenm=None, rename=False,
                mnemo=None, taxid=None, seqlen=60):
    '''
    The function writes, in fasta format, a sequences dictionary.
    It allows the user to print the fasta formated sequences into a file and
    rename the sequences according to an alphanumeric code. 

    Args:
        seqs (_type_): sequences dictionary, key = header, value = sequence
        ofilenm (_type_, optional): output file path. Defaults to None.
        rename (bool, optional): whether to rename the proteins. Defaults to False.
        mnemo (_type_, optional): mnemonic for the sequences. Defaults to None.
        taxid (_type_, optional): tax id for the sequences. Defaults to None.
        seqlen (int, optional): width of the fasta string. Defaults to 60.

    Returns:
        _type_: _description_
    '''                
    ostr = ''
    if rename:
        i = 0
        j = 0
        k = 0
    
    if rename and (taxid is None or mnemo is None):
        print('Error, you did not specify the mnemo or the taxid.')
        sys.exit(1)

    corrdict = {}
    for seq in seqs:
        if rename:
            # Defining the protein id and appending it as header to the
            # fasta outfile
            protid = alph[i] + alph[j] + str(k).zfill(4)
            header = '%s_%s_%s' % (taxid, mnemo, protid)
            ostr += ('>%s\n' % header)
            corrdict[seq] = {'header': header}

            if k < 1000:
                k += 1
            else:
                k = 1
                if j < len(alph) - 1:
                    j += 1
                else:
                    j = 0
                    if i < len(alph) - 1:
                        i += 1
                    else:
                        i = 0
        else:
            ostr += ('>%s\n' % seq)
        
        # Appending the sequence according to the length to the fasta outfile
        for l in range(0, len(seqs[seq]), seqlen):
            ostr += seqs[seq][l:l+seqlen] + '\n'
    
    # Writing the output file
    if ofilenm is not None:
        ofile = open(ofilenm, 'w')
        ofile.write(ostr)
        ofile.close()

    return {'fasta': ostr, 'correspondence': corrdict}


def clean_seq(seq):
    '''
    Cleaning the sequnce of unknown character in the amino acid code

    Args:
        seq (str): string containing the sequence

    Returns:
        str: cleaned seq
    '''

    oseq = seq.replace('J', 'X').replace('B', 'X').replace('O', 'X')
    oseq = oseq.replace('U', 'X').replace('Z', 'X').replace('*', '')

    return oseq


def basename(file, extension=None):
    if extension is None:
        return file.rsplit('/', 1)[1]
    else:
        return file.rsplit('/', 1)[1].replace(extension, '')

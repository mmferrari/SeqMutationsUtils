#!/usr/bin/env python3

import argparse
import hashlib
import re

"""
Script to generate CODON mutations of a given DNA sequence.


Copyright 2020 Margherita Maria Ferrari.


This file is part of SeqMutationsUtils.

SeqMutationsUtils is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SeqMutationsUtils is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SeqMutationsUtils.  If not, see <http://www.gnu.org/licenses/>.
"""


class SeqMutationsUtils:
    CODON_FULL_NAME_DICT = {
        'ala': ('GCT', 'GCC', 'GCA', 'GCG'),
        'arg': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
        'asn': ('AAT', 'AAC'),
        'asp': ('GAT', 'GAC'),
        'cys': ('TGT', 'TGC'),
        'gln': ('CAA', 'CAG'),
        'glu': ('GAA', 'GAG'),
        'gly': ('GGT', 'GGC', 'GGA', 'GGG'),
        'his': ('CAT', 'CAC'),
        'ile': ('ATT', 'ATC', 'ATA'),
        'leu': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
        'lys': ('AAA', 'AAG'),
        'met': ('ATG',),
        'phe': ('TTT', 'TTC'),
        'pro': ('CCT', 'CCC', 'CCA', 'CCG'),
        'ser': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
        'stop': ('TAA', 'TAG', 'TGA'),
        'thr': ('ACT', 'ACC', 'ACA', 'ACG'),
        'trp': ('TGG',),
        'tyr': ('TAT', 'TAC'),
        'val': ('GTT', 'GTC', 'GTA', 'GTG')
    }

    CODON_SHORT_NAME_DICT = {
        'A': ('GCT', 'GCC', 'GCA', 'GCG'),
        'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
        'N': ('AAT', 'AAC'),
        'D': ('GAT', 'GAC'),
        'C': ('TGT', 'TGC'),
        'Q': ('CAA', 'CAG'),
        'E': ('GAA', 'GAG'),
        'G': ('GGT', 'GGC', 'GGA', 'GGG'),
        'H': ('CAT', 'CAC'),
        'I': ('ATT', 'ATC', 'ATA'),
        'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
        'K': ('AAA', 'AAG'),
        'M': ('ATG',),
        'F': ('TTT', 'TTC'),
        'P': ('CCT', 'CCC', 'CCA', 'CCG'),
        'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
        'STOP': ('TAA', 'TAG', 'TGA'),
        'T': ('ACT', 'ACC', 'ACA', 'ACG'),
        'W': ('TGG',),
        'Y': ('TAT', 'TAC'),
        'V': ('GTT', 'GTC', 'GTA', 'GTG')
    }

    @classmethod
    def get_args(cls):
        parser = argparse.ArgumentParser(description='Mutations in substring')
        parser.add_argument('-i', '--input-file', metavar='INPUT_FILE', type=str, required=True,
                            help='File containing the coding DNA sequence', default=None)
        parser.add_argument('-p', '--initial-pos', metavar='INITIAL_POSITION', type=int, required=True,
                            help='Initial position of substring to modify starting from 1', default=0)
        parser.add_argument('-m', '--mutations', metavar='MUTATIONS', type=str, required=False,
                            help='List of mutations in the form of \'(ala,glu,...),(val,...),...\''
                                 ' -- NOTE: the "\'" must be included', default='')
        parser.add_argument('-n', '--num', metavar='NUM', type=int, required=True,
                            help='Number of consecutive codons', default=1)
        parser.add_argument('-o', '--output-prefix', metavar='OUTPUT_PREFIX', type=str, required=True,
                            help='Output prefix file name', default='output_')
        parser.add_argument('-a', '--all-codons', required=False, action='store_true',
                            help='Try all codons')
        parser.add_argument('-v', '--inv', required=False, action='store_true',
                            help='Also consider combinations with original codons without mutations')
        parser.add_argument('-x', '--single-output-file', required=False, action='store_true',
                            help='Create single fasta file')

        return parser.parse_args()

    @classmethod
    def __generate_mutations(cls, seq, mut, idx, inv=False, single_file=False, out_prefix='', out_suffix=1,
                             mut_hashes=None):
        if mut_hashes is None:
            mut_hashes = list()

        if len(mut) < 1:
            seq_hash = hashlib.sha3_512(seq.encode('utf8')).hexdigest()

            if seq_hash in mut_hashes:
                return out_suffix

            mut_hashes.append(seq_hash)
            file_name = out_prefix if single_file else (out_prefix + str(out_suffix))
            open_mode = 'a' if (single_file and out_suffix > 1) else 'w'

            with open(file_name + '.fa', open_mode) as fout:
                fout.write('>' + out_prefix + str(out_suffix) + ' range=' + out_prefix + str(out_suffix) + ':1-' +
                           str(len(seq)) + ' 5\'pad = 0 3\'pad=0 strand=+ repeatMasking=none\n')
                fout.write(seq)

                if single_file:
                    fout.write('\n\n')

            return out_suffix + 1

        if inv:
            out_suffix = cls.__generate_mutations(seq, mut[1:], idx + 3, inv, single_file, out_prefix, out_suffix,
                                                  mut_hashes)

        for m in mut[0]:
            new_codons = cls.CODON_FULL_NAME_DICT.get(m.lower(), cls.CODON_SHORT_NAME_DICT.get(m.upper(), list()))

            if len(new_codons) < 1:
                raise AssertionError('Unknown codon name: ' + m)

            for c in new_codons:
                seq = seq[:idx] + c + seq[idx + 3:]
                out_suffix = cls.__generate_mutations(seq, mut[1:], idx + 3, inv, single_file, out_prefix, out_suffix,
                                                      mut_hashes)

        return out_suffix

    @classmethod
    def generate_mutations(cls, input_file, mutations, start_pos, codon_num, inv=False, single_output=False,
                           output_prefix=''):
        if not input_file:
            raise AssertionError('You must specify an input file')

        with open(input_file, 'r') as fin:
            seq = fin.readline().strip().upper()

        if start_pos % 3 != 1:
            raise AssertionError('Initial position must be the beginning of a codon')

        start_pos -= 1

        if codon_num < 1:
            raise AssertionError('Number of codons must be positive')
        if start_pos + codon_num * 3 > len(seq):
            raise AssertionError('Number of codons is too large')
        if len(mutations) != codon_num:
            raise AssertionError('Number of mutations do not match the number of codons')

        cls.__generate_mutations(seq, mutations, start_pos, inv, single_output, output_prefix)


if __name__ == '__main__':
    pattern = re.compile(r'\([\w\s,]+\)')
    args = vars(SeqMutationsUtils.get_args())
    init_pos = args.get('initial_pos', 0)
    num = args.get('num', 0)
    in_file = args.get('input_file', None)
    try_all_codons = args.get('all_codons', False)
    mut = list()

    for i in pattern.findall(args.get('mutations', '')):
        tmp = list()
        for j in i.strip('()').split(','):
            if len(j.strip()) > 0:
                tmp.append(j.strip())
        if len(tmp) > 0:
            mut.append(list(tmp))

    if try_all_codons:
        mut.clear()
        for i in range(0, num):
            mut.append(list(SeqMutationsUtils.CODON_FULL_NAME_DICT.keys()))

    SeqMutationsUtils.generate_mutations(in_file, mut, args.get('initial_pos', 0), num, args.get('inv', False),
                                         args.get('single_output_file', False),
                                         args.get('output_prefix', ''))

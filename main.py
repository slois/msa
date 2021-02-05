# This is a MSA script
import json

import Bio
from Bio import AlignIO, pairwise2, Align, SeqIO
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from plotnine import *
from datetime import datetime
import re

id_parser = re.compile('^(hCoV19|hCoV_19){1}_([A-Za-z]+)_', re.IGNORECASE)

class Probe(object):
    def __init__(self, id, start, end, seq, reference=None):
        self.id = id
        self.start = start
        self.end = end
        self.sequence = seq

        if reference:
            self.status = self.check(reference=reference)

    def check(self, reference):
        return self.sequence == reference[self.start-1:self.end].seq.upper()

    def __contains__(self, position):
        return self.start <= position <= self.end

    def __repr__(self):
        return self.id

    def __str__(self):
        return self.id





def col_f(col, alphabet=['A', 'T', 'C', 'G', '-']):
    total = len(col)
    f = {letter: col.count(letter)/total for letter in alphabet}
    return f


def shannon_entropy(col):
    f = col_f(col)
    return -np.sum([freq * np.log2(freq) for letter, freq in f.items() if freq > 0.0])


def score_msa(msa):
    l = msa.get_alignment_length()
    bit_score = [shannon_entropy(msa[:, col]) for col in range(l)]
    return bit_score


def extract_info_from_id(record_id):
    info = record_id.split('|')

    if len(info) >= 3:
        identifier = info[0]
        epi_code = info[1]
        date = info[2]
        continent = info[3] if len(info) > 3 else None
        m = id_parser.match(identifier)
        region = m.group(2)

        fmt = "%Y_%m_%d" if '_' in date else "%Y%m%d"
        try:
            date = datetime.strptime(date, fmt)
            date = datetime.strftime(date, "%Y%m%d")
        except ValueError:
            date = None
    else:
        identifier = info[0]
        epi_code = date = continent = region = None

    annotation = {
        'id': identifier,
        'epi': epi_code,
        'region': region,
        'continent': continent,
        'date': date}
    return annotation


def search_primers(aligner, record, probes):
    for probe in probes:
        try:
            aln = aligner.align(record.seq, probe.sequence)[0]
            primer_location = aln.aligned[0]
            for loc in primer_location:
                exact_location = FeatureLocation(loc[0], loc[1] - 1)
                record.features.append(SeqFeature(location=exact_location, strand=1, type='probe', id=probe.id))
        except ValueError:
            print("Error")


def main(msa_input, fmt):
    msa = AlignIO.read(msa_input, format=fmt)

    # Load reference from GenBank record
    reference_sequence = SeqIO.read('data/MN908947_sequence.gb', 'genbank')

    # Load probe locations
    probes_data = json.load(open('data/IP4.coordinates.release2.json', mode='r'))
    probes = [Probe(**p, reference=reference_sequence) for p in probes_data]
    assert all([p.status for p in probes])

    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.25
    # TODO: put defined alphabet instead of adhoc alphabet
    aligner.alphabet = ['A', 'T', 'C', 'G', '-', 'N', 'W', 'Y', 'M', 'R', 'V', 'S', 'H']

    for record in msa:
        annot = extract_info_from_id(record.id)
        record.annotations = annot
        search_primers(aligner=aligner, record=record, probes=probes)


    print("Hola")



    #dat = pd.DataFrame({'pos': range(l), 'bit_score': bit_score})

    #plt = (ggplot(dat) + aes(x='pos', y='bit_score') + geom_line(size=0.4))






    for loc in primer_location:
        msa_primer = msa[:, loc[0]:loc[1]]
        score = score_msa(msa_primer)
        print(score)



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main('data/GISAID_100.maff.fa', 'fasta')



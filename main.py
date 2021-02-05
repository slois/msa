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


class ProbeFinder(object):
    def __init__(self, probes, open_gap_score=-0.5, extend_gap_score=-0.25, **kwargs):
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'local'
        self.aligner.open_gap_score = open_gap_score
        self.aligner.extend_gap_score = extend_gap_score
        # TODO: put defined alphabet instead of adhoc alphabet
        self.aligner.alphabet = ['A', 'T', 'C', 'G', '-', 'N', 'W', 'Y', 'M', 'R', 'V', 'S', 'H']
        self.probes = probes

    def search(self, record):
        features = list()
        for probe in self.probes:
            try:
                aln = self.aligner.align(record.seq, probe.sequence)[0]
                primer_location = aln.aligned[0]
                for loc in primer_location:
                    exact_location = FeatureLocation(loc[0], loc[1] - 1)
                    features.append(SeqFeature(location=exact_location, strand=1, type='probe', id=probe.id))
            except ValueError:
                print("Error")
        return features

def main(msa_input, fmt):
    msa = SeqIO.parse(msa_input, format=fmt)

    # Load reference from GenBank record
    reference_sequence = SeqIO.read('data/MN908947_sequence.gb', 'genbank')

    # Load probe locations
    probes_data = json.load(open('data/IP4.coordinates.release2.json', mode='r'))
    probes = [Probe(**p, reference=reference_sequence) for p in probes_data]
    assert all([p.status for p in probes])

    pf = ProbeFinder(probes=probes)
    aln = Align.PairwiseAligner()
    aln.mode = 'global'
    aln.open_gap_score = -5
    aln.extend_gap_score = -2
    aln.alphabet = ['A', 'T', 'C', 'G', '-', 'N', 'W', 'Y', 'M', 'R', 'V', 'S', 'H']

    for record in msa:
        annot = extract_info_from_id(record.id)
        record.annotations = annot

        # Ungap only if records are previously aligned
        record_aln = aln.align(reference_sequence.seq, record.seq.ungap('-'))[0]

        print(record_aln)

        feat = pf.search(record=record)
        record.features += feat

        print(record.annotations)
        for feat in record.features:
            print(feat.extract(record.seq))




if __name__ == '__main__':
    main('data/GISAID_100.maff.fa', 'fasta')



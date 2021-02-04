# This is a MSA script
import Bio
from Bio import AlignIO, pairwise2
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from plotnine import *

def col_f(col, alphabet=['A', 'T', 'C', 'G', '-']):
    total = len(col)
    f = {letter: col.count(letter)/total for letter in alphabet}
    return f


def shannon_entropy(col):
    f = col_f(col)
    return -np.sum([freq * np.log2(freq) for letter, freq in f.items() if freq > 0.0])


def main(msa_input, fmt):
    msa = AlignIO.read(msa_input, format=fmt)

    l = msa.get_alignment_length()
    bit_score = [shannon_entropy(msa[:, col]) for col in range(l)]

    dat = pd.DataFrame({'pos': range(l), 'bit_score': bit_score})

    plt = (ggplot(dat) + aes(x='pos', y='bit_score') + geom_line(size=0.4))

    primer = Seq('AGAGTTTGATCCTGGCTCAGGGTGA')

    xx = pairwise2.align.localxx(msa[1].seq.ungap('-'), primer, one_alignment_only=True)

    print(xx[0].start, xx[0].end)
    print(msa[1].seq.ungap('-')[1374:1470])

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main('data/16SRNA_Deino_87seq.aln', 'fasta')



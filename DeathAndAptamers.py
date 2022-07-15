import ExpressionLocation
import Gil

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

import shrna

def DcR3_shRNA(record):
    #SAM! Put your DNA sequence to target here
    #I think this is the right one, it says DcR3 is a synonym for this https://www.ncbi.nlm.nih.gov/nuccore/BC034349.1 (excluding polyA tail)
    DcR3_SEQ = "CTAGTTCATCGCGAGCGGCCGCACAACTTTGTACAAGAAAGTTTCGGTCAGGCACAGCAGGGTCCTGTGTCCGCGCTGAGCCGCGCTCTCCCTGCTCCAGCAAGGACCATGAGGGCGCTGGAGGGGCCAGGCCTGTCGCTGCTGTGCCTGGTGTTGGCGCTGCCTGCCCTGCTGCCGGTGCCGGCTGTACGCGGAGTGGCAGAAACACCCACCTACCCCTGGCGGGACGCAGAGACAGGGGAGCGGCTGGTGTGCGCCCAGTGCCCCCCAGGCACCTTTGTGCAGCGGCCGTGCCGCCGAGACAGCCCCACGACGTGTGGCCCGTGTCCACCGCGCCACTACACGCAGTTCTGGAACTACCTAGAGCGCTGCCGCTACTGCAACGTCCTCTGCGGGGAGCGTGAGGAGGAGGCACGGGCTTGCCACGCCACCCACAACCGTGCCTGCCGCTGCCGCACCGGCTTCTTCGCGCACGCTGGTTTCTGCTTGGAGCACGCATCGTGTCCACCTGGTGCCGGCGTGATTGCCCCGGGCACCCCCAGCCAGAACACGCAGTGCCAGCCGTGCCCCCCAGGCACCTTCTCAGCCAGCAGCTCCAGCTCAGAGCAGTGCCAGCCCCACCGCAACTGCACGGCCCTGGGCCTGGCCCTCAATGTGCCAGGCTCTTCCTCCCATGACACCCTGTGCACCAGCTGCACTGGCTTCCCCCTCAGCACCAGGGTACCAGGAGCTGAGGAGTGTGAGCGTGCCGTCATCGACTTTGTGGCTTTCCAGGACATCTCCATCAAGAGGCTGCAGCGGCTGCTGCAGGCCCTCGAGGCCCCGGAGGGCTGGGGTCCGACACCAAGGGCGGGCCGCGCGGCCTTGCAGCTGAAGCTGCGTCGGCGGCTCACGGAGCTCCTGGGGGCGCAGGACGGGGCGCTGCTGGTGCGGCTGCTGCAGGCGCTGCGCGTGGCCAGGATGCCCGGGCTGGAGCGGAGCGTCCGTGAGCGCTTCCTCCCTGTGCACTGATCCTGGCCCCCTCTTATTTATTCTACATCCTTGGCACCCCACTTGCACTGAAAGAGGCTTTTTTTTAAATAGAAGAAATGAGGTTTCTT"
    shRNA = shrna.Gen_shRNA(DcR3_SEQ)
    start = len(record.seq)
    record.seq += shRNA
    feature = SeqFeature(FeatureLocation(start=start, end=len(record.seq)), type="CDS", id='DcR3_shRNA')
    record.features.append(feature)

def silenced_DcR3(record):
    ExpressionLocation.LocalExpressionInCells(record, lambda record: (
        DcR3_shRNA(record)
    ))

def Tutu22_aptamer(record):
    start = len(record.seq)
    record.seq += 'taccagtgcgatgctcagtgccgtttcttctctttcgctttttttgcttttgagcatgctgacgcattcggttgac'
    feature = SeqFeature(FeatureLocation(start=start, end=len(record.seq)), type="CDS", id='Tutu22_aptamer')
    record.features.append(feature)

def Tutu22(record):
    ExpressionLocation.LocalExpressionInCells(record, lambda record: (
        Tutu22_aptamer(record)
    ))

def construct(record):
    silenced_DcR3(record)
    Tutu22(record)
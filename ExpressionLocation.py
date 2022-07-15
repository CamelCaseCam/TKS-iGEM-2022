from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from regex import E

import human
import Gil
import params

def Promoter_pro(record):
    start = len(record.seq)
    record.seq += 'tgaagtgtgtcgtgtggcgttgcgaaatgtcaaaggtggcgttatcatagaggattgttgtgtttccctaaggggtccttcaaagcgctttcccactcgccatcgaggactttcagggagcccacggataatgaggccgaaaccgccaagcctgaccagcgcggcagtccaagcccggaacacaataagcaagaaaaggtacacca'
    feature = SeqFeature(FeatureLocation(start=start, end=len(record.seq)), type='promoter', id='Promoter_pro')
    record.features.append(feature)

def Terminator_pro(record):
    start = len(record.seq)
    record.seq += "ccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata"
    feature = SeqFeature(FeatureLocation(start=start, end=len(record.seq)), type='terminator', id='Terminator_pro')
    record.features.append(feature)

def B_LongumPromoter(record, innercode):
    Promoter_pro(record)
    innercode(record)
    Terminator_pro(record)

def LocalExpressionInCells(record, innercode):
    B_LongumPromoter(record, lambda record, innerode=innercode: 
        (human.KozakSequence(record),
        Gil.AminoSequence(record, 'gggxxxxx', Gil.HumanCodons),
        innercode(record))
    )

def AddStr(record, str):
    record.seq += str

def ShineDalgarnoSeq(record):
    record.seq += 'AGAAAGGAG'

def LuxR(record):
    start = len(record.seq)
    if params.COMPILE_CLONED_SEQUENCES:
        record.seq += 'atgaaaaacataaatgccgacgacacatacagaataattaataaaattaaagcttgtagaagcaataatgatattaatcaatgcttatctgatatgactaaaatggtacattgtgaatattatttactcgcgatcatttatcctcattctatggttaaatctgatatttcaatcctagataattaccctaaaaaatggaggcaatattatgatgacgctaatttaataaaatatgatcctatagtagattattctaactccaatcattcaccaattaattggaatatatttgaaaacaatgctgtaaataaaaaatctccaaatgtaattaaagaagcgaaaacatcaggtcttatcactgggtttagtttccctattcatacggctaacaatggcttcggaatgcttagttttgcacattcagaaaaagacaactatatagatagtttatttttacatgcgtgtatgaacataccattaattgttccttctctagttgataattatcgaaaaataaatatagcaaataataaatcaaacaacgatttaaccaaaagagaaaaagaatgtttagcgtgggcatgcgaaggaaaaagctcttgggatatttcaaaaatattaggttgcagtgagcgtactgtcactttccatttaaccaatgcgcaaatgaaactcaatacaacaaaccgctgccaaagtatttctaaagcaattttaacaggagcaattgattgcccatactttaaaaattaataacactgatagtgctagtgtagatcac'
    feature = SeqFeature(FeatureLocation(start=start, end=len(record.seq)), type='CDS', id='LuxR')
    record.features.append(feature)

def AllowQuorumDetection(record):
    LocalExpressionInCells(record, lambda record:(
        ShineDalgarnoSeq(record),
        AddStr(record, 'AATAAACAA'),
        LuxR(record)
    ))

def LuxI(record):
    start = len(record.seq)
    record.seq += 'atgactataatgataaaaaaatcggattttttggcaattccatcggaggagtataaaggtattctaagtcttcgttatcaagtgtttaagcaaagacttgagtgggacttagttgtagaaaataaccttgaatcagatgagtatgataactcaaatgcagaatatatttatgcttgtgatgatactgaaaatgtaagtggatgctggcgtttattacctacaacaggtgattatatgctgaaaagtgtttttcctgaattgcttggtcaacagagtgctcccaaagatcctaatatagtcgaattaagtcgttttgctgtaggtaaaaatagctcaaagataaataactctgctagtgaaattacaatgaaactatttgaagctatatataaacacgctgttagtcaaggtattacagaatatgtaacagtaacatcaacagcaatagagcgatttttaaagcgtattaaagttccttgtcatcgtattggagacaaagaaattcatgtattaggtgatactaaatcggttgtattgtctatgcctattaatgaacagtttaaaaaagcagtcttaaatgctgcaaacgacgaaaactacgctttagtagcttaataactctgatagtgctagtgtagatctc'
    feature = SeqFeature(FeatureLocation(start=start, end=len(record.seq)), type='CDS', id='LuxI')
    record.features.append(feature)

def SendQuorumSignal(record):
    LocalExpressionInCells(record, lambda record:(
        ShineDalgarnoSeq(record),
        AddStr(record, 'AATAAACAA'),
        LuxI(record)
    ))

def Enhancer_Euk(record):
    record.seq += "gtggagtatttacggtaaactgcccacttggcagtacatcaagtgtatcatatgccaagtacgccccctattgacgtcaatgacggtaaatggcccgcctggcattatgcccagtacatgaccttatgggactttcctacttggcagtacatctacgtattagtcatcgctattaccatg"

def Terminator_Euk(record):
    record.seq += 'attccgataacttgtttattgcagcttataatggttacaaataaagcaatagcatcacaaatttcacaaataaagcatttttttcactgcattctagttgtggtttgtccaaactcatcaatgtatcttatcatgtctg'

def Bactofection(record, innercode):
    start = len(record.seq)
    Enhancer_Euk(record)
    ShineDalgarnoSeq(record)
    AddStr(record, 'AATAAACAA')
    Gil.AminoSequence(record, 'GGGxxxxx', Gil.BLongumCodons)
    innercode(record)
    Terminator_Euk(record)
    feature = SeqFeature(FeatureLocation(start=start, end=len(record.seq)), type="Promoter", id='Bactofection')
    record.features.append(feature)

def LocalExpressionInCells_SD(record, innercode):
    LocalExpressionInCells(record, lambda record, innerode=innercode: (
        ShineDalgarnoSeq(record),
        AddStr(record, 'AATAAACAA'),
        innerode(record)
    ))

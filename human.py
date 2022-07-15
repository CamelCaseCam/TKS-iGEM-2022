from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def KozakSequence(record):
    start = len(record.seq)
    record.seq += 'GACACCATGG'
    feature = SeqFeature(FeatureLocation(start=start, end=len(record.seq)), type='misc_feature', id='KozakSequence')
    record.features.append(feature)
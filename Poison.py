from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

import Gil
import ExpressionLocation

def Saporin(record):
    Gil.AminoSequence(record, "DPNLKYGGTDIAVIGPPSRDKFLRLNFQSSRGTVSLGLKRENLYVVAYLAMDNANVNRAYYFGTEITSAELTTLLPEATVANQKALEYTEDYQSIEKNAKITEGDKTRKELGLGINLLSTLMDAVNKKARVVKNEARFLLIAIQMTAEAARFRYIQNLVTKNFPNKFNSEDKVIQFQVNWSKISKAIYGDAKNGVFNKDYDFGFGKVRQVKDLQMGLLMYLGTTPNNAADRYRAEL", 
            Gil.HumanCodons)

def TrojanHorseSaporin(record):
    ExpressionLocation.Bactofection(record, lambda record: (
        Saporin(record)
    ))

def construct(record):
    TrojanHorseSaporin(record)
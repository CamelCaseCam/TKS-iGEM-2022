import ExpressionLocation
import Poison
import DeathAndAptamers
import SelfDestruct

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def main():
    # Create a record
    record = SeqRecord(
        Seq(""),
        name="TKS_iGEM_construct",
        description="TKS iGEM construct",
        annotations={"molecule_type": "DNA"}
    )
    
    # run all construct functions
    ExpressionLocation.AllowQuorumDetection(record)
    ExpressionLocation.SendQuorumSignal(record)

    Poison.construct(record)
    DeathAndAptamers.construct(record)
    SelfDestruct.construct(record)
    
    
    # Save as GenBank file
    output_file = open('example.gb', 'w')
    SeqIO.write(record, output_file, 'genbank')

main()
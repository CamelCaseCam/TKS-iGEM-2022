HumanCodons = {
    'g':['GGC', 'GGG', 'GGA', 'GGT'],
    'a':['GCC', 'GCT', 'GCA', 'GCG'],
    'v':['GTG', 'GTC', 'GTT', 'GTA'],
    'l':['CTG', 'CTC', 'CTT', 'TTG', 'TTA', 'CTA'],
    'i':['ATC', 'ATT', 'ATA'],
    'm':['ATG'],
    'p':['CCC', 'CCT', 'CCA', 'CCG'],
    'f':['TTC', 'TTT'],
    'w':['TGG'],
    's':['AGC', 'TCC', 'TCT', 'AGT', 'TCA', 'TCG'],
    't':['ACC', 'ACA', 'ACT', 'ACG'],
    'n':['AAC', 'AAT'],
    'q':['CAG', 'CAA'],
    'y':['TAC', 'TAT'],
    'c':['TGC', 'TGT'],
    'k':['AAG', 'AAA'],
    'r':['CGG', 'AGA', 'AGG', 'CGC', 'CGA', 'CGT'],
    'h':['CAC', 'CAT'],
    'd':['GAC', 'GAC'],
    'e':['GAG', 'GAA'],
    'x':['TGA', 'TAA', 'TAG'],
}

BLongumCodons = {
    'g':['GGC', 'GGT', 'GGG', 'GGA'],
    'a':['GCC', 'GCG', 'GCT', 'GCA'],
    'v':['GTG', 'GTC', 'GTT', 'GTA'],
    'l':['CUG', 'CUC', 'TTG', 'CTT', 'CTA', 'TTA'],
    'i':['ATC', 'ATT', 'ATA'],
    'm':['ATG'],
    'p':['CCG', 'CCC', 'CCT', 'CCA'],
    'f':['TTC', 'TTT'],
    'w':['TGG'],
    's':['TCC', 'AGC', 'TCG', 'TCT', 'TCA', 'AGU'],
    't':['ACC', 'ACG', 'ACT', 'ACA'],
    'n':['AAC', 'AAT'],
    'q':['CAG', 'CAA'],
    'y':['TAC', 'TAT'],
    'c':['TGC', 'TGT'],
    'k':['AAG', 'AAA'],
    'r':['CGC', 'CGG', 'CGT', 'AGG', 'CGC', 'AGA'],
    'h':['CAC', 'CAT'],
    'd':['GAC', 'GAT'],
    'e':['GAG', 'GAA'],
    'x':['TGA', 'TAA', 'TAG'],
}

def AminoSequence(record, aminos, target):
    for amino in aminos.lower():
        record.seq += target[amino][0]
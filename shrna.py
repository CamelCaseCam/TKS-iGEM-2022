import random

def IsGHG(chr1, chr3):
    return chr1 == "g" and chr3 == "g"

def GenRandBases(length: int) -> str:
    """
    Generate a random string of bases
    """
    return "".join(random.choice("acgt") for _ in range(length))

def ContainsPossibleGHG(region: str) -> bool:
    """
    Check if the region contains a possible GHG
    """
    for i in range(len(region) - 2):
        if IsGHG(region[i], region[i+2]):
            return True
    return False

def GetRegionWithoutStart(region: str, length) -> str:
    """
    Get the region without the start codon
    """
    #loop until you reach a region without a stop codon in it
    i = 0
    while True:
        if i + length >= len(region):
            raise Exception("Could not find a region without a stop codon")
        
        # Check if the region contains a start codon
        if "atg" in region[i:length]:
            i = region[i:length].index("atg")
        elif not ContainsPossibleGHG(region[i:length]):
            #advance by the length of the region
            i += length - 1
        else:
            return region[i:length]

        i += 1

def Arich(length):
    return "a" * length

complements = {
    "a": "u",
    "t": "a",
    "u": "a",
    "c": "g",
    "g": "c",
}
def GetReverseComplement(region: str) -> str:
    """
    Get the reverse complement of a region
    """
    return "".join(complements[base] for base in reversed(region.upper()))

def GetComplement(region: str) -> str:
    """
    Get the reverse complement of a region
    """
    return "".join(complements[base] for base in region)

def GenLoop(LoopRandSize: int) -> str:
    """
    Generate a random loop
    """
    return "ugug" + GenRandBases(LoopRandSize) + "a"


def ConvertToGHGMismatch(region, complement, idx) -> str:
    #First, check if idx + 1 is a G
    if region[idx + 1] == "g":
        region[idx + 1] = "u"
        return region, reversed(complement)
    if complement[idx + 1] == "g":
        complement[idx + 1] = "u"
        return region, reversed(complement)
    
    # Now, check if idx + 1 is a U
    if region[idx + 1] == "u" or region[idx + 1] == "t":
        region[idx + 1] = "a"
        return region, reversed(complement)

    # Now, check if idx + 1 is an A
    if region[idx + 1] == "a":
        region[idx + 1] = "u"
        return region, reversed(complement)


def GenStem(StemSize: int, Target: str) -> str:
    """
    Generate an shRNA stem
    """
    Stem = GetRegionWithoutStart(Target, StemSize)
    StemComplement = GetComplement(Stem)

    # now add the GHG mismatch motif
    i = 0
    while True:
        if i + 2 >= len(Stem):
            raise Exception("Could not create the GHG mismatch motif")
        
        if IsGHG(Stem[i], Stem[i+2]):
            Stem, StemComplement = ConvertToGHGMismatch(list(Stem), list(StemComplement), i)
            return Arich(12) + "ug" + "".join(Stem) + GenLoop(LOOP_LENGTH) + "".join(StemComplement) + "c" + Arich(5) + "c" + GenRandBases(2) + "c" + Arich(4)

        i += 1

    



# length is 24 for the overhangs
HAIRPIN_LENGTH = 35
LOOP_LENGTH = 10
FLANK_LENGTH = 12

def Gen_shRNA(target: str):
    # We want to generate an shRNA for the target that does not contain any start codons. Technically, this will be a pre-miRNA
    # Get the region to target
    target=list(target.replace("\n", "").replace("\r", "").replace("\t", "").replace(" ", "").lower())
    Hairpin = GenStem(HAIRPIN_LENGTH, target)

    return Hairpin
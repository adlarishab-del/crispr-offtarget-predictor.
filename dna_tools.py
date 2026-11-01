def main():
    pass
import pandas as pd

input_sequence = "ATGCCG"

pairs = {
    "A":"T",
    "C":"G",
    "G":"C",
    "T":"A"
}


def reverse_complement(seq):
    seq = seq.upper()
    result = ""
    for index in range(len(seq) - 1,-1,-1):
        base = seq[index]
        if base not in "ATGC":
            raise ValueError("Invalid Base")
        complement = pairs[base]
        result += complement
    return result

rc = reverse_complement(input_sequence)
print("The reverse complement of this seq is:" + rc)


def gc_content(seq):
    seq = seq.upper()
    if len(seq) == 0:
        return 0.0
    
    gc_count = 0
    for base in seq:
        if base in ["G","C"]:
            gc_count += 1
    return gc_count / len(seq)

gc = gc_content(input_sequence)

def find_pam_sites(seq):
    seq = seq.upper()
    sites = []
    for i in range(20, len(seq) - 2): # 20 is always used is because the guide RNA cas9 uses is always 20 letters long.
        chunk = seq[i:i+3]
        if chunk.endswith("GG"):
            sites.append(i)
    return sites


print(find_pam_sites("AATGGCC"))   # [2]
print(find_pam_sites("AGGTGG"))    # [0, 3]
print(find_pam_sites("aggtgg"))    # [0, 3]
print(find_pam_sites("AAAA"))      # []
print(find_pam_sites("A" * 20 + "TGG"))            # [20]
print(find_pam_sites("A" * 20 + "TGGAGG"))         # [20, 23]
print(find_pam_sites("TGG" + "A" * 20 + "TGG"))    # [23]  (the first TGG has no room)
print(find_pam_sites("AGGTGG"))                    # []    (too short)

if __name__ == '__main__':
    main()





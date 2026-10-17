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








print("The ratio of the GC content in the whole sequence is:", str(gc))
print(gc_content("ATGC"))    # 0.5
print(gc_content("AAAA"))    # 0.0
print(gc_content("GGCC"))    # 1.0
print(gc_content(""))        # no crash
print(gc_content("ATGCCG"))  # 0.6666...
print(gc_content(""))       # 0.0
print(gc_content("gcgc"))   # 1.0
print(reverse_complement("aaaa"))  # TTTT



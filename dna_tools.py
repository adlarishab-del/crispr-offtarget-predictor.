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

def find_guides(seq):
    seq = seq.upper()
    guides = []
    for i in find_pam_sites(seq):
        guide = (seq[i-20:i])
        guides.append(guide)
    return guides 

def find_all_guides(seq):
    forward =find_guides(seq)
    reverse = find_guides(reverse_complement(seq))
    return forward + reverse

def has_homopolymer(guide):
    guide = guide.upper()
    for run in ["AAAA","CCCC","GGGG","TTTT"]: # TLDR This just sees if guide in run has the AAAA,CCCC,GGGG or TTTT.
        if run in guide:
            return True 
    return False 

def score_guide(guide):
    score = 0
    if not has_homopolymer(guide): # TTTT can act as a stop signal, so it can get cut short. Long G runs also can cause the guide to turn out odd, so if the guide has no run of homopolymers, to cause a mishap, it adds one as it isn't a bad run.
        score += 1
    if 0.4 <= gc_content(guide) <= 0.7: # If the gc content is below 40%, the bind target is weak, and if its above 70% it binds too tightly. hence, we need it within this sweet spot.
        score +=1
    return score 
    


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





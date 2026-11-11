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

def rank_guides(seq):
    results = []
    for guide in find_all_guides(seq):
        results.append((guide, score_guide(guide))) # runs for each guide
    results.sort(key= lambda pair: pair[1], reverse = True) #after loop, runs once, lambda just is a easier way for me to return pair[1]
    return results 

test = "ACGT" * 5 + "TGG" + "A" * 20 + "AGG"

for guide, score in rank_guides(test):
    print(guide,score)

def count_mismatches(a,b): # Compares sequences "A" and "B", and counts the no. of mismatches between them
    if len(a) != len(b):
        raise ValueError("The length of a should be equal to the length of b.")    
    mismatches = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            mismatches += 1
    return mismatches 

def weighted_mismatches(a,b): # Index 0-19 works as a "target" as there are usually 20 bases next to PAM(NGG.) What this does is showos how far it is from the pam via weight.
    if len(a) != len(b):
        raise ValueError("The length of a should be equal to the length of b.") 
    penalty = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            if i >= 10:
               penalty += 2
            else:
                penalty += 1
    return penalty

def scan_for_matches(guide,genome,max_penalty):
    genome = genome.upper()
    hits = []
    for i in range(len(genome) - 22):
        window = genome[i:i+20]
        pam = genome[i+20:i+23]
        penalty = weighted_mismatches(guide,window)
        if pam[1:3] == "GG" and penalty <= max_penalty:
            hits.append((i,penalty))
    return hits 

def scan_both_strands(guide,genome,max_penalty):
    forward = scan_for_matches(guide,genome,max_penalty)
    reverse = scan_for_matches(reverse_complement(guide),genome,max_penalty)
    return forward + reverse 







guide = "A" * 20
genome = "C" * 5 + "A" * 20 + "C" * 5
g1 = "C" * 5 + "A" * 20 + "C" * 5
g2 = "C" * 5 + "A" * 20 + "TGG" + "C" * 5
g3 = "C" * 5 + "A" * 19 + "T" + "AGG" + "C" * 5

guide = "A" * 20
fwd = "C" * 5 + "A" * 20 + "TGG" + "C" * 5
rev = reverse_complement(fwd)

print(len(scan_both_strands(guide, fwd, 0)))   # 1
print(len(scan_both_strands(guide, rev, 0)))   # 1
print(scan_for_matches(guide, g1, 0))   # []  (perfect match but no PAM)
print(scan_for_matches(guide, g2, 0))   # [(5, 0)]
print(scan_for_matches(guide, g3, 2))   # [(5, 2)]
print(scan_for_matches(guide, g3, 1))   # []  (seed mismatch costs 2)

print(scan_for_matches(guide, genome, 0))   # [(5, 0)]
print(scan_for_matches(guide, genome, 1))   # [(4, 1), (5, 0)]
print(scan_for_matches(guide, "AAAA", 0))   # []

# print(count_mismatches("ACGT","AC"))              gi
print(count_mismatches("ACGT", "AGGT"))              # 1
print(count_mismatches("AAAA", "TTTT"))              # 4
print(count_mismatches("ACGTACGTAC", "ACGTTCGTAA"))  # 2

print(weighted_mismatches("A" * 20, "A" * 20))          # 0
print(weighted_mismatches("A" * 20, "T" + "A" * 19))    # 1
print(weighted_mismatches("A" * 20, "A" * 19 + "T"))    # 2
print(weighted_mismatches("A" * 20, "T" * 20))          # 30  (10 x 1 + 10 x 2)





if __name__ == '__main__':
    main()





pairs = {
    "A":"T",
    "C":"G",
    "G":"C",
    "T":"A"
}

GUIDE_LEN = 20

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

def gc_content(seq):
    seq = seq.upper()
    if len(seq) == 0:
        return 0.0
    
    gc_count = 0
    for base in seq:
        if base in ["G","C"]:
            gc_count += 1
    return gc_count / len(seq)

def find_pam_sites(seq):
    seq = seq.upper()
    sites = []
    for i in range(GUIDE_LEN, len(seq) - 2): # 20 is always used is because the guide RNA cas9 uses is always 20 letters long.
        chunk = seq[i:i+3]
        if chunk.endswith("GG"):
            sites.append(i)
    return sites

def find_guides(seq):
    seq = seq.upper()
    guides = []
    for i in find_pam_sites(seq):
        guide = (seq[i-GUIDE_LEN:i])
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
    for i in range(len(genome) - GUIDE_LEN):
        if genome[i+GUIDE_LEN+1:i+GUIDE_LEN + 3] != "GG": # PAM check
            continue
        penalty = weighted_mismatches(guide,genome[i:i+GUIDE_LEN])
        if penalty <= max_penalty:
            hits.append((i,penalty))
    return hits 

def scan_both_strands(guide,genome,max_penalty):
    forward = scan_for_matches(guide,genome,max_penalty)
    reverse = scan_for_matches(guide,reverse_complement(genome),max_penalty)
    return forward + reverse 

GENE_CHROM = "chr2"
GENE_START = 72916189

def window_to_genome(i):
    start = GENE_START + i - 1 
    end = start + GUIDE_LEN
    return (GENE_CHROM,start,end)

def reverse_window_to_genome(j,seq_len):
    forward_i = seq_len - j - GUIDE_LEN #  # j = position of a hit found in the FLIPPED copy of the gene (reverse_complement(seq)).
    return window_to_genome(forward_i)

def scan_for_genomic(guide,seq,max_penalty): # Scan both strands, return hits as (chrom, start, end, strand, penalty) in genome coordinates.
    results = []
    for i, penalty in scan_for_matches(guide,seq,max_penalty):
        chrom, start, end = window_to_genome(i)
        results.append((chrom,start,end,"+",penalty))
    for j, penalty in scan_for_matches(guide,reverse_complement(seq),max_penalty):
        chrom,start,end = reverse_window_to_genome(j,len(seq))
        results.append((chrom,start,end,"-",penalty))
    return results 

def load_peaks(path):
    peaks = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            peaks.append((parts[0],int(parts[1]),int(parts[2])))
    return peaks 
   # Reads a BED peaks file (one open-chromatin region per line: chrom, start, end)
    # returns a list of (chrom, start, end) tuples, with start/end as ints.

def overlaps_peak(chrom,start,end,peaks):
    for p_chrom,p_start,p_end in peaks:
        if chrom == p_chrom and start < p_end and end > p_start:
            return True
    return False
# Returns True if the window (chrom, start, end) overlaps any open-chromatin peak.
    # peaks = list of (chrom, start, end) tuples from load_peaks().
    # Two regions overlap when they're on the same chromosome AND each one starts
    # before the other ends:  window start < peak end  and  window end > peak start.

def annotate_hits(hits,peaks):
    annotations = []
    for chrom,start,end,strand,penalty in hits:
        is_open = overlaps_peak(chrom,start,end,peaks)
        annotations.append((chrom,start,end,strand,penalty,is_open))
    return annotations

def hit_risk(penalty,is_open):
    risk = 1/(1 + penalty)
    if not is_open:
        risk = risk*0.5  # closed chromatin: harder for Cas9 to reach, so halve it
    return risk

def offtarget_risk(guide,seq,peaks,max_penalty):
    hits = annotate_hits(scan_for_genomic(guide,seq,max_penalty),peaks)
    total = 0
    for chrom,start,end,strand,penalty,is_open in hits:
        if penalty > 0:
            total += hit_risk(penalty,is_open) # Sum hit_risk() over every off-target (penalty > 0) found on both strands: close matches in open chromatin add the most risk, closed sites are halved.
    return total 

def rank_by_risk(seq,peaks,max_penalty,top_n):
    ranked = rank_guides(seq)[:top_n]
    results = []
    for guide, score in ranked:
        results.append((guide,score,offtarget_risk(guide, seq, peaks, max_penalty)))
    results.sort(key=lambda r: (-r[1],r[2]))
    return results







        
        







def main():
    with open("gene.fna") as f:
        seq = "".join(line.strip() for line in f if not line.startswith(">"))
    peaks = [("chr2", 72916000, 72916200)]
    hits = [("chr2", 72916193, 72916213, "+", 0),
        ("chr2", 72920000, 72920020, "-", 2)]
    assert annotate_hits(hits, peaks) == [("chr2", 72916193, 72916213, "+", 0, True),
                                       ("chr2", 72920000, 72920020, "-", 2, False)]
    print("annotate_hits test passed")

    assert reverse_window_to_genome(0, 18703) == ("chr2", 72934871, 72934891)
    assert reverse_window_to_genome(18683, 18703) == ("chr2", 72916188, 72916208)

    peaks = load_peaks("toy_peaks.bed")
    assert overlaps_peak("chr2", 150, 170, peaks) == True
    assert overlaps_peak("chr2", 190, 250, peaks) == True
    assert overlaps_peak("chr2", 200, 300, peaks) == False
    assert overlaps_peak("chr2", 300, 400, peaks) == False
    assert overlaps_peak("chr4", 150, 170, peaks) == False
    print("overlaps_peak tests passed")

    assert load_peaks("toy_peaks.bed") == [("chr2", 100, 200), ("chr2", 500, 600), ("chr3", 50, 80)]
    print("load_peaks test passed")

    ranked = rank_guides(seq)
    input_sequence = "ATGCCG"

    site = "CCA" + "T" * 20
    assert len(scan_both_strands("A" * 20, site, 0)) == 1

    g1 = "C" * 5 + "A" * 20 + "C" * 5
    g2 = "C" * 5 + "A" * 20 + "TGG" + "C" * 5
    g3 = "C" * 5 + "A" * 19 + "T" + "AGG" + "C" * 5

    guide = "A" * 20
    fwd = "C" * 5 + "A" * 20 + "TGG" + "C" * 5
    rev = reverse_complement(fwd)
    print(scan_for_matches(guide, g1, 4))
    print(scan_for_matches(guide, g2, 4))
    print(scan_for_matches(guide, g3, 4))
    rc = reverse_complement(input_sequence)
    assert scan_for_matches(guide, g1, 4) == []
    assert scan_for_matches(guide, g2, 4) == [(5, 0)]
    assert scan_for_matches(guide, g3, 4) == [(5, 2)]
    print("scan_for_matches tests passed")
    assert scan_for_matches(guide, g3, 1) == []
    print(len(ranked))
    print(ranked[:5])
    print("The reverse complement of this seq is:" + rc)
    test = "ACGT" * 5 + "TGG" + "A" * 20 + "AGG"
    for guide, score in rank_guides(test):
        print(guide,score)
    print(len(seq))
    gc = gc_content(input_sequence)
    print(gc)

    fwd_site = "C" * 5 + "A" * GUIDE_LEN + "TGG" + "C" * 5
    assert scan_for_genomic("A" * GUIDE_LEN, fwd_site, 0) == [("chr2", 72916193, 72916213, "+", 0)]

    rev_site = "CCA" + "T" * GUIDE_LEN
    assert scan_for_genomic("A" * GUIDE_LEN, rev_site, 0) == [("chr2", 72916191, 72916211, "-", 0)]
    print("scan_genomic tests passed")

    assert hit_risk(0, True) == 1.0
    assert hit_risk(1, True) == 0.5
    assert hit_risk(1, False) == 0.25
    assert hit_risk(3, True) == 0.25
    assert hit_risk(3, False) == 0.125
    print("hit_risk tests passed")

    results = rank_by_risk(seq, [], 3, 2)
    for guide, score, risk in results:
     print(guide, score, round(risk, 3))

    for a, b in zip(results, results[1:]):
        assert a[1] > b[1] or (a[1] == b[1] and a[2] <= b[2])
    print("rank_by_risk ordering passed")



if __name__ == '__main__':
    main()





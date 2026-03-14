from collections import defaultdict

# =========================
# Parametry
# =========================
fasta_file = "C:/LM.txt"
k = 20  # długość szukanego motywu

# =========================
# Funkcje
# =========================
def read_fasta(path):
    seq = ""
    with open(path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip().upper()
    return seq

def reverse_complement(seq):
    complement = str.maketrans("ATCGN", "TAGCN")
    return seq.translate(complement)[::-1]

def find_repeats(seq, k):
    repeats = defaultdict(list)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        repeats[kmer].append(i + 1)  # pozycje od 1
    return {kmer: pos for kmer, pos in repeats.items() if len(pos) > 15}

# =========================
# Analiza
# =========================
sequence = read_fasta(fasta_file)
rev_comp_sequence = reverse_complement(sequence)

repeats_forward = find_repeats(sequence, k)
repeats_reverse = find_repeats(rev_comp_sequence, k)

# =========================
# Wyniki
# =========================
print("=== Powtórzenia w sekwencji 5'→3' ===")
for kmer, positions in repeats_forward.items():
    print(f"{kmer} -> pozycje: {positions}")

print("\n=== Powtórzenia w sekwencji komplementarnej (5'→3') ===")
for kmer, positions in repeats_reverse.items():
    print(f"{kmer} -> pozycje: {positions}")

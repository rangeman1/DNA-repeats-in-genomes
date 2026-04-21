import os
from collections import defaultdict
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext

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

def gc_content(seq):  # NEW
    gc = seq.count("G") + seq.count("C")
    return (gc / len(seq)) * 100 if len(seq) > 0 else 0

def find_repeats(seq, k, min_repeats, gc_min=None, gc_max=None):  # NEW
    repeats = defaultdict(list)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        repeats[kmer].append(i + 1)

    # filtrowanie po liczbie powtórzeń i GC
    filtered = {}
    for kmer, pos in repeats.items():
        if len(pos) >= min_repeats:
            gc = gc_content(kmer)
            if gc_min is not None and gc < gc_min:
                continue
            if gc_max is not None and gc > gc_max:
                continue
            filtered[kmer] = pos

    return filtered

# =========================
# Funkcja uruchamiająca analizę
# =========================
def run_analysis():
    path = file_path.get()
    
    if not os.path.exists(path):
        messagebox.showerror("Błąd", "Plik nie istnieje!")
        return
    
    try:
        k = int(k_entry.get())
        min_rep = int(min_entry.get())
        gc_min = float(gc_min_entry.get()) if gc_min_entry.get() else None  # NEW
        gc_max = float(gc_max_entry.get()) if gc_max_entry.get() else None  # NEW
    except ValueError:
        messagebox.showerror("Błąd", "Podaj poprawne wartości liczbowe!")
        return

    sequence = read_fasta(path)
    rev_comp_sequence = reverse_complement(sequence)

    repeats_forward = find_repeats(sequence, k, min_rep, gc_min, gc_max)  # NEW
    repeats_reverse = find_repeats(rev_comp_sequence, k, min_rep, gc_min, gc_max)  # NEW

    output.delete(1.0, tk.END)

    output.insert(tk.END, "=== Powtórzenia w sekwencji 5'→3' ===\n")
    for kmer, positions in repeats_forward.items():
        output.insert(tk.END, f"{kmer} (GC={gc_content(kmer):.1f}%) -> {positions}\n")

    output.insert(tk.END, "\n=== Powtórzenia w sekwencji komplementarnej (5'→3') ===\n")
    for kmer, positions in repeats_reverse.items():
        output.insert(tk.END, f"{kmer} (GC={gc_content(kmer):.1f}%) -> {positions}\n")

# =========================
# Wybór pliku
# =========================
def browse_file():
    filename = filedialog.askopenfilename(
        title="Wybierz plik FASTA",
        filetypes=[("Pliki tekstowe", "*.txt *.fasta"), ("Wszystkie pliki", "*.*")]
    )
    file_path.set(filename)

# =========================
# GUI
# =========================
root = tk.Tk()
root.title("Analiza powtórzeń w genomie")
root.geometry("700x550")

file_path = tk.StringVar()

# Wybór pliku
tk.Label(root, text="Plik FASTA:").pack()
tk.Entry(root, textvariable=file_path, width=60).pack()
tk.Button(root, text="Wybierz plik", command=browse_file).pack(pady=5)

# Parametry
tk.Label(root, text="Długość motywu (k):").pack()
k_entry = tk.Entry(root)
k_entry.pack()

tk.Label(root, text="Minimalna liczba powtórzeń:").pack()
min_entry = tk.Entry(root)
min_entry.pack()

# NEW — zakres GC
tk.Label(root, text="Minimalna zawartość GC (%) [opcjonalnie]:").pack()
gc_min_entry = tk.Entry(root)
gc_min_entry.pack()

tk.Label(root, text="Maksymalna zawartość GC (%) [opcjonalnie]:").pack()
gc_max_entry = tk.Entry(root)
gc_max_entry.pack()

# Przycisk start
tk.Button(root, text="Uruchom analizę", command=run_analysis, bg="lightgreen").pack(pady=10)

# Pole wyników
output = scrolledtext.ScrolledText(root, width=80, height=20)
output.pack(pady=10)

# Start GUI
root.mainloop()

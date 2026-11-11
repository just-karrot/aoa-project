import time
import os
import re
from pathlib import Path
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# 1. ALGORITHM IMPLEMENTATIONS
# ---------------------------------------------------------

def naive_search(text: str, pattern: str):
    n, m = len(text), len(pattern)
    matches = []
    comparisons = 0
    for i in range(n - m + 1):
        j = 0
        while j < m:
            comparisons += 1
            if text[i + j] != pattern[j]:
                break
            j += 1
        if j == m:
            matches.append(i)
    return matches, comparisons


def rabin_karp_search(text: str, pattern: str, q: int = 1_000_000_007):
    n, m = len(text), len(pattern)
    matches = []
    comparisons = 0
    if m > n:
        return matches, comparisons

    base = 4
    char_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    high = pow(base, m - 1, q)
    p_hash = 0
    t_hash = 0

    for i in range(m):
        p_hash = (p_hash * base + char_map.get(pattern[i], 0)) % q
        t_hash = (t_hash * base + char_map.get(text[i], 0)) % q

    for i in range(n - m + 1):
        if p_hash == t_hash:
            j = 0
            while j < m:
                comparisons += 1
                if text[i + j] != pattern[j]:
                    break
                j += 1
            if j == m:
                matches.append(i)

        if i < n - m:
            t_hash = ((t_hash - char_map.get(text[i], 0) * high) * base + char_map.get(text[i + m], 0)) % q

    return matches, comparisons


def kmp_search(text: str, pattern: str):
    n, m = len(text), len(pattern)
    matches = []
    comparisons = 0
    if m > n or m == 0:
        return matches, comparisons

    # Build failure table (preprocessing)
    table = [0] * m
    length = 0
    i = 1
    while i < m:
        if pattern[i] == pattern[length]:
            length += 1
            table[i] = length
            i += 1
        else:
            if length != 0:
                length = table[length - 1]
            else:
                table[i] = 0
                i += 1

    # Search
    i = 0  # text index
    j = 0  # pattern index
    while i < n:
        comparisons += 1
        if text[i] == pattern[j]:
            i += 1
            j += 1
        else:
            if j > 0:
                j = table[j - 1]
            else:
                i += 1

        if j == m:
            matches.append(i - m)
            j = table[m - 1]

    return matches, comparisons


def boyer_moore_search(text: str, pattern: str):
    n, m = len(text), len(pattern)
    matches = []
    comparisons = 0
    if m > n or m == 0:
        return matches, comparisons

    # Bad-character table
    bad_char = {ch: -1 for ch in "ACGT"}
    for idx, ch in enumerate(pattern):
        bad_char[ch] = idx

    # Good-suffix table
    border = [0] * (m + 1)
    shift = [0] * (m + 1)
    i = m
    j = m + 1
    border[i] = j
    while i > 0:
        while j <= m and pattern[i - 1] != pattern[j - 1]:
            if shift[j] == 0:
                shift[j] = j - i
            j = border[j]
        i -= 1
        j -= 1
        border[i] = j

    j = border[0]
    for i in range(m + 1):
        if shift[i] == 0:
            shift[i] = j
        if i == j:
            j = border[j]

    # Search
    s = 0
    while s <= n - m:
        j = m - 1
        while j >= 0:
            comparisons += 1
            if pattern[j] == text[s + j]:
                j -= 1
            else:
                break

        if j < 0:
            matches.append(s)
            s += max(shift[0], 1)
        else:
            bc_shift = j - bad_char.get(text[s + j], -1)
            gs_shift = shift[j + 1]
            s += max(bc_shift, gs_shift, 1)

    return matches, comparisons


# ---------------------------------------------------------
# 2. DATA LOADER
# ---------------------------------------------------------

def load_dna_data():
    """Finds and loads sequence from data/sample.txt or data/human.txt."""
    candidates = [
        Path("data/sample.txt"),
        Path("../data/sample.txt"),
        Path("data/human.txt"),
        Path("../data/human.txt"),
    ]
    for p in candidates:
        if p.exists():
            text = p.read_text(encoding="ascii", errors="replace")
            # Extract FASTA or TSV sequence
            lines = [l.strip() for l in text.splitlines() if l.strip() and not l.startswith(">") and not l.startswith("sequence")]
            seq = "".join(lines)
            seq = re.sub(r"[^ACGT]", "", seq.upper())
            print(f"Loaded DNA data from {p} (length: {len(seq)} bp)")
            return seq

    raise FileNotFoundError("Could not find sample.txt or human.txt in data/ folder.")


# ---------------------------------------------------------
# 3. BENCHMARK & EVALUATION
# ---------------------------------------------------------

def run_evaluation():
    dna_text = load_dna_data()

    # Pick test patterns: biological motifs + slices of different lengths
    test_cases = [
        ("TATA Box", "TATAAA"),
        ("Start Codon Region", "ATGGCC"),
        ("Short Motif (m=10)", dna_text[100:110]),
        ("Medium Motif (m=20)", dna_text[200:220]),
        ("Long Motif (m=40)", dna_text[300:340]),
        ("Very Long Motif (m=80)", dna_text[400:480]),
    ]

    algos = {
        "Naive": naive_search,
        "Rabin-Karp": rabin_karp_search,
        "KMP": kmp_search,
        "Boyer-Moore": boyer_moore_search,
    }

    results = {name: {"times": [], "comps": [], "matches": []} for name in algos}
    labels = []

    print("\n--- Running Evaluation ---")
    for label, pattern in test_cases:
        labels.append(f"{label}\n(m={len(pattern)})")
        print(f"\nPattern: '{label}' (len={len(pattern)})")

        for name, fn in algos.items():
            # CPU / High-resolution Wall-clock Timer
            t_start = time.perf_counter()
            matches, comps = fn(dna_text, pattern)
            t_end = time.perf_counter()

            elapsed_ms = (t_end - t_start) * 1000.0
            results[name]["times"].append(elapsed_ms)
            results[name]["comps"].append(comps)
            results[name]["matches"].append(len(matches))

            print(f"  {name:12s} | Time: {elapsed_ms:8.3f} ms | Comparisons: {comps:8d} | Matches: {len(matches)}")

    # ---------------------------------------------------------
    # 4. PLOT COMPARISON GRAPHS
    # ---------------------------------------------------------
    os.makedirs("results", exist_ok=True)
    colors = {"Naive": "#e74c3c", "Rabin-Karp": "#3498db", "KMP": "#2ecc71", "Boyer-Moore": "#f39c12"}

    # Graph 1: Execution Time Comparison
    plt.figure(figsize=(10, 5))
    x = range(len(labels))
    width = 0.2
    for idx, (name, data) in enumerate(results.items()):
        offset = (idx - 1.5) * width
        plt.bar([p + offset for p in x], data["times"], width=width, label=name, color=colors[name])
    plt.xticks(x, labels, fontsize=9)
    plt.ylabel("Execution Time (milliseconds)")
    plt.title("Algorithm Execution Time Comparison on Sample DNA")
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig("results/graph1_execution_time.png", dpi=150)
    plt.close()

    # Graph 2: Character Comparisons Comparison
    plt.figure(figsize=(10, 5))
    for idx, (name, data) in enumerate(results.items()):
        offset = (idx - 1.5) * width
        plt.bar([p + offset for p in x], data["comps"], width=width, label=name, color=colors[name])
    plt.xticks(x, labels, fontsize=9)
    plt.ylabel("Character Comparisons")
    plt.title("Character Comparisons Comparison on Sample DNA")
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig("results/graph2_comparisons.png", dpi=150)
    plt.close()

    # Graph 3: Comparisons vs Pattern Length (Trend line)
    plt.figure(figsize=(9, 5))
    m_lengths = [len(p) for _, p in test_cases]
    for name, data in results.items():
        plt.plot(m_lengths, data["comps"], marker='o', label=name, color=colors[name], linewidth=2)
    plt.xlabel("Pattern Length (m)")
    plt.ylabel("Comparisons Count")
    plt.title("Comparisons vs Pattern Length")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig("results/graph3_trend_pattern_length.png", dpi=150)
    plt.close()

    # Graph 4: Speedup Relative to Naive
    plt.figure(figsize=(10, 5))
    for name in ["Rabin-Karp", "KMP", "Boyer-Moore"]:
        speedups = [results["Naive"]["times"][i] / max(results[name]["times"][i], 1e-6) for i in range(len(labels))]
        plt.plot(labels, speedups, marker='s', label=f"{name} Speedup vs Naive", linewidth=2)
    plt.ylabel("Speedup Multiplier (x times faster)")
    plt.title("Speedup Relative to Naive Search")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig("results/graph4_speedup.png", dpi=150)
    plt.close()

    print("\n[Done] All 4 comparison graphs generated and saved in the 'results/' folder:")
    print("  - results/graph1_execution_time.png")
    print("  - results/graph2_comparisons.png")
    print("  - results/graph3_trend_pattern_length.png")
    print("  - results/graph4_speedup.png")


if __name__ == "__main__":
    run_evaluation()

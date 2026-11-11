# DNA Sequence Pattern Matching - Analysis of Algorithms

A comparative study and benchmark of four fundamental exact pattern matching algorithms applied to genomic DNA sequences:
1. **Naive (Brute-Force) Search**
2. **Rabin-Karp Rolling Hash Algorithm**
3. **Knuth-Morris-Pratt (KMP) Algorithm**
4. **Boyer-Moore Algorithm (Bad-Character + Good-Suffix Rules)**

---

## Project Structure

```
aoa-project/
├── data/
│   ├── sample.txt          # Real NCBI Homo sapiens Keratin FASTA sequence (~5.3 kb)
│   └── human.txt           # Additional genomic sequence dataset
├── results/                # Output comparison graphs
│   ├── graph1_execution_time.png       # CPU execution times per motif
│   ├── graph2_comparisons.png          # Total character comparisons per algorithm
│   ├── graph3_trend_pattern_length.png # Scaling curve across pattern length (m)
│   └── graph4_speedup.png              # Speedup relative to Naive baseline
├── run_evaluation.py       # Standalone execution and graphing script
├── requirements.txt        # Minimal dependencies (matplotlib, pytest)
├── README.md               # Quickstart guide & documentation
├── REPORT.md               # Detailed academic report (Markdown format)
└── REPORT.txt              # Detailed academic report (Plain text format)
```

---

## Algorithm Complexities

- $n$: Length of text sequence
- $m$: Length of motif/pattern
- $|\Sigma|$: Alphabet size ($|\Sigma| = 4$ for DNA: `A`, `C`, `G`, `T`)

| Algorithm | Preprocessing Time | Matching Time (Best) | Matching Time (Worst) | Space Complexity | Key Feature |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Naive** | $O(1)$ | $O(n)$ | $O(n \cdot m)$ | $O(1)$ | Simple baseline, shifts by 1 on mismatch |
| **Rabin-Karp** | $O(m)$ | $O(n + m)$ | $O(n \cdot m)$ | $O(1)$ | Polynomial rolling hash with prime modulus |
| **KMP** | $O(m)$ | $O(n)$ | $O(n + m)$ | $O(m)$ | Prefix failure table, monotonic text pointer |
| **Boyer-Moore** | $O(m + |\Sigma|)$ | $O(n / m)$ | $O(n \cdot m)$ | $O(m + |\Sigma|)$ | Right-to-left scan, bad-char & good-suffix shifts |

---

## Quickstart

### 1. Requirements & Setup

Make sure Python 3.10+ is installed. Activate your virtual environment and install dependencies:

```bash
pip install -r requirements.txt
```

### 2. Run the Evaluation

Execute the all-in-one benchmark script:

```bash
python run_evaluation.py
```

The script will:
- Load the genomic sequence from `data/sample.txt` (NCBI Keratin *Homo sapiens*).
- Benchmark all 4 algorithms across multiple motif lengths ($m = 3$ to $m = 80$).
- Output exact execution times (ms), character comparisons, and match indices to the terminal.
- Generate and save 4 comparison plots into the `results/` directory:
  1. `results/graph1_execution_time.png` (CPU execution time)
  2. `results/graph2_comparisons.png` (Character comparison count)
  3. `results/graph3_trend_pattern_length.png` (Performance vs motif length)
  4. `results/graph4_speedup.png` (Speedup relative to Naive baseline)

---

## Documentation

For full algorithmic derivations, invariants, and detailed empirical observations, refer to:
- **[`REPORT.md`](REPORT.md)** (Markdown format)
- **[`REPORT.txt`](REPORT.txt)** (Plain text format)

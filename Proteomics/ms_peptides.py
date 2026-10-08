#!/usr/bin/env python3
"""
Exact peptide→protein mapper for proteogenomics QC.

Inputs
------
1) Prodigal proteins: prodigal_translated.out           (FASTA)
2) MS peptides:      MS_peptides.fasta             (FASTA)
3) MS table:         MS_table.tsv                  (TSV: peptide,charge,intensity,m_z,rt)

Outputs
-------
- mapping_summary.tsv
- top10_proteins_by_peptide_count.tsv
- figures/
    peptide_length_hist.png
    charge_distribution.png
    coverage_heatmap.png         (coarse 50-bin heatmap of coverage across proteins)
    intensity_vs_mz.png

Run
---
python ms_peptides.py \
  --proteins ~/BIFS619/Bthuringiensis/03_prodigal/prodigal_translated.out \
  --peptides MS_peptides.fasta \
  --ms-table MS_table.tsv
"""
#!/usr/bin/env python3
import argparse, csv, os, sys
from collections import defaultdict, Counter

def fail(msg, code=1):
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(code)

def warn(msg):
    print(f"[WARN] {msg}", file=sys.stderr)

def info(msg):
    print(f"[INFO] {msg}", file=sys.stderr)

def read_fasta(path):
    if not os.path.exists(path):
        fail(f"FASTA not found: {path}")
    seqs, sid, buf = {}, None, []
    with open(path) as fh:
        print(f"Reading Fasta:\n{path}")
        for line in fh:
            line=line.rstrip("\n")
            if not line: 
                continue
            if line.startswith(">"):
                if sid is not None:
                    seqs[sid] = "".join(buf).replace("*","")
                sid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if sid is not None:
            seqs[sid] = "".join(buf).replace("*","")
    if not seqs:
        fail(f"No sequences parsed from FASTA: {path}")
    return seqs

# Accept common header variants
PEP_KEYS = ["peptide","Peptide","sequence","Sequence","peptide_sequence","PEPTIDE"]
CHARGE_KEYS = ["charge","z","Charge"]
INT_KEYS = ["intensity","Intensity","area","Area"]
MZ_KEYS = ["m_z","mz","m/z","MZ"]
RT_KEYS = ["rt","RT","retention_time","RetentionTime"]

def pick_key(d, candidates):
    for k in candidates:
        if k in d: return k
    return None

def read_ms_table(path):
    if not os.path.exists(path):
        fail(f"TSV not found: {path}")
    info_keys = {}
    rows = []
    with open(path, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        if r.fieldnames is None:
            fail("TSV appears to have no header row. Expect a tab-delimited header.")
        # map headers once
        hdr = {h:h for h in r.fieldnames}
        pk = pick_key(hdr, PEP_KEYS)
        if not pk: fail(f"Could not find a peptide column. Saw headers: {r.fieldnames}")
        ck = pick_key(hdr, CHARGE_KEYS)
        ik = pick_key(hdr, INT_KEYS)
        mk = pick_key(hdr, MZ_KEYS)
        rk = pick_key(hdr, RT_KEYS)
        info_keys = {"pep":pk, "charge":ck, "int":ik, "mz":mk, "rt":rk}

        for row in r:
            pep = (row.get(pk,"") or "").strip()
            if not pep: 
                continue
            d = {"peptide": pep}
            # charge
            if ck and row.get(ck,""):
                try: d["charge"] = int(row[ck])
                except: d["charge"] = None
            else:
                d["charge"] = None
            # intensity
            if ik and row.get(ik,""):
                try: d["intensity"] = float(row[ik])
                except: d["intensity"] = None
            else:
                d["intensity"] = None
            # m/z
            if mk and row.get(mk,""):
                try: d["m_z"] = float(row[mk])
                except: d["m_z"] = None
            else:
                d["m_z"] = None
            # rt
            if rk and row.get(rk,""):
                try: d["rt"] = float(row[rk])
                except: d["rt"] = None
            else:
                d["rt"] = None
            rows.append(d)
    if not rows:
        warn("No rows parsed from MS table; proceeding without intensity/mz/rt plots.")
    # keep first occurrence per peptide
    out = {}
    for d in rows:
        if d["peptide"] not in out:
            out[d["peptide"]] = d
    return out

def peptide_find_all(prot_seq, pep):
    starts, start = [], 0
    while True:
        idx = prot_seq.find(pep, start)
        if idx == -1: break
        starts.append(idx); start = idx + 1
    return starts

def main(args):
    print(args)
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    figdir = os.path.join(outdir, "ms_figures")
    os.makedirs(figdir, exist_ok=True)

    proteins = read_fasta(args.proteins)      # id -> AA seq
    peptides = read_fasta(args.peptides)      # id -> AA pep
    msinfo   = read_ms_table(args.ms_table)   # pep -> meta

    # Reverse index and per-protein coverage
    pep_map = defaultdict(list)               # pep_seq -> [(protein_id, start, end)]
    coverage = {pid: [0]*len(seq) for pid, seq in proteins.items()}

    # Map peptides
    for pep_id, pep_seq in peptides.items():
        if not pep_seq:
            continue
        for pid, pseq in proteins.items():
            starts = peptide_find_all(pseq, pep_seq)
            for s in starts:
                e = s + len(pep_seq)
                pep_map[pep_seq].append((pid, s, e))
                for i in range(s, e):
                    coverage[pid][i] = 1

    total_ms_peps = len(peptides)
    mapped_peps   = sum(1 for p in peptides.values() if pep_map.get(p))
    pct_mapped    = 100.0 * mapped_peps / total_ms_peps if total_ms_peps else 0.0

    proteins_with_pep = set(pid for hits in pep_map.values() for pid,_,_ in hits)

    # Distributions
    lengths = [len(seq) for seq in peptides.values() if seq]
    charges = [d.get("charge") for d in msinfo.values() if d.get("charge") is not None]

    # Per-protein stats
    from collections import defaultdict as dd
    pep_counts = Counter()
    cov_stats  = {}
    for pid, arr in coverage.items():
        covered = sum(arr); L = len(arr)
        cov_stats[pid] = (covered, L, (covered / L) if L else 0.0)

    for pep, hits in pep_map.items():
        for pid in set(pid for pid,_,_ in hits):
            pep_counts[pid] += 1

    top10 = pep_counts.most_common(10)

    # ---------- write tables ----------
    with open(os.path.join(outdir, "mapping_summary.tsv"), "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"total_ms_peptides\t{total_ms_peps}\n")
        out.write(f"mapped_peptides\t{mapped_peps}\n")
        out.write(f"percent_mapped\t{pct_mapped:.2f}\n")
        out.write(f"unique_proteins_supported\t{len(proteins_with_pep)}\n")
        unmapped = [p for p in peptides.values() if not pep_map.get(p)]
        out.write(f"unmapped_peptides\t{len(unmapped)}\n")

    with open(os.path.join(outdir, "top10_proteins_by_peptide_count.tsv"), "w") as out:
        out.write("protein_id\tpeptide_count\tcovered_aa\tprotein_len\tcoverage_fraction\n")
        for pid, cnt in top10:
            covered, L, frac = cov_stats[pid]
            out.write(f"{pid}\t{cnt}\t{covered}\t{L}\t{frac:.3f}\n")

    # ---------- plots (robust to empty data) ----------
    try:
        import matplotlib.pyplot as plt
        imported_matplotlib = True
    except Exception as e:
        warn(f"Matplotlib not available; skipping figures ({e})")
        imported_matplotlib = False

    if imported_matplotlib:
        # length hist
        if lengths:
            plt.figure()
            bins = range(max(5, min(lengths)), max(lengths)+2)
            plt.hist(lengths, bins=bins, edgecolor="black")
            plt.xlabel("Peptide length (aa)"); plt.ylabel("Count")
            plt.title("Peptide length distribution")
            plt.savefig(os.path.join(figdir, "peptide_length_hist.png"), bbox_inches="tight")
            plt.close()

        # charge hist
        if charges:
            plt.figure()
            bins = range(min(charges), max(charges)+2)
            plt.hist(charges, bins=bins, edgecolor="black")
            plt.xlabel("Charge state"); plt.ylabel("Count")
            plt.title("Observed charge distribution")
            plt.savefig(os.path.join(figdir, "charge_distribution.png"), bbox_inches="tight")
            plt.close()

        # coverage heatmap
        try:
            import numpy as np
            # choose proteins by coverage fraction
            best = sorted(proteins.keys(), key=lambda p: cov_stats[p][2], reverse=True)[:30]
            if best:
                mat = []
                for pid in best:
                    arr = coverage[pid]
                    L = len(arr)
                    if L == 0:
                        mat.append([0.0]*50)
                        continue
                    bins = []
                    step = max(1, L // 50)
                    for i in range(0, L, step):
                        chunk = arr[i:i+step]
                        bins.append(sum(chunk)/len(chunk))
                    # pad to 50 bins
                    while len(bins) < 50: bins.append(0.0)
                    mat.append(bins[:50])
                mat = np.array(mat)
                plt.figure(figsize=(8, max(4, len(best)*0.25)))
                plt.imshow(mat, aspect="auto", interpolation="nearest")
                plt.colorbar(label="Fraction covered")
                plt.yticks(range(len(best)), best, fontsize=6)
                plt.xlabel("Protein position (binned)")
                plt.title("Per-protein coverage heatmap (top by coverage)")
                plt.savefig(os.path.join(figdir, "coverage_heatmap.png"), bbox_inches="tight")
                plt.close()
        except Exception as e:
            warn(f"Skipping coverage heatmap ({e})")

        # intensity vs m/z
        xs, ys = [], []
        for d in msinfo.values():
            if d.get("m_z") is not None and d.get("intensity") is not None:
                xs.append(d["m_z"]); ys.append(d["intensity"])
        if xs and ys:
            plt.figure()
            plt.scatter(xs, ys, alpha=0.6)
            plt.xlabel("m/z"); plt.ylabel("Intensity")
            plt.title("Intensity vs m/z")
            plt.savefig(os.path.join(figdir, "intensity_vs_mz.png"), bbox_inches="tight")
            plt.close()

    # console preview
    info(f"% mapped: {pct_mapped:.2f} | unique proteins supported: {len(proteins_with_pep)}")
    # Show first 5 peptide → proteins
    shown = 0
    for pep, hits in pep_map.items():
        if hits:
            info(f"PEP {pep} → {len(set(h[0] for h in hits))} protein(s) (e.g., {hits[0][0]})")
            shown += 1
            if shown >= 5: break
    info(f"Wrote outputs to: {outdir}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--proteins", required=True, help="Prodigal protein FASTA (prodigal_translated.out)")
    ap.add_argument("--peptides", required=True, help="MS peptides FASTA")
    ap.add_argument("--ms-table", required=True, help="Peptide table TSV")
    ap.add_argument("--outdir", default="ms_results", help="Output directory (default: ms_results)")
    args = ap.parse_args()
    main(args)

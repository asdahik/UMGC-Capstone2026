#!/usr/bin/env python3
# usage: python3 03_build_ms_table.py PXD020259
import csv, re, sys
from collections import defaultdict
from pathlib import Path

PXD = sys.argv[1]
proj = Path("/home/adahik/BIOT6701/UMGC-Capstone2026/Proteomics")
psm_file = proj / "results" / f"{PXD}_PSM.tsv"
mgf_dir  = proj / "data" / PXD
out_file = proj / "results" / f"{PXD}_MS_table.tsv"
PROTON = 1.007276 # this is to calculate monoisotropic mass

def read_mgf(path):
    blocks, cur = [], None
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if line == "BEGIN IONS":
                cur = {}
            elif line == "END IONS":
                if cur is not None: blocks.append(cur)
                cur = None
            elif cur is not None and "=" in line and not line[0].isdigit():
                k, v = line.split("=", 1)
                if k == "PEPMASS":
                    f = v.split()
                    cur["mz"] = float(f[0])
                    cur["inten"] = float(f[1]) if len(f) > 1 else None
                elif k == "CHARGE":
                    m = re.search(r"\d+", v)
                    cur["charge"] = int(m.group()) if m else None
                elif k == "SCANS":
                    cur["scan"] = v
    return blocks

def locate(sid, blocks, scans, off):
    m = re.match(r"(?:index|query)=(\d+)$", sid or "")
    if m:
        i = int(m.group(1)) + off
        return blocks[i] if 0 <= i < len(blocks) else None
    m = re.match(r"scan=(\d+)$", sid or "")
    return scans.get(m.group(1)) if m else None

def matches(b, row):
    return (b is not None and b.get("mz") is not None
            and abs(b["mz"] - float(row["exp_mz"])) < 0.02
            and b.get("charge") == int(row["charge"]))

def mascot_score(s):
    m = re.search(r"[Ss]core[^=;]*=([\d.]+)", s or "")
    return float(m.group(1)) if m else 0.0

# group confident PSMs by MGF file
by_mgf = defaultdict(list)
with open(psm_file) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["passThreshold"].lower() == "true" and r["peptide"]:
            by_mgf[r["mgf"]].append(r)

best, tot = {}, defaultdict(int)       # (peptide, charge) -> (score, row, block)
for mgf, rows in sorted(by_mgf.items()):
    blocks = read_mgf(mgf_dir / mgf)
    scans = {b["scan"]: b for b in blocks if "scan" in b}
    # pick the index offset that matches precursor m/z best on a sample
    off = max((0, 1, -1), key=lambda o: sum(
        matches(locate(r["spectrumID"], blocks, scans, o), r) for r in rows[:300]))
    for r in rows:
        b = locate(r["spectrumID"], blocks, scans, off)
        if b is None:
            tot["unresolved"] += 1; continue
        if not matches(b, r):
            tot["mismatch"] += 1; continue
        tot["ok"] += 1
        key = (r["peptide"], r["charge"])
        sc = mascot_score(r["scores"])
        if key not in best or sc > best[key][0]:
            best[key] = (sc, r, b)
    print(mgf, "offset", off, file=sys.stderr)

with open(out_file, "w") as out:
    out.write("Peptide_ID\tSequence\tProtein_Header\tCharge\tMonoisotopic_Mass\tm/z\tIntensity\n")
    for n, ((pep, ch), (sc, r, b)) in enumerate(sorted(best.items()), 1):
        z, mz = int(ch), float(r["exp_mz"])
        inten = "NA" if b["inten"] is None else f'{b["inten"]:.0f}'
        prot = r["proteins"].split(";")[0]
        out.write(f"MS_{n:05d}\t{pep}\t{prot}\t{z}\t{mz*z - z*PROTON:.5f}\t{mz:.5f}\t{inten}\n")

print(dict(tot), len(best), "rows written", file=sys.stderr)
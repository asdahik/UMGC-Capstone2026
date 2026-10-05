#!/bin/env python3

import gzip, sys, csv, re
import xml.etree.ElementTree as ET
from pathlib import Path


def local(t): return t.rsplit("}", 1)[-1]

def parse(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    peptides, dbseq, evid = {}, {}, {}
    spectra_loc = {}
    with opener(path, "rb") as fh:
        for _, el in ET.iterparse(fh, events=("end",)):
            t = local(el.tag)
            if t == "Peptide":
                seq = next((c.text for c in el if local(c.tag) == "PeptideSequence"), "")
                mods = [f'{m.get("location")}:{m.get("monoisotopicMassDelta")}'
                        for m in el if local(m.tag) == "Modification"]
                peptides[el.get("id")] = (seq, ";".join(mods)); el.clear()
            elif t == "DBSequence":
                dbseq[el.get("id")] = el.get("accession"); el.clear()
            elif t == "PeptideEvidence":
                evid[el.get("id")] = dbseq.get(el.get("dBSequence_ref"), el.get("dBSequence_ref")); el.clear()
            elif t == "SpectraData":
                spectra_loc[el.get("id")] = Path(el.get("location", "")).name
            elif t == "SpectrumIdentificationResult":
                res_params = {c.get("name"): c.get("value") for c in el if local(c.tag) == "cvParam"}
                for item in el:
                    if local(item.tag) != "SpectrumIdentificationItem": continue
                    if item.get("rank") not in (None, "1"): continue
                    seq, mods = peptides.get(item.get("peptide_ref"), ("", ""))
                    prots = sorted({evid.get(r.get("peptideEvidence_ref"), "") for r in item
                                    if local(r.tag) == "PeptideEvidenceRef"})
                    ip = {c.get("name"): c.get("value") for c in item if local(c.tag) == "cvParam"}
                    yield {
                        "mzid": Path(path).name,
                        "mgf": spectra_loc.get(el.get("spectraData_ref"), ""),
                        "spectrumID": el.get("spectrumID"),
                        "result_params": ";".join(f"{k}={v}" for k, v in res_params.items()),
                        "charge": item.get("chargeState"),
                        "exp_mz": item.get("experimentalMassToCharge"),
                        "peptide": seq, "mods": mods,
                        "passThreshold": item.get("passThreshold"),
                        "scores": ";".join(f"{k}={v}" for k, v in ip.items()),
                        "proteins": ";".join(prots),
                    }
                el.clear()

out = csv.DictWriter(open(sys.argv[1], "w", newline=""), delimiter="\t",
    fieldnames=["mzid","mgf","spectrumID","result_params","charge","exp_mz","peptide",
                "mods","passThreshold","scores","proteins"])
out.writeheader()
for p in sys.argv[2:]:
    n = 0
    for row in parse(p): out.writerow(row); n += 1
    print(p, n, "PSMs", file=sys.stderr)
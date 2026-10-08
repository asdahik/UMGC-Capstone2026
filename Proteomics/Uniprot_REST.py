import json, time, urllib.parse, urllib.request, sys
from pathlib import Path

main_dir = Path.cwd()

# txt file containing list of ascension ids 
print(sys.argv)
asc_txt = str(sys.argv[1])
print(asc_txt)
# slice the PDX
pdx=asc_txt[asc_txt.index("PXD"):asc_txt.index("_Prot")]
print(pdx)
output_dir = sys.argv[2]



BASE = "https://rest.uniprot.org"

with open(main_dir/asc_txt, "r") as query:
    ids = query.read().split()
    # web request
    req = urllib.request.Request(f"{BASE}/idmapping/run", method="POST",
        data=urllib.parse.urlencode({"from": "UniProtKB_AC-ID", "to": "UniProtKB",
                                     "ids": ",".join(ids)}).encode())
    
    job = json.load(urllib.request.urlopen(req))["jobId"]

    while True:
        s = json.load(urllib.request.urlopen(f"{BASE}/idmapping/status/{job}"))
        if "jobStatus" not in s:break # finished
        time.sleep(3)

    # output an annotation csv file and an annotation fasta
    for fmt, fields, out in [  
        ("tsv", "&fields=accession,id,protein_name,gene_names,organism_name,length", f"{output_dir}/{pdx}_annotation.tsv"),
        ("fasta", "", f"{output_dir}/{pdx}_proteins.out")]:
        url = f"{BASE}/idmapping/uniprotkb/results/stream/{job}?format={fmt}{fields}"
        open(out, "wb").write(urllib.request.urlopen(url).read())
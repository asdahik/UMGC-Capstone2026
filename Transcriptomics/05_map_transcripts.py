import pandas as pd
import sys

res_dir=sys.argv[1]
ref_dir=sys.argv[2]

res = pd.read_csv(f"{res_dir}/DE_results.csv", index_col=0)
ann = pd.read_csv(f"{ref_dir}/potato_transcript_key.tsv", sep="\t", header=None,
                  names=["transcript", "annotation"], quoting=3)

# strip the .1 / .2 isoform suffix to get the locus ID
ann["gene"] = ann["transcript"].str.replace(r"\.\d+$", "", regex=True)

# one annotation per locus, preferring .1
ann = ann.sort_values("transcript").drop_duplicates("gene")

res["annotation"] = res.index.map(ann.set_index("gene")["annotation"])
res = res.rename(columns={"A": "SpudDBGeneID"})
res.to_csv(f"{res_dir}/deseq_results_annotated.csv")
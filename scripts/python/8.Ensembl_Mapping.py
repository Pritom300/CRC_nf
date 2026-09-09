import sys
import os
import json
import time
import requests
import pandas as pd

if len(sys.argv) < 3:
    print("Error: Missing arguments!")
    print("Usage: python3 8.Ensembl_Mapping.py <importance_csv> <output_dir>")
    sys.exit(1)

csv_path   = sys.argv[1]
output_dir = sys.argv[2]

json_path = os.path.join(output_dir, "target_gene_ensg_mapping.json")
missing_path = os.path.join(output_dir, "unresolved_genes.txt")

if not os.path.exists(csv_path):
    raise FileNotFoundError(
        f"CSV did not found: {csv_path}\n"
        "Colab"
    )


df = pd.read_csv(csv_path)

if "Gene" not in df.columns:
    raise ValueError(f"CSV has not gene columns: {list(df.columns)}")


genes = list(dict.fromkeys(
    df["Gene"]
    .dropna()
    .astype(str)
    .str.strip()
    .str.upper()
))

print(f"unique gene: {len(genes)}")

session = requests.Session()
session.headers.update({
    "Accept": "application/json",
    "User-Agent": "Google-Colab-Gene-Mapping"
})

mapping = {}
unresolved = []


def get_ensembl_id(gene):
    
    url = "https://mygene.info/v3/query"
    params = {
        "q": gene,
        "scopes": "symbol",
        "species": "human",
        "fields": "symbol,ensembl",
        "size": 10,
        "operator": "and"
    }

    for attempt in range(5 ):
        try:
            r = session.get(url, params=params, timeout=60)

            if r.status_code == 200:
                data = r.json()
                hits = data.get("hits", [])

                # Exact gene-symbol match 
                hits = sorted(
                    hits,
                    key=lambda x: str(x.get("symbol", "")).upper() != gene
                )

                for hit in hits:
                    ensembl = hit.get("ensembl")

                    if isinstance(ensembl, dict):
                        ensembl_id = ensembl.get("gene")
                        if isinstance(ensembl_id, str) and ensembl_id.startswith("ENSG"):
                            return ensembl_id

                    elif isinstance(ensembl, list):
                        for item in ensembl:
                            if isinstance(item, dict):
                                ensembl_id = item.get("gene")
                                if isinstance(ensembl_id, str) and ensembl_id.startswith("ENSG"):
                                    return ensembl_id

                return None

            if r.status_code in [429, 500, 502, 503, 504]:
                wait = min(2 ** attempt, 20)
                print(f"{gene}: HTTP {r.status_code}, retrying after {wait}s")
                time.sleep(wait)
                continue

            print(f"{gene}: HTTP {r.status_code}")
            return None

        except requests.RequestException as e:
            wait = min(2 ** attempt, 20)
            print(f"{gene}: connection error, retrying after {wait}s")
            time.sleep(wait)

    return None



for i, gene in enumerate(genes, start=1):
    print(f"[{i}/{len(genes)}] {gene}")

    ensembl_id = get_ensembl_id(gene)

    if ensembl_id:
        mapping[gene] = ensembl_id
    else:
        unresolved.append(gene)
        print(f"  ID did not find: {gene}")

    time.sleep(0.2)


with open(json_path, "w", encoding="utf-8") as f:
    json.dump(mapping, f, indent=4, ensure_ascii=False)
    f.write("\n")


with open(missing_path, "w", encoding="utf-8") as f:
    f.write("\n".join(unresolved))
    if unresolved:
        f.write("\n")

print("\nCompleted")
print(f"Total gene: {len(genes)}")
print(f"Mapping Found: {len(mapping)}")
print(f"Mapping did not found: {len(unresolved)}")
print(f"JSON : {json_path}")

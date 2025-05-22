# algaetax; db_queries.py
# Import libraries
import os
import pandas as pd
import time
import requests
# Import functions
from Bio import Entrez
from concurrent.futures import ThreadPoolExecutor, as_completed
from utils import is_database_enabled, load_config


# Load config
config = load_config()

# Set NCBI email and API key
Entrez.email = config["db_ncbi"].get("ncbi_email", "anonymous@ncbi.com")
api_key = config["db_ncbi"].get("api_key")
if api_key:
    Entrez.api_key = api_key

# Set path to PR2 database file
PR2_DB_PATH = config["database_path"]["db_pr2"]


# Load PR2 database from Excel
def load_pr2_db():
    if not os.path.exists(PR2_DB_PATH):
        print("ERROR: PR2 database file not found!")
        exit(1)
    try:
        df = pd.read_excel(PR2_DB_PATH, sheet_name=0)
        print(f"\n[at] PR2_DB loaded: {PR2_DB_PATH}\n")
        return df
    except Exception as e:
        print(f"ERROR loading PR2 data: {e}")
        exit(1)

# Search PR2 database for matching genus or species
def query_pr2(taxon, pr2_data):
    if pr2_data is None:
        return "No PR2 data loaded"

    t = taxon.replace(" ", "_")
    matches = pr2_data[
        (pr2_data["genus"] == taxon) | (pr2_data["genus"] == t) |
        (pr2_data["species"] == taxon) | (pr2_data["species"] == t)
    ]

    if not matches.empty:
        cols = ["domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species"]
        return matches[cols].to_dict(orient="records")
    else:
        return "Not found"


# Query NCBI for taxon IDs (multithreaded)
def fetch_taxon_id(taxon, delay):
    time.sleep(delay)
    try:
        handle = Entrez.esearch(db="taxonomy", term=taxon)
        record = Entrez.read(handle)
        return taxon, record["IdList"][0] if record["IdList"] else None
    except Exception as e:
        print(f"Error fetching ID for {taxon}: {e}")
        return taxon, None


# Parallel ID lookup using multiple threads
def get_taxon_ids_parallel(taxa_list):
    taxon_ids = {}
    delay = 0.5 if api_key else 0.4
    max_workers = 3 if api_key else 2

    print(f"{'-' * 20}")
    print(f"[NCBI] Start fetching taxon IDs (Delay: {delay}s, Threads: {max_workers})\n")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(fetch_taxon_id, taxon, delay): taxon for taxon in taxa_list}
        for i, future in enumerate(as_completed(futures), 1):
            taxon, tax_id = future.result()
            taxon_ids[taxon] = tax_id
            print(f"[{i}/{len(taxa_list)}] NCBI: {taxon}")

    return taxon_ids


# Query NCBI in batches to retrieve taxonomy info
def batch_query_ncbi(taxon_ids, batch_size=10):
    results = {}
    delay = 0.6 if api_key else 1.0
    batch_size = 40 if api_key else 10

    print(f"\n{'-' * 20}")
    print("[NCBI] Start API batch query...\n")
    taxon_names = list(taxon_ids.keys())

    for i in range(0, len(taxon_names), batch_size):
        batch_taxa = taxon_names[i:i + batch_size]
        batch_ids = [taxon_ids[t] for t in batch_taxa if taxon_ids[t]]

        if not batch_ids:
            continue

        print(f"[Process] Batch {i // batch_size + 1}")
        time.sleep(delay)

        try:
            handle = Entrez.efetch(db="taxonomy", id=",".join(batch_ids), retmode="xml")
            records = Entrez.read(handle)

            for record in records:
                tax_id = record["TaxId"]
                lineage = record.get("Lineage", "Not found")
                matched_taxon = next((name for name, tid in taxon_ids.items() if tid == tax_id), None)
                if matched_taxon:
                    results[matched_taxon] = {
                        "NCBI_ID": tax_id,
                        "Taxonomy": lineage
                    }

        except Exception as e:
            print(f"Error with batch {batch_ids}: {e}")

    return results


# Query AlgaeBase using API (genus or species)
def query_algaebase(taxon):
    api_key = config["db_algaebase"].get("api_key", "")
    genus_url = config["db_algaebase"].get("api_url", "").replace("/species", "/genus")
    species_url = config["db_algaebase"].get("api_url", "")

    result = {
        "ALGB_Status": "Not found",
        "ALGB_ID": "",
        "ALGB_Lineage": ""
    }

    parts = taxon.strip().split()
    if not parts:
        return result

    try:
        if len(parts) == 1:
            url = genus_url
            params = {"genus": f"[eq]{parts[0]}"}
        elif len(parts) == 2:
            url = species_url
            params = {
                "genus": f"[eq]{parts[0]}",
                "specificEpithet": f"[eq]{parts[1]}"
            }
        else:
            return result

        # Optional debug output (activate when needed)
        # prepared = requests.Request("GET", url, headers={"abapikey": api_key}, params=params).prepare()
        # print(f"Request: {prepared.url}\n")

        response = requests.get(url, headers={"abapikey": api_key}, params=params)
        if not response.ok:
            return result  # skip silently if 404 or other error

        data = response.json()

        if data.get("result"):
            entry = data["result"][0]
            result["ALGB_Status"] = "Found"
            result["ALGB_ID"] = str(entry.get("dwc:scientificNameID", ""))

            lineage_parts = []
            for rank in ["dwc:class", "dwc:order", "dwc:family", "dwc:genus", "dwc:scientificName"]:
                val = entry.get(rank)
                if val:
                    lineage_parts.append(val.strip())
            result["ALGB_Lineage"] = "; ".join(lineage_parts)

    except Exception as e:
        return result
        # print(f"AlgaeBase error for {taxon}: {e}")

    return result

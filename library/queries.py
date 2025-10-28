# algaetax; file: queries.py

# Import libraries
import os, time, requests
import pandas as pd
from Bio import Entrez
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import functions from library.utils
from library.utils import ( 
    load_config, 
    extract_valid_taxon
    )



# Load configuration
config = load_config()

# Configure NCBI Entrez API
Entrez.email = config["db_ncbi"].get("ncbi_email", "anonymous@ncbi.com")
api_key = config["db_ncbi"].get("api_key")
if api_key:
    Entrez.api_key = api_key

# Path to PR2 file
PR2_DB_PATH = config["database_path"]["db_pr2"]



# Load PR2 database from Excel file.
def load_pr2_db():
    if not os.path.exists(PR2_DB_PATH):
        print("[ERROR] PR2 database file not found!")
        exit(1)
    try:
        df = pd.read_excel(PR2_DB_PATH, sheet_name=0)
        print(f"\n[algaetax] PR2_DB loaded: \n{PR2_DB_PATH}\n")
        return df
    except Exception as e:
        print(f"[ERROR] Failed to load PR2 data: {e}")
        exit(1)


# Search PR2 for matching genus or species.
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


# Query NCBI for a taxon ID.
def fetch_taxon_id(taxon, delay):
    time.sleep(delay)
    try:
        handle = Entrez.esearch(db="taxonomy", term=taxon)
        record = Entrez.read(handle)
        return taxon, record["IdList"][0] if record["IdList"] else None
    except Exception as e:
        print(f"[ERROR] Failed to fetch ID for {taxon}: {e}")
        return taxon, None


# Run fetch_taxon_id() for multiple taxa in parallel.
def get_taxon_ids_parallel(taxa_list):
    taxon_ids = {}
    delay = 0.5 if api_key else 0.4
    max_workers = 3 if api_key else 2

    print(f"\n\n{'-' * 20}")
    print(f"[NCBI] Start fetching taxon IDs \nDelay: {delay}s; Threads: {max_workers}\n")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(fetch_taxon_id, taxon, delay): taxon for taxon in taxa_list}
        for i, future in enumerate(as_completed(futures), 1):
            taxon, tax_id = future.result()
            taxon_ids[taxon] = tax_id
            print(f"[{i}/{len(taxa_list)}] NCBI: {taxon}")
    return taxon_ids


# Fetch full taxonomy from NCBI using taxon IDs in batches.
def batch_query_ncbi(taxon_ids, batch_size=10):
    results = {}
    delay = 0.6 if api_key else 1.0
    batch_size = 40 if api_key else 10

    print("\n[NCBI] Start API batch query...")
    taxon_names = list(taxon_ids.keys())

    for i in range(0, len(taxon_names), batch_size):
        batch_taxa = taxon_names[i:i + batch_size]
        batch_ids = [taxon_ids[t] for t in batch_taxa if taxon_ids[t]]
        if not batch_ids:
            continue

        print(f"Process: Batch {i // batch_size + 1}")
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


# Query AlgaeBase API for genus/species info and full lineage.
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
        if len(parts) == 2:
            params = {
                "genus": f"[eq]{parts[0]}",
                "specificEpithet": f"[eq]{parts[1]}"
            }

            responseSP = requests.get(species_url, headers={"abapikey": api_key}, params=params)
            if not responseSP.ok:
                return result
            dataSP = responseSP.json()
            if dataSP.get("result"):
                entrySP = dataSP["result"][0]
                result["ALGB_Name_status"] = entrySP.get("dwc:taxonomicStatus")
                if result["ALGB_Name_status"] == "currently accepted taxonomically":
                    result["ALGB_CurrentName"] = parts[0] + " " + parts[1]
                else:
                    result["ALGB_CurrentName"] = extract_valid_taxon(entrySP.get("dwc:acceptedNameUsage"))
        elif len(parts) == 1:
            pass  # valid genus-only query; handled below
        else:
            return result
        
        params = {"genus": f"[eq]{parts[0]}"}
        response = requests.get(genus_url, headers={"abapikey": api_key}, params=params)  
        if not response.ok:
            return result

        data = response.json()
        if data.get("result"):
            entry = data["result"][0]
            result["ALGB_Status"] = "Found"
            result["ALGB_ID"] = str(entry.get("dwc:scientificNameID", ""))
            result["ALGB_Empire"] = str(entry.get("dwc:empire", ""))
            result["ALGB_Kingdom"] = str(entry.get("dwc:kingdom", ""))
            result["ALGB_Phylum"] = str(entry.get("dwc:phylum", ""))
            result["ALGB_Subphylum"] = str(entry.get("dwc:subphylum", ""))
            result["ALGB_Class"] = str(entry.get("dwc:class", ""))
            result["ALGB_Order"] = str(entry.get("dwc:order", ""))
            result["ALGB_Family"] = str(entry.get("dwc:family", ""))
            result["ALGB_Genus"] = str(entry.get("dwc:genus", ""))
            result["ALGB_scientificName"] = str(entry.get("dwc:scientificName", ""))
            lineage_parts = []
            for rank in [
                "dwc:empire", "dwc:kingdom", "dwc:phylum", "dwc:subphylum", 
                "dwc:class", "dwc:order", "dwc:family", "dwc:genus", "dwc:scientificName"
            ]:
                val = entry.get(rank)
                if val:
                    lineage_parts.append(val.strip())
            result["ALGB_Lineage"] = "; ".join(lineage_parts)

    except Exception as e:
        return result
    return result

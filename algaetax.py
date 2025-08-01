# algaetax; MAIN_file: algaetax.py
# Main script to run taxonomic queries

# Import libraries
import os, sys, time, csv
import shutil, requests
import pandas as pd
from collections import Counter
from pathlib import Path

# Import metadata from library.metadata
from library.metadata import NAME, VERSION, SUPPORTED_DATABASES

# Import functions from library.utils
from library.utils import (
    load_config,
    load_blacklist,
    load_taxa_from_excel,
    is_database_enabled,
    extract_valid_taxon,
    safe_save_path
)

# Import functions from library.queries
from library.queries import (
    load_pr2_db,
    batch_query_ncbi,
    query_pr2,
    query_algaebase,
    get_taxon_ids_parallel,
)



# Load config and blacklist
config = load_config()
blacklist = set(load_blacklist(config))

# Set output directory and output file name
OUTPUT_DIR = Path(config["general"]["output_dir"])
OUTPUT_FILE_CSV = "algaetax_taxa_results.csv"
output_file = OUTPUT_DIR/OUTPUT_FILE_CSV

# Check if output file already exists
output_file = safe_save_path(output_file)

# Start runtime timer
start_time = time.time()

# Backup config.yaml if enabled
if config["general"].get("backup_config", False):
    config_source_path = Path("config.yaml")
    backup_dir = OUTPUT_DIR / "backups"
    config_target_path = backup_dir / "config.yaml"

    try:
        backup_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(config_source_path, config_target_path)
    except Exception as e:
        print(f"[WARN] Could not copy config.yaml: {e}")



# FUNCTION: STARTUP_INFO
# Print metadata and active database info
def startup_info(config):
    print(f"\n\n{NAME} {VERSION}\n")
    print("Metadata for Databases:")
    if is_database_enabled("NCBI", config):
        email = config["db_ncbi"].get("ncbi_email")
        api_key = config["db_ncbi"].get("api_key")
        print(f"[NCBI_DB] Email: {email or '[WARN] Not set in config'}")
        print(f"[NCBI_DB] API Key: {'Set in config' if api_key else '[WARN] Not set in config'}")

    if is_database_enabled("PR2", config):
        path = config["database_path"].get("db_pr2", "[WARN] Not set in config")
        print(f"\n[PR2_DB] Path: \n{path}")

    if is_database_enabled("ALGB", config):
        api_key = config["db_algaebase"].get("api_key", "").strip()
        api_url = config["db_algaebase"].get("api_url", "").strip()

        # Print info about API key and URL
        print(f"\n[ALGB_DB] API URL: \n{api_url or '[WARN] Not set in config'}")
        print(f"[ALGB_DB] API Key: {'Set in config' if api_key else '[WARN] Not set in config'}")

        # Case 1: No API key set: warn user and allow disabling AlgaeBase
        if not api_key:
            print("\n\n[WARN] No AlgaeBase API key set in config.yaml.")
            choice = input("Continue without AlgaeBase DB? (y/yes or n/no): ").strip().lower()
            if choice not in ["y", "yes"]:
                print("\n[INFO] algaetax stopped by user (no API key).")
                sys.exit(0)
            else:
                # Disable AlgaeBase for this run to avoid invalid API calls later
                config["database"]["ALGB"] = False
                print("\n\n[INFO] AlgaeBase DB disabled for this run (no API key).")
        else:
            # Case 2: API key is set: test if it works with a common genus
            try:
                print("\n[algaetax] Checking AlgaeBase API key...")
                test_params = {"genus": "[eq]Chlorella"}
                r = requests.get(api_url.replace("/species", "/genus"),
                                headers={"abapikey": api_key},
                                params=test_params,
                                timeout=10)

                # If request fails or returns no results: key likely invalid
                if not r.ok or not r.json().get("result"):
                    print("\n\n[ERROR] AlgaeBase API key seems invalid or expired.")
                    choice = input("Continue anyway without AlgaeBase DB? (y/yes or n/no): ").strip().lower()
                    if choice not in ["y", "yes"]:
                        print("\n[INFO] algaetax stopped by user (invalid API key).")
                        sys.exit(0)
                    else:
                        # Disable AlgaeBase for this run to avoid invalid API calls later
                        config["database"]["ALGB"] = False
                        print("\n[INFO] AlgaeBase DB disabled for this run (invalid API key).")
            except Exception as e:
                print(f"\n\n[ERROR] Could not verify AlgaeBase API key: {e}")
                sys.exit(1)

    # List all active databases
    print("\nEnabled Databases:")
    active = [db for db, enabled in config["database"].items() if enabled]
    if active:
        for db in active:
            name = SUPPORTED_DATABASES.get(db, "Unknown")
            print(f"{db}: {name}")
    else:
        print("[WARN] No databases enabled. Check config.yaml.")
        sys.exit(1)



# FUNCTION: QUERY_AND_SAVE
# Run query workflow and write results
def query_and_save(taxa_list, output_file, config):
    pr2_data = load_pr2_db() if is_database_enabled("PR2", config) else None
    taxon_ids = get_taxon_ids_parallel(taxa_list) if is_database_enabled("NCBI", config) else {}
    ncbi_results = batch_query_ncbi(taxon_ids) if is_database_enabled("NCBI", config) else {}

    if is_database_enabled("PR2", config):
        print(f"\n\n{'-' * 20}")
        print("[PR2] Start LOCAL DB query...")
    if is_database_enabled("ALGB", config):
        print(f"\n\n{'-' * 20}")
        print("[ALGB] Start API query...\n")

    results = []
    pr2_counter = 0
    algb_counter = 0

    for taxon in taxa_list:
        taxon_result = {}
        
        # Fetch results for each active database
        if is_database_enabled("NCBI", config):
            taxon_result["NCBI"] = ncbi_results.get(taxon, "Not found")
        
        if is_database_enabled("PR2", config):
            pr2_counter += 1
            pr2_result = query_pr2(taxon, pr2_data)
            taxon_result["PR2"] = pr2_result
        
        if is_database_enabled("ALGB", config):
            algb_counter += 1
            print(f"[{algb_counter}/{len(taxa_list)}] ALGB: {taxon}")
            algae_result = query_algaebase(taxon)
            taxon_result["ALGB"] = algae_result

        results.append((taxon, taxon_result))

    # Save results
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    write_results_to_csv(results, output_file, config)

    # Runtime summary
    elapsed = time.time() - start_time
    minutes = int(elapsed // 60)
    seconds = int(elapsed % 60)
    
    print(f"\n\n{'-' * 20}")
    print(f"[algaetax] All queries completed.")
    print(f"Runtime: {minutes} min {seconds} sec; {elapsed:.2f} seconds")
    print(f"\n[algaetax] Results saved to: \n{output_file}\n\n")


# FUNCTION: WRITE_RESULTS_TO_CSV
# Format and write query results to CSV file
def write_results_to_csv(results_list, output_file, config):
    rows = []

    ncbi_enabled = is_database_enabled("NCBI", config)
    pr2_enabled = is_database_enabled("PR2", config)
    algae_enabled = is_database_enabled("ALGB", config)

    for taxon, result in results_list:
        row = {"Taxon": taxon}

        # Process NCBI results
        if ncbi_enabled:
            ncbi = result.get("NCBI")
            if isinstance(ncbi, dict):
                row["NCBI_Status"] = "Found"
                row["NCBI_ID"] = str(ncbi.get("NCBI_ID", ""))
                row["NCBI_Lineage"] = ncbi.get("Taxonomy", "")
            else:
                row["NCBI_Status"] = "Not found"

        # Process PR2 results
        if pr2_enabled:
            pr2 = result.get("PR2")
            if isinstance(pr2, list) and pr2:
                row["PR2_Status"] = "Found"
                # Build lineage string from PR2 fields
                entry = pr2[0]
                lineage_parts = []
                for key in ["domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus"]:
                    val = entry.get(key, "").strip()
                    if val:
                        lineage_parts.append(val)
                row["PR2_Lineage"] = "; ".join(lineage_parts)
                row["PR2_domain"] = entry.get("domain", "")
                row["PR2_supergroup"] = entry.get("supergroup", "")
                row["PR2_division"] = entry.get("division", "")
                row["PR2_subdivision"] = entry.get("subdivision", "")
                row["PR2_class"] = entry.get("class", "")
                row["PR2_order"] = entry.get("order", "")
                row["PR2_family"] = entry.get("family", "")
                row["PR2_genus"] = entry.get("genus", "")
            else:
                row["PR2_Status"] = "Not found"

        # Process AlgaeBase result
        if algae_enabled:
            algae = result.get("ALGB")
            if isinstance(algae, dict):
                row["ALGB_Status"] = algae.get("ALGB_Status", "")
                row["ALGB_Name_status"] = algae.get("ALGB_Name_status", "")
                row["ALGB_CurrentName"] = algae.get("ALGB_CurrentName", "")
                row["ALGB_ID"] = algae.get("ALGB_ID", "")
                row["ALGB_Empire"] = algae.get("ALGB_Empire", "")
                row["ALGB_Kingdom"] = algae.get("ALGB_Kingdom", "")
                row["ALGB_Phylum"] = algae.get("ALGB_Phylum", "")
                row["ALGB_Subphylum"] = algae.get("ALGB_Subphylum", "")
                row["ALGB_Class"] = algae.get("ALGB_Class", "")
                row["ALGB_Order"] = algae.get("ALGB_Order", "")
                row["ALGB_Family"] = algae.get("ALGB_Family", "")
                row["ALGB_Genus"] = algae.get("ALGB_Genus", "")
                row["ALGB_scientificName"] = algae.get("ALGB_scientificName", "")
                row["ALGB_Lineage"] = algae.get("ALGB_Lineage", "")
            else:
                row["ALGB_Status"] = "Not found"

        rows.append(row)

    df = pd.DataFrame(rows)

    # Set column order based on enabled databases
    columns = ["Taxon"]
    if ncbi_enabled:
        columns += ["NCBI_Status", "NCBI_ID", "NCBI_Lineage"]
    if pr2_enabled:
        columns += ["PR2_Status", "PR2_domain", "PR2_supergroup", "PR2_division",
                    "PR2_subdivision", "PR2_class", "PR2_order", "PR2_family", "PR2_genus", "PR2_Lineage"]
    if algae_enabled:
        columns += ["ALGB_Status", "ALGB_Name_status", "ALGB_CurrentName", "ALGB_ID",
                    "ALGB_Empire", "ALGB_Kingdom", "ALGB_Phylum", "ALGB_Subphylum", "ALGB_Class",
                    "ALGB_Order", "ALGB_Family", "ALGB_Genus", "ALGB_scientificName", "ALGB_Lineage"]
    df = df[columns]

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    df.to_csv(output_file, index=False, encoding="utf-8", sep=";", quoting=csv.QUOTE_MINIMAL)

    # Optional: Export taxa where all ACTIVE databases returned "Not found"
    if config["general"].get("export_not_found_taxa_list", False):
        status_cols = [col for col in df.columns if "Status" in col]
        enabled_cols = [col for col in status_cols if col.split("_")[0] in config["database"] and config["database"][col.split("_")[0]]]

        if enabled_cols:
            all_not_found_df = df[(df[enabled_cols] == "Not found").all(axis=1)]
            all_not_found_df = all_not_found_df[["Taxon"] + enabled_cols]

            not_found_path = os.path.join(os.path.dirname(output_file), "taxa_not_found.csv")

            if not all_not_found_df.empty:
                all_not_found_df = all_not_found_df[["Taxon"] + enabled_cols]
                all_not_found_df.to_csv(not_found_path, index=False, encoding="utf-8", sep=";")
                print(f"\n\n\n{'-' * 20}")
                print(f"[algaetax] 'Not found' taxa saved: \n{not_found_path}")
            else:
                empty_df = pd.DataFrame(columns=["Taxon"] + enabled_cols)
                empty_df.to_csv(not_found_path, index=False, encoding="utf-8", sep=";")
                print(f"\n[algaetax] No taxa missing in all DBs. Empty file: {not_found_path}")



# RUN_MAIN_PROCESS
# Main execution: load data, clean/filter, and query taxonomic databases
if __name__ == "__main__":
    # Print metadata and check active databases
    startup_info(config)
    print(f"\n\n\n{'-' * 20}")
    print("[algaetax] RUNNING QUERY...\n")

    # Load taxa list from Excel file (raw data)
    input_file = config["general"]["input_data"]
    taxa_col_number = config["general"]["taxa_column_number"]
    df, taxa_list = load_taxa_from_excel(input_file, taxa_col_number, output_file, config)

    # Clean and filter taxa based on blacklist
    valid_taxa = []
    invalid_taxa = []
    filter_log = []

    for taxon in taxa_list:
        cleaned = extract_valid_taxon(taxon, blacklist) # Remove unwanted terms
        if cleaned:
            valid_taxa.append(cleaned)
            if cleaned != taxon.strip():
                # Track cleaned values
                filter_log.append((taxon, cleaned)) # Keep original cleaned log
        else:
            invalid_taxa.append(taxon) # Entirely excluded

    # Check for duplicate taxa in original input
    taxa_index = taxa_col_number - 1
    raw_column = df.iloc[:, taxa_index].dropna().astype(str).str.strip()
    counter = Counter(raw_column)
    duplicates = {k: v for k, v in counter.items() if v > 1}


    if duplicates:
        print("\n[WARN] Duplicate taxa found:")
        for taxon, count in sorted(duplicates.items(), key=lambda x: x[1], reverse=True):
            print(f"-{taxon}: {count}×")
    else:
        print("\n[algaetax] No duplicate taxa found.")

    # Summary of loaded and cleaned taxa
    print(f"\n[algaetax] TAXA loaded: {len(taxa_list)}")
    print(f"Preview: {taxa_list[:2]} ...")
    print(f"Filtered: {valid_taxa[:2]} ...")


    # Save filter log (original with cleaned names)
    if filter_log:
        filepath = os.path.join(OUTPUT_DIR, "filtered_taxa_log.csv")
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        with open(OUTPUT_DIR / "filtered_taxa_log.csv", "w", encoding="utf-8") as f:
            f.write("Taxa,Filtered\n")
            for orig, filt in filter_log:
                f.write(f"{orig},{filt}\n")
        print(f"\n[algaetax] Writing filtered taxa to file: \n{filepath}")



    # Query all enabled databases and save results
    query_and_save(valid_taxa, output_file, config)

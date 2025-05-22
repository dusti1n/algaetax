# algaetax; algaetax.py
# Import libraries
import os, sys, time, csv
import pandas as pd
# Import functions
from collections import Counter
from pathlib import Path
from utils import load_config, is_database_enabled, load_taxa_from_excel, col_letter_to_index, extract_valid_taxon
from db_queries import batch_query_ncbi, get_taxon_ids_parallel, load_pr2_db, query_pr2, query_algaebase
from metadata import NAME, VERSION, SUPPORTED_DATABASES


# Start execution timer
start_time = time.time()

# Load config
config = load_config()

# Define output directory and filename
OUTPUT_DIR = Path(config["general"]["output_dir"])
OUTPUT_FILE_CSV = "algaetax_taxa_results.csv"
output_file = OUTPUT_DIR / OUTPUT_FILE_CSV


# Print algaetax metadata
def startup_info():
    print(f"\n{NAME} {VERSION}\n")

    print("algaetax Metadata:")
    if is_database_enabled("NCBI"):
        email = config["db_ncbi"].get("ncbi_email", "Not set")
        api_key = config["db_ncbi"].get("api_key")
        print(f"  NCBI_Email   : {email}")
        print(f"  NCBI_API_Key : {'API Key set' if api_key else 'No API Key'}")

    if is_database_enabled("PR2"):
        path = config["database_path"].get("db_pr2", "Not set")
        print(f"  PR2_DB_Path  : {path}")

    if is_database_enabled("ALGB"):
        api_key = config["db_algaebase"].get("api_key")
        api_url = config["db_algaebase"].get("api_url")
        print(f"  ALGB_API_Key : {'API Key set' if api_key else '[WARN] No API Key'}")
        print(f"  ALGB_API_URL : {api_url or '[WARN] No API URL'}")

    print("\nActive databases:")
    active = [db for db, enabled in config["database"].items() if enabled]
    if active:
        for db in active:
            name = SUPPORTED_DATABASES.get(db, "Unknown")
            print(f"  {db} ({name})")
    else:
        print("  None enabled. Check config.yaml.")
        sys.exit(1)
    print(f"{'-' * 20}\n")

# Query available databases and save results
def query_and_save(taxa_list, output_file):
    pr2_data = load_pr2_db() if is_database_enabled("PR2") else None
    taxon_ids = get_taxon_ids_parallel(taxa_list) if is_database_enabled("NCBI") else {}
    ncbi_results = batch_query_ncbi(taxon_ids) if is_database_enabled("NCBI") else {}

    # Print headers for PR2 and ALGB queries
    if is_database_enabled("PR2"):
        print(f"\n{'-' * 20}")
        print("[PR2] Start LOCAL query...")
    if is_database_enabled("ALGB"):
        print(f"\n{'-' * 20}")
        print("[ALGB] Start API query...\n")

    results = []
    pr2_counter = 0
    algb_counter = 0

    for taxon in taxa_list:
        taxon_result = {}
        
        if is_database_enabled("NCBI"):
            taxon_result["NCBI"] = ncbi_results.get(taxon, "Not found")
        
        if is_database_enabled("PR2"):
            pr2_counter += 1
            # print(f"[{pr2_counter}/{len(taxa_list)}] PR2: {taxon}")
            pr2_result = query_pr2(taxon, pr2_data)
            taxon_result["PR2"] = pr2_result
        
        if is_database_enabled("ALGB"):
            algb_counter += 1
            print(f"[{algb_counter}/{len(taxa_list)}] ALGB: {taxon}")
            algae_result = query_algaebase(taxon)
            taxon_result["ALGB"] = algae_result

        results.append((taxon, taxon_result))

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    write_results_to_csv(results, output_file)

    elapsed = time.time() - start_time
    minutes = int(elapsed // 60)
    seconds = int(elapsed % 60)
    
    print(f"\n{'-' * 20}")
    print(f"[at] All queries completed.")
    print(f"[at] Results saved to: {output_file}")
    print(f"[at] Total runtime: {minutes} min {seconds} sec ({elapsed:.2f} seconds)\n")


# Format and write query results to CSV file
def write_results_to_csv(results_list, output_file):
    rows = []

    ncbi_enabled = is_database_enabled("NCBI")
    pr2_enabled = is_database_enabled("PR2")
    algae_enabled = is_database_enabled("ALGB")


    for taxon, result in results_list:
        row = {"Taxon": taxon}

        # Process NCBI result
        if ncbi_enabled:
            ncbi = result.get("NCBI")
            if isinstance(ncbi, dict):
                row["NCBI_Status"] = "Found"
                row["NCBI_ID"] = str(ncbi.get("NCBI_ID", ""))
                row["NCBI_Lineage"] = ncbi.get("Taxonomy", "")
            else:
                row["NCBI_Status"] = "Not found"

        # Process PR2 result
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
            else:
                row["PR2_Status"] = "Not found"

        # Process AlgaeBase result
        if algae_enabled:
            algae = result.get("ALGB")
            if isinstance(algae, dict):
                row["ALGB_Status"] = algae.get("ALGB_Status", "")
                row["ALGB_ID"] = algae.get("ALGB_ID", "")
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
        columns += ["PR2_Status", "PR2_Lineage"]
    if algae_enabled:
        columns += ["ALGB_Status", "ALGB_ID", "ALGB_Lineage"]
    df = df[columns]

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    df.to_csv(output_file, index=False, encoding="utf-8", sep=";", quoting=csv.QUOTE_MINIMAL)
    print(f"\n[at] CSV file saved: {output_file}")

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
                print(f"\n[at] Exported list: Not found taxa (active DBs) - {not_found_path}")
            else:
                empty_df = pd.DataFrame(columns=["Taxon"] + enabled_cols)
                empty_df.to_csv(not_found_path, index=False, encoding="utf-8", sep=";")
                print(f"\n[at] No taxa with all active databases 'Not found'. Empty file created: {not_found_path}")

# Run main process
if __name__ == "__main__":
    startup_info()
    print("RUNNING algaetax ...\n")

    # Load taxa list from Excel file
    input_file = config["general"]["input_data"]
    taxa_col_letter = config["general"]["taxa_column_letter"]
    df, taxa_list = load_taxa_from_excel(input_file, taxa_col_letter, output_file)


    # Clean and filter taxa
    valid_taxa = []
    invalid_taxa = []
    filter_log = []

    for taxon in taxa_list:
        cleaned = extract_valid_taxon(taxon)
        if cleaned:
            valid_taxa.append(cleaned)
            if cleaned != taxon.strip():
                # Track cleaned values
                filter_log.append((taxon, cleaned))
        else:
            invalid_taxa.append(taxon)


    # Check for duplicates
    taxa_index = col_letter_to_index(taxa_col_letter)
    raw_column = df.iloc[:, taxa_index].dropna().astype(str).str.strip()
    counter = Counter(raw_column)
    duplicates = {k: v for k, v in counter.items() if v > 1}


    if duplicates:
        print("\n[at] Duplicate taxa found:")
        for taxon, count in sorted(duplicates.items(), key=lambda x: x[1], reverse=True):
            print(f"- {taxon}: {count}×")
    else:
        print("\n[at] No duplicate taxa found.")

    print(f"\n[at] TAXA loaded: {len(taxa_list)}")
    print(f"[at] Input preview: {taxa_list[:2]} ...")
    print(f"\n[at] Filtered input: {valid_taxa[:2]} ...")


    # Save cleaned taxa log if needed
    if filter_log:
        filepath = os.path.join(OUTPUT_DIR, "filtered_taxa_log.csv")
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        with open(OUTPUT_DIR / "filtered_taxa_log.csv", "w", encoding="utf-8") as f:
            f.write("Taxa,Filtered\n")
            for orig, filt in filter_log:
                f.write(f"{orig},{filt}\n")
        print(f"[at] Saved filtered taxa to: {filepath}")


    # Run taxonomic queries
    query_and_save(valid_taxa, output_file)

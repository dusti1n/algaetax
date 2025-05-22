# algaetax; utils.py
# Import libraries
import yaml
import pandas as pd
import shutil
import os
import string
import re

# Load YAML config
def load_config():
    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)

config = load_config()


# Check if a database is enabled in config
def is_database_enabled(db_name):
    return config["database"].get(db_name, False)


# Convert Excel column letter (e.g. 'A') to index (e.g. 0)
def col_letter_to_index(letter):
    return string.ascii_uppercase.index(letter.upper())


# Load taxa column from Excel and make backup
def load_taxa_from_excel(filepath, taxa_col_letter, output_path):
    header_setting = config["general"].get("header_row", None)
    backup_file = config["general"].get("backup_file", True)

    # Interpret header config
    if isinstance(header_setting, int):
        header = header_setting
    elif header_setting is None or str(header_setting).lower() == "none":
        header = None
    else:
        raise ValueError("Invalid 'header_row' value in config.yaml")

    # Create backup of original Excel file
    if backup_file:
        output_dir = os.path.dirname(output_path)
        os.makedirs(output_dir, exist_ok=True)
        backup_path = os.path.join(output_dir, os.path.basename(filepath).replace(".xlsx", "_backup.xlsx"))
        if not os.path.exists(backup_path):
            shutil.copy(filepath, backup_path)
            print(f"[at] Backup file created: {backup_path}")
        else:
            print("[at] Backup file already exists.")

    try:
        df = pd.read_excel(filepath, header=header)
        col_index = col_letter_to_index(taxa_col_letter)
        taxa_series = df.iloc[:, col_index].dropna().astype(str).str.strip()

        # Remove title row if no header and first cell looks like "Taxa"
        if header is None and isinstance(taxa_series.iloc[0], str) and "taxa" in taxa_series.iloc[0].lower():
            print("Note: Detected header 'Taxa'. It will be removed.")
            taxa_series = taxa_series.iloc[1:]

        return df, taxa_series.tolist()
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        raise


# Extract valid taxon name (e.g. remove 'sp.', 'cf.', etc.)
def extract_valid_taxon(name):
    # Check if input is a string; return None if not
    if not isinstance(name, str):
        return None

    # Remove leading/trailing whitespace
    name = name.strip()

    # If multiple names are listed with a slash, take the first
    if "/" in name:
        name = name.split("/")[0].strip()

    # Lowercase blacklist from config (e.g. sp, cf, aff, gr)
    blacklist = set(word.lower().strip(".") for word in config.get("filter", {}).get("blacklist", []))

    # Remove measurement annotations like "35-40u", "10µm", "20 um", etc.
    name = re.sub(r"\b\d+[\d\-\sxumµ]*", "", name, flags=re.IGNORECASE)

    # Remove unwanted characters (everything except letters and space)
    name = re.sub(r"[^A-Za-z\s]", " ", name)

    # Normalize multiple spaces
    name = re.sub(r"\s+", " ", name).strip()

    # Extract only words that consist of letters (no digits/symbols)
    words = re.findall(r"\b[A-Za-z]+\b", name)

    if not words:
        return None

    # First word is genus
    genus = words[0]

    # Check if there is a valid species and not in blacklist
    if len(words) > 1:
        species = words[1]
        if species.lower() in blacklist:
            return genus
        return f"{genus} {species}"

    return genus

# algaetax; file: utils.py

# Import libraries
import yaml, os, re, shutil
import sys
import pandas as pd
import argparse


# Cache for loaded configuration and its file path
_CONFIG = None
_CONFIG_PATH = None

# Return path of the loaded configuration file
def get_config_path():
    return _CONFIG_PATH

# Load YAML config file (only once, then cached)
def load_config():
    global _CONFIG, _CONFIG_PATH
    if _CONFIG is not None:
        return _CONFIG

    # Parse optional --configfile argument
    parser = argparse.ArgumentParser(description="Run algaetax workflow.")
    parser.add_argument("--configfile", type=str, default="config.yaml",
                        help="Path to configuration file (default: config.yaml)")
    args, _ = parser.parse_known_args()

    config_path = args.configfile
    if not os.path.exists(config_path):
        print(f"[ERROR] Config file not found: {config_path}")
        sys.exit(1)

    # Read and load YAML config
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    print(f"[algaetax] Using configuration: {config_path}\n")

    # Cache config and path
    _CONFIG = config
    _CONFIG_PATH = config_path
    return _CONFIG


# Check if a database is enabled in config
def is_database_enabled(db_name, config):
    return config["database"].get(db_name, False)


# Read taxa from Excel file and optionally create a backup
def load_taxa_from_excel(filepath, taxa_col_number, output_path, config):
    header_setting = config["general"].get("header_row", 1)
    backup_file = config["general"].get("backup_input", True)
    id_col_number = config["general"].get("id_column_number", False)

    # Determine Excel header row
    if header_setting is False:
        header = None  # No header
    elif isinstance(header_setting, int):
        if header_setting < 1:
            raise ValueError("header_row must be >= 1 or false")
        header = header_setting - 1  # Convert to 0-based
    else:
        raise ValueError("Invalid value for 'header_row' in config.yaml")

    # Create backup of input Excel file (if enabled)
    if backup_file:
        output_dir = os.path.dirname(output_path)
        backup_dir = os.path.join(output_dir, "backups")
        os.makedirs(backup_dir, exist_ok=True)

        original_filename = os.path.basename(filepath)
        backup_filename = original_filename.replace(".xlsx", "_backup.xlsx")
        backup_path = os.path.join(backup_dir, backup_filename)

        if not os.path.exists(backup_path):
            shutil.copy(filepath, backup_path)
            print(f"[algaetax] Backup file created: \n{backup_path}")
        else:
            print("[WARN] Backup file already exists.")
            
    # Read Excel and extract taxa column
    try:
        df = pd.read_excel(filepath, header=header)
        col_index = taxa_col_number - 1
        taxa_series = df.iloc[:, col_index].dropna().astype(str).str.strip()

        # Remove header-like first row if no official header
        if header is None and isinstance(taxa_series.iloc[0], str) and "taxa" in taxa_series.iloc[0].lower():
            print("Note: Detected header 'Taxa'. It will be removed.")
            taxa_series = taxa_series.iloc[1:]

        # id_col_number; Optional extraction of ID column
        id_series = None
        id_column_name = None

        if id_col_number and isinstance(id_col_number, int):
            try:
                id_index = id_col_number - 1
                id_series = df.iloc[:, id_index].astype(str).str.strip()

                # Detect actual column name
                if header_setting:
                    id_column_name = df.columns[id_index]
                else:
                    id_column_name = "ID"  # Fallback if no header

                # Debug: show what column was read
                col_name = df.columns[id_index] if header_setting else f"Column {id_col_number}"
                print(f"\n[algaetax] ID column detected: '{col_name}' (index {id_col_number})")
                # print(f"[algaetax] Example values: {id_series.head(3).tolist()}")

                # If the column is completely empty, warn user
                if id_series.dropna().empty:
                    print(f"[WARN] ID column (col {id_col_number}) is empty or contains only NaN values.")
                    id_series = None

            except Exception as e:
                print(f"[ERROR] Could not read ID column (col {id_col_number}): {e}")
                id_series = None
        else:
            print("[INFO] No ID column defined")

        # Return with debug output
        if id_series is not None:
            print(f"[algaetax] Loaded ID column ({len(id_series)} entries)")
            return df, taxa_series.tolist(), id_series.tolist(), id_column_name  # Include name
        else:
            return df, taxa_series.tolist(), None, None  # Return name placeholder

    except Exception as e:
        print(f"Error reading Excel file: {e}")
        raise


# Clean and extract valid taxon name (removes noise like 'sp.', digits)
def extract_valid_taxon(name, blacklist):
    # Skip invalid entries (e.g., None, numbers)
    if not isinstance(name, str):
        return None

    # Trim whitespace at beginning and end
    s = name.strip()

    # Remove anything after a slash (e.g., "Genus species / comment")
    if "/" in s:
        s = s.split("/")[0].strip()

    # Normalize unicode dashes (– — − etc.) to a standard hyphen
    s = re.sub(r"[\u2010\u2011\u2012\u2013\u2014\u2212]", "-", s)

    # Remove measurements like "10 µm", "5-10um", etc.
    s = re.sub(r"\b\d+[\d\-\sxumµ]*", "", s, flags=re.IGNORECASE)

    # Keep only letters, spaces, and hyphens; remove all other characters
    s = re.sub(r"[^A-Za-z\s\-]", " ", s)

    # Collapse multiple spaces into one
    s = re.sub(r"\s+", " ", s).strip()

    # Extract tokens, allowing internal hyphens (e.g., "flos-aquae")
    tokens = re.findall(r"\b[A-Za-z]+(?:-[A-Za-z]+)*\b", s)
    if not tokens:
        return None

    # Capitalize the genus and lowercase the species
    genus = tokens[0].lower().capitalize()

    # If a species name is present, check and clean it
    if len(tokens) > 1:
        species = tokens[1].lower()

        # Skip uncertain species names (e.g., "sp.", "cf.") but keep genus
        if species in blacklist:
            return genus

        # Return clean "Genus species"
        return f"{genus} {species}"

    # If only one valid word found, return genus only
    return genus


# Load blacklist from file and (optionally) back it up
def load_blacklist(config):
    blacklist_file = config["filter"].get("blacklist_file")
    if not blacklist_file:
        raise ValueError("Missing 'blacklist_file' entry in config.yaml under 'filter'")

    if not os.path.exists(blacklist_file):
        raise FileNotFoundError(f"Blacklist file not found: {blacklist_file}")

    # Backup if configured
    if config["filter"].get("backup_blacklist", False):
        output_dir = config["general"].get("output_dir", "results")
        backup_dir = os.path.join(output_dir, "backups")
        os.makedirs(backup_dir, exist_ok=True)

        backup_path = os.path.join(backup_dir, os.path.basename(blacklist_file))
        try:
            shutil.copy(blacklist_file, backup_path)
        except Exception as e:
            print(f"[WARN] Could not back up blacklist: {e}")

    # Return cleaned list of blacklist terms
    with open(blacklist_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # Return cleaned list of terms, ignoring empty lines and comments
    return [line.strip().lower() for line in lines if line.strip() and not line.strip().startswith("#")]


def safe_save_path(output_path):
    """
    Checks if the output file already exists at the very beginning of the workflow.
    If it exists, the user is asked whether they want to overwrite it.
    Accepts 'y'/'yes' to overwrite, 'n'/'no' to cancel the workflow.
    """
    if os.path.exists(output_path):
        print(f"\n[WARNING] The file '{output_path}' already exists.")
        choice = input("Do you want to overwrite it? (y/yes or n/no): ").strip().lower()

        if choice in ["y", "yes"]:
            # Clear terminal for a clean workflow start
            os.system('cls' if os.name == 'nt' else 'clear')
            print("\n[INFO] Overwriting existing file...")
            return output_path
        elif choice in ["n", "no"]:
            print("\n[INFO] AlgaeTax was stopped because the output file already exists.")
            sys.exit(0)
        else:
            print("\n[ERROR] Invalid input. Please enter 'y'/'yes' or 'n'/'no'.")
            sys.exit(1)
    return output_path

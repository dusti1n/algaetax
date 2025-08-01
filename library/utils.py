# algaetax; file: utils.py

# Import libraries
import yaml, os, re, shutil
import sys
import pandas as pd



# Load YAML config file
def load_config():
    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)


# Check if a database is enabled in config
def is_database_enabled(db_name, config):
    return config["database"].get(db_name, False)


# Read taxa from Excel file and optionally create a backup
def load_taxa_from_excel(filepath, taxa_col_number, output_path, config):
    header_setting = config["general"].get("header_row", 1)
    backup_file = config["general"].get("backup_input", True)

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
        backup_dir = os.path.join(output_dir, "backups")  # 🔧 Unterordner 'backups'
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

        return df, taxa_series.tolist()
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        raise


# Clean and extract valid taxon name (removes noise like 'sp.', digits)
def extract_valid_taxon(name, blacklist):
    if not isinstance(name, str):
        return None
    name = name.strip()
    if "/" in name:
        name = name.split("/")[0].strip()

    # Remove measurements like "10 µm"
    name = re.sub(r"\b\d+[\d\-\sxumµ]*", "", name, flags=re.IGNORECASE)

    # Remove non-letter characters
    name = re.sub(r"[^A-Za-z\s]", " ", name)
    name = re.sub(r"\s+", " ", name).strip()

    words = re.findall(r"\b[A-Za-z]+\b", name)
    if not words:
        return None

    genus = words[0]
    if len(words) > 1:
        species = words[1]
        if species.lower() in blacklist:
            return genus
        return f"{genus} {species}"

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
            print("\n[INFO] Overwriting existing file...")
            return output_path
        elif choice in ["n", "no"]:
            print("\n[INFO] AlgaeTax was stopped because the output file already exists.")
            sys.exit(0)
        else:
            print("\n[ERROR] Invalid input. Please enter 'y'/'yes' or 'n'/'no'.")
            sys.exit(1)

    return output_path

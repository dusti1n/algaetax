# algaetax; file: metadata.py

# Basic tool information
NAME = "algaetax"
VERSION = "v1.5.0"
INSTITUTION = "University of Duisburg-Essen (UDE)"
AUTHOR = "Dustin Finke, Prof. Dr. Bank Beszteri"
LICENSE = "MIT"

# Description
DESCRIPTION = (
    "algaetax is a bioinformatics tool for querying and extracting taxonomic data from various databases. "
    "It automates taxonomic lookups and structures the retrieved data into a standardized format for further analysis."
)

# Supported file formats
SUPPORTED_OUTPUT_FILE = [".csv"]
SUPPORTED_INPUT_FORMATS = [".xlsx"]

# Supported taxonomic databases
SUPPORTED_DATABASES = {
    "NCBI": "National Center for Biotechnology Information",
    "PR2": "Protist Ribosomal Reference",
    "ALGB": "AlgaeBase; Global Database of Algae"
}

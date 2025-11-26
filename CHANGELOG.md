#### algaetax; bioinformatics

# Changelog

This document records all notable changes to the **algaetax** project.  
**Note**: Use the **dev** branch for the latest updates and features.

---

### [2025-11-26]
- Update `README.md`
- Add acknowledgements section to `README.md`

---

### [2025-11-20]
- **Update algaetax to version v1.7.0**  
- Update `README.md` with clearer API key information  
- Add troubleshooting notes for missing or invalid AlgaeBase API keys  
- Add Skiplist support to skip defined taxonomy terms  
- Keep Skiplist taxa in results but mark them as skipped  
- Add “Skipped” status to the output  
- Save skipped taxa to `skipped_taxa.txt`  
- Fix ALGB counter to count only queried taxa  
- Fix NCBI counter to count only queried taxa  
- Add optional Skiplist debug output  
- Clean up blacklist to avoid conflicts with the Skiplist

---

### [2025-10-28]
- Update `.gitignore` to include additional project files and folders  
- Update `README.md` with clear setup, configuration, and usage guidance  
- Extend PR2 output by adding a `species` column with optional values  
- Improve taxon name parsing to correctly handle hyphenated entries  
- Add `id_column_number` option in `config.yaml` for including an ID column  
- Enable inclusion of optional input columns like sequence or sample IDs  
- Clean and refactor code for better structure and readability  
- Update CLI output for improved terminal readability  
- Remove `metadata.py`; merged functionality into main workflow  
- Add `config_presets/` folder for custom configuration files  
- Add `--configfile <path>/<config.yaml>` CLI parameter for specifying custom configs  
- Update `load_config()` to cache configuration for better efficiency and consistency  

---

### [2025-08-01]
- Update `CHANGELOG.md`  
- Refactor codebase  
  - Clean up code and add meaningful comments in various files  
  - Restructure entire root directory  
  - Improve efficiency and readability  
  - Improve modular code structure for Python scripts  
- Create `library/` folder for important Python files to improve root directory organization  
- Update `algaetax.yaml` environment configuration  
- Extract blacklist into standalone `blacklist.txt`  
  - Improve word list visibility  
  - Provide user-specific customization  
- Update `config.yaml` with improved structure and values  
- Adjust `taxa_column_number` parameter to accept numeric values for the target taxa column  
- Add backup parameters: `backup_config`, `backup_blacklist`  
- Revise `header_row` parameter to allow `false` and make row indexing start at 1 instead of 0  
- Improve lineage retrieval from AlgaeBase for cleaner output  
- Improve analysis progress display  
- Add safety check to prompt when output file already exists  
- Add AlgaeBase API key validation with safety confirmation  

---

### [2025-05-22]
- Add blacklist-based taxa filter option in `config.yaml`  
- Add `backup_file` option to `config.yaml`  
- Add `export_not_found_taxa_list` option to `config.yaml`  
- Add logic to export taxa not found in all active databases  
- Update `algaetax.py` and `db_queries.py` for refactoring and new logic  

---

### [2025-05-18]
- Update PR2 database to version `5.1.0`  
- Fix bug in taxon validation  
- Fix bug by using the PR2 database  

---

### [2025-05-08]
- Add dependencies to `algaetax.yaml` environment file  
- Add `taxa_column_letter` option to config file  
- Add `header_row` option to config file  
- Add new file: `metadata.py`  
- Update `.gitignore` file  
- Include AlgaeBase database API request  
- Include dynamic input loading for `.xlsx` files  
- Add logfile for filtered taxon names  
- Revise files: `algaetax.py`, `db_queries.py`, `utils.py`  
- Add support for CSV output file  
- Enable automatic backup for input file  

---

### [2025-04-22]
- Update files: `config.yaml`, `.gitignore`  
- Integrate taxonomy PR2 database (v5.0.0)  
- Modify data structure in the root directory  
- Add NCBI API key option  
- Process NCBI data in batch mode  

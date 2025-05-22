#### algaetax; Bioinformatics
# Changelog
This document records all notable changes to algaetax project.

## [2025-05-22]
- Add blacklist-based taxa filter option in config.yaml
- Add backup_file option to config.yaml
- Add export_not_found_taxa_list option to config.yaml
- Add logic to export taxa not found in all active databases
- Update algaetax.py and db_queries.py for refactoring and new logic

## [2025-05-18]
- Update PR2 database to version 5.1.0
- Fix bug in taxon validation
- Fix bug by using the PR2 database

## [2025-05-08]
- Add dependencies to algaetax.yaml ENV  
- Add taxa_column_letter option; config file  
- Add header_row option; config file
- Add file: metadata.py
- Update gitignore file  
- Include AlgaeBase database; API request  
- Include dynamic input load with (.xlsx) file
- Include logfile for filtered taxon names   
- Revise files: algaetax.py, db_queries.py, utils.py
- Support CSV output file   
- Automatic backup for input file  

## [2025-04-22]
- Update files: config.yaml, .gitignore  
- Integrate taxonomy PR2 database (v5.0.0)  
- Modify data structure in the root directory  
- Add NCBI API key option  
- Process NCBI data with batch 
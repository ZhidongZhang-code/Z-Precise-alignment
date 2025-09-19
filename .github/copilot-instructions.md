# GitHub Copilot Instructions for Z-Precise-alignment

This document provides guidance for AI coding agents to quickly understand the architecture, workflows, and project-specific conventions of this repository.

## Project Overview

- **Purpose:** Implements sequence alignment and precise species extraction pipelines.
- **Major Components:**
  - Species checking and report generation (e.g., `src/species_checker_format.py`).
  - Sequence processing modules like BLAST and HMM handling (e.g., `src/blast_module.py`, `src/hmm_module.py`).
  - Protein and CDS extraction (e.g., `src/excat_site_cds.py`, `src/excat_site_protein.py`).

## Architecture and Data Flow

- **Input:** FASTA files with sequence annotations, and configuration via YAML (e.g., `cfg/config.yaml`).
- **Processing:**
  - Parsing FASTA headers using `_**_` as a delimiter to extract GCF identifiers and protein IDs.
  - Annotating sequences based on a CSV file of GCF annotations.
  - Checking for protein presence using simple string matching within input files.
- **Output:** CSV reports (tab-delimited) with detailed annotations and homology check results.

## Developer Workflows

- **Running Scripts:**
  - Execute main scripts directly from the command line. For example:
    ```cmd
    python src/species_checker_format.py --extractprotein input.fasta --site_protein_output site.out --config cfg/config.yaml --log_file log.txt --formatoutput report.tsv
    ```
- **Logging:**
  - Use `setup_logging` from `src/log_config.py` for detailed runtime logs; inspect log files to diagnose issues.
- **Configuration:**
  - YAML configuration files (e.g., in `cfg/`) define paths and parameters for annotations and file inputs.

## Project Conventions & Patterns

- **Data Parsing:**
  - FASTA headers are split using the `_**_` delimiter. The first part (after removing `>` and splitting by `_`) is used as the GCF identifier; the second part is taken as the protein ID.
- **Report Generation:**
  - DataFrames are constructed with a specific column order (e.g., `GCF_ID`, `strain_Name`, `taxon`, `Species`, `Genus`, `Family`, `Order`, `Class`, `Phylum`, `homology`, `homology_with_site`).
  - CSV outputs are written with tab delimiters to ensure proper formatting.
- **Error Handling:**
  - Exceptions in data extraction and file operations are logged, not necessarily raised, to allow for iterative processing.

## Integration and Dependencies

- **External Libraries:**
  - Utilizes `pandas` for data handling and `yaml` for configuration parsing.
- **File Interactions:**
  - Multiple scripts read from common data sources (e.g., annotation CSVs) and write reports in a structured format.

## Testing & Debugging

- **Local Execution:**
  - Modules can be executed individually to test specific functionalities; command line arguments ensure reusability and modularity.
- **Log Analysis:**
  - Review log messages from `setup_logging` to troubleshoot failures during report generation and data extraction.

## Additional Resources

- **Code Examples:** Refer to `src/species_checker_format.py` and related modules for specific implementation patterns.
- **Feedback:** This is a living document. Please provide feedback on any unclear or incomplete sections to iteratively improve these instructions.

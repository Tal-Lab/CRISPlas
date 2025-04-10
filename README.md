
#### **Overview**
The `CRISPlas` pipeline is a Python-based tool designed to analyze plasmid sequences against a spacers database using BLAST and perform downstream taxonomic lineage extraction, unique organism identification, and host range analysis. The pipeline consists of four main steps:
1. **BLAST Search**: Matches plasmid sequences against a spacers database.
2. **Lineage Extraction**: Extracts taxonomic lineage information for matched sequences.
3. **Unique Organism Identification**: Processes lineage data to identify unique plasmid-organism combinations.
4. **Host Range Analysis**: Determines the first identical taxonomic level for each plasmid.

---

#### **Requirements**
- **Python 3.7+**
- **Dependencies**:
  - `pandas`
  - `BioPython`
- **BLAST+ Tools**:
  - Ensure `blastn` is installed and accessible in your system's PATH.
  - Create a BLAST database for the spacers using `makeblastdb`.

---

#### **Installation**
1. Install Python dependencies:
   ```bash
   pip install pandas biopython
   ```
2. Install BLAST+ tools:
   - Download and install from [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
   - Add the BLAST+ binaries to your system's PATH.

3. Create a BLAST database for the spacers:
   ```bash
   makeblastdb -in spacers_db.fasta -dbtype nucl -out Spacers_db
   ```

---

#### **Usage**
1. **Prepare Input Files**:
   - **Plasmid FASTA File**: A file containing plasmid sequences (e.g., `plasmids.fasta`).
   - **Spacers Database**: A BLAST-compatible database created from spacer sequences (e.g., `Spacers_db`).

2. **Run the Pipeline**:
   Update the file paths in the `main` function:
   ```python
   plasmid_fasta = r"\path\to\your\working\directory\plasmids.fasta"
   spacers_db = r"\path\to\your\working\directory\Spacers_db"
   main(plasmid_fasta, spacers_db)
   ```

   Execute the script:
   ```bash
   python CRISPlas.py
   ```

---

#### **Pipeline Steps**
1. **BLAST Search**:
   - Matches plasmid sequences against the spacers database using the following parameters:
     - `-evalue 1e-5`: E-value threshold for significant matches.
     - `-max_target_seqs 1`: Limits results to the top hit per query.
     - Optional: `-perc_identity 100`: Filters matches to only those with 100% identity.
   - Output: `blast_results.txt`.

2. **Lineage Extraction**:
   - Reads BLAST results and fetches taxonomic lineage for matched sequences using the NCBI Entrez API.
   - Output: `source_target_lineage.csv`.

3. **Unique Organism Identification**:
   - Processes lineage data to identify unique plasmid-organism combinations.
   - Output: `uniques.csv`.

4. **Host Range Analysis**:
   - Analyzes the processed lineage data to determine the first identical taxonomic level for each plasmid.
   - Output: `Host_Range.csv`.

---

#### **Output Files**
1. **`blast_results.txt`**:
   - Tab-delimited BLAST results containing query and subject sequence IDs.

2. **`source_target_lineage.csv`**:
   - CSV file containing plasmid-target pairs and their taxonomic lineage.

3. **`uniques.csv`**:
   - CSV file with unique plasmid-organism combinations.

4. **`Host_Range.csv`**:
   - CSV file showing the first identical taxonomic level for each plasmid.

---

#### **Code Structure**
1. **`run_blast`**:
   - Runs the BLAST search and saves results to `blast_results.txt`.

2. **`lineage_extraction`**:
   - Extracts taxonomic lineage for matched sequences using the NCBI Entrez API.

3. **`unique_csv`**:
   - Processes lineage data to identify unique plasmid-organism combinations.

4. **`host_range_csv`**:
   - Analyzes the processed lineage data to determine host range.

5. **`main`**:
   - Orchestrates the pipeline by calling the above functions in sequence.

---

#### **Example Output**
- **`Host_Range.csv`**:
   | Plasmid       | First_Identical_Column |
   |---------------|------------------------|
   | plasmid1      | Genus                 |
   | plasmid2      | Family                |

---

#### **Troubleshooting**
1. **BLAST Command Not Found**:
   - Ensure `blastn` is installed and added to your system's PATH.

2. **Empty BLAST Results**:
   - Verify that the plasmid sequences and spacers database are compatible.
   - Check the e-value threshold and percentage identity settings.

3. **NCBI API Errors**:
   - Ensure `Entrez.email` is set to a valid email address.
   - Add delays between API requests if rate-limiting occurs:
     ```python
     time.sleep(0.5)
     ```

4. **FileNotFoundError**:
   - Ensure all input file paths are correct and accessible.

---

#### **License**
This work is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.
You are free to:

Share: Copy and redistribute the material in any medium or format.
Adapt: Remix, transform, and build upon the material for any purpose, even commercially.
Attribution:
You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.

For more details, see the Creative Commons License (https://creativecommons.org/licenses/by/4.0/).
---

#### **Contact**
For questions or issues, please contact [shay.tal@ocean.org.il].

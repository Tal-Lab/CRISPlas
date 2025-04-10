import os
import subprocess
import pandas as pd
from Bio import Entrez, SeqIO

def run_blast(plasmid_fasta, spacers_db):
    # Step 1: Run BLAST
    print("Running BLAST...")
    blast_output_path = os.path.join(os.path.dirname(plasmid_fasta), "blast_results.txt")
    
    # Construct the BLAST command
    blast_command = [
        "blastn",
        "-query", plasmid_fasta,
        "-db", spacers_db,
        "-out", blast_output_path,
        "-outfmt", "6",  # Output format - you can change this as needed. As long as the first 2 columns are qseqid and sseqid, not further modification is needed.
        "-evalue", "1e-5",
        "-max_target_seqs", "1"  # Limit to 1 target sequence per query, which is enough for the purpose of host prediction. This can be adjusted based on your needs.
        #"-perc_identity", "100"  # Use this to only include matches with 100% identity.

    ]
    
    # Run the BLAST command
    try:
        subprocess.run(blast_command, check=True)
        print(f"BLAST completed. Results saved to {blast_output_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error running BLAST: {e}")
        raise
    
    return blast_output_path

def lineage_extraction(blast_results_path):
    # Step 2: Run Lineage Extraction
    print("Running Lineage Extraction...")
    lineage_output_path = os.path.join(os.path.dirname(blast_results_path), "source_target_lineage.csv")
    Entrez.email = "your_email@example.com"  # Replace with your email

    # Read the BLAST results (only the first two columns: qseqid and sseqid)
    blast_results = pd.read_csv(blast_results_path, sep="\t", header=None, usecols=[0, 1], names=["source_accession", "target_accession"])
    results = []

    for index, row in blast_results.iterrows():
        target_accession = row["target_accession"].split('_')[0]
        try:
            handle = Entrez.efetch(db="nucleotide", id=target_accession, rettype="gb", retmode="text")
            seq_record = SeqIO.read(handle, "genbank")
            handle.close()
            organism = seq_record.annotations.get("organism", "")
            taxonomy = seq_record.annotations.get("taxonomy", [])
            result = [row["source_accession"], target_accession, organism] + taxonomy
            results.append(result)
        except Exception as e:
            print(f"Error fetching lineage for {target_accession}: {e}")
            results.append([row["source_accession"], target_accession, "Error"])

    expanded_df = pd.DataFrame(results).fillna('')
    expanded_df.to_csv(lineage_output_path, index=False)
    print(f"Lineage extraction completed. Results saved to {lineage_output_path}")
    return lineage_output_path

def unique_csv(lineage_csv):
    # Step 3: Process Unique CSV
    print("Processing Unique CSV...")
    df = pd.read_csv(lineage_csv, usecols=range(10))  # Read only the first 10 columns
    df.columns = ["Plasmid", "HOST_ID", "Organism", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
    df["Organism"] = df["Organism"].apply(lambda x: " ".join(str(x).split()[:2]))
    df_combined = df.groupby(["Plasmid", "Organism"], as_index=False).first()
    
    # Reorder columns to match the desired format
    column_order = ["Plasmid", "Organism", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Domain"]
    df_combined = df_combined[column_order]

    unique_csv_path = os.path.join(os.path.dirname(lineage_csv), "uniques.csv")
    df_combined.to_csv(unique_csv_path, index=False)
    print(f"Unique CSV processing completed. Results saved to {unique_csv_path}")
    return unique_csv_path

def host_range_csv(unique_csv_path):
    # Step 4: Process Host Range CSV
    print("Processing Host Range CSV...")
    df = pd.read_csv(unique_csv_path)
    df.columns = ["Plasmid", "Organism", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Domain"]
    columns_to_check = ["Organism", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Domain"]
    plasmid_first_identical = {}

    for plasmid, group in df.groupby("Plasmid"):
        for column in columns_to_check:
            if group[column].nunique() == 1:
                plasmid_first_identical[plasmid] = column
                break

    result_df = pd.DataFrame(list(plasmid_first_identical.items()), columns=["Plasmid", "First_Identical_Column"])
    host_range_csv_path = os.path.join(os.path.dirname(unique_csv_path), "Host_Range.csv")
    result_df.to_csv(host_range_csv_path, index=False)
    print(f"Host Range CSV processing completed. Results saved to {host_range_csv_path}")
    return host_range_csv_path

def main(plasmid_fasta, spacers_db):
    blast_results_path = run_blast(plasmid_fasta, spacers_db)
    lineage_csv_path = lineage_extraction(blast_results_path)
    unique_csv_path = unique_csv(lineage_csv_path)
    host_range_csv_path = host_range_csv(unique_csv_path)
    print("All steps completed successfully!")

# Example usage
plasmid_fasta = r"\path\to\your\working\directory\plasmids.fasta"  # Replace with your plasmid FASTA file path
spacers_db = r"\path\to\your\working\directory\Spacers_db"  # Replace with your spacers database path
main(plasmid_fasta, spacers_db)
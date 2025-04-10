import pandas as pd
from Bio import Entrez
from Bio import SeqIO

# Function to fetch taxonomic lineage for an accession number
def fetch_lineage(accession):
    Entrez.email = "YourMail@server.org "  # Replace with your email
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    seq_record = SeqIO.read(handle, "genbank")
    handle.close()
    organism = seq_record.annotations.get("organism", "")  # Extract organism information
    taxonomy = seq_record.annotations.get("taxonomy", [])  # Extract taxonomy information
    
    return organism, taxonomy

# Read the BLASTn results (first two columns only)
blast_results = pd.read_csv(r"\path\to\your\working\directory\blast_results.txt", sep="\t", usecols=[0, 1], header=None, names=["source_accession", "target_accession"])

# Fetch taxonomic lineage for each target accession number
results = []  # Initialize an empty list to store the results
for index, row in blast_results.iterrows():
    target_accession = row["target_accession"].split('_')[0]  # Trim after the first underscore
    try:
        organism, taxonomy = fetch_lineage(target_accession)
        result = [row[“source_accession”], target_accession, organism] + taxonomy
        print(result)  # Print the value for each round
        results.append(result)
        expanded_df = pd.DataFrame(results).fillna('')  # Fill missing values with empty strings

        # Save the final DataFrame to a CSV file
        expanded_df.to_csv(r"\path\to\your\working\directory\source_target_lineage.csv", index=False)
    except Exception as e:
        print(f"Error fetching lineage for {target_accession}: {e}")
        results.append([row["source_accession"], target_accession, {}])  # Append empty dict if error occurs

print("New table with taxonomic levels created successfully!")

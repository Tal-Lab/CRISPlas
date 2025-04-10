import pandas as pd

# Load the CSV file without a header row
df = pd.read_csv(r"\path\to\your\working\directory\source_target_lineage.csv")

# Extract the first two words from column "2"
df["2"] = df["2"].apply(lambda x: " ".join(str(x).split()[:2]))

# Group by Plasmid accession number and species
df_combined = df.groupby(["0", "2"], as_index=False) [["9", "8", "7", "6", "5", "4", "3"]].first()

# Save to a new CSV file
df_combined.to_csv(r"\path\to\your\working\directory\unique_CSV.csv", index=False)

print("Processed table saved as unique_CSV.csv")

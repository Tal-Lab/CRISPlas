import pandas as pd

# Load the uniques file
file_path = r"\path\to\your\working\directory\uniques.csv"  # Update this with the correct file path
df = pd.read_csv(file_path)

# Rename columns
df.columns = ["Plasmid", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Domain"]

# Define column order (starting from Organism to Domain)
columns_to_check = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Domain"]

# Dictionary to store results
plasmid_first_identical = {}

# Loop through each Plasmid group
for plasmid, group in df.groupby("Plasmid"):
    for column in columns_to_check:
        if group[column].nunique() == 1:  # Check if all values are identical
            plasmid_first_identical[plasmid] = column
            break  # Stop at the first match

# Convert results to a DataFrame and save to Excel
result_df = pd.DataFrame(list(plasmid_first_identical.items()), columns=["Plasmid", "First_Identical_Column"])
result_df.to_excel(r"\path\to\your\working\directory\Host_Range.xlsx", index=False)

print(“Host_Range.xlsx table was created successfully!") 

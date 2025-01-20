import numpy as np
import glob
import shutil
import pandas as pd
import os

def load_deepfri_predictions(path):
    '''
    Loads deepFRI predictions from multiple CSV files into a single DataFrame, assuming all predictions are in the same folder.
    input: path do deepFRI predictions folder.
    output: DataFrame with protein name, go-term IDs, go-term names, deepFRI scores and aspect.
    '''
    all_data = []
    file_pattern = "*.csv"
    protein_names = set()

    for file in glob.glob(os.path.join(path, file_pattern)):
        # Extract protein name and aspect from file name
        file_basename = os.path.basename(file)
        parts = file_basename.split("_")
        if len(parts) == 3:
            protein_name = parts[0]
            aspect = parts[1]
            # Check for duplicate protein names
            if protein_name in protein_names:
                print(f"Duplicate protein name found: {protein_name}")
            else:
                protein_names.add(protein_name)
        else:
            print("Non-standard file name: ", file_basename)
            print("Expected format: proteinname_aspect_predictions.csv")
        
        # Load the CSV file, skipping the first row
        data = pd.read_csv(file, skiprows=1)

        # Add protein name and aspect to the DataFrame
        data["protein"] = protein_name
        data["aspect"] = aspect

        # Select and rename columns
        data = data[["protein", "GO_term/EC_number", "GO_term/EC_number name", "Score", "aspect"]]
        data.columns = ["protein", "go_term_id", "go_term_name", "score", "aspect"]

        all_data.append(data)

    # Filter out all-NA DataFrames
    all_data = [df for df in all_data if not df.isna().all().all()]

    return pd.concat(all_data, ignore_index=True)


def extract_file_subset(source_path, destination_path, n_files, file_pattern="*"):
    '''
    Extract a subset of files from a folder to another folder.
    input:
        source_path: path to the folder with files
        destination_path: path to the folder where the files will be copied
        n_files: number of files to be copied
        file_pattern (optional): pattern to filter files, e.g., "*.txt"
    '''
    os.makedirs(destination_path, exist_ok=True)
    files = glob.glob(os.path.join(source_path, file_pattern))
    for file in files[:n_files]:
        shutil.copy(file, destination_path)
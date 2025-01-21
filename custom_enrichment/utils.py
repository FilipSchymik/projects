import numpy as np
import glob
import shutil
import pandas as pd
import os
from scipy.stats import fisher_exact 
from statsmodels.stats.multitest import multipletests
import scipy.stats as stats

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


def Fischer_enrich(background_df, sample_df, score_threshold=0.0, fdr_threshold=1):
    
    '''
    Fischer's exact test for go-term enrichment analysis.

    Assumptions:
    - Each protein appears only once in each dataset.
    - No duplicated protein-GO term predictions.

    Args:
        background_df: DataFrame with background protein-GO term associations.
        sample_df: DataFrame with sample protein-GO term associations.
        score_threshold: Minimum score to consider for enrichment analysis (default: 0.0).
        fdr_threshold: False discovery rate threshold for significance (default: 1).

    Returns:
        DataFrame with enrichment results, including GO term ID, GO term name, p-value,
        adjusted p-value (FDR), and aspect.
    '''

    # Filter by score
    background_df = background_df[background_df["score"] >= score_threshold]
    sample_df = sample_df[sample_df["score"] >= score_threshold]

    # Create a dictionary for faster GO term name lookup
    go_term_name_dict = {}
    for _, row in background_df.iterrows():
        go_term_name_dict[row["go_term_id"]] = row["go_term_name"]
    for _, row in sample_df.iterrows():
        go_term_name_dict[row["go_term_id"]] = row["go_term_name"]

    results = []

    for aspect in ["BP", "MF", "CC"]:
        # Subset dataframes for the current aspect
        aspect_background_df = background_df[background_df["aspect"] == aspect]
        aspect_sample_df = sample_df[sample_df["aspect"] == aspect]

        # Get unique GO terms for the current aspect
        aspect_all_go_terms = set(aspect_background_df["go_term_id"]).union(set(aspect_sample_df["go_term_id"]))

        # Vectorized calculations for the contingency table
        go_term_counts_sample = aspect_sample_df.groupby("go_term_id")["protein"].count()
        go_term_counts_background = aspect_background_df.groupby("go_term_id")["protein"].count()

        total_go_terms_sample = go_term_counts_sample.sum()
        total_go_terms_background = go_term_counts_background.sum()

        for go_term in aspect_all_go_terms:
            # Contingency table using vectorized operations
            go_term_count_sample = go_term_counts_sample.get(go_term, 0)
            go_term_count_background = go_term_counts_background.get(go_term, 0)

            other_go_term_count_sample = total_go_terms_sample - go_term_count_sample
            other_go_term_count_background = total_go_terms_background - go_term_count_background

            # Perform Fisher's exact test (one-sided for enrichment)
            _, p_value = stats.fisher_exact(
                [[go_term_count_sample, other_go_term_count_sample],
                 [go_term_count_background, other_go_term_count_background]],
                alternative="greater"
            )

            results.append({
                "go_term_id": go_term,
                "go_term_name": go_term_name_dict.get(go_term, "Unknown"),
                "p_value": p_value,
                "aspect": aspect
            })

    # Create DataFrame from results
    results_df = pd.DataFrame(results)

    # Multiple testing correction (FDR)
    _, results_df["adjusted_p_value"], _, _ = multipletests(results_df["p_value"], method="fdr_bh")

    # Filter by FDR threshold
    results_df = results_df[results_df["adjusted_p_value"] <= fdr_threshold].reindex(
        columns=["go_term_id", "go_term_name", "p_value", "adjusted_p_value","aspect"])

    return results_df

def run_deepfri(pdb_folder, script_folder, output_folder, threads=1, conda_env_name="deepfri"):
    '''
    Run deepFRI predictions for all PDB files in a folder.
    input:
        pdb_folder: path to the folder with PDB files.
        script_folder: path to the folder with deepFRI predicty.py script.
        output_folder: path to the folder where the predictions will be saved.
        threads: number of threads to use (default: 1).
        conda_env_name: name of the conda environment with deepFRI installed (default: "deepfri").
    '''
    
    output_prefix = output_folder / pdb_folder.stem

    deepfri_command = [
    "conda",
    "run",
    "-n",
    conda_env_name,
    "python",
    f"{str(deepfri_script)}/predict.py",
    "--pdb_dir",
    str(pdb_file),
    "-ont",
    "mf", "bp", "cc",
    "--output_fn_prefix",
    str(output_prefix)
]
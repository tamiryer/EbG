import os
import warnings
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.linear_model import BayesianRidge
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

warnings.filterwarnings('ignore')

# Set category and directories
category = "Phylum"
phylum_directories = [
    ("NEON", f"/home/projects/zeevid/Analyses/2023-Tamir/NEON/Soil/all_samples_5k_freq/{category}"),
    ("Geotraces", f"/home/projects/zeevid/Analyses/2023-Tamir/Geotraces/all_samples_5k_freq/{category}"),
    ("Tara", f"/home/projects/zeevid/Analyses/2023-Tamir/Tara/PRJEB1787_prokaryote_DNA/all_samples_5k_freq/{category}")
]

general_samples_directories = [
    ("Tara", f"/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/fastp/Tara/Tara_combined/4mers/kmers_paired_tnf"),
    ("Geotraces", f"/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/fastp/GEOTRACES/mapped_tnf/4mers/kmers_paired_tnf"),
    ("NEON", f"/home/projects/zeevid/Analyses/2023-EbG/EnvironmentByGenome2024/GIT_IGNORED/fastp/NEON/Soil/mapped_tnf/4mers/kmers_paired_tnf_runID")
]

# Initialize a dictionary to store the count of each Phylum
phylum_counts = {}
phylum_project_counts = {}  # To store the distribution of projects for each phylum

# Collect Phyla and update counts
for project_name, tsv_directory in phylum_directories:
    for file_name in os.listdir(tsv_directory):
        if not file_name.endswith(".tsv"):
            continue
        tsv_df = pd.read_csv(os.path.join(tsv_directory, file_name), sep='\t')
        phyla = tsv_df.iloc[:, 0].tolist()
        for phylum in phyla:
            if phylum not in phylum_counts:
                phylum_counts[phylum] = 0
                phylum_project_counts[phylum] = {}
            phylum_counts[phylum] += 1
            phylum_project_counts[phylum][project_name] = phylum_project_counts[phylum].get(project_name, 0) + 1

# Filter Phyla
selected_phyla = [phylum for phylum, count in phylum_counts.items() if count >= 100]

# Initialize lists to store results for Phyla and General samples
results = {
    'mae_list_phyla': [], 'r2_list_phyla': [], 'pearson_list_phyla': [],
    'spearman_list_phyla': [], 'mae_list_general': [], 'r2_list_general': [],
    'pearson_list_general': [], 'spearman_list_general': [],
}

csv_file_path = "/home/projects/zeevid/tamirye/Data/United_all_Ocean-with_fasta_no_dup.csv"
csv_df = pd.read_csv(csv_file_path)

# Colors for different directories
dir_colors = ['blue', 'red', 'limegreen']
color_legend = ['NEON', 'Geotraces', 'Tara Oceans']

column_names = None  # Placeholder for feature names

# Function to perform Bayesian Ridge regression and return performance metrics
def perform_regression(X, y):
    model = BayesianRidge()
    y_pred = np.zeros(len(y))
    
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    for train, test in kf.split(X, y):
        X_train, X_test = X.iloc[train], X.iloc[test]
        y_train, y_test = y.iloc[train], y.iloc[test]
        model.fit(X_train, y_train)
        y_pred[test] = model.predict(X_test)
    
    mae = mean_absolute_error(y.values.flatten(), y_pred)
    r2 = r2_score(y.values.flatten(), y_pred)
    pearson_corr, _ = pearsonr(y.values.flatten(), y_pred.flatten())
    spearman_corr, _ = spearmanr(y.values.flatten(), y_pred.flatten())
    
    return mae, r2, pearson_corr, spearman_corr

# Dictionary to store all data for CSV
all_data = {key: [] for key in ['phylum', 'project_distribution', 'mae_phyla', 'r2_phyla', 'pearson_phyla', 'spearman_phyla', 
                                'bootstrap_mae', 'bootstrap_r2', 'bootstrap_pearson', 'bootstrap_spearman', 
                                'outliers_mae', 'outliers_r2', 'outliers_pearson', 'outliers_spearman']}

for phylum in selected_phyla:
    kmer_freq = []
    actual_temp = []

    # Collect data for the current phylum from specific directories
    for project_name, tsv_directory in phylum_directories:
        for file_name in os.listdir(tsv_directory):
            if not file_name.endswith(".tsv"):
                continue
            sample_name = os.path.splitext(file_name)[0]
            csv_row = csv_df[csv_df['RunID'] == sample_name]
            if csv_row.empty:
                continue
            temperature = csv_row.iloc[0, 4]
            tsv_df = pd.read_csv(os.path.join(tsv_directory, file_name), sep='\t')
            if column_names is None:
                column_names = tsv_df.columns.tolist()
            phylum_df = tsv_df[tsv_df.iloc[:, 0] == phylum]
            if phylum_df.empty:
                continue
            tetramer_freq = phylum_df.iloc[:, 1:137].values.flatten().tolist()
            kmer_freq.append(tetramer_freq)
            actual_temp.append(temperature)

    # Convert to DataFrames
    kmer_freq = pd.DataFrame(kmer_freq, columns=column_names[1:137])
    actual_temp = pd.DataFrame(actual_temp, columns=["Temperature"])
    valid_indices = ~actual_temp["Temperature"].isna()
    kmer_freq = kmer_freq[valid_indices]
    actual_temp = actual_temp[valid_indices]

    if actual_temp.empty: continue # Skip if no valid samples

    X, y = kmer_freq, actual_temp
    mae, r2, pearson_corr, spearman_corr = perform_regression(X, y)
    bootstrap_mae, bootstrap_r2, bootstrap_pearson, bootstrap_spearman = [], [], [], []

    for _ in range(1000):
        general_kmer_freq, general_temp = [], []
        for project_name, general_dir in general_samples_directories:
            if project_name not in phylum_project_counts[phylum]:
                continue
            num_samples_project = phylum_project_counts[phylum][project_name]
            available_files = [f for f in os.listdir(general_dir) if f.endswith(".csv")]
            sampled_files = np.random.choice(available_files, num_samples_project, replace=False)
            for file_name in sampled_files:
                general_df = pd.read_csv(os.path.join(general_dir, file_name))
                if '10000' not in general_df.columns:
                    continue
                sample_name = os.path.splitext(file_name)[0]
                csv_row = csv_df[csv_df['RunID'] == sample_name]
                if csv_row.empty:
                    continue
                temperature = csv_row.iloc[0, 4]
                tetramer_freq = general_df['10000'].values.flatten().tolist()
                general_kmer_freq.append(tetramer_freq)
                general_temp.append(temperature)

        general_kmer_freq = pd.DataFrame(general_kmer_freq)
        general_temp = pd.DataFrame(general_temp, columns=["Temperature"])
        valid_indices = ~general_temp["Temperature"].isna()
        general_kmer_freq = general_kmer_freq[valid_indices]
        general_temp = general_temp[valid_indices]

        if general_temp.empty: continue
        mae, r2, pearson_corr, spearman_corr = perform_regression(general_kmer_freq, general_temp)
        bootstrap_mae.append(mae)
        bootstrap_r2.append(r2)
        bootstrap_pearson.append(pearson_corr)
        bootstrap_spearman.append(spearman_corr)

    # Calculate outliers
    outliers_mae = [x for x in bootstrap_mae if np.abs(x - np.mean(bootstrap_mae)) > 2 * np.std(bootstrap_mae)]
    outliers_r2 = [x for x in bootstrap_r2 if np.abs(x - np.mean(bootstrap_r2)) > 2 * np.std(bootstrap_r2)]
    outliers_pearson = [x for x in bootstrap_pearson if np.abs(x - np.mean(bootstrap_pearson)) > 2 * np.std(bootstrap_pearson)]
    outliers_spearman = [x for x in bootstrap_spearman if np.abs(x - np.mean(bootstrap_spearman)) > 2 * np.std(bootstrap_spearman)]

    # Update results
    for key, value in zip(
        ['phylum', 'project_distribution', 'mae_phyla', 'r2_phyla', 'pearson_phyla', 'spearman_phyla', 
         'bootstrap_mae', 'bootstrap_r2', 'bootstrap_pearson', 'bootstrap_spearman', 
         'outliers_mae', 'outliers_r2', 'outliers_pearson', 'outliers_spearman'], 
        [phylum, phylum_project_counts[phylum], mae, r2, pearson_corr, spearman_corr, 
         bootstrap_mae, bootstrap_r2, bootstrap_pearson, bootstrap_spearman, 
         outliers_mae, outliers_r2, outliers_pearson, outliers_spearman]):
        all_data[key].append(value)

# Finally, save all results to a DataFrame
results_df = pd.DataFrame(all_data)
results_df.to_csv("/path/to/output/results.csv", index=False)
import scanpy as sc
import memento
import numpy as np
import pandas as pd
import os
import sys 
import math

# packages to manage input arguments
from argparse import ArgumentParser
import yaml

# plotting libraries (not needed, right?)
from matplotlib import pyplot as plt 
import seaborn as sns 
import scipy.stats as stats

# ----------------------------------------------------------------------- #
# Main
# ----------------------------------------------------------------------- #

def main():

    # initialize argument parser
    parser = ArgumentParser()

    # add command line arguments
    parser.add_argument('--config', help='the configuration .yml file containing all the settings for running the differential analysis snakemake', type=str, default='config.yml')
    parser.add_argument('--table', help='the .tsv file containing information about the different comparisons to run', type=str, default='config.tsv')
    parser.add_argument('--comparison', help='the name of the current comparison to perform the analysis on', type=str)
    parser.add_argument('--cluster', help='the name of the current cluster to perform the analysis on', type=str)
    parser.add_argument('--test', help='the name of the current DEG test, corresponding to a named entry inside config["tests"]', type=str)
    parser.add_argument('--outdir', help='path to output directory', type=str)

    # ----------------------------------------------------------------------- #
    # Part 0: parse arguments, set variables, perform checks
    # ----------------------------------------------------------------------- #

    # maybe I should make a function to validate meta-data?

    args = parser.parse_args() 
    config_path = args.config
    config_df_path = args.table
    cur_comparison = args.comparison
    cur_cluster = args.cluster
    cur_test = args.test
    outdir = args.outdir

    # load the config file:
    with open(config_path, "r") as config_file:
        config = yaml.safe_load(config_file)

    # load the comparison df 
    config_df = pd.read_table(config_df_path)

    # subset the table to get the options:
    cur_config = config_df[
        (config_df.comparison == cur_comparison) &
        (config_df.cluster == cur_cluster)
    ]

    # check cur_config 
    if cur_config.shape[0] == 0:
        raise ValueError(f"config info missing after subsetting for cluster: {cur_cluster} and comparison: {cur_comparison}. Check that the arguments --comparison and --cluster have matching entries in the config.tsv file.")

    # Load AnnData object
    if not os.path.exists(config["adata_path"]):
        raise FileNotFoundError(f"AnnData file '{adata_filepath}' not found.")
    adata = sc.read_h5ad(config["adata_path"])
                
    # cluster column
    group_var = cur_config["groupby"].values[0]
    cur_cluster = str(cur_cluster)
    ensure_categorical_groups(adata, group_var)
    validate_obs(adata, group_var, [cur_cluster])

    # condition info for the DEG comparison
    condition = str(cur_config["condition"].values[0])
    g1_name = str(cur_config["group1"].values[0])
    g2_name = str(cur_config["group2"].values[0])
    ensure_categorical_groups(adata, condition)
    validate_obs(adata, condition, [g1_name, g2_name])

    # replicate info
    replicate_col = cur_config["replicate"].values[0]
    ensure_categorical_groups(adata, replicate_col)
    validate_obs(adata, replicate_col)

    # additional subsetting
    if cur_config['subset_cols'].isnull().values[0] or cur_config['subset_groups'].isnull().values[0]:
        # Use empty lists for easy iteration later
        subset_cols = [] 
        subset_groups = [] 
    else:
        subset_cols = cur_config["subset_cols"].astype(str).values[0].split(';')
        subset_groups_str = cur_config["subset_groups"].astype(str).values[0].split(';')

        if len(subset_cols) != len(subset_groups_str):
            raise ValueError("Mismatched number of subset columns and subset groups.")

        subset_groups = [] 
        for cur_col, cur_groups_str in zip(subset_cols, subset_groups_str):
            
            # Convert each group value individually
            groups = [group_str.strip() for group_str in cur_groups_str.split(',')]
            ensure_categorical_groups(adata, cur_col)
            validate_obs(adata, cur_col, groups)
            subset_groups.append(groups)

    # additional exclusion
    if cur_config['exclude_cols'].isnull().values[0] or cur_config['exclude_groups'].isnull().values[0]:
        # Use empty lists for easy iteration later
        exclude_cols = [] 
        exclude_groups = [] 
    else:
        exclude_cols = cur_config["exclude_cols"].astype(str).values[0].split(';')
        exclude_groups_str = cur_config["exclude_groups"].astype(str).values[0].split(';')

        if len(exclude_cols) != len(exclude_groups_str):
            raise ValueError("Mismatched number of exclude columns and exclude groups.")

        exclude_groups = [] 
        for cur_col, cur_groups_str in zip(exclude_cols, exclude_groups_str):
            
            # Convert each group value individually
            groups = [group_str.strip() for group_str in cur_groups_str.split(',')]
            ensure_categorical_groups(adata, cur_col)
            validate_obs(adata, cur_col, groups)
            exclude_groups.append(groups)

    # additional memento settings
    test_options = config["tests"][cur_test]["options"]
    cell_thresh = test_options["cell_thresh"]
    filter_thresh = test_options["filter_thresh"]
    capture_rate = test_options["capture_rate"]
    num_boot = test_options["num_boot"]
    num_cpus = test_options["num_cpus"]
    verbose = test_options["verbose"]
    min_cells_group = config['min_cells_group']

    if not isinstance(cell_thresh, (int, float)) or cell_thresh <= 0:
        raise ValueError(f"Invalid cell_thresh value: {cell_thresh}. Must be a positive number.")
    
    if not isinstance(filter_thresh, (int, float)) or filter_thresh < 0:
        raise ValueError(f"Invalid filter_thresh value: {filter_thresh}. Must be non-negative.")

    if not isinstance(capture_rate, (int, float)) or not (0 <= capture_rate <= 1):
        raise ValueError(f"Invalid capture_rate value: {capture_rate}. Must be between 0 and 1.")

    if not isinstance(num_boot, int) or num_boot <= 0:
        raise ValueError(f"Invalid num_boot value: {num_boot}. Must be a positive integer.")

    if not isinstance(num_cpus, int) or num_cpus <= 0:
        raise ValueError(f"Invalid num_cpus value: {num_cpus}. Must be a positive integer.")

    if not isinstance(verbose, bool):
        raise ValueError(f"Invalid verbose value: {verbose}. Must be True or False.")

    print("All checks passed, proceeding with analysis.")

    # ----------------------------------------------------------------------- #
    # Part 1: Subset the data by cluster of interest & other variables
    # ----------------------------------------------------------------------- #

    # subset by the current cluster
    adata_test = adata[adata.obs[group_var] == cur_cluster].copy() 

    # subset by g1_name and g2_name 
    adata_test = adata_test[
        adata_test.obs[condition].isin([g1_name, g2_name])
    ].copy()
    print(adata_test.shape)

    # subset by additional groups
    if subset_cols: 
        for cur_subset, cur_groups in zip(subset_cols, subset_groups):
            
            print(cur_subset)
            
            adata_test = adata_test[
                adata_test.obs[cur_subset].isin(cur_groups) 
            ].copy()
            
            print(adata_test.shape)
    else:
        print("No additional subsetting applied.")

    # exclude additional groups:
    if exclude_cols: 
        for cur_exclude, cur_groups in zip(exclude_cols, exclude_groups):
            
            print(cur_exclude)
            
            adata_test = adata_test[
                ~ adata_test.obs[cur_exclude].isin(cur_groups)
            ].copy()
            
            print(adata_test.shape)
    else:
        print("No additional exclusions applied.")

    # finally, remove replicates with very few cells:
    replicates_remove = cur_config['remove_replicates'].astype(str).values[0].split(';')
    adata_test = adata_test[~adata_test.obs[replicate_col].isin(replicates_remove)].copy()

    # ----------------------------------------------------------------------- #
    # Part 1.1: Check if we have sufficient cells to run DEG
    # ----------------------------------------------------------------------- #

    # Check the counts of the two groups
    counts = adata_test.obs[condition].value_counts()

    # Get counts for the specified groups, using .get() to safely handle missing groups (count will be 0)
    count_g1 = counts.get(g1_name, 0)
    count_g2 = counts.get(g2_name, 0)

    can_run_deg = (count_g1 >= min_cells_group) and (count_g2 >= min_cells_group)

    if not can_run_deg:
        print(f"⚠️ Skipping DEG for {cur_cluster}: Insufficient cells.")
        print(f"Group 1 ({g1_name}): {count_g1} cells. Group 2 ({g2_name}): {count_g2} cells. Required: {min_cells_group}")
        empty_df = pd.DataFrame(
            columns=[
                'gene', 'cluster', 'group1', 'group2', 'de_pval', 'de_fdr', 'log2_de',
                'de_se', 'dv_pval', 'dv_fdr', 'log2_dv', 'dv_se', 'pct_group1', 'pct_group2'
            ] 
        )
        empty_df.to_csv(
            "{}/{}_DEGs.csv".format(outdir, cur_cluster),
            index=False
        )
        sys.exit(0)

    # ----------------------------------------------------------------------- #
    # Part 2: Set up data for running memento
    # ----------------------------------------------------------------------- #

    # set the .X to the raw counts
    # TODO: make the layer name an option in the config?
    adata_test.X = adata_test.layers['counts'].copy()

    # set up the binary test column based on g1_name and g2_name
    adata_test.obs['test'] = adata_test.obs[condition].astype(str).apply(lambda x: 1 if x == g1_name else 0)
    treatment_col = 'test'

    # additional setup for memento
    adata_test.obs['capture_rate'] = capture_rate
    memento.setup_memento(adata_test, q_column='capture_rate')
    memento.create_groups(adata_test, label_columns=[treatment_col, replicate_col])
    memento.compute_1d_moments(
        adata_test, 
        min_perc_group=0.9
    )

    # format sample meta
    sample_meta = memento.get_groups(adata_test)
    sample_meta[replicate_col] = sample_meta[replicate_col].astype('category')
    sample_meta.head()

    # format covariate DF 
    # the documentation is unclear for this step so I am not sure that
    # we are doing this correctly !!!
    treatment_df = sample_meta[[treatment_col]]
    cov_df = pd.get_dummies(sample_meta[replicate_col].astype('category'))
    cov_df.head()

    # ----------------------------------------------------------------------- #
    # Part 3: Run the memento DEG test
    # ----------------------------------------------------------------------- #

    memento.ht_1d_moments(
        adata_test, 
        treatment=treatment_df,
        covariate=cov_df,
        num_boot=num_boot, 
        verbose=verbose,
        num_cpus=num_cpus)
    deg_df = memento.get_1d_ht_result(adata_test)

    # calculate FDR:
    deg_df['de_fdr'] = stats.false_discovery_control(deg_df['de_pval'])
    deg_df['dv_fdr'] = stats.false_discovery_control(deg_df['dv_pval'])

    # convert the effect size from ln to log2
    deg_df['log2_de'] = deg_df['de_coef'] / np.log(2)
    deg_df['log2_dv'] = deg_df['dv_coef'] / np.log(2)

    # add info about g1, g2, cluster, etc
    deg_df['cluster'] = cur_cluster
    deg_df['group1'] = g1_name
    deg_df['group2'] = g2_name

    # drop unused columns
    cols_drop = ['tx', 'de_coef', 'dv_coef']
    deg_df = deg_df.drop(cols_drop, axis=1)

    # re-order the columns:
    cols_order = [
        "gene", "cluster", "group1", "group2",
        "de_pval", "de_fdr", "log2_de", "de_se",
        "dv_pval", "dv_fdr", "log2_dv", "dv_se"
    ]
    deg_df = deg_df[cols_order]

    # add % expressed for each gene in group1 & group2 similar 
    deg_df = percent_expressed(adata_test, deg_df, condition)

    # sort values by the effect size:
    deg_df = deg_df.sort_values(by=["log2_de"], ascending=False)

    # save the result as a .csv file
    deg_df.to_csv(
        "{}/{}_DEGs.csv".format(outdir, cur_cluster),
        index=False
    )

# ----------------------------------------------------------------------- #
# Helper functions
# ----------------------------------------------------------------------- #

# validates metadata in the anndata object
def validate_obs(adata, col, values=[]):
    
    if col not in adata.obs.columns:
        raise ValueError(f"Column '{col}' not found in adata.obs.")

    # check that the selected values are valid in this column
    valid = adata.obs[col].unique()
    valid_str = adata.obs[col].astype(str).unique()
    
    for value in values:
        if (value not in valid) and (value not in valid_str):
            print("Valid values:")
            print(valid)
            raise ValueError(f"Value '{value}' not found in adata.obs['{col}'].")

# calculate the percentage of cells in group1 & group2 expressing at least count for each gene
def percent_expressed(adata, deg_df, condition):
    
    # Binary matrix indicating whether each gene is expressed in each cell
    binary_matrix = adata.layers['counts'] > 0 # Adjust if data is log-transformed or sparse
    
    # Initialize new columns
    pct_group1 = []
    pct_group2 = []
    
    # Iterate over each cluster in the marker gene DataFrame
    cur_groups = [deg_df['group1'].values[0], deg_df['group2'].values[0]]
    
    # Create a mask for the cluster of interest
    cluster_mask = adata.obs[condition] == cur_groups[0]
    
    # Calculate percentages
    pct_in = np.asarray(binary_matrix[cluster_mask, :].mean(axis=0)).flatten() * 100
    pct_out = np.asarray(binary_matrix[~cluster_mask, :].mean(axis=0)).flatten() * 100
    
    # Map percentages back to marker genes
    for gene in deg_df.gene:
        gene_idx = adata.var_names.get_loc(gene)
        pct_group1.append(pct_in[gene_idx])
        pct_group2.append(pct_out[gene_idx])
    
    # Add the calculated percentages to the DataFrame
    deg_df["pct_group1"] = pct_group1
    deg_df["pct_group2"] = pct_group2

    return deg_df

def ensure_categorical_groups(adata, column_name):
    """
    Checks if an adata.obs column is numeric (float or int) and converts 
    it to a string/categorical type, first casting to integer to eliminate 
    decimal places (e.g., 2.0 -> "2").

    Args:
        adata (anndata.AnnData): The AnnData object.
        column_name (str): The name of the column in adata.obs to check and convert.
    """
    obs_series = adata.obs[column_name]
    target_dtype = obs_series.dtype

    if pd.api.types.is_numeric_dtype(target_dtype):
        print(f"Converting numeric column '{column_name}' to categorical strings...")
        
        try:
            temp_int = obs_series.astype(int) 
            adata.obs[column_name] = temp_int.astype(str)
            adata.obs[column_name] = adata.obs[column_name].astype('category')
            print(f"Successfully converted '{column_name}' to category/string.")
            
        except TypeError as e:
            print(f"Warning: Column '{column_name}' appears numeric but contains non-integer values or NaNs that prevent int conversion. Keeping as original type for now.")
            print(f"Original error: {e}")
            
    else:
        if not isinstance(target_dtype, pd.CategoricalDtype):
             adata.obs[column_name] = adata.obs[column_name].astype('category')
        print(f"Column '{column_name}' is already categorical/string. No conversion performed.")

# ----------------------------------------------------------------------- #
# execute the main
# ----------------------------------------------------------------------- #

if __name__ == '__main__':
    exit(main())

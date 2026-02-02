
```bash

# install pretty-jupyter
# pip install pretty-jupyter
# mamba install -c conda-forge -c bioconda snakemake
# mamba install snakemake-executor-plugin-cluster-generic

# launch an interactive node
interactive 8 256G 12:00:00 genD

# to launch:
module load GCC
conda activate compact

# launching a jupyter notebook (for testing & developing)
jupyter lab --no-browser --ip=$(hostname -s) --port=8366
ssh -L 8367:cnd09:8367 cnag6

# python jupyter notebook
module load GCC
conda activate scanpy-081124
jupyter lab --no-browser --ip=$(hostname -s) --port=8355
ssh -L 8367:cnd09:8367 cnag6

# set up the directory for the snakemake pipeline
cd ~/analysis/SERPENTINE/July_2025/hdWGCNA/TumorEpi/
# cd ~/analysis/SERPENTINE/July_2025/hdWGCNA/Tcells/

mkdir logs
mkdir results
mkdir notebooks

python -m ipykernel install --user --name myenv --display-name "Python (myenv)"


# dry-run
snakemake --dry-run

# unlock (if necessary)
snakemake --unlock

# snakemake \
#     --jobs 50 \
#     --latency-wait 240 \
#     --cluster "sbatch --time=4:00:00 --mem=128G --cpus-per-task=8 \
#     --output=logs/slurm-%j.out \
#     --error=logs/slurm-%j.err" \
#     --keep-going


# run the snakemake pipeline!!!
# this is snakemake V8, but I should think about downgrading to V7
snakemake \
    --jobs 50 \
    --latency-wait 60 \
    --keep-going \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "
    sbatch \
        --time=5:00:00 \
        --mem={resources.mem_mb} \
        --cpus-per-task={threads} \
        --output=logs/slurm-%j.out \
        --error=logs/slurm-%j.err"


#-----------------------------------------------------------#
# testing below
#-----------------------------------------------------------#

cd ~/analysis/SERPENTINE/July_2025/hdWGCNA/


cur_net="Liver_T1"
out_nb="test/$cur_net.ipynb"
echo $out_nb


papermill run_hdWGCNA.ipynb $out_nb \
    -p config_file config.yml \
    -p wgcna_name $cur_net



# export XDG_RUNTIME_DIR=$HOME/.quarto_runtime
# mkdir -p $XDG_RUNTIME_DIR

# quarto render run_hdWGCNA.ipynb \
#     --execute \
#     --to html \
#     --to ipynb \
#     -P config_file=config.yml \
#     -P wgcna_name=$cur_net


# convert it to an .html

jupyter nbconvert --to html $out_nb


jupyter nbconvert --to html --template pj $out_nb

jupyter nbconvert --to html --template pj $out_nb --TagRemovePreprocessor.remove_all_outputs_tags='["hide_output"]' --TagRemovePreprocessor.enabled=True


# need to do this in order to format the output
jupyter nbconvert run_hdWGCNA.ipynb --to html  --TagRemovePreprocessor.remove_all_outputs_tags='["hide_output"]' --TagRemovePreprocessor.enabled=True


```
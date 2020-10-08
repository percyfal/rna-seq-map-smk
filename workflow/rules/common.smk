import re
from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")

##############################
## Wildcard constraints
##############################
wildcard_constraints:
    sample = "|".join(samples.index),
    unit = "|".join(units["unit"])


##############################
## get functions
##############################
def get_final_output(wildcards):
    """Retrieve a list of the workflow final outputs"""
    gmap = expand("data/interim/gmap/{db}-{sample}.gmap.psl",
                  db=list(config["ref"].keys()),
                  sample=samples.index)
    return gmap

def get_sample(wildcards):
    return [re.sub("(.gz|.zip)$", "", x) for x in units.loc[wildcards.sample]["fasta"].tolist()]



def get_reference(wildcards):
    # In principle could use ensembl wrapper if db not present:
    # file:// would be better
    return config["ref"][wildcards.genome]["db"]

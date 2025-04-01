"""
Snakemake file for flood-wealth
"""

import os.path

import numpy as np
import requests
import geopandas as gpd


configfile: "config/config.yaml"

##### load rules #####
include: "rules/download/admin_boundaries.smk"
include: "rules/download/jrc_flood.smk"
include: "rules/download/relative_wealth_index.smk"
include: "rules/download/ghs_pop.smk"

include: "rules/prepare/ghs_pop.smk"
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
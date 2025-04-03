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
include: "rules/download/giri_flood.smk"
include: "rules/download/wri_flood.smk"
include: "rules/download/relative_wealth_index.smk"
include: "rules/download/ghs_pop.smk"

include: "rules/prepare/ghs_pop.smk"
include: "rules/prepare/jrc_flood.smk"
include: "rules/prepare/giri_flood.smk"
include: "rules/prepare/wri_flood.smk"
include: "rules/prepare/relative_wealth_index.smk"

include: "rules/analyze/flood_risk.smk"
include: "rules/analyze/social_flood_analysis.smk"

include: "rules/plot/figures.smk"
include: "rules/plot/maps.smk"

"""
Script for extracting country event maps for specific global flood database events.
"""

import logging
import sys
import glob
import os
from pathlib import Path
import subprocess

import numpy as np
import ast
import rasterio
import json
import yaml

if __name__ == "__main__":
    try:
        raw_flood_file: str = snakemake.input["raw_flood_file"]
        json_file: str = snakemake.input["json_file"]
        output_dir: str = snakemake.output["flood_event_dir"]
        event_id: str = snakemake.wildcards["event_id"]
    except:
        raise ValueError("Must be run via snakemake.")

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

# Load config file and pull country list
current_file = Path(__file__).resolve()
config_path = current_file.parents[2] / "config" / "config.yaml"
# Load config file
with open(config_path, "r") as file:
    config = yaml.safe_load(file)
# Pull ISO list and gfd code mapping from config
valid_countries = config.get("iso_codes", [])
gfd_code_mapping = config.get("gfd_code_mapping", {})

logging.info(f"Preparing country data for DFO event {event_id}.")

# Read in the json file to get country codes
logging.info("Reading json file for flood event information.")
with open(json_file, "r") as json_in:
    data = json.load(json_in)

raw_codes = data.get("gfd_country_code", "[]")

try:
    gfd_codes = ast.literal_eval(raw_codes)
except (ValueError, SyntaxError):
    gfd_codes = []

country_list = []
invalid_gfd_codes = []
seen = set()

for code in gfd_codes:
    iso3 = gfd_code_mapping.get(code)

    if iso3 is None:
        invalid_gfd_codes.append(code)
        continue

    if iso3 not in valid_countries:
        continue

    if iso3 in seen:
        continue

    seen.add(iso3)
    country_list.append(iso3)

logging.info(f"Working on the following countries: {country_list}. Checking validitiy...")
# Comparing lists to check validity
invalid_countries = [country for country in country_list if country not in valid_countries]
valid_countries = [country for country in country_list if country in valid_countries]
# Create output folder (if it doesn't already exist)
os.makedirs(output_dir, exist_ok=True)
# Write valid and invalid countries to json file
with open(os.path.join(output_dir, "countries.json"), 'w') as json_file:
    json.dump({"valid": valid_countries, "invalid": invalid_countries}, json_file, indent=4)

logging.info(f"{len(valid_countries)} valid countries found out of {len(country_list)}.")

logging.info("Extracting country-specific flood event data.")
for country in valid_countries:
    logging.info(f"Working on {country}...")
    # Load boundary and pop path (pop will be used for clipping) WARNING HARD CODING PATHS HERE
    boundary_path = current_file.parents[2] / "data" / "inputs" / "boundaries" / f"{country}" / f"geobounds_{country}.geojson"
    pop_path = current_file.parents[2] / "data" / "inputs" / "analysis" / "countries"/ f"{country}" / f"{country}_ghs-pop.tif"
    if not os.path.exists(boundary_path):
        sys.exit(f"Boundary file for {country} not found. Consider running all_boundaries snakemake rule.")
    if not os.path.exists(pop_path):
        sys.exit(f"No pop file found for {country}. Consider running clip_ghs_pop snakemake rule.")
    
    # Get the bounding box of the population file
    results = subprocess.run(['gdalinfo', "-json", pop_path], capture_output=True, text=True, check=True)
    info = json.loads(results.stdout)
    cc = info.get('cornerCoordinates', {})
    te_values = [cc['upperLeft'][0], cc['upperRight'][1], cc['lowerLeft'][0], cc['lowerRight'][1]]
    te_args = list(map(str, te_values))
    # Build the gdalwarp command
    gdal_cmd = [
        'gdalwarp',
        "-cutline", boundary_path,
        "-crop_to_cutline",
        "-tr", "0.00083333333333333", "0.00083333333333333",
        "-tap",
        "-te_srs", "EPSG:4326",
        "-te", *te_args,
        "-of", "GTiff",
        "-co", "COMPRESS=LZW",
        "-co", "BIGTIFF=YES",
        raw_flood_file,
        os.path.join(output_dir, f"{country}_{event_id}.tif")
    ]

    logging.info("Running gdalwarp command...")
    subprocess.run(gdal_cmd, check=True)

logging.info("Done.")


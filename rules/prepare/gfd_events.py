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
import geopandas as gpd
import ast
import rasterio
from rasterio.features import geometry_mask
from shapely.geometry import box
import json
import yaml

if __name__ == "__main__":
    try:
        raw_flood_file: str = snakemake.input["raw_flood_file"]
        adm0_file: str = snakemake.input["adm0_file"]
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

logging.info(f"Preparing country data for DFO event {event_id}.")

# Get flood raster bounds and identify overlapping countries
logging.info("Reading flood raster bounds and overlaying with global ADM0...")
with rasterio.open(raw_flood_file) as src:
    flood_bounds = src.bounds
    flood_crs = src.crs
    flood_transform = src.transform
    flood_shape = (src.height, src.width)
    
    # Create bounding box polygon
    bbox = box(flood_bounds.left, flood_bounds.bottom, flood_bounds.right, flood_bounds.top)
    
    # Load ADM0 boundaries that intersect with flood bbox (for efficiency)
    logging.info("Loading ADM0 boundaries intersecting flood extent...")
    adm0 = gpd.read_file(adm0_file, bbox=bbox)
    
    if adm0.crs != flood_crs:
        adm0 = adm0.to_crs(flood_crs)
    
    # Read the flood data to check actual pixel overlap
    flood_data = src.read(1)
    nodata = src.nodata
    
    # Create mask of valid flood pixels
    if nodata is not None:
        valid_flood_mask = (flood_data != nodata) & (flood_data > 0)
    else:
        valid_flood_mask = flood_data > 0

# Check which countries actually have flood pixels
logging.info("Checking which countries have actual flood pixel overlap...")
country_list = []

for idx, row in adm0.iterrows():
    geom = row.geometry
    # ISO3 column name in geoboundaries file
    iso3 = row.get('shapeGroup')
    
    if iso3 is None:
        continue
    
    # Create a mask for this country's geometry
    try:
        country_mask = geometry_mask(
            [geom],
            out_shape=flood_shape,
            transform=flood_transform,
            invert=True  # True inside geometry
        )
        
        # Check if any valid flood pixels fall within this country
        overlap = np.any(valid_flood_mask & country_mask)
        
        if overlap:
            country_list.append(iso3)
            logging.info(f"  Found flood overlap in {iso3}")
            
    except Exception as e:
        logging.warning(f"Could not process geometry for {iso3}: {e}")
        continue

logging.info(f"Found {len(country_list)} countries with flood pixel overlap: {country_list}")
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


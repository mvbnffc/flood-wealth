"""
Script for creating a merged raster of the Google historical inundation datasets

Raw data comes in km^2 GeoJSON tiles classfied as high, medium, and low risk:
    High: wet at least 5% of the time
    Medium: wet at least 1% of the time
    Low: wet at least 0.5% of the time

This script merges the GeoJSONs into one global (km^2 resolution) GeoTiff with cell values
corresponding to the probability of flooding in each grid cell (high=0.05, medium=0.01, low=0.005)
"""

import logging
import sys
import glob
import os

import numpy as np
import rasterio
import zipfile

if __name__ == "__main__":
    try:
        input_path: str = snakemake.input["raw_gfd_folder"]
        output_path: str = snakemake.output["merge_gfd_folder"]
    except:
        raise ValueError("Must be run via snakemake.")

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Preparing the Global Flood Database maps for a global merge.")

#### Define functions for the analysis ####
def extract_files(zip_path):
    # Create unzip folder (if it doesn't already exist)
    os.makedirs(os.path.join(input_path, "unzipped"), exist_ok=True)

    with zipfile.ZipFile(zip_path, 'r') as zf:
        # List all files that end with .tif
        tif_files = [f for f in zf.namelist() if f.lower().endswith('.tif')]
        json_files = [f for f in zf.namelist() if f.lower().endswith('.json')]

        if not tif_files:
            raise ValueError(f"No TIFF file found in {zip_path}")
        if not json_files:
            raise ValueError(f"No JSON file found in {zip_path}")

        tif_file = tif_files[0]
        json_file = json_files[0]
        output_tif = os.path.join(os.path.join(input_path, "unzipped"), os.path.basename(tif_file))
        output_json = os.path.join(os.path.join(input_path, "unzipped"), os.path.basename(json_file))

        # Extract the TIFF file's content and write it to output file
        with open(output_tif, 'wb') as tif_out:
            tif_out.write(zf.read(tif_file))
        # Extract the TIFF file's content and write it to output file
        with open(output_json, 'wb') as json_out:
            json_out.write(zf.read(json_file))


def process_raster(file, output_folder):
    with rasterio.open(file) as src:
        raster = src.read(1)  # Read the first band

        # Replace no-data values and NaNs with zero
        nodata = src.nodata
        if nodata is not None:
            raster[raster == nodata] = 0
        raster = np.nan_to_num(raster)  # Converts NaN values to zero

        output_file = os.path.join(output_folder, os.path.basename(file))

        # Save the modified raster
        with rasterio.open(output_file, 'w', driver='GTiff', height=raster.shape[0], width=raster.shape[1], count=1, dtype=raster.dtype, crs=src.crs, transform=src.transform) as dst:
            dst.write(raster, 1)

logging.info("Creating output directory (if it doesn't already exist)")
os.makedirs(output_path, exist_ok=True)

logging.info("Extract the rasters (and JSON files) from the zipped folder.")
zipped_files = glob.glob(os.path.join(input_path, "*.zip"))
for zipped_file in zipped_files:
    extract_files(zipped_file)

logging.info("Loop over TIF files and converting NaNs - skipping those for specified in config folder.")
# Finding GFD IDs to skip (those associated with dambreak or coastal flooding)
skip_files = []
with open("config/gfd_ignore.txt", "r") as f:
    for line in f.readlines():
        skip_files.append(line.strip())
tif_files = glob.glob(os.path.join(input_path, "unzipped", "*.tif"))
json_files = glob.glob(os.path.join(input_path, "unzipped", "*.json"))

for file in tif_files:
    if not any(skip_id in os.path.basename(file) for skip_id in skip_files):
        process_raster(file, output_path)

logging.info("Copying relevant JSON files.")
# Create json folder
os.makedirs(os.path.join(output_path, "json"), exist_ok=True)
# Also want to extract the JSON files.
for file in json_files:
    if not any(skip_id in os.path.basename(file) for skip_id in skip_files):
        with open(file, 'r') as json_in:
            json_out_path = os.path.join(output_path, "json", os.path.basename(file))
            with open(json_out_path, 'w') as json_out:
                json_out.write(json_in.read())

logging.info("Done.")


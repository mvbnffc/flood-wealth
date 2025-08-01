"""
This script prepares the river network dataset by clipping it to the country boundary
and then filtering out rivers that are smaller than a certain threshold (UDA 500 km^2).
"""

import logging
import sys
import glob
import os

import geopandas as gpd

if __name__ == "__main__":
    try:
        river_path: str = snakemake.input["river_network_file"]
        boundary_path: str = snakemake.input["boundary_file"]
        output_path: str = snakemake.output["river_file"]
        country: str = snakemake.wildcards["ISO3"]
    except:
        raise ValueError("Must be run via snakemake.")

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Preparing the river network dataset for {country}.")

logging.info("Creating output directory (if it doesn't already exist)")
out_dir = os.path.dirname(output_path)
os.makedirs(out_dir, exist_ok=True)

logging.info("Load files.")
rivers = gpd.read_file(river_path)
boundary = gpd.read_file(boundary_path)
# Make sure both layers are in same CRS
rivers = rivers.to_crs(boundary.crs)

logging.info("Clipping the river network to the country boundary.")
rivers_clipped = gpd.clip(rivers, boundary)

logging.info("Filtering the river network to 500 km^2 catchment area.")
catchment_area = 500 # WARNING HAD CODING THIS
rivers_filtered = rivers_clipped[rivers_clipped['UPLAND_SKM'] >= catchment_area]

logging.info("Write the new river dataset to file")
rivers_filtered.to_file(output_path, driver="GPKG")

logging.info("Done.")


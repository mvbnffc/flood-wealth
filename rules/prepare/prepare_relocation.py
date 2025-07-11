"""
This script creates an (adapted) relocation protection GeoTiff to use in the risk analysis.
This GeoTiff is created by specifying a (minimum) degree of urbanization where relocation will take place.
Those elgible for relocation are those exposed to flooding greater than 1 m at a 10-year return period.
We use 10-year RP as this is the lowest RP that is common across datasets. Will only relocate people if the 
baseline level of flood protection is below 10-year RP.

To integrate relocation within the risk analysis we burn in a 1000-year (the highest possible RP) protection level
for the grid cells that are eligible for relocation.

Outputs are a relocation protection GeoTiff and a text file with the number of people that will be relocated.
"""

import logging
import sys
import glob
import os

import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.mask import mask
import geopandas as gpd
from collections import Counter
from pyproj import Geod
from shapely.geometry import LineString, MultiLineString

if __name__ == "__main__":
    try:
        flopros_path: str = snakemake.input["flopros_path"]
        flood_path: str = snakemake.input["flood_path"]
        pop_path: str = snakemake.input["pop_path"]
        urbanization_path: str = snakemake.input["urbanization_path"]
        output_path: str = snakemake.output["flood_protection"]
        cost_output_path: str = snakemake.output["people_relocated"]
        country: str = snakemake.wildcards["ISO3"]
        urban_class: str = snakemake.wildcards["urban_class"]
        model: str = snakemake.wildcards["model"]
    except:
        raise ValueError("Must be run via snakemake.")

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Creating a relocation protection layer for {country} using {model} model. Also estimating # of people relocated. Urbanization class: {urban_class}.")

logging.info("Creating output directories (if they doesn't already exist)")
out_dir = os.path.dirname(output_path)
cost_out_dir = os.path.dirname(cost_output_path)
os.makedirs(out_dir, exist_ok=True)
os.makedirs(cost_out_dir, exist_ok=True)

logging.info("Load files.")
with rasterio.open(flopros_path) as flopros_src, rasterio.open(flood_path) as flood_src, \
     rasterio.open(pop_path) as pop_src, rasterio.open(urbanization_path) as urban_src:
    flopros = flopros_src.read(1)
    flood = flood_src.read(1)
    pop = pop_src.read(1)
    urban = urban_src.read(1)
    profile = flood_src.profile.copy()

logging.info("Building a mask for the cells to relocate.")
mask = (
    (urban <= int(urban_class)) &
    (flopros < 10) & # Only relocate if flood protection is below 10-year RP
    (flood >= 1) # Only relocate if flood depth is greater than (or equal to) 1 m
)

logging.info("Calculating the number of people to relocate.")
people_to_relocate = int(np.nansum(pop[mask]))
logging.info(f"Number of people to relocate: {people_to_relocate}")

logging.info("Creating the relocation protection layer.")
new_protection = flopros.copy()
new_protection[mask] = 1000  # Set the relocation protection level to 1000-year RP for the cells to relocate

logging.info("Writing the relocation protection layer to a GeoTiff and the number of people relocated to a text file.")
profile.update(
    dtype=rasterio.float32,
    count=1,
    compress="lzw"
)
with rasterio.open(output_path, "w", **profile) as dst:
    dst.write(new_protection.astype("float32"), 1)

with open(cost_output_path, "w") as f:
    f.write(f"{people_to_relocate} people relocated.\n")

logging.info("Done.")



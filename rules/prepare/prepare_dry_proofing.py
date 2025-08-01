"""
This script creates an (adapted) dry proofing GeoTiff to use in the risk analysis.
This GeoTiff is created by passing a flood map and a residential buidling area dataset.
Buildings elgible for dry proofing are those exposed to a 500-year flood extent and those outside
the 1 m flood depth extent for the 10-year RP flood.
We use 10-year and 500-year RP as these are the lowest and highest RP, respectively, across the GFM datasets. 

Approach is based off the following paper: Mortensen et al (2024) https://nhess.copernicus.org/articles/24/1381/2024/

Outputs are a masked built up area dataset and a text file with the area of buildings to be retrofitted.
"""

import logging
import sys
import glob
import os

import numpy as np
import rasterio

if __name__ == "__main__":
    try:

        rp10_path: str = snakemake.input["rp10_path"]
        rp500_path: str = snakemake.input["rp500_path"]
        res_path: str = snakemake.input["res_path"]
        output_path: str = snakemake.output["dry_proofed_buildings"]
        cost_output_path: str = snakemake.output["area_protected"]
        country: str = snakemake.wildcards["ISO3"]
        model: str = snakemake.wildcards["model"]
    except:
        raise ValueError("Must be run via snakemake.")

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Creating a dry proofing protection layer for {country} using {model} model. Also estimating area of buildings to be retrofitted.")

logging.info("Creating output directories (if they doesn't already exist)")
out_dir = os.path.dirname(output_path)
cost_out_dir = os.path.dirname(cost_output_path)
os.makedirs(out_dir, exist_ok=True)
os.makedirs(cost_out_dir, exist_ok=True)

logging.info("Load files.")
with rasterio.open(rp10_path) as rp10_src, rasterio.open(rp500_path) as rp500_src, rasterio.open(res_path) as res_src:
    rp10 = rp10_src.read(1)
    rp500 = rp500_src.read(1)
    res = res_src.read(1)
    profile = rp10_src.profile.copy()

logging.info("Building a mask for the cells to dry proof.")
mask = (
    (rp10 < 1) & # Only dry_proof if exposure is outside the > 1 m RP10 flood extent
    (rp500 > 0) # Only dry_proof if building is within the 500-year RP flood extent
)

logging.info("Calculating the area of buildings to dry-proof.")
area_to_dryproof = int(np.nansum(res[mask]))
logging.info(f"Area of buildings to dry-proof: {area_to_dryproof} m^2.")

logging.info("Creating the relocation protection layer.")
dry_proof_layer = res.copy()
dry_proof_layer[~mask] = 0  # Set the areas that won't be dry-proofed to 0

logging.info("Writing the dry-proofing layer to a GeoTiff and the area of buildings to be protected to a text file.")
profile.update(
    dtype=rasterio.float32,
    count=1,
    compress="lzw"
)
with rasterio.open(output_path, "w", **profile) as dst:
    dst.write(dry_proof_layer.astype("float32"), 1)

with open(cost_output_path, "w") as f:
    f.write(f"{area_to_dryproof} m^2 dry-proofed.\n")

logging.info("Done.")



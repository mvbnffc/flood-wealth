"""
Distribute IMF capital stock data across an infrastructure density layer.
This is used to estimate the infrastructure capital stock at risk from flooding.
"""

import logging
import os

import argparse
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling



if __name__ == "__main__":

    try:
        infra_path: str = snakemake.input["infrastructure_raster"]
        capstock_path: str = snakemake.input["capital_stock"]
        output_path: str = snakemake.output["infra_value"]
        country: str = snakemake.wildcards["ISO3"]
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Assigning capital stock values to the the OSM infrastructure layer for {country}.")

logging.info("Reading 2019 public capital stock value from IMF spreadsheet.")
df = pd.read_excel(capstock_path, sheet_name="Dataset")
df['isocode'] = df['isocode'].astype(str).str.upper()  # Ensure ISO codes are uppercase
# Find the row for the given country
rows = df[(df['isocode'] == country) & (df['year'] == 2019)] # want latest year, which for this dataset is 2019
if rows.empty:
    raise ValueError(f"No capital stock data found for {country} in 2019.")
if len(rows) > 1:
    raise ValueError(f"Multiple rows found for {country} in 2019. Please check the dataset.")
capital_stock = rows.iloc[0]['kgov_rppp'] # public capital stock in 2019 in billions of USD
logging.info(f"Public capital stock for {country} in 2019: {capital_stock} billion USD.")

logging.info("Reading infrastructure raster.")
with rasterio.open(infra_path) as src:
    meta = src.meta.copy()
    meta.update({
        "count": 1,
        "dtype": "float32",
        "compress": "lzw",
        "BIGTIFF": "YES",
        "nodata": 0.0,
    })
    height, width = src.height, src.width
    weights_ma = src.read(1, masked=True)

# Check for valid pixels to weight by (not masked by no-data)
valid_mask = ~weights_ma.mask
weights = np.where(valid_mask, np.asarray(weights_ma.filled(0.0), dtype=np.float64), 0.0)

# Disallow negative weights
neg_count = (weights < 0).sum()
if neg_count > 0:
    logging.warning(f"Found {neg_count} negative weights in the infrastructure raster. Setting them to 0.")
    weights[weights < 0] = 0.0

logging.info("Allocate capital stock value to infrastructure density layer.")
w_sum = float(weights.sum())
if w_sum <= 0.0:
    # Fallback: uniform allocatino over valid pixels
    n_valid = int(valid_mask.sum())
    if n_valid == 0:
        raise ValueError("No valid pixels found in the infrastructure raster.")
    logging.warning("No positive weights found in the infrastructure raster. Allocating capital stock uniformly over valid pixels.")
    alloc = np.zeros_like(weights, dtype=np.float64)
    alloc[valid_mask] = capital_stock * 1e9 / n_valid  # convert billion USD to USD
else:
    alloc = (capital_stock * 1e9 * weights) / w_sum  # convert billion USD to USD

logging.info("Writing the capital stock allocation to a GeoTIFF.")
os.makedirs(os.path.dirname(output_path), exist_ok=True)
with rasterio.open(output_path, "w", **meta) as dst:
    dst.write(alloc.astype("float32"), 1)
    dst.write_mask(valid_mask)  # Write the valid mask to the output file

logging.info("Sanity check")
sum_out = float(alloc.sum())
logging.info(f"Check: sum(output) = {sum_out:.0f} (billions USD) vs IMF total {(capital_stock*1e9):.0f}.")

logging.info("Done.")
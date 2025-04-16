"""
This script plots the concentration curve for the administration level of interest.
"""

import logging
import json
import os

import rasterio
from rasterio.features import geometry_mask
import pandas as pd
import geopandas as gpd
import numpy as np
import shapely
import matplotlib.pyplot as plt
from tqdm import tqdm

if __name__ == "__main__":

    try:
        rwi_path: str = snakemake.input["rwi_file"]
        pop_path: str = snakemake.input["pop_file"]
        mask_path : str = snakemake.input["mask_file"]
        flood_path: str = snakemake.input["flood_file"]
        output_path: str = snakemake.output["figure_file"]
        country: str = snakemake.wildcards.ISO3
        event_id: str = snakemake.wildcards.event_id
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Plotting Concentration Curves for event DFO_{event_id} in {country}.")

logging.info("Creating Concentration Curve Figure Directory")
os.makedirs(os.path.dirname(output_path), exist_ok=True)

logging.info("Reading raster data.")
with rasterio.open(rwi_path) as rwi_src, rasterio.open(pop_path) as pop_src, \
    rasterio.open(mask_path) as mask_src, rasterio.open(flood_path) as flood_src:
    rwi = rwi_src.read(1)
    pop = pop_src.read(1)
    mask = mask_src.read(1)
    flood = flood_src.read(1)
    affine = flood_src.transform

logging.info("Preparing rasters for plotting (applying mask, etc.)") 
rwi[rwi==-999] = np.nan # convert -999 in RWI dataset to NaN    
water_mask = np.where(mask>50, np.nan, 1) # WARNING WE ARE HARD CODING PERM_WATER > 50% mask here
mask_areas = (
    ~np.isnan(pop) &
    ~np.isnan(rwi) &
    ~np.isnan(flood) &
    ~np.isnan(water_mask)
    )
# Flatten data
pop_flat = pop[mask_areas]
rwi_flat = rwi[mask_areas]
flood_flat = flood[mask_areas]
# Mask out zero-populatoin cells
valid = pop_flat > 0
pop_flat = pop_flat[valid]
rwi_flat = rwi_flat[valid]
flood_flat = flood_flat[valid]

logging.info("Preparing dataframes for plotting")
df = pd.DataFrame({
    'pop': pop_flat,
    'rwi': rwi_flat,
    'flood': flood_flat,
})

logging.info("Plotting Concentration Curves")
# Define both functions 
def compute_concentration_curve(data):
    # Sort by RWI ascending
    data = data.sort_values(by="rwi", ascending=True).copy()
    data['cum_pop'] = data['pop'].cumsum()
    data['flood_pop'] = data['flood'] * data['pop']  # flood burden
    data['cum_flood_pop'] = data['flood_pop'].cumsum()
    total_pop = data['pop'].sum()
    total_flood_pop = data['flood_pop'].sum()
    data['frac_pop'] = data['cum_pop'] / total_pop
    data['frac_flood'] = data['cum_flood_pop'] / total_flood_pop
    # Insert 0 at the beginning for a proper starting point at (0,0)
    x_vals = np.insert(data['frac_pop'].values, 0, 0)
    y_vals = np.insert(data['frac_flood'].values, 0, 0)
    return x_vals, y_vals

def plot_concentration_curve(x, y, title, output_path):
    plt.plot(x, y, label="Concentration Curve")
    plt.plot([0,1], [0,1], "--", color='gray', label='Equality Line')
    plt.xlabel('Cumulative Population (Wealth Rank)')
    plt.ylabel('Cumulative Flood Risk Exposure')
    plt.title(title)
    plt.legend()
    plt.savefig(output_path)
    plt.close()

x, y = compute_concentration_curve(df)

plot_concentration_curve(x, y, f"DFO_{event_id} Concentration Curve - {country}", output_path)

logging.info("Done.")
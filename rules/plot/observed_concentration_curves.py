"""
This script plots the observed concentration curve for the administration level of interest.
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
        admin_path: str = snakemake.input["admin_areas"]
        rwi_path: str = snakemake.input["rwi_file"]
        pop_path: str = snakemake.input["pop_file"]
        mask_path : str = snakemake.input["mask_file"]
        risk_path: str = snakemake.input["risk_file"]
        output_path: str = snakemake.output["figure_directory"]
        administrative_level: int = snakemake.wildcards.ADMIN_SLUG
        country: str = snakemake.wildcards.ISO3
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

admin_level = int(administrative_level.replace("ADM-", ""))
logging.info(f"Plotting Observed Concentration Curves for {country} at Admin Level {admin_level}.")

if not os.path.exists(output_path):
    logging.info("Creating Concentration Curve Figure Directory")
    os.makedirs(output_path)

logging.info("Reading raster data.")
with rasterio.open(rwi_path) as rwi_src, rasterio.open(pop_path) as pop_src, rasterio.open(risk_path) as risk_src, \
    rasterio.open(mask_path) as mask_src:
    rwi = rwi_src.read(1)
    pop = pop_src.read(1)
    mask = mask_src.read(1)
    risk = risk_src.read(1)
    affine = risk_src.transform
rwi[rwi==-999] = np.nan # convert -999 in RWI dataset to NaN
water_mask = np.where(mask>50, np.nan, 1) # WARNING WE ARE HARD CODING PERM_WATER > 50% mask here

logging.info("Reading level {admin_level} admin boundaries")

logging.info(f"Reading level {administrative_level} admin boundaries")
layer_name = f"ADM_ADM_{admin_level}"
admin_areas: gpd.GeoDataFrame = gpd.read_file(admin_path, layer=layer_name)
area_unique_id_col = f"GID_{admin_level}"
# retain names of larger, encompassing adminstrative units
contextual_name_cols = [f"GID_{i}" for i in range(0, admin_level)]
admin_areas = admin_areas[[area_unique_id_col, *contextual_name_cols, "geometry"]]
logging.info(f"There are {len(admin_areas)} admin areas to analyze.")

logging.info("Looping over admin regions and plotting concentration curves")
results = [] # List for collecting results
 # Loop over each admin region
for idx, region in tqdm(admin_areas.iterrows()):
    # Get the geometry for the current admin region
    geom = region["geometry"].__geo_interface__
    # Create a mask from the geometry
    mask_array = geometry_mask([geom],
                                transform=affine,
                                invert=True,
                                out_shape=rwi.shape)
    # Use the mask to clip each raster by setting values outside the region to nan
    rwi_clip = np.where(mask_array, rwi, np.nan)
    pop_clip = np.where(mask_array, pop, np.nan)
    risk_clip = np.where(mask_array, risk, np.nan)
    water_mask_clip = np.where(mask_array, water_mask, np.nan)
    # Mask out areas where not all rasters are valid
    mask = (
        ~np.isnan(pop_clip) &
        ~np.isnan(rwi_clip) &
        ~np.isnan(risk_clip) &
        ~np.isnan(water_mask_clip)
    )
    # Flatten data
    pop_flat = pop_clip[mask]
    rwi_flat = rwi_clip[mask]
    risk_flat = risk_clip[mask]
    # Mask out zero-populatoin cells
    valid = pop_flat > 0
    pop_flat = pop_flat[valid]
    rwi_flat = rwi_flat[valid]
    risk_flat = risk_flat[valid]

    logging.info("Preparing dataframes for plotting")
    df = pd.DataFrame({
        'pop': pop_flat,
        'rwi': rwi_flat,
        'flood': risk_flat,
    })
    
    # Pull region code for plotting (title, name, etc.)
    region_code = region[area_unique_id_col]
    region_code = region_code.replace(".", "-")
    region_code = region_code.replace("_", "-")

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

    def plot_concentration_curve(x, y, title, admin_output_dir):
        plt.plot(x, y, label="Concentration Curve")
        plt.plot([0,1], [0,1], "--", color='gray', label='Equality Line')
        plt.xlabel('Cumulative Population (Wealth Rank)')
        plt.ylabel('Cumulative Flood Risk Exposure')
        plt.title(title)
        plt.legend()
        plt.savefig(os.path.join(admin_output_dir, "%s.png" % title))
        plt.close()

    x, y = compute_concentration_curve(df)
    plot_concentration_curve(x, y, f"{region_code}observed_concentration_curve", output_path)

logging.info("Done.")
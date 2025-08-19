"""
This script calculates the capital stock losses for a country and sums them per admin region.
"""

import logging

import rasterio
from rasterio.features import geometry_mask
import pandas as pd
import geopandas as gpd
import numpy as np
import shapely
from tqdm import tqdm

if __name__ == "__main__":

    try:
        admin_path: str = snakemake.input["admin_areas"]
        res_risk_path: str = snakemake.input["res_risk_file"]
        nres_risk_path: str = snakemake.input["nres_risk_file"]
        infr_risk_path: str = snakemake.input["infr_risk_file"]
        res_capstock_path: str = snakemake.input["res_capstock_file"]
        nres_capstock_path: str = snakemake.input["nres_capstock_file"]
        infr_capstock_path: str = snakemake.input["infr_capstock_file"]
        mask_path: str = snakemake.input["mask_file"]
        output_path: str = snakemake.output["regional_losses"]
        administrative_level: int = snakemake.wildcards.ADMIN_SLUG
        model: str = snakemake.wildcards.MODEL
        country: str = snakemake.wildcards.ISO3
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

# Update notation for GADM
admin_level = int(administrative_level.replace("ADM", ""))

logging.info(f"Calculating capital stock losses for {country} at admin level {admin_level}.")

logging.info("Reading raster data.")
with rasterio.open(res_risk_path) as res_risk_src, rasterio.open(nres_risk_path) as nres_risk_src, \
     rasterio.open(infr_risk_path) as infr_risk_src, rasterio.open(res_capstock_path) as res_capstock_src, \
     rasterio.open(nres_capstock_path) as nres_capstock_src, rasterio.open(infr_capstock_path) as infr_capstock_src, \
     rasterio.open(mask_path) as mask_src:
    res_risk = res_risk_src.read(1)
    nres_risk = nres_risk_src.read(1)
    infr_risk = infr_risk_src.read(1)
    res_capstock = res_capstock_src.read(1)
    nres_capstock = nres_capstock_src.read(1)
    infr_capstock = infr_capstock_src.read(1)
    water_mask = mask_src.read(1)
    affine = res_risk_src.transform
    
# Create the water mask
water_mask = np.where(water_mask>50, np.nan, 1) # WARNING WE ARE HARD CODING PERM_WATER > 50% mask here

logging.info(f"Reading level {administrative_level} admin boundaries")
layer_name = f"ADM{admin_level}"
admin_areas: gpd.GeoDataFrame = gpd.read_file(admin_path, layer=layer_name)
area_unique_id_col = "shapeName"
admin_areas = admin_areas[[area_unique_id_col, "geometry"]]
logging.info(f"There are {len(admin_areas)} admin areas to analyze.")

logging.info("Looping over admin regions and calculating capital stock losses")
results = [] # List for collecting results
 # Loop over each admin region
for idx, region in tqdm(admin_areas.iterrows()):
    # Get the geometry for the current admin region
    geom = region["geometry"].__geo_interface__
    # Create a mask from the geometry
    mask_array = geometry_mask([geom],
                                transform=affine,
                                invert=True,
                                out_shape=res_risk.shape)
    # Use the mask to clip each raster by setting values outside the region to nan
    res_risk_clip = np.where(mask_array, res_risk, np.nan)
    nres_risk_clip = np.where(mask_array, nres_risk, np.nan)
    infr_risk_clip = np.where(mask_array, infr_risk, np.nan)
    res_capstock_clip = np.where(mask_array, res_capstock, np.nan)
    nres_capstock_clip = np.where(mask_array, nres_capstock, np.nan)
    infr_capstock_clip = np.where(mask_array, infr_capstock, np.nan)
    water_mask_clip = np.where(mask_array, water_mask, np.nan)
    
    # Mask out areas where not all rasters are valid
    mask = (
        ~np.isnan(res_risk_clip) &
        ~np.isnan(nres_risk_clip) &
        ~np.isnan(infr_risk_clip) &
        ~np.isnan(res_capstock_clip) &
        ~np.isnan(nres_capstock_clip) &
        ~np.isnan(infr_capstock_clip) &
        ~np.isnan(water_mask_clip)
    )
    # Flatten data
    res_risk_flat = res_risk_clip[mask]
    nres_risk_flat = nres_risk_clip[mask]
    infr_risk_flat = infr_risk_clip[mask]
    res_capstock_flat = res_capstock_clip[mask]
    nres_capstock_flat = nres_capstock_clip[mask]
    infr_capstock_flat = infr_capstock_clip[mask]

    # Calculate sectoral capital stock losses for the region
    res_losses = np.nansum(res_risk_flat * res_capstock_flat)
    nres_losses = np.nansum(nres_risk_flat * nres_capstock_flat)
    infr_losses = np.nansum(infr_risk_flat * infr_capstock_flat)

    # Append risk metrics to results list
    results.append({
         area_unique_id_col: region[area_unique_id_col],
         "res_losses": res_losses,
         "nres_losses": nres_losses,
         "infr_losses": infr_losses,
         "total_losses": res_losses + nres_losses + infr_losses,
         "geometry": region["geometry"]
    })

logging.info("Writing reults to GeoPackage.")
results_gdf = gpd.GeoDataFrame(results, geometry="geometry")
results_gdf.to_file(output_path, driver="GPKG")

logging.info("Done.")
"""
This script calculates the capital stock losses for a country and sums them per admin region.
"""

import logging

import rasterio
from rasterio.features import geometry_mask, rasterize
import pandas as pd
from scipy import ndimage
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
    
logging.info("Pre-computing masks.")
water_mask = np.where(water_mask > 50, False, True)  # Boolean mask instead of NaN

# Pre-compute global validity mask (areas where all rasters have valid data)
global_valid_mask = (
    ~np.isnan(res_risk) &
    ~np.isnan(nres_risk) &
    ~np.isnan(infr_risk) &
    ~np.isnan(res_capstock) &
    ~np.isnan(nres_capstock) &
    ~np.isnan(infr_capstock) &
    water_mask
)

logging.info(f"Reading level {administrative_level} admin boundaries")
layer_name = f"ADM{admin_level}"
admin_areas: gpd.GeoDataFrame = gpd.read_file(admin_path, layer=layer_name)
if layer_name == "ADM0":
    area_unique_id_col = "shapeName"
else:
    area_unique_id_col = "shapeID"
admin_areas = admin_areas[[area_unique_id_col, "shapeName", "geometry"]]
logging.info(f"There are {len(admin_areas)} admin areas to analyze.")

# OPTIMIZATION: vectorize geometry masking
# Create a dictionary mapping region index to geometry
geom_dict = {idx: geom for idx, geom in enumerate(admin_areas.geometry)}

# Create a single raster where each pixel contains the region ID it belongs to
region_ids = rasterize(
    [(geom, idx) for idx, geom in geom_dict.items()],
    out_shape=res_risk.shape,
    transform=affine,
    fill=-1,  # -1 for pixels not in any region
    dtype=np.int32
)

logging.info("Precompute national risk maps.")
res_loss_arr = res_risk * res_capstock
nres_loss_arr = nres_risk * nres_capstock
infr_loss_arr = infr_risk * infr_capstock
# Set invalid areas to 0 for faster summing
res_loss_arr[~global_valid_mask] = 0
nres_loss_arr[~global_valid_mask] = 0
infr_loss_arr[~global_valid_mask] = 0

logging.info("Calculating capital stock losses across the admin regions.")
res_losses = ndimage.sum_labels(res_loss_arr, labels=region_ids, index=np.arange(len(admin_areas)))
nres_losses = ndimage.sum_labels(nres_loss_arr, labels=region_ids, index=np.arange(len(admin_areas)))
infr_losses = ndimage.sum_labels(infr_loss_arr, labels=region_ids, index=np.arange(len(admin_areas)))

logging.info("Writing reults to GeoPackage.")
results_gdf = admin_areas.copy()
results_gdf["res_losses"] = res_losses
results_gdf["nres_losses"] = nres_losses
results_gdf["infr_losses"] = infr_losses
results_gdf["total_losses"] = res_losses + nres_losses + infr_losses
results_gdf.to_file(output_path, driver="GPKG")

logging.info("Done.")
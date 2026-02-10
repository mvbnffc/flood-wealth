"""
This script calculates flood risk statistics at the administrative level of interest
Outputs results to a CSV
"""

import logging

import rasterio
from rasterio.features import geometry_mask, rasterize
import numpy as np
from tqdm import tqdm
import pandas as pd
import geopandas as gpd

if __name__ == "__main__":

    try:
        admin_path: str = snakemake.input["admin_areas"]
        pop_path: str = snakemake.input["pop_file"]
        risk_path: str = snakemake.input["risk_file"]
        mask_path: str = snakemake.input["mask_file"]
        quintile_path: str = snakemake.input["quintile_file"]
        output_path: str = snakemake.output["risk_stats"]
        administrative_level: int = snakemake.wildcards.ADMIN_SLUG
        model = snakemake.wildcards.MODEL
        vuln = snakemake.wildcards.VULN
        country: str = snakemake.wildcards.ISO3
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

# Update notation for GADM
admin_level = int(administrative_level.replace("ADM", ""))
layer_name = f"ADM{admin_level}"

logging.info(f"Calculating flood risk stats in {country} for model {model}, vulnerability {vuln}.")

logging.info("Load the raster data.")
with rasterio.open(pop_path) as pop_ds, rasterio.open(quintile_path) as q_ds, rasterio.open(risk_path) as risk_ds, rasterio.open(mask_path) as mask_ds:
    pop_arr = pop_ds.read(1).astype("float32")
    q_arr = q_ds.read(1).astype("float32")
    risk_arr = risk_ds.read(1).astype("float32")
    water_mask = mask_ds.read(1).astype("float32")
    affine = pop_ds.transform

logging.info("Pre-computing masks.")
water_mask = np.where(water_mask > 50, False, True)  # Boolean mask instead of NaN

# Pre-compute global validity mask (areas where all rasters have valid data)
global_valid_mask = (
    ~np.isnan(pop_arr) &
    ~np.isnan(q_arr) &
    ~np.isnan(risk_arr) &
    water_mask
)

logging.info(f"Reading level {administrative_level} admin boundaries")
layer_name = f"ADM{admin_level}"
admin_areas: gpd.GeoDataFrame = gpd.read_file(admin_path, layer=layer_name)
if layer_name == "ADM0":
    area_unique_id_col = "shapeName"
    admin_areas = admin_areas[[area_unique_id_col, "geometry"]]
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
    out_shape=pop_arr.shape,
    transform=affine,
    fill=-1,  # -1 for pixels not in any region
    dtype=np.int32
)

logging.info("Precompute natioanl risk maps")
total_risk = pop_arr * risk_arr
q1_risk = np.where(q_arr == 1, total_risk, 0)
q2_risk = np.where(q_arr == 2, total_risk, 0)
q3_risk = np.where(q_arr == 3, total_risk, 0)
q4_risk = np.where(q_arr == 4, total_risk, 0)
q5_risk = np.where(q_arr == 5, total_risk, 0)
# Set invalid areas to zero for faster summing
total_risk[~global_valid_mask] = 0
q1_risk[~global_valid_mask] = 0
q2_risk[~global_valid_mask] = 0
q3_risk[~global_valid_mask] = 0
q4_risk[~global_valid_mask] = 0
q5_risk[~global_valid_mask] = 0

logging.info("Calculating flood risk stats for each admin region")
# Flatten arrays
flat_region_ids = region_ids.flatten()
flat_total_risk = total_risk.flatten()
flat_q1_risk = q1_risk.flatten()
flat_q2_risk = q2_risk.flatten()
flat_q3_risk = q3_risk.flatten()
flat_q4_risk = q4_risk.flatten()
flat_q5_risk = q5_risk.flatten()
# Filter valid regions
valid_mask = flat_region_ids >= 0
flat_region_ids = flat_region_ids[valid_mask]
flat_total_risk = flat_total_risk[valid_mask]
flat_q1_risk = flat_q1_risk[valid_mask]
flat_q2_risk = flat_q2_risk[valid_mask]
flat_q3_risk = flat_q3_risk[valid_mask]
flat_q4_risk = flat_q4_risk[valid_mask]
flat_q5_risk = flat_q5_risk[valid_mask]
# Sum risk by region using numpy's bincount
total_risk_by_region = np.bincount(flat_region_ids, weights=flat_total_risk, minlength=len(admin_areas))
q1_risk_by_region = np.bincount(flat_region_ids, weights=flat_q1_risk, minlength=len(admin_areas))
q2_risk_by_region = np.bincount(flat_region_ids, weights=flat_q2_risk, minlength=len(admin_areas))
q3_risk_by_region = np.bincount(flat_region_ids, weights=flat_q3_risk, minlength=len(admin_areas))
q4_risk_by_region = np.bincount(flat_region_ids, weights=flat_q4_risk, minlength=len(admin_areas))
q5_risk_by_region = np.bincount(flat_region_ids, weights=flat_q5_risk, minlength=len(admin_areas))

logging.info("Write results to CSV")
results_df = pd.DataFrame({
    area_unique_id_col: admin_areas[area_unique_id_col],
    "shapeName": admin_areas["shapeName"],
    "total_risk": total_risk_by_region,
    "q1_risk": q1_risk_by_region,
    "q2_risk": q2_risk_by_region,
    "q3_risk": q3_risk_by_region,
    "q4_risk": q4_risk_by_region,
    "q5_risk": q5_risk_by_region,
})

logging.info("Done.")
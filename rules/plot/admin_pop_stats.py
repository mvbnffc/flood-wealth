"""
This script calculates population statistics at the administrative level of interest
Outputs results to a CSV
"""

import logging

import rasterio
from rasterio.features import geometry_mask
import numpy as np
from tqdm import tqdm
import pandas as pd
import geopandas as gpd

if __name__ == "__main__":

    try:
        admin_path: str = snakemake.input["admin_areas"]
        pop_path: str = snakemake.input["pop_file"]
        quintile_path: str = snakemake.input["quintile_file"]
        output_path: str = snakemake.output["pop_stats"]
        administrative_level: int = snakemake.wildcards.ADMIN_SLUG
        country: str = snakemake.wildcards.ISO3
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

# Update notation for GADM
admin_level = int(administrative_level.replace("ADM", ""))
layer_name = f"ADM{admin_level}"

logging.info(f"Calculating population stats in {country}.")

logging.info("Load the raster data.")
with rasterio.open(pop_path) as pop_ds, rasterio.open(quintile_path) as q_ds:
    pop_arr = pop_ds.read(1).astype("float32")
    q_arr = q_ds.read(1).astype("float32")

    # Replace nodata with nan (for BOTH rasters)
    if pop_ds.nodata is not None:
        pop_arr[pop_arr == pop_ds.nodata] = np.nan
    if q_ds.nodata is not None:
        q_arr[q_arr == q_ds.nodata] = np.nan

    affine = pop_ds.transform
    out_shape = pop_arr.shape

logging.info(f"Reading level {administrative_level} admin boundaries")
layer_name = f"ADM{admin_level}"
admin_areas: gpd.GeoDataFrame = gpd.read_file(admin_path, layer=layer_name)
if layer_name == "ADM0":  
    # 🔧 Ensure one feature per country
    admin_areas = admin_areas.dissolve(by="shapeName", as_index=False)
    area_unique_id_col = "shapeName"
else:
    area_unique_id_col = "shapeID"
    admin_areas = admin_areas[[area_unique_id_col, "shapeName", "geometry"]].copy()
logging.info(f"There are {len(admin_areas)} admin areas to analyze.")

logging.info("Looping over admin regions and calculating population stats")
results = [] # List for collecting results
 # Loop over each admin region
for _, region in tqdm(admin_areas.iterrows(), total=len(admin_areas)):
    geom = region["geometry"].__geo_interface__

    mask_array = geometry_mask(
        [geom],
        transform=affine,
        invert=True,          # True inside polygon
        out_shape=out_shape
    )

    pop_clip = np.where(mask_array, pop_arr, np.nan)
    quintile_clip = np.where(mask_array, q_arr, np.nan)

    mask = (~np.isnan(pop_clip)) & (~np.isnan(quintile_clip))

    pop_values = pop_clip[mask]
    quintile_values = quintile_clip[mask]

    results.append({
        area_unique_id_col: region[area_unique_id_col],
        "shapeName": region["shapeName"],
        "total_pop": float(np.nansum(pop_values)),
        "q1_pop": float(np.nansum(pop_values[quintile_values == 1])),
        "q2_pop": float(np.nansum(pop_values[quintile_values == 2])),
        "q3_pop": float(np.nansum(pop_values[quintile_values == 3])),
        "q4_pop": float(np.nansum(pop_values[quintile_values == 4])),
        "q5_pop": float(np.nansum(pop_values[quintile_values == 5])),
    })

logging.info("Write to CSV")
results_df = pd.DataFrame(results)
results_df.to_csv(output_path, index=False)

logging.info("Done.")
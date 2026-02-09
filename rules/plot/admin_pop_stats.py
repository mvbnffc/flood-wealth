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

logging.info(f"Calculating population stats in {country}.")

logging.info("Load the raster data.")
pop = rasterio.open(pop_path)
quintiles = rasterio.open(quintile_path)
# Get raster info
profile = pop.meta.copy()
profile.update(dtype=rasterio.float32, compress='lzw', nodata=0, count=1)
affine = pop.transform
nodata = pop.nodata

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
for idx, region in tqdm(admin_areas.iterrows()):
    # Get the geometry for the current admin region
    geom = region["geometry"].__geo_interface__
    # Create a mask from the geometry
    mask_array = geometry_mask([geom],
                                transform=affine,
                                invert=True,
                                out_shape=pop.shape)
    # Use the mask to clip each raster by setting values outside the region to nan
    pop_clip = np.where(mask_array, pop, np.nan)
    quintile_clip = np.where(mask_array, quintiles, np.nan)
    # Mask out areas where not all rasters are valid
    mask = (
        ~np.isnan(pop_clip) &
        ~np.isnan(quintile_clip)
    )
    # Flatten data
    pop_values = pop_clip[mask]
    quintile_values = quintile_clip[mask]
    # Calculate pop per quintile
    q1_pop = pop_values[quintile_values==1].sum()
    q2_pop = pop_values[quintile_values==2].sum()
    q3_pop = pop_values[quintile_values==3].sum()
    q4_pop = pop_values[quintile_values==4].sum()
    q5_pop = pop_values[quintile_values==5].sum()

    results.append({
        area_unique_id_col: region[area_unique_id_col],
        "shapeName": region["shapeName"],
        "total_pop": pop_values.sum(),
        "q1_pop": q1_pop,
        "q2_pop": q2_pop,
        "q3_pop": q3_pop,
        "q4_pop": q4_pop,
        "q5_pop": q5_pop
    })

logging.info("Write to CSV")
results_df = pd.DataFrame(results)
results_df.to_csv(output_path, index=False)

logging.info("Done.")
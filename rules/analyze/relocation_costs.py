"""
This script calculates the cost of relocating people from high-risk urban areas to lower-risk areas.
Those elgible for relocation are those exposed to flooding greater than 1 m at a 10-year return period.
We use 10-year RP as this is the lowest RP that is common across datasets. Will only relocate people if the 
baseline level of flood protection is below 10-year RP.

Outputs are a relocation costs at the admin region.
"""

import logging
import sys
import glob
import os

import numpy as np
import rasterio
from tqdm import tqdm
from rasterio.features import geometry_mask, rasterize
from rasterio.mask import mask
import geopandas as gpd
from collections import Counter
from pyproj import Geod
from shapely.geometry import LineString, MultiLineString

if __name__ == "__main__":
    try:
        admin_path: str = snakemake.input["admin_areas"]
        flopros_path: str = snakemake.input["flopros_path"]
        flood_path: str = snakemake.input["flood_path"]
        pop_path: str = snakemake.input["pop_path"]
        res_area: str = snakemake.input["res_area"]
        res_capstock: str = snakemake.input["res_capstock"]
        urbanization_path: str = snakemake.input["urbanization_path"]
        output_path: str = snakemake.output["relocation_costs"]
        country: str = snakemake.wildcards["ISO3"]
        urban_class: str = snakemake.wildcards["urban_class"]
        model: str = snakemake.wildcards["model"]
        administrative_level: int = snakemake.wildcards.ADMIN_SLUG
    except:
        raise ValueError("Must be run via snakemake.")

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

# Update notation for GADM
admin_level = int(administrative_level.replace("ADM", ""))

logging.info(f"Calculating the relocation costs at Admin Level {admin_level} for {country} using {model} model. Reloaction threshold for Urbanization class: {urban_class}.")

logging.info("Creating output directories (if they doesn't already exist)")
out_dir = os.path.dirname(output_path)
os.makedirs(out_dir, exist_ok=True)

logging.info("Load files.")
with rasterio.open(flopros_path) as flopros_src, rasterio.open(flood_path) as flood_src, \
     rasterio.open(pop_path) as pop_src, rasterio.open(urbanization_path) as urban_src, \
     rasterio.open(res_area) as res_area_src, rasterio.open(res_capstock) as res_capstock_src:
    flopros = flopros_src.read(1)
    flood = flood_src.read(1)
    pop = pop_src.read(1)
    urban = urban_src.read(1)
    res_area = res_area_src.read(1)
    res_capstock = res_capstock_src.read(1)
    profile = flood_src.profile.copy()
    affine = flood_src.transform

logging.info("Building a mask for the cells to relocate.")
global_relocation_mask = (
    (urban <= int(urban_class)) &
    (flopros < 10) & # Only relocate if flood protection is below 10-year RP
    (flood >= 1) # Only relocate if flood depth is greater than (or equal to) 1 m
)

# Pre-compute global validity mask (areas where all rasters have valid data)
global_valid_mask = (
    ~np.isnan(flood) &
    ~np.isnan(pop) &
    ~np.isnan(urban) &
    ~np.isnan(res_capstock) &
    ~np.isnan(res_area)
)

logging.info(f"Reading level {administrative_level} admin boundaries")
layer_name = f"ADM{admin_level}"
admin_areas: gpd.GeoDataFrame = gpd.read_file(admin_path, layer=layer_name)
area_unique_id_col = "shapeName"
admin_areas = admin_areas[[area_unique_id_col, "geometry"]]
logging.info(f"There are {len(admin_areas)} admin areas to analyze.")

# OPTIMIZATION: vectorize geometry masking
# Create a dictionary mapping region index to geometry
geom_dict = {idx: geom for idx, geom in enumerate(admin_areas.geometry)}
# Create a single raster where each pixel contains the region ID it belongs to
region_ids = rasterize(
    [(geom, idx) for idx, geom in geom_dict.items()],
    out_shape=flood.shape,
    transform=affine,
    fill=-1,  # -1 for pixels not in any region
    dtype=np.int32
)

logging.info("Looping over admin regions and calculating relocation_costs")
results = [] # List for collecting results
 # Loop over each admin region
for idx, region in tqdm(admin_areas.iterrows()):
    # Create boolean mask for this specific region
    region_mask = (region_ids == idx) & global_valid_mask & global_relocation_mask
    
    if not np.any(region_mask):  # Skip if no valid pixels in this region
        results.append({
            area_unique_id_col: region[area_unique_id_col],
            "people_relocated": 0.0,
            "area_relocated": 0.0,
            "capstock_relocated": 0.0,
            "geometry": region["geometry"]
        })
        continue

    # Calculate sectoral capital stock losses for the region
    people_relocated = int(np.nansum(pop[region_mask]))
    area_relocated = int(np.nansum(res_area[region_mask]))
    capstock_relocated = int(np.nansum(res_capstock[region_mask]))
    
    # Append risk metrics to results list
    results.append({
         area_unique_id_col: region[area_unique_id_col],
         "people_relocated": people_relocated,
         "area_relocated": area_relocated,
         "capstock_relocated": capstock_relocated,
         "geometry": region["geometry"]
    })

logging.info("Writing reults to GeoPackage.")
results_gdf = gpd.GeoDataFrame(results, geometry="geometry")
results_gdf.to_file(output_path, driver="GPKG")

logging.info("Done.")
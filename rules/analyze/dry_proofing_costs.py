"""
This script calculates the cost of dry-proofing at the admin level.
Buildings elgible for dry proofing are those exposed to a 500-year flood extent and those outside
the 1 m flood depth extent for the 10-year RP flood.
We use 10-year and 500-year RP as these are the lowest and highest RP, respectively, across the GFM datasets. 

Approach is based off the following paper: Mortensen et al (2024) https://nhess.copernicus.org/articles/24/1381/2024/

Outputs are the dry-proofing costs at the admin region.
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
        rp10_path: str = snakemake.input["rp10_path"]
        rp500_path: str = snakemake.input["rp500_path"]
        res_path: str = snakemake.input["res_path"]
        cost_path: str = snakemake.input["res_unit_cost"]
        output_path: str = snakemake.output["dry_proofing_costs"]
        country: str = snakemake.wildcards["ISO3"]
        model: str = snakemake.wildcards["model"]
        administrative_level: int = snakemake.wildcards.ADMIN_SLUG
    except:
        raise ValueError("Must be run via snakemake.")

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

# Update notation for GADM
admin_level = int(administrative_level.replace("ADM", ""))

logging.info(f"Calculating the dry-proofing costs at admin Level {admin_level} for {country} using {model} model.")

logging.info("Creating output directories (if they don't already exist)")
out_dir = os.path.dirname(output_path)
os.makedirs(out_dir, exist_ok=True)

logging.info("Load files.")
with rasterio.open(rp10_path) as rp10_src, rasterio.open(rp500_path) as rp500_src, \
     rasterio.open(res_path) as res_area_src, rasterio.open(cost_path) as cost_src:
    rp10 = rp10_src.read(1)
    rp500 = rp500_src.read(1)
    res_area = res_area_src.read(1)
    cost = cost_src.read(1)
    profile = rp10_src.profile.copy()
    affine = rp10_src.transform

# Pre-compute global dry-proofing mask
global_dry_proofing_mask = (
    (rp10 < 1) & # Only dry_proof if exposure is outside the > 1 m RP10 flood extent
    (rp500 > 0) # Only dry_proof if building is within the 500-year RP flood extent
)

# Pre-compute global validity mask (areas where all rasters have valid data)
global_valid_mask = (
    ~np.isnan(rp10) &
    ~np.isnan(rp500) &
    ~np.isnan(res_area) &
    ~np.isnan(cost)
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
    out_shape=rp500.shape,
    transform=affine,
    fill=-1,  # -1 for pixels not in any region
    dtype=np.int32
)

logging.info("Looping over admin regions and calculating dry-proofing costs.")
results = [] # List for collecting results
 # Loop over each admin region
for idx, region in tqdm(admin_areas.iterrows()):
    # Create boolean mask for this specific region
    region_mask = (region_ids == idx) & global_valid_mask & global_dry_proofing_mask
    
    if not np.any(region_mask):  # Skip if no valid pixels in this region
        results.append({
            area_unique_id_col: region[area_unique_id_col],
            "area_dry-proofed": 0.0,
            "geometry": region["geometry"]
        })
        continue

    # Calculate sectoral capital stock losses for the region
    area_protected = int(np.nansum(res_area[region_mask]))
    # Calculate average unit cost for the region
    avg_unit_cost = np.nanmean(cost[region_mask])
    
    # Append risk metrics to results list
    results.append({
         area_unique_id_col: region[area_unique_id_col],
         "area_dry-proofed": area_protected,
         "average_unit_cost": avg_unit_cost,
         "geometry": region["geometry"]
    })

logging.info("Writing reults to GeoPackage.")
results_gdf = gpd.GeoDataFrame(results, geometry="geometry")
results_gdf.to_file(output_path, driver="GPKG")

logging.info("Done.")
"""
Use rasterio zonalstatistics to summarize BEM data within each admin unit.
"""

import logging
import sys
import glob
import os
import gc
from tqdm import tqdm
import rasterio
import exactextract as ee
import pandas as pd
import geopandas as gpd
import warnings
import shapely
from shapely.validation import make_valid  # Shapely ≥2.0

# Suppress GDAL warnings
warnings.filterwarnings('ignore', category=FutureWarning, module='osgeo')
warnings.filterwarnings('ignore', message='.*gdal.UseExceptions.*')

if __name__ == "__main__":
    try:
        adm2_path: str = snakemake.input["adm2"]
        adm1_path: str = snakemake.input["adm1"]
        adm0_path: str = snakemake.input["adm0"]
        res_path: str = snakemake.input["res_raster"]
        nres_path: str = snakemake.input["nres_raster"]
        adm2_output_path: str = snakemake.output["adm2"]
        adm1_output_path: str = snakemake.output["adm1"]
        adm0_output_path: str = snakemake.output["adm0"]
    except:
        raise ValueError("Must be run via snakemake.")

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

# ----------------------------- Common functions -----------------------------

def clean_admin_geoms(gdf):
    '''
    Function to fix invalid geometries
    '''
    # Fix invalids (fallback to buffer(0) if make_valid unavailable)
    if hasattr(shapely, "make_valid"):
        gdf["geometry"] = gdf["geometry"].apply(
            lambda g: shapely.make_valid(g) if g is not None else None
        )
    else:
        gdf["geometry"] = gdf["geometry"].buffer(0)

    # Drop Nones/empties/zero-area
    gdf = gdf[ gdf["geometry"].notna() ]
    gdf = gdf[ ~gdf.geometry.is_empty ]
    # If CRS is projected, also drop zero area; if geographic, area may be 0—skip this line.
    # gdf = gdf[gdf.geometry.area > 0]

    # Normalize to polygons only (drop GeometryCollections with no surface)
    gdf = gdf[gdf.geometry.geom_type.isin(["Polygon","MultiPolygon"])]

    # Split MultiPolygons; keep attributes
    gdf = gdf.explode(index_parts=False, ignore_index=True)

    # Add a stable ID to help debugging
    if "shapeID" not in gdf.columns:
        gdf["shapeID"] = gdf.index.astype("int64")

    return gdf

def summarize_bem(adm_path: str, bem_res_raster_path: str, bem_nres_raster_path: str, output_path: str, ADM_level: str):
    logging.info(f"Summarizing BEM data for ADM level {ADM_level}")
    gdf = gpd.read_file(adm_path)
    logging.info(f"Loaded {len(gdf)} polygons")

    # Open the rasters ONCE and reuse the handles in exact_extract
    with rasterio.open(bem_res_raster_path) as res_src, rasterio.open(bem_nres_raster_path) as nres_src:
        logging.info(f"Raster dimensions: {res_src.width} x {res_src.height}")
        logging.info(f"Raster CRS: {res_src.crs}")

        # Align CRS
        if gdf.crs != res_src.crs:
            logging.info("Reprojecting admin boundaries to match raster CRS.")
            gdf = gdf.to_crs(res_src.crs)

        # Clean geometries
        logging.info("Cleaning geometries.")
        gdf = clean_admin_geoms(gdf)
        logging.info(f"{len(gdf)} valid polygons after cleaning geometries.")

        chunk_size = 200  # can be larger now
        all_results = []
        logging.info(f"Processing {len(gdf)} polygons in chunks of {chunk_size}")

        for i in tqdm(range(0, len(gdf), chunk_size), desc=f"Processing {ADM_level}"):
            chunk = gdf.iloc[i:i + chunk_size].copy()

            try:
                # Pass the OPEN DATASET HANDLES, not the file paths
                res_stats  = ee.exact_extract(res_src,  chunk, ['sum'], include_geom=False)
                nres_stats = ee.exact_extract(nres_src, chunk, ['sum'], include_geom=False)

                res_df  = pd.DataFrame([f['properties'] for f in res_stats]).add_prefix('res_')
                nres_df = pd.DataFrame([f['properties'] for f in nres_stats]).add_prefix('nres_')

                chunk_results = pd.concat(
                    [chunk.reset_index(drop=True), res_df.reset_index(drop=True), nres_df.reset_index(drop=True)],
                    axis=1
                )
                all_results.append(chunk_results)

                del res_stats, nres_stats, res_df, nres_df, chunk_results, chunk
                gc.collect()

            except Exception as e:
                logging.error(f"Error processing chunk {(i // chunk_size) + 1}: {e}")
                logging.info("Attempting to free memory and continue...")
                try: del chunk
                except: pass
                gc.collect()
                continue

    if not all_results:
        logging.error("No chunks processed successfully!")
        return

    logging.info("Combining all chunks.")
    results_gdf = pd.concat(all_results, ignore_index=True)

    logging.info("Saving summarized data to CSV.")
    columns_to_keep = ['shapeName', 'res_sum', 'nres_sum']
    results_gdf.to_csv(output_path, index=False,
                       columns=[c for c in columns_to_keep if c in results_gdf.columns])



# -----------------------------------------------------------------------

# Run the analysis for each admin level
summarize_bem(adm2_path, res_path, nres_path, adm2_output_path, 'ADM2')
summarize_bem(adm1_path, res_path, nres_path, adm1_output_path, 'ADM1')
summarize_bem(adm0_path, res_path, nres_path, adm0_output_path, 'ADM0')

logging.info("Done.")

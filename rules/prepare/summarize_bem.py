"""
Use rasterio zonalstatistics to summarize BEM data within each admin unit.
"""

import logging
import sys
import glob
import os

import rasterio
import exactextract as ee
import pandas as pd
import geopandas as gpd

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
def summarize_bem(adm_path: str, bem_res_raster_path: str, bem_nres_raster_path: str, output_path: str, ADM_level: str):
    logging.info(f"Summarizing BEM data for ADM level {ADM_level}")
    logging.info("Loading the admin layer.")
    gdf = gpd.read_file(adm_path)
    
    with rasterio.open(bem_res_raster_path) as res_src:
        # Check that the shapefile and rasters have the same CRS
        if gdf.crs != res_src.crs:
            logging.info("Reprojecting admin boundaries to match raster CRS.")
            gdf = gdf.to_crs(res_src.crs)

    # Process in chunks to avoid memory issues
    chunk_size = 100  # Adjust this if still having memory issues
    all_results = []
    
    logging.info(f"Processing {len(gdf)} polygons in chunks of {chunk_size}")
    
    for i in range(0, len(gdf), chunk_size):
        chunk = gdf.iloc[i:i + chunk_size]
        logging.info(f"Processing chunk {i//chunk_size + 1}/{(len(gdf)-1)//chunk_size + 1}")
        
        logging.info("Calculating zonal statistics for residential BEM.")
        res_stats = ee.exact_extract(
            bem_res_raster_path,
            chunk, 
            ['sum'],
            include_geom=False
        )
        
        logging.info("Calculating zonal statistics for non-residential BEM.")
        nres_stats = ee.exact_extract(
            bem_nres_raster_path,
            chunk, 
            ['sum'],
            include_geom=False
        )
        
        # Convert stats to DataFrame
        res_df = pd.DataFrame(res_stats).add_prefix('res_')
        nres_df = pd.DataFrame(nres_stats).add_prefix('nres_')
        # Combined with chunk
        chunk_results = chunk.copy()
        chunk_results = pd.concat([chunk_results.reset_index(drop=True), res_df.reset_index(drop=True), nres_df.reset_index(drop=True)], axis=1)
        
        all_results.append(chunk_results)
    
    logging.info("Combining all chunks.")
    results_gdf = pd.concat(all_results, ignore_index=True)

    logging.info(f"Saving summarized data to CSV.")
    # Only keep relevant columns
    columns_to_keep = ['shapeName', 'res_sum', 'nres_sum']
    results_gdf[columns_to_keep].to_csv(output_path, index=False)

# -----------------------------------------------------------------------

# Run the analysis for each admin level
summarize_bem(adm2_path, res_path, nres_path, adm2_output_path, 'ADM2')
summarize_bem(adm1_path, res_path, nres_path, adm1_output_path, 'ADM1')
summarize_bem(adm0_path, res_path, nres_path, adm0_output_path, 'ADM0')

logging.info("Done.")

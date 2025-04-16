"""
This script calculates inequality metrics (concentration index and quantile ratio) 
at a given administrative level

The GADM administrative layers we are interested in:
    - National (level 0)
    - State/province/equivalent (level 1)
    - County/district/equivalent (level 2)
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
        rwi_path: str = snakemake.input["rwi_file"]
        pop_path: str = snakemake.input["pop_file"]
        risk_path: str = snakemake.input["risk_file"]
        output_path: str = snakemake.output["regional_CI"]
        administrative_level: int = snakemake.wildcards.ADMIN_SLUG
        model: str = snakemake.wildcards.MODEL
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

# Update notation for GADM
admin_level = int(administrative_level.replace("ADM-", ""))

logging.info(f"Calculating concentration indices at admin level {admin_level}.")

logging.info("Reading raster data.")
with rasterio.open(rwi_path) as rwi_src, rasterio.open(pop_path) as pop_src, rasterio.open(risk_path) as risk_src:
    rwi = rwi_src.read(1)
    pop = pop_src.read(1)
    risk = risk_src.read(1)
    affine = risk_src.transform 
rwi[rwi==-999] = np.nan # convert -999 in RWI dataset to NaN

logging.info(f"Reading level {administrative_level} admin boundaries")
layer_name = f"ADM_ADM_{admin_level}"
admin_areas: gpd.GeoDataFrame = gpd.read_file(admin_path, layer=layer_name)
area_unique_id_col = f"GID_{admin_level}"
# retain names of larger, encompassing adminstrative units
contextual_name_cols = [f"GID_{i}" for i in range(0, admin_level)]
admin_areas = admin_areas[[area_unique_id_col, *contextual_name_cols, "geometry"]]
logging.info(f"There are {len(admin_areas)} admin areas to analyze.")

logging.info("Looping over admin regions and calculating concentration indices")
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
    # Mask out areas where not all rasters are valid
    mask = (
        ~np.isnan(pop_clip) &
        ~np.isnan(rwi_clip) &
        ~np.isnan(risk_clip)
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

    # Going to calculate concentration indices for total, urban, and rural populations (un)protected
    df = pd.DataFrame({
        'pop': pop_flat,
        'rwi': rwi_flat,
        'flood': risk_flat,
    })

    # Define function
    def calculate_CI(df):
        # Sort dataframe by wealth
        df = df.sort_values(by="rwi", ascending=True)
        # Calculate cumulative population rank (to represent distribution of people)
        df['cum_pop'] = df['pop'].cumsum()
        # Calculate total pop of sample
        total_pop = df['pop'].sum()
        if total_pop == 0:
            return np.nan
        # Calculate fractional rank of each row
        df['rank'] = (df['cum_pop'] - 0.5*df['pop']) / total_pop
        try:
            # Calcualte weighted mean of flood risk
            weighted_mean_flood = np.average(df['flood'], weights=df['pop'])
        except ZeroDivisionError:
            return np.nan
        if weighted_mean_flood == 0:
            return np.nan
        # Calculate weighted sum of (flood * rank * pop)
        sum_xR = (df['flood'] * df['rank'] * df['pop']).sum()
        # Calculate Concentration Index
        CI = (2 * sum_xR) / (df['pop'].sum() * weighted_mean_flood) - 1
        return CI

    CI = calculate_CI(df)
    
    
    def calculate_quantile_ratio(df, quantile=0.2):
                # Sort the DataFrame by RWI (ascending)
                df_sorted = df.sort_values(by='rwi', ascending=True).copy()
            
                total_pop = df_sorted['pop'].sum()
                if total_pop == 0:
                    return np.nan
                
                # Calculate cumulative population
                df_sorted['cum_pop'] = df_sorted['pop'].cumsum()
                # Get bottom quantile: cells that add up to the first quantile share of the population
                bottom_df = df_sorted[df_sorted['cum_pop'] <= quantile * total_pop]
                # Get top quantile: cells that add up to the top quantile share of the population
                top_df = df_sorted[df_sorted['cum_pop'] >= (1 - quantile) * total_pop]
                
                try:
                    # Compute population-weighted average flood exposure
                    bottom_weighted_avg = np.average(bottom_df["flood"], weights=bottom_df["pop"])
                    top_weighted_avg = np.average(top_df["flood"], weights=top_df["pop"])
                except ZeroDivisionError:
                    return np.nan
                    
                # Return the quantile ratio. (Make sure you don’t divide by zero.)
                return top_weighted_avg / bottom_weighted_avg if bottom_weighted_avg != 0 else np.nan
    
    QR = calculate_quantile_ratio(df, quantile=0.2)

    # Calculate the number of cells where population and rwi overlaps
    total_pop = np.nansum(pop_clip)
    pop_rwi = np.nansum(np.where(~np.isnan(rwi_clip), pop_clip, 0))

    results.append({
        area_unique_id_col: region[area_unique_id_col],
        "CI": CI,
        "QR": QR,
        "Population": total_pop,
        "Population Coverage (%)": (pop_rwi/total_pop)*100,
        "rwi_count": np.count_nonzero(rwi_flat),
        "geometry": region["geometry"]
    })

logging.info("Writing reults to GeoPackage.")
results_gdf = gpd.GeoDataFrame(results, geometry="geometry")
results_gdf.to_file(output_path, driver="GPKG")

logging.info("Done.")
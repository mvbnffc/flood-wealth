"""
This script calculates inequality metrics (concentration index and quantile ratio)
and flood risk metrics at a given administrative level.

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
        social_path: str = snakemake.input["social_file"]
        pop_path: str = snakemake.input["pop_file"]
        mask_path: str = snakemake.input["mask_file"]
        risk_path: str = snakemake.input["risk_file"]
        output_path: str = snakemake.output["regional_CI"]
        administrative_level: int = snakemake.wildcards.ADMIN_SLUG
        model: str = snakemake.wildcards.MODEL
        social_name: str = snakemake.wildcards.SOCIAL
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

# Update notation for GADM
admin_level = int(administrative_level.replace("ADM-", ""))

logging.info(f"Calculating concentration indices at admin level {admin_level}.")

logging.info("Reading raster data.")
with rasterio.open(social_path) as social_src, rasterio.open(pop_path) as pop_src, \
     rasterio.open(mask_path) as mask_src, rasterio.open(risk_path) as risk_src:
    social = social_src.read(1)
    pop = pop_src.read(1)
    water_mask = mask_src.read(1)
    risk = risk_src.read(1)
    affine = risk_src.transform 

if social_name == "rwi":
    social[social==-999] = np.nan # convert -999 in RWI dataset to NaN
# Create the water mask
water_mask = np.where(water_mask>50, np.nan, 1) # WARNING WE ARE HARD CODING PERM_WATER > 50% mask here

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
                                out_shape=social.shape)
    # Use the mask to clip each raster by setting values outside the region to nan
    social_clip = np.where(mask_array, social, np.nan)
    pop_clip = np.where(mask_array, pop, np.nan)
    risk_clip = np.where(mask_array, risk, np.nan)
    water_mask_clip = np.where(mask_array, water_mask, np.nan)
    
    # Mask out areas where not all rasters are valid
    mask = (
        ~np.isnan(pop_clip) &
        ~np.isnan(social_clip) &
        ~np.isnan(risk_clip) &
        ~np.isnan(water_mask_clip)
    )
    # Flatten data
    pop_flat = pop_clip[mask]
    social_flat = social_clip[mask]
    risk_flat = risk_clip[mask]
    # Mask out zero-populatoin cells
    valid = pop_flat > 0
    pop_flat = pop_flat[valid]
    social_flat = social_flat[valid]
    risk_flat = risk_flat[valid]

    # Calculate total flood risk (pop * risk) for the region
    total_flood_risk = np.nansum(pop_flat * risk_flat)

    # Prepare dataframe for metric calculation
    df = pd.DataFrame({
        'pop': pop_flat,
        'social': social_flat,
        'flood': risk_flat,
    })

    # Define function
    def calculate_CI(df):
        # Sort dataframe by wealth
        df = df.sort_values(by="social", ascending=True)
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
                df_sorted = df.sort_values(by='social', ascending=True).copy()
            
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

    def calculate_flood_risk_per_quantile(df, quantile=0.2):
        """
        Calculate flood risk (pop * risk) for all quantiles
        NOTE: we have hardcoded quintiles here, but this can be adjusted
        """
        # Sort the DataFrame by social indicator (ascending)
        df_sorted = df.sort_values(by='social', ascending=True).copy()
        
        total_pop = df_sorted['pop'].sum()
        if total_pop == 0:
            return np.nan, np.nan
        
        # Calculate cumulative population
        df_sorted['cum_pop'] = df_sorted['pop'].cumsum()
        
        # Get first quantile (cells that add up to the first quantile share of the population)
        q1_df = df_sorted[df_sorted['cum_pop'] <= quantile * total_pop]
        # Get second quantile (cells that add up to the second quantile share of the population)
        q2_df = df_sorted[(df_sorted['cum_pop'] > quantile * total_pop) & (df_sorted['cum_pop'] <= 2 * quantile * total_pop)]
        # Get third quantile (cells that add up to the third quantile share of the population)
        q3_df = df_sorted[(df_sorted['cum_pop'] > 2 * quantile * total_pop) & (df_sorted['cum_pop'] <= 3 * quantile * total_pop)]
        # Get fourth quantile (cells that add up to the fourth quantile share of the population)
        q4_df = df_sorted[(df_sorted['cum_pop'] > 3 * quantile * total_pop) & (df_sorted['cum_pop'] <= 4 * quantile * total_pop)]
        # Get top quantile: cells that add up to the top quantile share of the population
        q5_df = df_sorted[df_sorted['cum_pop'] >= (1 - quantile) * total_pop]
        
        # Calculate flood risk (pop * risk) for each quantile
        q1_flood_risk = np.sum(q1_df['pop'] * q1_df['flood'])
        q2_flood_risk = np.sum(q2_df['pop'] * q2_df['flood'])
        q3_flood_risk = np.sum(q3_df['pop'] * q3_df['flood'])
        q4_flood_risk = np.sum(q4_df['pop'] * q4_df['flood'])
        q5_flood_risk = np.sum(q5_df['pop'] * q5_df['flood'])
        
        return q1_flood_risk, q2_flood_risk, q3_flood_risk, q4_flood_risk, q5_flood_risk
    
    Q1_risk, Q2_risk, Q3_risk, Q4_risk, Q5_risk = calculate_flood_risk_per_quantile(df, quantile=0.2)

    # Calculate the number of cells where population and rwi overlaps
    total_pop = np.nansum(pop_clip)
    pop_social = np.nansum(np.where(~np.isnan(social_clip), pop_clip, 0))

    results.append({
        area_unique_id_col: region[area_unique_id_col],
        "CI": CI,
        "QR": QR,
        "Population": total_pop,
        "Population Coverage (%)": (pop_social/total_pop)*100,
        "Total Flood Risk": total_flood_risk,
        "Q1 Flood Risk": Q1_risk,
        "Q2 Flood Risk": Q2_risk,
        "Q3 Flood Risk": Q3_risk,
        "Q4 Flood Risk": Q4_risk,
        "Q5 Flood Risk": Q5_risk,
        "geometry": region["geometry"]
    })

logging.info("Writing reults to GeoPackage.")
results_gdf = gpd.GeoDataFrame(results, geometry="geometry")
results_gdf.to_file(output_path, driver="GPKG")

logging.info("Done.")
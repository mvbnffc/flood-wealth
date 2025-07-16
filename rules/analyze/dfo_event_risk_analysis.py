"""
Runs a flood risk assessment for a given DFO flood event.
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
        rwi_path_list = snakemake.input["rwi_file"]
        pop_path_list = snakemake.input["pop_file"]
        mask_path_list = snakemake.input["mask_file"]
        flood_path_list = snakemake.input["flood_file"]
        output_path: str = snakemake.output["results"]
        event_id: str = snakemake.wildcards.event_id
        iso3_list = snakemake.params.iso3_list
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

# # Update notation for GADM
# admin_level = int(administrative_level.replace("ADM-", ""))

logging.info(f"Calculating risk metrics for DFO event {event_id} for {len(iso3_list)} countries: {iso3_list}.")

results = [] # List for collecting results
logging.info("Analyzing one country at a time.")
for iso3 in iso3_list:
    logging.info(f"Working on {iso3}.")

    # Get the paths for the current country
    rwi_path = [path for path in rwi_path_list if iso3 in path][0]
    pop_path = [path for path in pop_path_list if iso3 in path][0]
    mask_path = [path for path in mask_path_list if iso3 in path][0]
    flood_path = [path for path in flood_path_list if iso3 in path][0]

    logging.info("Reading raster data.")
    with rasterio.open(rwi_path) as rwi_src, rasterio.open(pop_path) as pop_src, \
          rasterio.open(mask_path) as mask_src, rasterio.open(flood_path) as flood_src:
        rwi = rwi_src.read(1)
        pop = pop_src.read(1)
        mask = mask_src.read(1)
        flood = flood_src.read(1)
        affine = flood_src.transform 
    
    rwi[rwi==-999] = np.nan # convert -999 in RWI dataset to NaN    

    water_mask = np.where(mask>50, np.nan, 1) # WARNING WE ARE HARD CODING PERM_WATER > 50% mask here
    # Remove all areas where one of the datasets is NaN
    mask_areas = (
        ~np.isnan(pop) &
        ~np.isnan(rwi) &
        ~np.isnan(flood) &
        ~np.isnan(water_mask)
    )
    # Flatten data
    pop_flat = pop[mask_areas]
    rwi_flat = rwi[mask_areas]
    flood_flat = flood[mask_areas]
    # Mask out zero-populatoin cells
    valid = pop_flat > 0
    pop_flat = pop_flat[valid]
    rwi_flat = rwi_flat[valid]
    flood_flat = flood_flat[valid]

    # Create dataframe for analysis
    df = pd.DataFrame({
        'pop': pop_flat,
        'rwi': rwi_flat,
        'flood': flood_flat,
    })

    def calculate_flood_exposure(df):
        # Calculate the sum of people exposed to flooding (flood > 0)
        exposed_population = df.loc[df['flood'] > 0, 'pop'].sum()
        return exposed_population
    
    FE = calculate_flood_exposure(df)

    def calculate_bottom_quintile_exposure(df):
        # Sort the DataFrame by RWI (ascending order)
        df_sorted = df.sort_values(by='rwi', ascending=True).copy()
        
        # Compute total population
        total_pop = df_sorted['pop'].sum()
        if total_pop == 0:
            return np.nan

        # Compute cumulative population
        df_sorted['cum_pop'] = df_sorted['pop'].cumsum()
        
        # Identify cells that make up the bottom 20% of the population
        bottom_df = df_sorted[df_sorted['cum_pop'] <= 0.2 * total_pop]
        
        # Sum the population in these cells, but only for cells where flood is > 0
        bottom_exposed = bottom_df.loc[bottom_df['flood'] > 0, 'pop'].sum()
        
        return bottom_exposed
    
    BQFE = calculate_bottom_quintile_exposure(df)

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

    def calculate_flood_exposure_per_quantile(df, quantile=0.2):
        """
        Calculate flood exposure for all quantiles
        NOTE: we have hardcoded quintiles here, but this can be adjusted
        """
        # Sort the DataFrame by social indicator (ascending)
        df_sorted = df.sort_values(by='social', ascending=True).copy()
        
        total_pop = df_sorted['pop'].sum()
        if total_pop == 0:
            logging.warning("Total pop is ZERO - returning NaN risk.")
            return np.nan, np.nan, np.nan, np.nan, np.nan
        
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
        q1_flood_exp = np.sum(q1_df['pop'] * q1_df['flood'])
        q2_flood_exp = np.sum(q2_df['pop'] * q2_df['flood'])
        q3_flood_exp = np.sum(q3_df['pop'] * q3_df['flood'])
        q4_flood_exp = np.sum(q4_df['pop'] * q4_df['flood'])
        q5_flood_exp = np.sum(q5_df['pop'] * q5_df['flood'])
        
        return q1_flood_exp, q2_flood_exp, q3_flood_exp, q4_flood_exp, q5_flood_exp
    
    Q1_exp, Q2_exp, Q3_exp, Q4_exp, Q5_exp = calculate_flood_exposure_per_quantile(df, quantile=0.2)

    results.append({
        "ISO3": iso3,
        "FE": FE,
        "BQFE": BQFE,
        "CI": CI,
        "QR": QR,
        "Q1_exp": Q1_exp,
        "Q2_exp": Q2_exp,
        "Q3_exp": Q3_exp,
        "Q4_exp": Q4_exp,
        "Q5_exp": Q5_exp
    })


logging.info("Writing reults to GeoPackage.")
results_df = pd.DataFrame(results)
results_df.to_csv(output_path)

logging.info("Done.")
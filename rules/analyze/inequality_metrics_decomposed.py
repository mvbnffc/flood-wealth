"""
This script calculates inequality metrics (concentration index and quantile ratio)
and flood risk metrics at a given administrative level. It also decomposes the metric by urbanization.

Note: now using geoboundaries rather than GADM for admin boundaries.
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
        urban_path: str = snakemake.input["urban_file"]
        risk_path: str = snakemake.input["risk_file"]
        output_path: str = snakemake.output["regional_CI"]
        administrative_level: int = snakemake.wildcards.ADMIN_SLUG
        model: str = snakemake.wildcards.MODEL
        social_name: str = snakemake.wildcards.SOCIAL
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

# Update notation for GADM
admin_level = int(administrative_level.replace("ADM", ""))

logging.info(f"Calculating concentration indices at admin level {admin_level}.")

logging.info("Reading raster data.")
with rasterio.open(social_path) as social_src, rasterio.open(pop_path) as pop_src, \
     rasterio.open(mask_path) as mask_src, rasterio.open(urban_path) as urban_src, \
    rasterio.open(risk_path) as risk_src:
    social = social_src.read(1)
    pop = pop_src.read(1)
    water_mask = mask_src.read(1)
    urban = urban_src.read(1).astype('float32') # convert to float to avoid errors
    risk = risk_src.read(1)
    affine = risk_src.transform 

if social_name == "rwi":
    social[social==-999] = np.nan # convert -999 in RWI dataset to NaN
urban[urban==10] = 0 # convert 10 in urban dataset (water class) to zero
urban[urban==-200] = 0 # convert -200 in urban dataset (no data) to NaN
# Create the water mask
water_mask = np.where(water_mask>50, np.nan, 1) # WARNING WE ARE HARD CODING PERM_WATER > 50% mask here

logging.info(f"Reading level {administrative_level} admin boundaries")
layer_name = f"ADM{admin_level}"
admin_areas: gpd.GeoDataFrame = gpd.read_file(admin_path, layer=layer_name)
if layer_name == "ADM0":
    area_unique_id_col = "shapeName"
else:
    area_unique_id_col = "shapeID"
admin_areas = admin_areas[[area_unique_id_col, "geometry"]]
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
    urban_clip = np.where(mask_array, urban, np.nan)
    risk_clip = np.where(mask_array, risk, np.nan)
    water_mask_clip = np.where(mask_array, water_mask, np.nan)
    
    # Mask out areas where not all rasters are valid
    mask = (
        ~np.isnan(pop_clip) &
        ~np.isnan(social_clip) &
        ~np.isnan(risk_clip) &
        ~np.isnan(urban_clip) &
        ~np.isnan(water_mask_clip)
    )
    # Flatten data
    pop_flat = pop_clip[mask]
    social_flat = social_clip[mask]
    urban_flat = urban_clip[mask]
    risk_flat = risk_clip[mask]
    # Mask out zero-populatoin cells
    valid = pop_flat > 0
    pop_flat = pop_flat[valid]
    urban_flat = urban_flat[valid]
    social_flat = social_flat[valid]
    risk_flat = risk_flat[valid]

    # Calculate total flood risk (pop * risk) for the region
    total_flood_risk = np.nansum(pop_flat * risk_flat)

    # Prepare dataframe for metric calculation
    df = pd.DataFrame({
        'pop': pop_flat,
        'social': social_flat,
        'flood': risk_flat,
        'urban': urban_flat
    })

    # Define function
    def calculate_CI(df):
        '''
        Calculates the concentration index (CI) for the given dataframe.
        Also Uses Clarke (2002) component decomposition of the concentration index.
        This is essentially the weighted average of the urban class concentration indices.
        Returns the overall CI and a dictionary with a tuple containing the (1) the CI contribution
        of each urban class, (2) the CI of the urban class (3) flood risk share of each urban class.
        (1, 2, 3)
        '''
        # Sort dataframe by wealth
        df = df.sort_values(by="social", ascending=True)
        # Calculate cumulative population rank (to represent distribution of people)
        df['cum_pop'] = df['pop'].cumsum()
        # # Calculate total pop of sample
        total_pop = df['pop'].sum()
        if total_pop == 0:
            return np.nan, {}
        # Calculate fractional rank of each row
        df['rank'] = (df['cum_pop'] - 0.5*df['pop']) / total_pop
        try:
            # Calcualte weighted mean of flood risk
            weighted_mean_flood = np.average(df['flood'], weights=df['pop'])
        except ZeroDivisionError:
            return np.nan, {}
        if weighted_mean_flood == 0:
            return np.nan, {}
        # Calculate weighted sum of (flood * rank * pop)
        sum_xR = (df['flood'] * df['rank'] * df['pop']).sum()
        # Calculate Concentration Index
        CI = (2 * sum_xR) / (df['pop'].sum() * weighted_mean_flood) - 1

        # Calculate contributions of each urban class
        contrib = {}
        for g, sub in df.groupby("urban"):
            if sub.empty:
                continue
            sub_pop = sub['pop'].sum()
            sub_weighted_mean_flood = np.average(sub['flood'], weights=sub['pop'])
            sub_sum_xR = (sub['flood'] * sub['rank'] * sub['pop']).sum()
            sub_CI = (2 * sub_sum_xR) / (sub_pop * sub_weighted_mean_flood) - 1
            sub_share = (sub_pop * sub_weighted_mean_flood) / (total_pop * weighted_mean_flood)
            contrib[g] = (sub_CI * sub_share, sub_CI, sub_share)

        # Ensure all urban classes are represented in the output
        valid_codes = [11, 12, 13, 21, 22, 23, 30]
        for c in valid_codes:
            if c not in contrib:
                contrib[c] = (np.nan, np.nan, np.nan)   # (net, CI_k, share)

        return CI, contrib


    def print_ci(ci_total, contrib_dict):
        df = (pd.DataFrame.from_dict(contrib_dict, orient="index",
                        columns=["Contribution", "CI", "Risk Share"])
                .rename_axis("Urban class")
                .sort_index())

        # Share of the overall CI (helpful when CI ≠ 1)
        df["Risk Share %"] = 100 * df["Risk Share"]
        df.drop(columns="Risk Share", inplace=True)

        # Share of the overall CI
        df["% of total CI"] = np.where(
            ci_total == 0, np.nan, 100 * df["Contribution"] / ci_total
        )

        # ------------------------------------------------------------------
        # Console output
        # ------------------------------------------------------------------
        print(f"\n\033[1mOverall concentration index (CI):\033[0m {ci_total:+.6f}\n")
        print(df.to_markdown(floatfmt=".4f"))

    CI, contrib = calculate_CI(df)
    print_ci(CI, contrib)

    # Calculate the number of cells where population and rwi overlaps
    total_pop = np.nansum(pop_clip)
    pop_social = np.nansum(np.where(~np.isnan(social_clip), pop_clip, 0))

    def safe(comp, code, idx):
        # Function to safely get values from the contribution dictionary
        return comp.get(code, (np.nan, np.nan, np.nan))[idx]

    results.append({
        area_unique_id_col: region[area_unique_id_col],
        "shapeName": region["shapeName"],
        "CI": CI,
        "Population": total_pop,
        "Population Coverage (%)": (pop_social/total_pop)*100,
        "Total Flood Risk": total_flood_risk,
        "DUC11 CI": safe(contrib, 11, 1),
        "DUC12 CI": safe(contrib, 12, 1),
        "DUC13 CI": safe(contrib, 13, 1),
        "DUC21 CI": safe(contrib, 21, 1),
        "DUC22 CI": safe(contrib, 22, 1),
        "DUC23 CI": safe(contrib, 23, 1),
        "DUC30 CI": safe(contrib, 30, 1),
        "DUC11 Contribution": safe(contrib, 11, 0),
        "DUC12 Contribution": safe(contrib, 12, 0),
        "DUC13 Contribution": safe(contrib, 13, 0),
        "DUC21 Contribution": safe(contrib, 21, 0),
        "DUC22 Contribution": safe(contrib, 22, 0),
        "DUC23 Contribution": safe(contrib, 23, 0),
        "DUC30 Contribution": safe(contrib, 30, 0),
        "DUC11 Risk Share": safe(contrib, 11, 2),
        "DUC12 Risk Share": safe(contrib, 12, 2),
        "DUC13 Risk Share": safe(contrib, 13, 2),
        "DUC21 Risk Share": safe(contrib, 21, 2),
        "DUC22 Risk Share": safe(contrib, 22, 2),
        "DUC23 Risk Share": safe(contrib, 23, 2),
        "DUC30 Risk Share": safe(contrib, 30, 2),
        "geometry": region["geometry"]
    })

logging.info("Writing reults to GeoPackage.")
results_gdf = gpd.GeoDataFrame(results, geometry="geometry")
results_gdf.to_file(output_path, driver="GPKG")

logging.info("Done.")
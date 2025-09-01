"""
This script calculates inequality metrics (concentration index and quantile ratio)
and flood risk metrics at a given administrative level. It also decomposes the metric by admin regions.

Note: now using geoboundaries rather than GADM for admin boundaries.
"""

import logging

import rasterio
from rasterio.features import geometry_mask, rasterize
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

logging.info(f"Calculating National CI and spatial decomoposition at admin level {admin_level}.")

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
    out_shape = risk.shape

if social_name == "rwi":
    social[social==-999] = np.nan # convert -999 in RWI dataset to NaN
urban[urban==10] = np.nan # convert 10 in urban dataset (water class) to NaN
urban[urban==-200] = np.nan # convert -200 in urban dataset (no data) to NaN
# Create the water mask
water_mask = np.where(water_mask>50, np.nan, 1) # WARNING WE ARE HARD CODING PERM_WATER > 50% mask here

logging.info(f"Reading level {administrative_level} admin boundaries.")
layer_name = f"ADM{admin_level}"
admin_areas: gpd.GeoDataFrame = gpd.read_file(admin_path, layer=layer_name)
area_unique_id_col = "shapeName"
admin_areas = admin_areas[[area_unique_id_col, "geometry"]].reset_index(drop=True)

logging.info(f"Rasterizing {len(admin_areas)} admin areas to a grid...")
# Map each admin feature to an integer ID (0 = no admin)
admin_areas['region_id'] = np.arange(1, len(admin_areas)+1, dtype=np.int32)
shapes = [(geom, rid) for geom, rid in zip(admin_areas.geometry, admin_areas['region_id'])]
region_id_arr = rasterize(
    shapes=shapes,
    out_shape=out_shape,
    transform=affine,
    fill=0,
    dtype=np.int32,
)

logging.info("Build nationa dataframe for metric calculation")
valid_mask = (
    (pop > 0) &
    np.isfinite(pop) &
    np.isfinite(social) &
    np.isfinite(risk) &
    np.isfinite(urban) &
    np.isfinite(water_mask) &
    (region_id_arr > 0)  # Only consider areas that fall within an admin region
)

pop_flat = pop[valid_mask]
social_flat = social[valid_mask]
urban_flat = urban[valid_mask].astype(np.int32)  # Ensure urban codes are integers
risk_flat = risk[valid_mask]
region_id_flat = region_id_arr[valid_mask].astype(np.int32)

df_national = pd.DataFrame({
    'pop': pop_flat,
    'social': social_flat,
    'flood': risk_flat,
    'urban': urban_flat,
    'region_id': region_id_flat
})


# Define some helper functions
def calculate_decomposed_CI(df, group_col=None, expected_groups=None):
    """
    Compute national CI using national ranks, and optionally decompose into subgroups
    using Clarke-style decomposition. If group_col is provided, subgroup CI uses the
    *national* rank variable (not within-group rank).

    Returns:
        ci_total (float),
        contrib_dict: {group_value: (net_contribution, subgroup_CI, risk_share)}
    """
    if df.empty or df['pop'].sum() == 0:
        return np.nan, {}
    
    # Sort dataframe by wealth
    df = df.sort_values(by="social", ascending=True).copy()
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
    # Calculate NATIONAL Concentration Index
    CI = (2 * sum_xR) / (df['pop'].sum() * weighted_mean_flood) - 1

    # Calculate contributions of each group
    contrib = {}
    if group_col is not None:
        for g, sub in df.groupby(group_col, sort=False):
            if sub.empty:
                continue
            sub_pop = sub['pop'].sum()
            sub_weighted_mean_flood = np.average(sub['flood'], weights=sub['pop'])
            sub_sum_xR = (sub['flood'] * sub['rank'] * sub['pop']).sum()
            sub_CI = (2 * sub_sum_xR) / (sub_pop * sub_weighted_mean_flood) - 1
            sub_share = (sub_pop * sub_weighted_mean_flood) / (total_pop * weighted_mean_flood)
            contrib[int(g)] = (sub_CI * sub_share, sub_CI, sub_share)

        # Ensure all urban classes are represented in the output
        if expected_groups is not None:
            for g in expected_groups:
                if int(g) not in contrib:
                    contrib[int(g)] = (np.nan, np.nan, np.nan)

    return CI, contrib

def region_standalone_ci(df_region):
    """
    CI computed *within* one region only (ranks recomputed inside region).
    """
    if df_region.empty or df_region["pop"].sum() == 0:
        return np.nan
    d = df_region.sort_values(by="social", ascending=True).copy()
    d["cum_pop"] = d["pop"].cumsum()
    total_pop = d["pop"].sum()
    d["rank"] = (d["cum_pop"] - 0.5 * d["pop"]) / total_pop
    wm = np.average(d["flood"], weights=d["pop"])
    if wm == 0:
        return np.nan
    sum_xR = (d["flood"] * d["rank"] * d["pop"]).sum()
    return (2 * sum_xR) / (d["pop"].sum() * wm) - 1

logging.info("Calculating national CI and decomposition...")
all_region_ids = sorted(admin_areas['region_id'].tolist())
CI, contrib_region = calculate_decomposed_CI(
    df_national, group_col="region_id", expected_groups=all_region_ids
)

logging.info("Computing standalone CI for each region...")
standalone_region_ci = (
    df_national.groupby("region_id")
    .apply(region_standalone_ci)
    .rename("CI_region_only")
    .to_dict()
)

def safe_get(d, key, idx):
    """d[key] -> (contrib, sub_ci, share); returns np.nan if missing."""
    try:
        return d.get(int(key), (np.nan, np.nan, np.nan))[idx]
    except Exception:
        return np.nan
    
logging.info("Combining regional metrics and writing output..." )
admin_areas["Nat_Contrib"] = admin_areas["region_id"].apply(lambda x: safe_get(contrib_region, x, 0))
admin_areas["Nat_SubCI"] = admin_areas["region_id"].apply(lambda x: safe_get(contrib_region, x, 1))
admin_areas["Nat_RiskShare"] = admin_areas["region_id"].apply(lambda x: safe_get(contrib_region, x, 2))
admin_areas["CI_region_only"] = admin_areas["region_id"].apply(lambda x: standalone_region_ci.get(x, np.nan))   
# Add a % of total CI column
admin_areas["% of total CI"] = np.where(
    np.isfinite(CI) & (CI != 0),
    100 * admin_areas["Nat_Contrib"] / CI,
    np.nan
)
# Add national CI back to column (for reference)
admin_areas["Nat_CI"] = CI

logging.info("Writing results to GeoPackage.")
admin_areas.drop(columns=["region_id"], inplace=True)
admin_areas.to_file(output_path, driver="GPKG")

logging.info(f"Overall NATIONAL CI: {CI:+.4f}")
logging.info(f"Sum of regional contributions: {np.nansum([v[0] for v in contrib_region.values()]):+.4f}")
logging.info("Done.")
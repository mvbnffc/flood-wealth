"""
This script prepares the GHS-DUC data for further analysis by mergine the urbanisation classfication CSV file
with the country's GADM file. The GADM level to use is extracted from the info spreadsheet.
"""

import logging
import sys
import glob
import os

import pandas as pd
import geopandas as gpd

if __name__ == "__main__":
    try:
        info_path: str = snakemake.input["duc_info"]
        adm_path: str = snakemake.input["adm_file"]
        output_path: str = snakemake.output["urban_file"]
        country: str = snakemake.wildcards["ISO3"]
    except:
        raise ValueError("Must be run via snakemake.")

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Preparing the GHS-DUC dataset by merging with GADM for {country}.")

logging.info("Creating output directory (if it doesn't already exist)")
out_dir = os.path.dirname(output_path)
os.makedirs(out_dir, exist_ok=True)

logging.info("Loading the Excel sheet into a dataframe.")
# Name of the sheet we are interested in 
year = 2020 # WARNING hard coding this....
sheet_name = f"Country DEGURBA {year}"
info_df = pd.read_excel(info_path, sheet_name=sheet_name, engine='openpyxl')

logging.info("Extracting the GADM level to use from the info dataframe.")
matched = info_df.loc[info_df['GADM ISO'] == country, 'Selected GADM Level']
if matched.empty:
    raise KeyError(f"No GADM level found for {country} in sheet {sheet_name}.")
gadm_level = int(matched.iloc[0])
logging.info(f"GADM level to use for {country} is {gadm_level}.")

logging.info("Load relevant GHS-DUC CSV file based on GADM level.")
csv_path = os.path.join(os.path.dirname(info_path), f"GHS_DUC_GLOBE_R2023A_V2_0_GADM41_{year}_level{gadm_level}.csv")
duc_df = pd.read_csv(csv_path)
duc_df = duc_df[duc_df['GID_0GHSL'] == country]

logging.info("Loading the GADM file")
layer_name = f"ADM_ADM_{gadm_level}"
gadm_gdf = gpd.read_file(adm_path, layer=layer_name)

logging.info("Mergeing the GHS-DUC dataframe with the GADM GeoDataFrame.")
merged_gdf = gadm_gdf.merge(duc_df, on=f"GID_{gadm_level}", how="left")

logging.info("Saving the merged GeoDataFrame to a new file.")
merged_gdf.to_file(output_path, driver='GPKG')

logging.info("Done.")


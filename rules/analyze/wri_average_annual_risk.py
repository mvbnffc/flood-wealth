"""
Given relative risk maps for GIRI flood return periods - calculate the average annual relative risk
for each flooded grid cell.
"""

import logging

import rasterio
import numpy as np

if __name__ == "__main__":

    try:
        RP2_path: str = snakemake.input["flood_rp_2"]
        RP5_path: str = snakemake.input["flood_rp_5"]
        RP10_path: str = snakemake.input["flood_rp_10"]
        RP25_path: str = snakemake.input["flood_rp_25"]
        RP50_path: str = snakemake.input["flood_rp_50"]
        RP100_path: str = snakemake.input["flood_rp_100"]
        RP250_path: str = snakemake.input["flood_rp_250"]
        RP500_path: str = snakemake.input["flood_rp_500"]
        RP1000_path: str = snakemake.input["flood_rp_1000"]
        aar_output_path: str = snakemake.output["flood_aar"]
        vuln_dataset: str = snakemake.wildcards["VULN_CURVE"]
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Calculating WRI average annual relative risk using {vuln_dataset} vulnerability curve.")

logging.info("Reading raster data.")
raster_paths = [RP2_path, RP5_path, RP10_path, RP25_path, RP50_path, RP100_path, RP250_path, RP500_path, RP1000_path]
flood_maps = [] # going to load rasters into this list
for path in raster_paths:
    with rasterio.open(path) as src:
        flood_maps.append(src.read(1))
        transform = src.transform
        meta = src.meta

logging.info("Calculating annual average risk - unprotected")
RPs = np.array([2, 5, 10, 25, 50, 100, 200, 500, 1000]) # define return peridos
aep = 1 / RPs # convert to annual exceedance probability
flood_maps = np.array(flood_maps) # convert to numpy array for calculation
aar = np.trapz(flood_maps[::-1], x=aep[::-1], axis=0) # debug: need to reverse list to avoid negative AAR values

# Update for compression
meta.update(
    compress='lzw',
    tiled=True
)

logging.info("Writing output rasters.")
with rasterio.open(aar_output_path, "w", **meta) as dst:
    dst.write(aar.astype(np.float32), 1)

logging.info("Done.")
"""
Given relative risk maps for GIRI flood return periods - calculate the average annual relative risk
for each flooded grid cell. This is for the flood protection adaptation scenario.
"""

import logging

import rasterio
import numpy as np

if __name__ == "__main__":

    try:
        RP5_path: str = snakemake.input["flood_rp_5"]
        RP10_path: str = snakemake.input["flood_rp_10"]
        RP25_path: str = snakemake.input["flood_rp_25"]
        RP50_path: str = snakemake.input["flood_rp_50"]
        RP100_path: str = snakemake.input["flood_rp_100"]
        RP200_path: str = snakemake.input["flood_rp_200"]
        RP500_path: str = snakemake.input["flood_rp_500"]
        RP1000_path: str = snakemake.input["flood_rp_1000"]
        flopros_path: str = snakemake.input["flood_protection"]
        aar_output_path: str = snakemake.output["flood_aar_adapted"]
        rp_protection: str = snakemake.wildcards["RP"]
        urban_class: str = snakemake.wildcards["urban_class"]
        vuln_dataset: str = snakemake.wildcards["VULN_CURVE"]
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Calculating (adapted) GIRI average annual relative risk using {vuln_dataset} vulnerability curve." \
             "Flood protection adaptation scenario. Return period for flood protection: {rp_protection}, Urbanization class: {urban_class}.")

logging.info("Reading raster data.")
raster_paths = [RP5_path, RP10_path, RP25_path, RP50_path, RP100_path, RP200_path, RP500_path, RP1000_path]
flood_maps = [] # going to load rasters into this list
for path in raster_paths:
    with rasterio.open(path) as src:
        flood_maps.append(src.read(1))
        transform = src.transform
        meta = src.meta
with rasterio.open(flopros_path) as src:
    flopros = src.read(1)

# Create an empty array for bankful (2-year) flood
RP2 = np.zeros_like(flood_maps[0])
flood_maps.insert(0, RP2) # insert this array at beginning of list.
# Convert to numpy array
flood_maps = np.array(flood_maps)

### Function for Loss-Probabiltiy Curve
def integrate_truncated_risk(risk_curve, T, RPs):
    """
    Integrate the risk curve above a protection threshold T.
    
    Parameters:
      risk_curve : 1D numpy array of risk values corresponding to each return period in RPs.
      T          : Protection threshold for the cell (e.g., 7, 10, 100, etc.)
      RPs        : 1D numpy array of discrete return periods.
    
    Returns:
      Integrated risk value after truncating the risk curve at T.
    """
    # If the protection threshold is less than or equal to the smallest RP,
    # no truncation is applied.
    if T <= RPs[0]:
        new_RPs = RPs
        new_risk = risk_curve
    # If the protection threshold is above the highest RP,
    # assume full protection (i.e. no risk).
    elif T >= RPs[-1]:
        return 0.0
    else:
        # Find the first index where the discrete RP is >= T.
        idx = np.searchsorted(RPs, T, side='left')
        # If T exactly matches one of the RPs, use that risk value.
        if T == RPs[idx]:
            risk_T = risk_curve[idx]
        else:
            # Linearly interpolate between the two adjacent RPs.
            risk_T = risk_curve[idx-1] + (risk_curve[idx] - risk_curve[idx-1]) * (T - RPs[idx-1]) / (RPs[idx] - RPs[idx-1])
        # Construct a new risk curve starting at T.
        new_RPs = np.concatenate(([T], RPs[idx:]))
        new_risk = np.concatenate(([risk_T], risk_curve[idx:]))
    
    # Calculate the annual exceedance probabilities for the new return periods.
    new_aep = 1 / new_RPs
    # Because aep decreases with increasing RP, reverse arrays to get increasing x values.
    sorted_aep = new_aep[::-1]
    sorted_risk = new_risk[::-1]
    
    # Compute the integrated risk using the trapezoidal rule.
    return np.trapezoid(sorted_risk, x=sorted_aep)

logging.info("Reshaping flood maps")
RPs = np.array([2, 5, 10, 25, 50, 100, 200, 500, 1000]) # define return periods
rows, cols = flopros.shape # for writing back to normal shape later
n_rps = len(RPs)
risk_flat = flood_maps.reshape(n_rps, -1) # each column is now a cell's risk curve
flopros_flat = flopros.flatten()
aar_protected_flat = np.empty(flopros_flat.shape) # array to store integrated risk for each cell.

logging.info("Calculating average annual risk (adapted) using vectorization")
for i in range(flopros_flat.size):
    aar_protected_flat[i] = integrate_truncated_risk(risk_flat[:, i], flopros_flat[i], RPs)
aar_protected = aar_protected_flat.reshape(rows, cols) # reshape to original dimensions

# Update for compression
meta.update(
    compress='lzw',
    tiled=True
)

logging.info("Writing output rasters.")
with rasterio.open(aar_output_path, "w", **meta) as dst:
    dst.write(aar_protected.astype(np.float32), 1)

logging.info("Done.")
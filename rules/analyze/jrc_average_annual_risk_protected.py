"""
Given relative risk maps for JRC flood return periods - calculate the average annual relative risk
for each flooded grid cell.
"""

import logging
import sys

import numba
import numpy as np
import rasterio

if __name__ == "__main__":

    try:
        RP10_path: str = snakemake.input["flood_rp_10"]
        RP20_path: str = snakemake.input["flood_rp_20"]
        RP50_path: str = snakemake.input["flood_rp_50"]
        RP75_path: str = snakemake.input["flood_rp_75"]
        RP100_path: str = snakemake.input["flood_rp_100"]
        RP200_path: str = snakemake.input["flood_rp_200"]
        RP500_path: str = snakemake.input["flood_rp_500"]
        flopros_path: str = snakemake.input["flopros"]
        aar_output_path: str = snakemake.output["flood_aar_protected"]
        vuln_dataset: str = snakemake.wildcards["VULN_CURVE"]
    except NameError:
        try:
            RP10_path: str = sys.argv[1]
            RP20_path: str = sys.argv[2]
            RP50_path: str = sys.argv[3]
            RP75_path: str = sys.argv[4]
            RP100_path: str = sys.argv[5]
            RP200_path: str = sys.argv[6]
            RP500_path: str = sys.argv[7]
            flopros_path: str = sys.argv[8]
            aar_output_path: str = sys.argv[9]
            vuln_dataset: str = sys.argv[10]
        except IndexError:
            print("""\
Example usage: python rules/analyze/jrc_average_annual_risk_protected.py \\
    data/results/flood_risk/countries/CRI/CRI_jrc-flood-risk_RP10_V-JRC.tif \\
    data/results/flood_risk/countries/CRI/CRI_jrc-flood-risk_RP20_V-JRC.tif \\
    data/results/flood_risk/countries/CRI/CRI_jrc-flood-risk_RP50_V-JRC.tif \\
    data/results/flood_risk/countries/CRI/CRI_jrc-flood-risk_RP75_V-JRC.tif \\
    data/results/flood_risk/countries/CRI/CRI_jrc-flood-risk_RP100_V-JRC.tif \\
    data/results/flood_risk/countries/CRI/CRI_jrc-flood-risk_RP200_V-JRC.tif \\
    data/results/flood_risk/countries/CRI/CRI_jrc-flood-risk_RP500_V-JRC.tif \\
    data/inputs/analysis/countries/CRI/CRI_flopros.tif \\
    data/results/flood_risk/countries/CRI/CRI_jrc-flood-risk_protected_AAR_V-JRC.tif \\
    JRC
""")
            raise ValueError("Must be run via snakemake or with all parameters supplied.")
    
logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info(f"Calculating JRC average annual relative risk using {vuln_dataset} vulnerability curve.")

logging.info("Reading raster data.")
raster_paths = [RP500_path, RP200_path, RP100_path, RP75_path, RP50_path, RP20_path, RP10_path]
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
flood_maps.append(RP2) # insert this array at beginning of list.
# Convert to numpy array
flood_maps = np.array(flood_maps)

### Function for Loss-Probabiltiy Curve
@numba.jit
def integrate_truncated_risk(risk_curve, T, RPs, RPs_r, aep):
    """
    Integrate the risk curve above a protection threshold T.
    
    Parameters:
      risk_curve : 1D numpy array of risk values corresponding to each return period in RPs.
      T          : Protection threshold for the cell (e.g., 7, 10, 100, etc.)
      RPs        : 1D numpy array of discrete return periods.
    
    Returns:
      Integrated risk value after truncating the risk curve at T.
    """
    # If the protection threshold is above the highest RP,
    # assume full protection (i.e. no risk).
    if T >= RPs[0]:
        return 0.0
    # Or if the highest RP risk value is zero, assume no risk
    if risk_curve[0] <= 0:
        return 0.0
    
    # If the protection threshold is less than or equal to the smallest RP,
    # no truncation is applied.
    if T <= RPs[-1]:
        protected_risk = risk_curve
        protected_aep = aep
    else:
        # Find the first index where the discrete RP is >= T.
        ridx = np.searchsorted(RPs_r, T, side='left')
        idx = (RPs.size - ridx) 
        # If T exactly matches one of the RPs, use that risk value.
        if T == RPs[idx]:
            risk_T = risk_curve[idx]
        else:
            # Linearly interpolate between the two adjacent RPs.
            risk_T = risk_curve[idx-1] + (risk_curve[idx] - risk_curve[idx-1]) * (T - RPs[idx-1]) / (RPs[idx] - RPs[idx-1])
        # Construct a new risk curve starting at T.
        protected_risk = risk_curve[:idx+1].copy()
        protected_risk[idx] = risk_T
        # Calculate the annual exceedance probabilities for the new return periods.
        protected_aep = aep[:idx+1].copy()
        protected_aep[idx] = 1 / T

    # Compute the integrated risk using the trapezoidal rule.
    return np.trapezoid(protected_risk, x=protected_aep, dx=None)

@numba.guvectorize([(numba.float32[:,:], numba.float32[:], numba.float32[:], numba.float32[:])], '(m,n),(n),(m)->(n)')
def protected_risk(risk, protection, rps, out):
    rps_r = rps[::-1].copy()
    aep = np.pow(rps, -1)
    for i in range(protection.shape[0]):
        out[i] = integrate_truncated_risk(risk[:, i], protection[i], rps, rps_r, aep)

logging.info("Reshaping flood maps")
RPs = np.array([500, 200, 100, 75, 50, 20, 10, 2], dtype="float32") # define return periods
rows, cols = flopros.shape # for writing back to normal shape later
n_rps = len(RPs)
risk_flat = flood_maps.reshape(n_rps, -1) # each column is now a cell's risk curve
flopros_flat = flopros.flatten()
aar_protected_flat = np.empty(flopros_flat.shape, dtype="float32") # array to store integrated risk for each cell.

logging.info("Calculating average annual risk (protected) using vectorization")
protected_risk(risk_flat, flopros_flat, RPs, aar_protected_flat)
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
"""
Download Google Historical Inundation Data

https://console.cloud.google.com/storage/browser/flood-forecasting/inundation_history

The inundation models are ML-based models, trained based on past flood events.

The past events are Satellite-based flood inundation maps: The synthetic aperture radar
ground range detected (SAR GRD) data from the Sentinel-1 satellite constellation are 
used to determine flood inundation maps at known timepoints and locations . At any AOI 
(area of interest), a SAR image is available once every several days, from which an 
inundation map was inferred using a binary classifier (Torres, 2012). Every pixel within 
a SAR image is classified as wet/dry via a Gaussian mixture based classification algorithm. 
In order to calibrate and evaluate the classification algorithm, we have collected a dataset 
of Sentinel 2 multispectral images of flood events that coincide with the SAR image dates 
and locations. Reference Sentinel-2 flood maps were created by calculating per-pixel 
Normalized Difference Water Index (NDWI=(B3-B8)/(B3 + B8), B3 and B8 are green and near 
infrared bands, respectively) and applied a threshold of 0.

These flood images are the ground truth for the ML model, while the inputs are the gauge 
value, the latest available satellite image, and flood history of the region.

The inundation models are trained and validated based on historical flood events, where 
flood inundation extent maps from satellite data, along with the corresponding gauge water 
stage measurements, are available. Similar to the stage forecast models, a 1-year leave 
out cross validation scheme is used for training and validation.
"""

import os
import requests
from urllib.parse import urljoin
import logging

if __name__ == "__main__":
    try:
        raw_folder: str = snakemake.output["raw_folder"]
    except NameError:
        raise ValueError("Must be run via snakemake.")
    
# Google Cloud bucket info
bucket = "flood-forecasting"
prefix = "inundation_history"

logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

logging.info("Creating output directory (if it doesn't already exist)")
os.makedirs(raw_folder, exist_ok=True)

logging.info("Defining functions")
def list_all_files(bucket, prefix):
    """
    List all objects in the given bucket that start with the provided prefix.
    Returns a list of object names.
    """
    base_api_url = f"https://www.googleapis.com/storage/v1/b/{bucket}/o"
    params = {"prefix": prefix}
    all_files = []

    while True:
        print(f"Fetching: {base_api_url} with params {params}")
        response = requests.get(base_api_url, params=params)
        if response.status_code != 200:
            print(f"Error fetching listing: {response.status_code}")
            break

        data = response.json()
        items = data.get("items", [])
        for item in items:
            name = item.get("name")
            if name:
                all_files.append(name)

        # Check for a next page token; if present, update params for next request
        next_page_token = data.get("nextPageToken")
        if next_page_token:
            params["pageToken"] = next_page_token
        else:
            break

    return all_files

def download_file(bucket, file_name, local_dir):
    """
    Downloads a single file from the bucket and saves it to the local directory,
    preserving the file's folder structure.
    """
    # Skip directory markers
    if file_name.endswith("/"):
        print(f"Skipping directory marker: {file_name}")
        return
    
    base_download_url = f"https://storage.googleapis.com/{bucket}/"
    file_url = urljoin(base_download_url, file_name)
    local_path = os.path.join(local_dir, file_name)
    os.makedirs(os.path.dirname(local_path), exist_ok=True)

    print(f"Downloading: {file_name} -> {local_path}")
    response = requests.get(file_url)
    if response.status_code == 200:
        with open(local_path, "wb") as f:
            f.write(response.content)
    else:
        print(f"Failed to download {file_name}. Status code: {response.status_code}")

logging.info("Collecting bucket file information.")
all_files = list_all_files(bucket, prefix)
logging.info(f"Total files in bucket: {len(all_files)}")

logging.info("Downloading all files...")
for file_name in all_files:
    download_file(bucket, file_name, raw_folder)

logging.info("Done.")


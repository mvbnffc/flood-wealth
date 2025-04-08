"""
Prepare the GFD Flood data. This involves:
    - 1 preparing data for merge (burn files to global extent)
    - 2 merging into a single global flood frequency file
    - 3 clipping and resampling for country of interest
"""

rule prepare_gfd_merge:
    """
    This rule extracts the relevant GFD raster bands (flood occurance and permanent water)
    and converts NaN values to zero - which prepares the files for a global merge.
    It also filters the relevant GFD maps and removes any caused by dam outburts 
    and storm surge flooding (as we are focussed on inland)
    """
    input:
        raw_gfd_folder="data/inputs/gfd/raw/"
    output:
        merge_gfd_folder=directory("data/inputs/gfd/merged/prep/")
    script:
        "./prepare_gfd_merge.py"

rule merge_gfd:
    """
    This rule merges GFD rasters into one global file
    """
    input:
        merge_gfd_folder="data/inputs/gfd/merged/prep"
    output:
        merge_gfd_file="data/inputs/gfd/merged/gfd.tif"
    script:
        "./merge_gfd.py"
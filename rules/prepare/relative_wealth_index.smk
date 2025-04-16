"""
Prepare (clip and resample) the RWI data
"""

rule clip_and_resample_rwi:
    """
    Resample the and clip the rwi (using population data as reference)
    """
    input:
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
        rwi_file="data/inputs/rwi/rwi.tif",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
    output:
        rwi_resampled = "data/inputs/analysis/countries/{ISO3}/{ISO3}_rwi.tif"
    script:
        "./rwi_prep.py"
""" 
Test with
snakemake -c1 data/inputs/analysis/KEN/KEN_rwi.tif
"""
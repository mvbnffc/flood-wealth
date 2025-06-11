"""
Prepare (clip and resample) the GDP data
"""

rule clip_and_resample_gdp:
    """
    Resample clip the GDP data (using population data as reference)
    """
    input:
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
        gdp_file="data/inputs/gdp/2019gdp.tif",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
    output:
        gdp_resampled = "data/inputs/analysis/countries/{ISO3}/{ISO3}_gdp.tif"
    script:
        "./gdp_prep.py"
""" 
Test with
snakemake -c1 data/inputs/analysis/countries/KEN/KEN_gdp.tif
"""
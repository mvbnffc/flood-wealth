# Going to create a rasterized infrastructure layer based on OSM road data

rule rasterize_osm_infra:
    """
    This rule rasterizes an infrastructure layer based on OSM data.
    Population dataset is used as reference for gridded layer
    """
    input:
        osm_folder="data/inputs/analysis/countries/{ISO3}/{ISO3}_osm/",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif"
    output:
        infrastructure_raster="data/inputs/analysis/countries/{ISO3}/{ISO3}_infra.tif"
    script:
        "./rasterize_infra.py"

"""
Test with
snakemake -c1 data/inputs/analysis/countries/RWA/RWA_infra.tif"
"""
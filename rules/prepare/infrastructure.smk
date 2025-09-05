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

rule giri_infrastructure_valuation:
    """
    This rule takes the GIRI infra capital stock data and disaggregates
    it across the rasterized infrastructure layer for a country.
    """
    input:
        infrastructure_raster="data/inputs/analysis/countries/{ISO3}/{ISO3}_infra.tif",
        capital_stock="config/giri_infra_data.csv"
    output:
        infra_value="data/inputs/analysis/countries/{ISO3}/{ISO3}_inf_capstock.tif"
    script:
        "./distribute_giri_infr_capital.py"

"""
Test with
snakemake -c1 data/inputs/analysis/countries/RWA/RWA_inf_capstock.tif"
"""

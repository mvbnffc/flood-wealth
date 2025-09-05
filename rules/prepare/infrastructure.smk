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
        infrastructure_raster="data/inputs/analysis/countries/{ISO3}/{ISO3}_infra_raw.tif"
    script:
        "./rasterize_infra.py"

"""
Test with
snakemake -c1 data/inputs/analysis/countries/RWA/RWA_infra.tif"
"""

rule clip_rasterized_infra:
    """
    This rule clips the rasterized infrastructure layer.
    """
    input:
        raw_infrastructure_raster="data/inputs/analysis/countries/{ISO3}/{ISO3}_infra_raw.tif",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif"
    output:
        infrastructure_raster="data/inputs/analysis/countries/{ISO3}/{ISO3}_infra.tif"
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.infrastructure_raster})
        
        # Clip raster using GeoJSON geometry
        gdalwarp \
            -cutline {input.boundary_file} \
            -crop_to_cutline \
            -tr 0.00083333333333333 0.00083333333333333 \
            -tap \
            -te_srs EPSG:4326 \
            -te $(gdalinfo -json {input.pop_file} | jq -r '.cornerCoordinates | [.upperLeft[0], .lowerLeft[1], .lowerRight[0], .upperRight[1]] | join(" ")') \
            -of GTiff \
            -co BIGTIFF=YES \
            -co compress=lzw \
            {input.raw_infrastructure_raster} \
            {output.infrastructure_raster}
        """

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

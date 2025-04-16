"""
Merge, rasterize, and clip the Google historical inundation data
"""

rule merge_and_rasterize_google_inundation:
    input:
        raw_folder = "data/inputs/google_inun/raw/",
        pop_path = "data/inputs/ghs-pop/GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"
    output:
        merged_file="data/inputs/google_inun/merged/google_historical_inun.tif"
    script:
        "./rasterize_google_inun.py"

rule clip_google_inundation:
    """
    Clip GFD flood raster to country boundary.
    """
    input:
        raw_flood_file="data/inputs/google_inun/merged/google_historical_inun.tif",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
    output:
        trimmed_flood_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_google-flood.tif",
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.trimmed_flood_file})
        
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
            {input.raw_flood_file} \
            {output.trimmed_flood_file}
        """
""" 
Test with
snakemake -c1 data/inputs/analysis/KEN/KEN_google-flood.tif
"""
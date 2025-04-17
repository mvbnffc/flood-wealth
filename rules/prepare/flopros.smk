"""
Rasterize and clip FLOPROS
"""

rule rasterize_flopros:
    """
    We are going to rasterize the flopros layer (using the merged protection value column).
    """
    input:
        flopros="data/inputs/flopros/Scussolini_etal_Suppl_info/FLOPROS_shp_V1/FLOPROS_shp_V1.shp",
        pop_file="data/inputs/ghs-pop/GHS_POP_E2020_GLOBE_R2023A_4326_3ss_V1_0.tif"
    output:
        rasterized_flopros="data/inputs/flopros/flopros.tif"
    script:
        "./rasterize_flopros.py"
"""
Test with 
snakemake -c1 data/inputs/flopros/flopros.tif
"""

rule clip_flopros:
    """
    Clip FLOPROS raster to country boundary.
    """
    input:
        raw_flopros="data/inputs/flopros/flopros.tif",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
    output:
        trimmed_flopros_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_flopros.tif",
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.trimmed_flopros_file})
        
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
            {input.raw_flopros} \
            {output.trimmed_flopros_file}
        """
""" 
Test with
snakemake -c1 data/inputs/analysis/countries/KEN/KEN_flopros.tif
"""

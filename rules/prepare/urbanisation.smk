"""
Prepare urbanisation dataset (combine GADM with CSV) and then rasterize
Also clip, reproject, and resample GHS-MOD dataset
"""

rule prepare_urbanisation:
    """
    Merge the smallest GADM level with the GHS DUC CSV values on urbanisation
    """
    input:
        duc_info="data/inputs/ghs-duc/GHS_DUC_GLOBE_R2023A_V2_0.xlsx",
        adm_file="data/inputs/boundaries/{ISO3}/gadm_{ISO3}.gpkg",
    output:
        urban_file= "data/inputs/analysis/countries/{ISO3}/{ISO3}_urbanization.gpkg"
    script:
        "./prepare_urban.py"

rule rasterize_urbanization:
    """
    Rasterize the urbanization layer
    """
    input:
        urbanization="data/inputs/analysis/countries/{ISO3}/{ISO3}_urbanization.gpkg",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
    output:
        rasterized="data/inputs/analysis/countries/{ISO3}/{ISO3}_urbanization.tif",
    script:
        "./rasterize_urban.py"

rule clip_ghs_mod:
    """
    Clip GHS-MOD dataset to country boundary.
    Also reproject from native Mollweide projection to WGS84
    Also resample from 1 km resolution to 0.000833 degree res (~90 m) using pop dataset as reference
    """
    input:
        ghs_mod_file="data/inputs/ghs-mod/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000_V2_0.tif",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
    output:
        trimmed_ghs_mod="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-mod.tif",
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.trimmed_ghs_mod})
        
        # Clip raster using GeoJSON geometry
        gdalwarp \
            -cutline {input.boundary_file} \
            -crop_to_cutline \
            -tr 0.00083333333333333 0.00083333333333333 \
            -tap \
            -s_srs ESRI:54009 \
            -t_srs EPSG:4326 \
            -te_srs EPSG:4326 \
            -te $(gdalinfo -json {input.pop_file} | jq -r '.cornerCoordinates | [.upperLeft[0], .lowerLeft[1], .lowerRight[0], .upperRight[1]] | join(" ")') \
            -of GTiff \
            -co BIGTIFF=YES \
            -co compress=lzw \
            {input.ghs_mod_file} \
            {output.trimmed_ghs_mod}
        """
""" 
Test with
snakemake -c1 data/inputs/analysis/countries/KEN/KEN_ghs-mod.tif
"""
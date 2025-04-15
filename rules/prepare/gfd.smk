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

rule clip_gfd:
    """
    Clip GFD flood raster to country boundary.
    """
    input:
        raw_flood_file="data/inputs/gfd/merged/gfd.tif",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
        pop_file="data/inputs/analysis/{ISO3}/{ISO3}_ghs-pop.tif",
    output:
        trimmed_flood_file="data/inputs/analysis/{ISO3}/{ISO3}_gfd-flood.tif",
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
snakemake -c1 data/inputs/analysis/KEN/KEN_gfd-flood.tif
"""
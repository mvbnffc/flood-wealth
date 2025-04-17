"""
Prepare (merge and clip) permanent water body data
"""

rule surface_water_vrt:
    """
    Creates a virtaul raster of all surface water tiles.
    """
    input:
        raw_folder="data/inputs/jrc-perm/occurrence/"
    output:
        surface_water_vrt="data/inputs/jrc-perm/merged/surface_water.vrt"
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.surface_water_vrt})
        
        # Build a virtual raster (VRT) from all raw .tif files in the input folder.
        gdalbuildvrt {output.surface_water_vrt} {input.raw_folder}/*.tif
    """
""" 
Test with
snakemake -c1 data/inputs/jrc-perm/merged/surface_water.vrt
"""

rule clip_surface_water:
    """
    Clips and resamples surface water dataset. Will use average resampling.
    """
    input:
        surface_water_vrt="data/inputs/jrc-perm/merged/surface_water.vrt",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
    output:
        surface_water_clipped="data/inputs/analysis/countries/{ISO3}/{ISO3}_surface_water.tif"
    shell:
        """
        gdalwarp \
        -cutline {input.boundary_file} \
        -crop_to_cutline \
        -tap \
        -tr 0.00083333333333333 0.00083333333333333 \
        -te_srs EPSG:4326 \
        -te $(gdalinfo -json {input.pop_file} | jq -r '.cornerCoordinates | [.upperLeft[0], .lowerLeft[1], .lowerRight[0], .upperRight[1]] | join(" ")') \
        -of GTiff \
        -r average \
        -co BIGTIFF=YES \
        -co compress=lzw \
        {input.surface_water_vrt} \
        {output.surface_water_clipped}
        """
""" 
Test with
snakemake -c1 data/inputs/analysis/KEN/KEN_surface_water.tif
"""
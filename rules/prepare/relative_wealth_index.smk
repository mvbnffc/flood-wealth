"""
Prepare (clip and resample) the RWI data
"""

rule clip_and_resample_rwi:
    """
    Resample the and clip the rwi (using population data as reference)
    """
    input:
        pop_file="data/inputs/analysis/{ISO3}/{ISO3}_ghs-pop_{RESOLUTION}.tif",
        rwi_file="data/inputs/rwi/rwi.tif",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
    output:
        rwi_resampled = "data/inputs/analysis/{ISO3}/{ISO3}_rwi_{RESOLUTION}.tif"
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.rwi_resampled})
        

        # Pull reference pixel resolution
        tr=$(gdalinfo {input.pop_file} | grep "Pixel Size" | head -n1 | sed -e 's/.*(//' -e 's/).*//' | awk -F, '{{print $1, $2}}')
        
        # Clip and resample rwi
        gdalwarp \
            -cutline {input.boundary_file} \
            -crop_to_cutline \
            -dstalpha \
            -t_srs EPSG:4326 \
            -r nearest \
            -tr $tr \
            -co compress=lzw \
            {input.rwi_file} \
            {output.rwi_resampled}
        """
""" 
Test with
snakemake -c1 data/inputs/analysis/KEN/KEN_rwi_3ss.tif
"""
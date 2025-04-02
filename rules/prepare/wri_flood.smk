"""
Prepare (clip) WRI Flood Data
"""

rule clip_wri_flood:
    """
    Clip WRI flood raster to country boundary. 
    """
    input:
        raw_flood_file=lambda wc: (
            f"data/inputs/flood/WRI/inunriver_historical_inunriver_historical_000000000WATCH_1980_rp{'%05d' % int(wc.RP)}.tif"
        ),
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson"
    output:
        trimmed_flood_file="data/inputs/analysis/{ISO3}/{ISO3}_wri-flood_RP{RP}.tif"
    wildcard_constraints:
        RP="2|5|10|25|50|100|250|500|1000"
    shell:
        """
        set -ex

        mkdir --parents $(dirname {output.trimmed_flood_file})
        
        # Clip raster using GeoJSON geometry
        gdalwarp \
            -cutline {input.boundary_file} \
            -crop_to_cutline \
            -dstalpha \
            -of GTiff \
            -co compress=lzw \
            {input.raw_flood_file} \
            {output.trimmed_flood_file}
        """
""" 
Test with
snakemake -c1 data/inputs/analysis/KEN/KEN_wri-flood_RP10.tif
"""
"""
Download WRI Global River Flood Data

References
---------
WRI: https://www.wri.org/data/aqueduct-floods-hazard-maps
"""

rule download_wri_flood:
    output:
        flood_file="data/inputs/flood/WRI/inunriver_historical_inunriver_historical_000000000WATCH_1980_rp0{RP}.tif"
    wildcard_constraints:
        RP="0002|0005|0010|0025|0050|0100|0250|0500|1000"
    shell:
        """
        mkdir -p $(dirname {output.flood_file})
        wget -nc https://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/inunriver_historical_000000000WATCH_1980_rp0{wildcards.RP}.tif -O {output.flood_file}
        """
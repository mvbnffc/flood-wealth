"""
Download JRC Global River Flood Data

References
---------
JRC: https://data.jrc.ec.europa.eu/dataset/jrc-floods-floodmapgl_rp50y-tif
"""

rule all_jrc_flood:
    input:
        tiffs=expand("data/inputs/flood/JRC/merged/jrc_global_flood_RP{RP}.tif", RP=[10, 20, 50, 75, 100, 200, 500])

rule download_jrc_flood:
    """
    Downloads JRC global flood maps (in tiles) and merges into a single global file
    """
    output:
        RP10="data/inputs/flood/JRC/merged/jrc_global_flood_RP10.tif",
        RP20="data/inputs/flood/JRC/merged/jrc_global_flood_RP20.tif",
        RP50="data/inputs/flood/JRC/merged/jrc_global_flood_RP50.tif",
        RP75="data/inputs/flood/JRC/merged/jrc_global_flood_RP75.tif",
        RP100="data/inputs/flood/JRC/merged/jrc_global_flood_RP100.tif",
        RP200="data/inputs/flood/JRC/merged/jrc_global_flood_RP200.tif",
        RP500="data/inputs/flood/JRC/merged/jrc_global_flood_RP500.tif"
    script:
        "./download_jrc.py"
"""
Test with
snakemake -c1 data/inputs/flood/JRC/merged/jrc_global_flood_RP10.tif
"""
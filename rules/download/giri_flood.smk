"""
Download GIRI Global River Flood Data

References
---------
GIRI: https://giri.unepgrid.ch/
"""

rule download_giri_flood:
    output:
        flood_file="data/inputs/flood/GIRI/global_pc_h{RP}glob.tif"
    wildcard_constraints:
        RP="5|10|25|50|100|200|500|1000"
    shell:
        """
        mkdir -p $(dirname {output.flood_file})
        wget -nc https://hazards-data.unepgrid.ch/global_pc_h{wildcards.RP}glob.tif -O {output.flood_file}
        """

rule convert_giri_flood:
    """
    convert flood depths from cm to m
    """
    input:
        flood_file="data/inputs/flood/GIRI/global_pc_h{RP}glob.tif"
    output:
        converted_flood_file="data/inputs/flood/GIRI/global_pc_h{RP}glob_m.tif"
    wildcard_constraints:
        RP="5|10|25|50|100|200|500|1000"
    threads: 4  # Adjust based on your system capabilities
    resources:
        mem_mb=8000  # Adjust based on your system memory
    shell:
        """
        gdal_calc.py -A {input.flood_file} --outfile={output.converted_flood_file} \
        --calc="A/100" --NoDataValue=-9999 \
        --co="TILED=YES" --co="BLOCKSIZE=512" \
        --co="COMPRESS=LZW" --co="BIGTIFF=YES" \
        --overwrite --quiet \
        --CPU {threads}
        """

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

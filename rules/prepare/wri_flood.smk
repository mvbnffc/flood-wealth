"""
Prepare (clip) WRI Flood Data
"""

rule clip_wri_flood:
    """
    Clip and resample WRI flood map (using population dataset as reference)
    Note: we are using rwi_prep.py script as it works the same.
    """
    input:
        pop_file="data/inputs/analysis/{ISO3}/{ISO3}_ghs-pop.tif",
        rwi_file=lambda wc: (
            f"data/inputs/flood/WRI/inunriver_historical_inunriver_historical_000000000WATCH_1980_rp{'%05d' % int(wc.RP)}.tif"
        ),
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
    wildcard_constraints:
        RP="2|5|10|25|50|100|250|500|1000"
    output:
        rwi_resampled = "data/inputs/analysis/{ISO3}/{ISO3}_wri-flood_RP{RP}.tif"
    script:
        "./rwi_prep.py"
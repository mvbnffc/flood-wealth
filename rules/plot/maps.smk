"""
Rule book for map plotting.
"""

rule flood_risk_exposure:
    """
    This rule creates a flood risk exposure GeoTiff by multiplying AAR with the population map
    """
    input:
        pop_file="data/inputs/analysis/{ISO3}/{ISO3}_ghs-pop.tif",
        risk_file="data/results/flood_risk/{ISO3}/{ISO3}_{MODEL}-flood-risk_AAR_V-{VULN_CURVE}.tif",
    output:
        risk_exposure="data/results/social_flood/{ISO3}/map_layers/{ISO3}_{MODEL}-flood-risk-exposure_V-{VULN_CURVE}.tif",
    wildcard_constraints:
        MODEL="giri|jrc|wri",
        VULN_CURVE="BER|JRC|EXP",
    script:
        "./flood_risk_exposure.py"
"""
Test with
snakemake -c1 data/results/social_flood/RWA/map_layers/RWA_jrc-flood-risk-exposure_V-JRC.tif
"""

rule pop_wealth_quntiles:
    """
    This rule returns a population map classified by wealth quintile (values 1-5 from bottom to top quantile)
    """
    input:
        pop_file="data/inputs/analysis/{ISO3}/{ISO3}_ghs-pop.tif",
        rwi_file="data/inputs/analysis/{ISO3}/{ISO3}_rwi.tif"
    output:
        wealth_quintiles="data/results/social_flood/{ISO3}/map_layers/{ISO3}_wealth_quintiles.tif",
    script:
        "./pop_wealth_quintiles.py"
"""
Test with
snakemake -c1 data/results/social_flood/RWA/map_layers/RWA_wealth_quintiles.tif
"""
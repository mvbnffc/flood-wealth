"""
Rulebook for the social flood analyses (concentration curves and inequality metrics)
"""

rule inequality_metrics:
    """
    This rule calcualtes two inequality metrics at the specified administrative level.
    Inequality metrics:
        - Concentration Index (CI) - understand the inequality of flood risk across the wealth distribution
        - Quantile Ratio (QR) - understand the tail inequality (20:80)
    """
    input:
        admin_areas = "data/inputs/boundaries/{ISO3}/gadm_{ISO3}.gpkg",
        rwi_file="data/inputs/analysis/{ISO3}/{ISO3}_rwi.tif",
        pop_file="data/inputs/analysis/{ISO3}/{ISO3}_ghs-pop.tif",
        risk_file="data/results/flood_risk/{ISO3}/{ISO3}_{MODEL}-flood-risk_AAR_V-{VULN_CURVE}.tif",
    output:
        regional_CI = "data/results/social_flood/{ISO3}/inequality_metrics/{ISO3}_{ADMIN_SLUG}_metrics_{MODEL}-flood_V-{VULN_CURVE}.gpkg",
    wildcard_constraints:
        MODEL="giri|jrc|wri",
        VULN_CURVE="BER|JRC|EXP",
        ADMIN_SLUG="ADM-0|ADM-1|ADM-2"
    script:
        "./inequality_metrics.py"
"""
Test with
snakemake -c1 data/results/social_flood/KEN/inequality_metrics/KEN_ADM-0_metrics_jrc-flood_V-JRC.gpkg
"""

# Run for all ISO3 codes, models, and admins

configfile: "config/config.yaml"
ADMINS = ["ADM-0", "ADM-1"]
MODELS = ["jrc", "wri", "giri"]

rule metrics_for_all_countries:
    input:
        expand("data/results/social_flood/{ISO3}/inequality_metrics/{ISO3}_{ADM}_metrics_{MODEL}-flood_V-JRC.gpkg",
                ISO3=config['iso_codes'], ADM=ADMINS, MODEL=MODELS)
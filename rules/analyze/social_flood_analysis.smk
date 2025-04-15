"""
Rulebook for the social flood analyses (concentration curves and inequality metrics)
"""

import json

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
        risk_file="data/results/flood_risk/{ISO3}/{ISO3}_{MODEL}-flood-risk_{TYPE}_V-{VULN_CURVE}.tif",
    output:
        regional_CI = "data/results/social_flood/{ISO3}/inequality_metrics/{ISO3}_{ADMIN_SLUG}_metrics_{MODEL}-flood_{TYPE}_V-{VULN_CURVE}.gpkg",
    wildcard_constraints:
        MODEL="giri|jrc|wri",
        TYPE="AAR|RP100",
        VULN_CURVE="BER|JRC|EXP",
        ADMIN_SLUG="ADM-0|ADM-1|ADM-2"
    script:
        "./inequality_metrics.py"
"""
Test with
snakemake -c1 data/results/social_flood/KEN/inequality_metrics/KEN_ADM-0_metrics_jrc-flood_V-JRC.gpkg
"""

rule inequality_metrics_observed:
    """
    This rule calcualtes two inequality metrics at the specified administrative level. For observed flooding datasets
    Inequality metrics:
        - Concentration Index (CI) - understand the inequality of flood risk across the wealth distribution
        - Quantile Ratio (QR) - understand the tail inequality (20:80)
    """
    input:
        admin_areas = "data/inputs/boundaries/{ISO3}/gadm_{ISO3}.gpkg",
        rwi_file="data/inputs/analysis/{ISO3}/{ISO3}_rwi.tif",
        pop_file="data/inputs/analysis/{ISO3}/{ISO3}_ghs-pop.tif",
        risk_file="data/inputs/analysis/{ISO3}/{ISO3}_{MODEL}-flood.tif",
    output:
        regional_CI = "data/results/social_flood/{ISO3}/inequality_metrics/{ISO3}_{ADMIN_SLUG}_metrics_{MODEL}-flood.gpkg",
    wildcard_constraints:
        MODEL="gfd|google",
        ADMIN_SLUG="ADM-0|ADM-1|ADM-2"
    script:
        "./inequality_metrics.py"
"""
Test with
snakemake -c1 data/results/social_flood/KEN/inequality_metrics/KEN_ADM-0_metrics_gfd-flood.gpkg
"""

def get_event_iso3s(wildcards):
    """Return a list of valid ISO3 codes for an event by reading its properties file."""
    props_path = f"data/inputs/analysis/events/DFO_{wildcards.event_id}/countries.json"
    with open(props_path, "r") as f:
        props = json.load(f)
    return props['valid']

rule dfo_event_analysis:
    """
    This rule carries out a risk assessment for individual dfo events.
    Reporting various metrics and returning a CSV file of results.
    """
    input:
        rwi_file = lambda wc: expand("data/inputs/analysis/{ISO3}/{ISO3}_rwi.tif", ISO3=get_event_iso3s(wc)),
        pop_file = lambda wc: expand("data/inputs/analysis/{ISO3}/{ISO3}_ghs-pop.tif", ISO3=get_event_iso3s(wc)),
        mask_file = lambda wc: expand("data/inputs/analysis/{ISO3}/{ISO3}_surface_water.tif", ISO3=get_event_iso3s(wc)),
        flood_file = lambda wc: expand("data/inputs/analysis/events/DFO_{event_id}/{ISO3}_{event_id}.tif", ISO3=get_event_iso3s(wc), event_id=wc.event_id)
    output:
        results = "data/results/social_flood/events/DFO_{event_id}/DFO_{event_id}_results.csv"
    params:
        iso3_list = lambda wc: get_event_iso3s(wc)
    script:
        "./dfo_event_risk_analysis.py"

"""
Test with
snakemake -c1 data/results/social_flood/events/DFO_1595/DFO_1595_results.csv
"""


############################## BULK ANALYSIS ###########################################

# Run modelled metrics for all ISO3 codes, models, and admins

configfile: "config/config.yaml"
ADMINS = ["ADM-0"]
MODELS = ["jrc", "wri", "giri"]
TYPES = ["RP100", "AAR"]
VULN_CURVES = ["JRC", "EXP", "BER"]

rule metrics_for_all_countries:
    input:
        expand("data/results/social_flood/{ISO3}/inequality_metrics/{ISO3}_{ADM}_metrics_{MODEL}-flood_{TYPE}_V-{VULN_CURVE}.gpkg",
                ISO3=config['iso_codes'], ADM=ADMINS, MODEL=MODELS, TYPE=TYPES, VULN_CURVE=VULN_CURVES)


# Run observed modelled metrics for all ISO3 codes
rule observed_metrics_for_all_countries:
    input:
        expand("data/results/social_flood/{ISO3}/inequality_metrics/{ISO3}_ADM-0_metrics_google-flood.gpkg", ISO3=config['iso_codes'])
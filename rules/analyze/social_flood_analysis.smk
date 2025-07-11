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
        social_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_{SOCIAL}.tif",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
        mask_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_surface_water.tif",
        risk_file="data/results/flood_risk/countries/{ISO3}/{ISO3}_{MODEL}-flood-risk_{TYPE}_V-{VULN_CURVE}.tif",
    output:
        regional_CI = "data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_{ADMIN_SLUG}_metrics_{MODEL}-flood_{TYPE}_V-{VULN_CURVE}_S-{SOCIAL}.gpkg",
    wildcard_constraints:
        MODEL="giri|jrc|wri",
        TYPE="AAR|RP100",
        SOCIAL="rwi|gdp",
        VULN_CURVE="BER|JRC|EXP",
        ADMIN_SLUG="ADM-0|ADM-1|ADM-2"
    script:
        "./inequality_metrics.py"
"""
Test with
snakemake -c1 data/results/social_flood/countries/KEN/inequality_metrics/KEN_ADM-0_metrics_jrc-flood_AAR_V-JRC_S-rwi.gpkg
"""

rule inequality_metrics_protected:
    """
    This rule calcualtes two inequality metrics at the specified administrative level. FLOPROS protection is ON.
    Inequality metrics:
        - Concentration Index (CI) - understand the inequality of flood risk across the wealth distribution
        - Quantile Ratio (QR) - understand the tail inequality (20:80)
    """
    input:
        admin_areas = "data/inputs/boundaries/{ISO3}/gadm_{ISO3}.gpkg",
        social_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_{SOCIAL}.tif",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
        mask_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_surface_water.tif",
        risk_file="data/results/flood_risk/countries/{ISO3}/{ISO3}_{MODEL}-flood-risk_protected_AAR_V-{VULN_CURVE}.tif",
    output:
        regional_CI = "data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_{ADMIN_SLUG}_metrics_{MODEL}-flood_protected_AAR_V-{VULN_CURVE}_S-{SOCIAL}.gpkg",
    wildcard_constraints:
        MODEL="giri|jrc|wri",
        VULN_CURVE="BER|JRC|EXP",
        SOCIAL="rwi|gdp",
        ADMIN_SLUG="ADM-0|ADM-1|ADM-2"
    script:
        "./inequality_metrics.py"
"""
Test with
snakemake -c1 data/results/social_flood/countries/KEN/inequality_metrics/KEN_ADM-0_metrics_jrc-flood_protected_AAR_V-JRC_S-rwi.gpkg
"""

rule inequality_metrics_flood_protection:
    """
    This rule calcualtes two inequality metrics at the specified administrative level for the flood protection adaptation scenario.
    Adaptation parameters required as input are RP protection and min level of urbanization to protect.
    Inequality metrics:
        - Concentration Index (CI) - understand the inequality of flood risk across the wealth distribution
        - Quantile Ratio (QR) - understand the tail inequality (20:80)
    """
    input:
        admin_areas = "data/inputs/boundaries/{ISO3}/gadm_{ISO3}.gpkg",
        social_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_{SOCIAL}.tif",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
        mask_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_surface_water.tif",
        risk_file="data/results/flood_risk/countries/{ISO3}/{ISO3}_{MODEL}-flood-risk_adapted_AAR_V-{VULN_CURVE}_fp_rp{RP}_duc{urban_class}.tif",
    output:
        regional_CI = "data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_{ADMIN_SLUG}_metrics_{MODEL}-flood_adapted_AAR_V-{VULN_CURVE}_S-{SOCIAL}_fp_rp{RP}_duc{urban_class}.gpkg",
    wildcard_constraints:
        MODEL="giri|jrc|wri",
        VULN_CURVE="BER|JRC|EXP",
        SOCIAL="rwi|gdp",
        urban_class="11|12|13|21|22|23|30",
        ADMIN_SLUG="ADM-0|ADM-1|ADM-2"
    script:
        "./inequality_metrics.py"
"""
Test with
snakemake -c1 data/results/social_flood/countries/KEN/inequality_metrics/KEN_ADM-0_metrics_jrc-flood_adapted_AAR_V-JRC_S-rwi_fp_rp100_duc23.gpkg
"""

rule inequality_metrics_relocation:
    """
    This rule calcualtes two inequality metrics at the specified administrative level for the relocation adaptation scenario.
    Adaptation parameter required as input is the max level of urbanization to relocate people from.
    Inequality metrics:
        - Concentration Index (CI) - understand the inequality of flood risk across the wealth distribution
        - Quantile Ratio (QR) - understand the tail inequality (20:80)
    """
    input:
        admin_areas = "data/inputs/boundaries/{ISO3}/gadm_{ISO3}.gpkg",
        social_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_{SOCIAL}.tif",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
        mask_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_surface_water.tif",
        risk_file="data/results/flood_risk/countries/{ISO3}/{ISO3}_{MODEL}-flood-risk_adapted_AAR_V-{VULN_CURVE}_rl_duc{urban_class}.tif",
    output:
        regional_CI = "data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_{ADMIN_SLUG}_metrics_{MODEL}-flood_adapted_AAR_V-{VULN_CURVE}_S-{SOCIAL}_rl_duc{urban_class}.gpkg",
    wildcard_constraints:
        MODEL="giri|jrc|wri",
        VULN_CURVE="BER|JRC|EXP",
        SOCIAL="rwi|gdp",
        urban_class="11|12|13|21|22|23|30",
        ADMIN_SLUG="ADM-0|ADM-1|ADM-2"
    script:
        "./inequality_metrics.py"
"""
Test with
snakemake -c1 data/results/social_flood/countries/KEN/inequality_metrics/KEN_ADM-0_metrics_jrc-flood_adapted_AAR_V-JRC_S-rwi_rl_duc23.gpkg
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
        social_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_{SOCIAL}.tif",
        pop_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif",
        mask_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_surface_water.tif",
        risk_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_{MODEL}-flood.tif",
    output:
        regional_CI = "data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_{ADMIN_SLUG}_metrics_{MODEL}-flood_S-{SOCIAL}.gpkg",
    wildcard_constraints:
        MODEL="gfd|google",
        SOCIAL="rwi|gdp",
        ADMIN_SLUG="ADM-0|ADM-1|ADM-2"
    script:
        "./inequality_metrics.py"
"""
Test with
snakemake -c1 data/results/social_flood/countries/KEN/inequality_metrics/KEN_ADM-0_metrics_gfd-flood_S-rwi.gpkg
"""

def get_event_iso3s(wildcards):
    """Return a list of valid ISO3 codes for an event by reading its properties file."""
    props_path = f"data/inputs/analysis/events/DFO_{wildcards.event_id}/countries.json"
    with open(props_path, "r") as f:
        props = json.load(f)
    # TEMP: PC DEBUG
    problem_isos = config['problem_iso_codes'] # DEBUG code
    return [iso for iso in props['valid'] if iso not in problem_isos] # DEBUG code
    # return props['valid'] # original

rule dfo_event_analysis:
    """
    This rule carries out a risk assessment for individual dfo events.
    Reporting various metrics and returning a CSV file of results.
    """
    input:
        rwi_file = lambda wc: expand("data/inputs/analysis/countries/{ISO3}/{ISO3}_rwi.tif", ISO3=get_event_iso3s(wc)),
        pop_file = lambda wc: expand("data/inputs/analysis/countries/{ISO3}/{ISO3}_ghs-pop.tif", ISO3=get_event_iso3s(wc)),
        mask_file = lambda wc: expand("data/inputs/analysis/countries/{ISO3}/{ISO3}_surface_water.tif", ISO3=get_event_iso3s(wc)),
        flood_file = lambda wc: expand("data/inputs/analysis/events/DFO_{event_id}/{ISO3}_{event_id}.tif", ISO3=get_event_iso3s(wc), event_id=wc.event_id),
        country_json="data/inputs/analysis/events/DFO_{event_id}/countries.json"
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
MODELS = ["jrc", "giri", "wri"]
TYPES = ["AAR"]
VULN_CURVES = ["JRC"]
SOCIALS = ['rwi']
RPs = [100]
DUC_protection = [21, 22, 23, 30]
DUC_relocation = [11, 12, 13]

rule metrics_for_all_countries:
    input:
        expand("data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_{ADM}_metrics_{MODEL}-flood_{TYPE}_V-{VULN_CURVE}_S-{SOCIAL}.gpkg",
                ISO3=config['iso_codes'], ADM=ADMINS, MODEL=MODELS, TYPE=TYPES, VULN_CURVE=VULN_CURVES, SOCIAL=SOCIALS)

rule protected_metrics_for_all_countries:
    input:
        expand("data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_{ADM}_metrics_{MODEL}-flood_protected_AAR_V-{VULN_CURVE}_S-{SOCIAL}.gpkg",
                ISO3=config['iso_codes'], ADM=ADMINS, MODEL=MODELS, VULN_CURVE=VULN_CURVES, SOCIAL=SOCIALS)

rule adapted_flood_protection_metrics_for_all_countries:
    input:
        expand("data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_{ADM}_metrics_{MODEL}-flood_adapted_AAR_V-{VULN_CURVE}_S-{SOCIAL}_fp_rp{RP}_duc{DUC}.gpkg",
                ISO3=config['iso_codes'], ADM=ADMINS, MODEL=MODELS, VULN_CURVE=VULN_CURVES, SOCIAL=SOCIALS, RP=RPs, DUC=DUC_protection)

rule adapted_relocation_metrics_for_all_countries:
    input:
        expand("data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_{ADM}_metrics_{MODEL}-flood_adapted_AAR_V-{VULN_CURVE}_S-{SOCIAL}_rl_duc{DUC}.gpkg",
                ISO3=config['iso_codes'], ADM=ADMINS, MODEL=MODELS, VULN_CURVE=VULN_CURVES, SOCIAL=SOCIALS, DUC=DUC_protection)


# Run observed modelled metrics for all ISO3 codes
rule observed_metrics_for_all_countries:
    input:
        expand("data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_ADM-0_metrics_gfd-flood_S-{SOCIAL}.gpkg", ISO3=config['iso_codes'], SOCIAL=SOCIALS)

# Run metrics on all DFO flood events
# Find all events in the prep folder
events = glob_wildcards("data/inputs/gfd/prep/DFO_{event_id}.tif").event_id
rule metrics_all_gfd_events:
    input:
        expand("data/results/social_flood/events/DFO_{event_id}/DFO_{event_id}_results.csv", event_id=events)
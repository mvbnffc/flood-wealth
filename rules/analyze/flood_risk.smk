"""
Run flood risk analysis
"""

rule relative_flood_risk:
    """
    This rule calculates the relative risk per grid cell according to a vulnerability curve.
    There are three approahces for calulation risk, using three different vulnerability curves:
        V-JRC: JRC Global residential curve (average of Asia, Africa, and SAmerica)
        V-BER: The Bernhofen et al (2023) curve
        V-EXP: Exposure approach (i.e. 100% damage at any flood depth)
    """
    input:
        flood_file="data/inputs/analysis/countries/{ISO3}/{ISO3}_{MODEL}-flood_RP{RP}.tif"
    output:
        risk_file="data/results/flood_risk/countries/{ISO3}/{ISO3}_{MODEL}-flood-risk_RP{RP}_V-{VULN_CURVE}.tif"
    wildcard_constraints:
        VULN_CURVE="BER|JRC|EXP",
        MODEL="giri|jrc|wri"
    script:
        "./relative_flood_risk.py"
""" 
Test with
snakemake -c1 data/results/flood_risk/countries/KEN/KEN_giri-flood-risk_RP10_V-JRC.tif
"""

rule giri_average_annual_risk:
    """
    This rule calculates one layer (average annual relative risk) given a set of 
    GIRI return period flood maps. 
    """
    input:
        flood_rp_5="data/results/flood_risk/countries/{ISO3}/{ISO3}_giri-flood-risk_RP5_V-{VULN_CURVE}.tif",
        flood_rp_10="data/results/flood_risk/countries/{ISO3}/{ISO3}_giri-flood-risk_RP10_V-{VULN_CURVE}.tif",
        flood_rp_25="data/results/flood_risk/countries/{ISO3}/{ISO3}_giri-flood-risk_RP25_V-{VULN_CURVE}.tif",
        flood_rp_50="data/results/flood_risk/countries/{ISO3}/{ISO3}_giri-flood-risk_RP50_V-{VULN_CURVE}.tif",
        flood_rp_100="data/results/flood_risk/countries/{ISO3}/{ISO3}_giri-flood-risk_RP100_V-{VULN_CURVE}.tif",
        flood_rp_200="data/results/flood_risk/countries/{ISO3}/{ISO3}_giri-flood-risk_RP200_V-{VULN_CURVE}.tif",
        flood_rp_500="data/results/flood_risk/countries/{ISO3}/{ISO3}_giri-flood-risk_RP500_V-{VULN_CURVE}.tif",
        flood_rp_1000="data/results/flood_risk/countries/{ISO3}/{ISO3}_giri-flood-risk_RP1000_V-{VULN_CURVE}.tif"
    output:
        flood_aar="data/results/flood_risk/countries/{ISO3}/{ISO3}_giri-flood-risk_AAR_V-{VULN_CURVE}.tif"
    wildcard_constraints:
        VULN_CURVE="BER|JRC|EXP"
    script:
        "./giri_average_annual_risk.py"
"""
Test with
snakemake -c1 data/results/flood_risk/countries/KEN/KEN_giri-flood-risk_AAR_V-JRC.tif
"""

rule jrc_average_annual_risk:
    """
    This rule calculates one layer (average annual relative risk) given a set of 
    JRC return period flood maps. 
    """
    input:
        flood_rp_10="data/results/flood_risk/countries/{ISO3}/{ISO3}_jrc-flood-risk_RP10_V-{VULN_CURVE}.tif",
        flood_rp_20="data/results/flood_risk/countries/{ISO3}/{ISO3}_jrc-flood-risk_RP20_V-{VULN_CURVE}.tif",
        flood_rp_50="data/results/flood_risk/countries/{ISO3}/{ISO3}_jrc-flood-risk_RP50_V-{VULN_CURVE}.tif",
        flood_rp_75="data/results/flood_risk/countries/{ISO3}/{ISO3}_jrc-flood-risk_RP75_V-{VULN_CURVE}.tif",
        flood_rp_100="data/results/flood_risk/countries/{ISO3}/{ISO3}_jrc-flood-risk_RP100_V-{VULN_CURVE}.tif",
        flood_rp_200="data/results/flood_risk/countries/{ISO3}/{ISO3}_jrc-flood-risk_RP200_V-{VULN_CURVE}.tif",
        flood_rp_500="data/results/flood_risk/countries/{ISO3}/{ISO3}_jrc-flood-risk_RP500_V-{VULN_CURVE}.tif"
    output:
        flood_aar="data/results/flood_risk/countries/{ISO3}/{ISO3}_jrc-flood-risk_AAR_V-{VULN_CURVE}.tif"
    wildcard_constraints:
        VULN_CURVE="BER|JRC|EXP"
    script:
        "./jrc_average_annual_risk.py"
"""
Test with
snakemake -c1 data/results/flood_risk/countries/KEN/KEN_jrc-flood-risk_AAR_V-JRC.tif
"""

rule wri_average_annual_risk:
    """
    This rule calculates one layer (average annual relative risk) given a set of 
    WRI return period flood maps. 
    """
    input:
        flood_rp_2="data/results/flood_risk/countries/{ISO3}/{ISO3}_wri-flood-risk_RP2_V-{VULN_CURVE}.tif",
        flood_rp_5="data/results/flood_risk/countries/{ISO3}/{ISO3}_wri-flood-risk_RP5_V-{VULN_CURVE}.tif",
        flood_rp_10="data/results/flood_risk/countries/{ISO3}/{ISO3}_wri-flood-risk_RP10_V-{VULN_CURVE}.tif",
        flood_rp_25="data/results/flood_risk/countries/{ISO3}/{ISO3}_wri-flood-risk_RP25_V-{VULN_CURVE}.tif",
        flood_rp_50="data/results/flood_risk/countries/{ISO3}/{ISO3}_wri-flood-risk_RP50_V-{VULN_CURVE}.tif",
        flood_rp_100="data/results/flood_risk/countries/{ISO3}/{ISO3}_wri-flood-risk_RP100_V-{VULN_CURVE}.tif",
        flood_rp_250="data/results/flood_risk/countries/{ISO3}/{ISO3}_wri-flood-risk_RP250_V-{VULN_CURVE}.tif",
        flood_rp_500="data/results/flood_risk/countries/{ISO3}/{ISO3}_wri-flood-risk_RP500_V-{VULN_CURVE}.tif",
        flood_rp_1000="data/results/flood_risk/countries/{ISO3}/{ISO3}_wri-flood-risk_RP1000_V-{VULN_CURVE}.tif"
    output:
        flood_aar="data/results/flood_risk/countries/{ISO3}/{ISO3}_wri-flood-risk_AAR_V-{VULN_CURVE}.tif"
    wildcard_constraints:
        VULN_CURVE="BER|JRC|EXP"
    script:
        "./wri_average_annual_risk.py"
"""
Test with
snakemake -c1 data/results/flood_risk/countries/KEN/KEN_wri-flood-risk_AAR_V-JRC.tif
"""
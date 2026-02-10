"""
Run bulk analytics for each section in the academic paper
"""

configfile: "config/config.yaml"

"""
Section 1: Observed Flooding Analysis
"""

# Run observed metrics for all countries 
rule observed_metrics_for_all_countries:
    input:
        expand("data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_ADM0_metrics_gfd-flood_S-rwi.gpkg", ISO3=config['iso_codes'])

# Run decomposed observed metrics for all countries
rule observed_metrics_decomposed_for_all_countries:
    input:
        expand("data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_ADM0_decomposed_metrics_gfd-flood_S-rwi.gpkg", ISO3=config['iso_codes'])

# Run individual DFO event CI analysis

# Find all events in the prep folder
events = glob_wildcards("data/inputs/gfd/prep/DFO_{event_id}.tif").event_id

# Before running below rule run clip_gfd_event rule for all events in the prep gfd folder

# Run metrics analysis for all DFO events
rule metrics_all_gfd_events:
    input:
        expand("data/results/social_flood/events/DFO_{event_id}/DFO_{event_id}_results.csv", event_id=events)

"""
Section 2: Modelled Flooding Analysis
"""

MODELS = ['jrc', 'wri', 'giri']

# Run country level flood model CI metrics
rule flood_model_metrics_ADM0_all_countries:
    input:
        expand("data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_ADM0_metrics_{MODEL}-flood_protected_AAR_V-JRC_S-rwi.gpkg",
            ISO3=config['iso_codes'], MODEL=MODELS)

# Run country level admin 1 decomposed CI metrics
rule flood_model_admin_CI_decomposed:
    input:
        expand("data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_ADM1_admin-decomposed_metrics_{MODEL}-flood_protected_AAR_V-JRC_S-rwi.gpkg",
             ISO3=config['iso_codes'], MODEL=MODELS)

"""
Section 3: the outputs from Section 1 and 2 should be enough
"""

"""
Section 4: Flood Adaptation Analysis
"""

MODELS = ['jrc', 'wri', 'giri']
RPs = [100] # in paper we use 100
fp_urban = [30] # in paper we use 30 (cities)
rl_urban = [13] # in paper we use 13 (densest rural)

rule bulk_flood_risk_and_adaptation_analysis:
    input:
        expand("data/results/flood_risk/summary/countries/{ISO3}/{ISO3}_ADM0_metrics_{MODEL}-flood_AALs_baseline_capstock.gpkg",
                ISO3=config['iso_codes'], MODEL=MODELS),
        expand("data/results/flood_risk/summary/countries/{ISO3}/{ISO3}_ADM0_metrics_{MODEL}-flood_AALs_adapted_fp_rp{RP}_duc{urban}_capstock.gpkg",
                ISO3=config['iso_codes'], MODEL=MODELS, RP=RPs, urban=fp_urban),
        expand("data/results/flood_risk/summary/countries/{ISO3}/{ISO3}_ADM0_metrics_{MODEL}-flood_AALs_adapted_rl_duc{urban}_capstock.gpkg",
                ISO3=config['iso_codes'], MODEL=MODELS, urban=rl_urban),
        expand("data/results/flood_risk/summary/countries/{ISO3}/{ISO3}_ADM0_metrics_{MODEL}-flood_AALs_adapted_dp_capstock.gpkg",
                ISO3=config['iso_codes'], MODEL=MODELS),
        expand("data/results/adaptation/costs/countries/{ISO3}/{ISO3}_adaptation-cost_fp_rp{RP}_duc{urban}_ADM0.gpkg",
                ISO3=config['iso_codes'], RP=RPs, urban=fp_urban),
        expand("data/results/adaptation/costs/countries/{ISO3}/{ISO3}_adaptation-cost_rl_m-{MODEL}_duc{urban}_ADM0.gpkg",
                ISO3=config['iso_codes'], MODEL=MODELS, urban=rl_urban),
        expand("data/results/adaptation/costs/countries/{ISO3}/{ISO3}_adaptation-cost_dp_m-{MODEL}_ADM0.gpkg",
                ISO3=config['iso_codes'], MODEL=MODELS)

rule bulk_social_metrics_adaptation:
    input:
        expand("data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_ADM0_metrics_{MODEL}-flood_adapted_AAR_V-JRC_S-rwi_fp_rp{RP}_duc{urban}.gpkg",
            ISO3=config['iso_codes'], MODEL=MODELS, RP=RPs, urban=fp_urban),
        expand("data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_ADM0_metrics_{MODEL}-flood_adapted_AAR_V-JRC_S-rwi_dp.gpkg",
            ISO3=config['iso_codes'], MODEL=MODELS),
        expand("data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_ADM0_metrics_{MODEL}-flood_adapted_AAR_V-JRC_S-rwi_rl_duc{urban}.gpkg",
            ISO3=config['iso_codes'], MODEL=MODELS, urban=rl_urban),    

"""
Section Supplementary: Sensitivity Analysis 
"""

MODELS = ['jrc', 'wri', 'giri']
METRICS = ['RP100', 'AAR']
VULN_CURVES = ['JRC', 'EXP', 'BER']
SOCIALS = ['gdp', 'rwi']

rule flood_model_sensitivity_run:
    input:
        expand("data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_ADM0_metrics_{MODEL}-flood_{METRIC}_V-{VULN_CURVE}_S-{SOCIAL}.gpkg",
            ISO3=config['iso_codes'], MODEL=MODELS, METRIC=METRICS, VULN_CURVE=VULN_CURVES, SOCIAL=SOCIALS),
        expand("data/results/social_flood/countries/{ISO3}/inequality_metrics/{ISO3}_ADM0_metrics_{MODEL}-flood_protected_AAR_V-{VULN_CURVE}_S-{SOCIAL}.gpkg",
                    ISO3=config['iso_codes'], MODEL=MODELS, METRIC=METRICS, VULN_CURVE=VULN_CURVES, SOCIAL=SOCIALS)
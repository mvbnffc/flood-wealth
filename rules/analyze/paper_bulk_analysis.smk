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
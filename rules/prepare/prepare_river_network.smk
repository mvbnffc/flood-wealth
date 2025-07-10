"""
Prepare river network data by clipping to country and applying upstream drainage area threshold (500 km^2)
"""

rule prepare_river_network:
    """
    Load global river network, clip by boundary file, and adjust network to align with flood maps (500 km^2 upstream drainage area)
    """
    input:
        river_network_file="data/inputs/rivers/hydroRIVERS/HydroRIVERS_v10.shp",
        boundary_file="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
    output:
        river_file= "data/inputs/analysis/countries/{ISO3}/{ISO3}_river_network.gpkg"
    script:
        "./prepare_rivers.py"
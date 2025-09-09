# Download and extract openstreetmap data

#
# Download OpenStreetMap Planet file, create country file and sector extracts
#

rule download_osm:
    output:
        pbf=protected("data/inputs/openstreetmap/planet-250901.osm.pbf"),
    shell:
        """
        mkdir -p data/inputs/openstreetmap
        cd data/inputs/openstreetmap
        aws s3 cp --no-sign-request s3://osm-planet-eu-central-1/planet/pbf/2025/planet-250901.osm.pbf .
        aws s3 cp --no-sign-request s3://osm-planet-eu-central-1/planet/pbf/2025/planet-250901.osm.pbf.md5 .
        md5sum --check planet-250901.osm.pbf.md5
        """

rule extract_osm_data:
    input:
        pbf="data/inputs/openstreetmap/planet-250901.osm.pbf",
        json="data/inputs/boundaries/{ISO3}/geobounds_{ISO3}.geojson",
    output:
        pbf="data/inputs/analysis/countries/{ISO3}/{ISO3}_openstreetmap.osm.pbf",
    shell:
        """
        osmium extract \
            --polygon {input.json} \
            --set-bounds \
            --strategy=complete_ways \
            --output={output.pbf} \
            {input.pbf}
        """

rule filter_osm_data:
    input:
        pbf="data/inputs/analysis/countries/{ISO3}/{ISO3}_openstreetmap.osm.pbf",
    output:
        pbf="data/inputs/analysis/countries/{ISO3}/{ISO3}_osm_roads.osm.pbf",
    shell:
        """
        osmium tags-filter \
            --expressions=config/osm-roads.txt \
            --output={output.pbf} \
            {input.pbf}
        """

rule convert_osm_data:
    input:
        pbf="data/inputs/analysis/countries/{ISO3}/{ISO3}_osm_roads.osm.pbf",
    output:
        gpkg="data/inputs/analysis/countries/{ISO3}/{ISO3}_osm_roads.gpkg",
    shell:
        """
        OSM_CONFIG_FILE=config/osm-roads.conf.ini ogr2ogr -f GPKG -overwrite {output.gpkg} {input.pbf}
        """
# Download and extract openstreetmap data

rule download_osm_zip:
    """
    Download the Geofabrik OSM Shapefile ZIP for a given ISO3.
    Expects lines like: RWA, https://download.geofabrik.de/africa/rwanda-latest-free.shp.zip
    """
    output:
        "data/inputs/analysis/countries/{ISO3}/{ISO3}_osm.zip"
    shell:
        r"""
        set -euo pipefail
        out="{output}"
        iso="{wildcards.ISO3}"
        mkdir -p "$(dirname "$out")"

        # Grab URL from config/osm_links.txt (case-insensitive ISO3, trims spaces)
        url="$(awk -F',' -v iso="$iso" '
            BEGIN{{IGNORECASE=1}}
            toupper($1)==toupper(iso) {{
                gsub(/^[ \t]+|[ \t]+$/, "", $2);
                print $2;
                exit
            }}' config/osm_links.txt)"

        if [ -z "$url" ]; then
          echo "No URL found for ISO3=$iso in config/osm_links.txt" >&2
          exit 1
        fi

        # Skip if already downloaded; otherwise fetch
        if [ ! -f "$out" ]; then
          wget -q -O "$out" "$url"
        fi
        """
"""
Test with
snakemake -c1 data/inputs/analysis/countries/RWA/RWA_osm.zip
"""

rule extract_osm_zip:
    """
    Unzip the downloaded OSM Shapefile bundle into a clean folder.
    """
    input:
        zipfile = "data/inputs/analysis/countries/{ISO3}/{ISO3}_osm.zip"
    output:
        directory("data/inputs/analysis/countries/{ISO3}/{ISO3}_osm")
    shell:
        r"""
        mkdir -p {output}
        unzip -oq {input.zipfile} -d {output}
        """

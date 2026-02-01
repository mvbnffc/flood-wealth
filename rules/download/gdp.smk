"""
Download gridded GDP data

References
---------
GDP: https://www.nature.com/articles/s41597-022-01322-5
"""

rule download_gdp:
    output:
        "data/inputs/gdp/2019gdp.tif"
    params:
        api_url="https://api.figshare.com/v2/file/download/31456837",
        top_zip="Real_GDP.zip",
        nested_zip="updated real GDP/2019.zip",
        nested_member="2019GDP.tif"
    shell:
        r"""
        set -euo pipefail

        output_dir=$(dirname {output})
        mkdir -p "$output_dir"

        outer="$output_dir/{params.top_zip}"
        tmp_nested="$output_dir/tmp_2019.zip"

        # Download outer zip via Figshare API (bypasses AWS WAF)
        curl -L -o "$outer" "{params.api_url}"

        # Basic sanity check
        file "$outer" | grep -qi zip

        # Extract nested zip
        unzip -p "$outer" "{params.nested_zip}" > "$tmp_nested"

        # Extract final GeoTIFF
        unzip -p "$tmp_nested" "{params.nested_member}" > "{output}"

        rm -f "$tmp_nested"
        """

"""
Download Meta Relative Wealth Index

References
---------
RWI: https://data.humdata.org/dataset/relative-wealth-index
"""

rule download_rwi:
    """
    Meta Relative Wealth Index

    Description:
    https://dataforgood.facebook.com/dfg/tools/relative-wealth-index

    Catalogue page: https://data.humdata.org/dataset/relative-wealth-index


    Citation:

    Microestimates of wealth for all low- and middle-income countries Guanghua
    Chi, Han Fang, Sourav Chatterjee, Joshua E. Blumenstock Proceedings of the
    National Academy of Sciences Jan 2022, 119 (3) e2113658119; DOI:
    10.1073/pnas.2113658119

    License: CC-BY-NC 4.0
    """
    output:
        raw_folder = directory("data/inputs/rwi/raw/"),
    run:
        import os
        from pathlib import Path

        output_dir = Path(output.raw_folder)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Loop over rwi file links and add to list
        with open("config/rwi_links.txt") as f:
            urls = [line.strip() for line in f if line.strip()]
        
        for url in urls:
            filename = url.split("/")[-1]
            target_path = output_dir / filename

            # Only download if file doesn't exist
            if not target_path.exists():
                shell(f"wget -nc {url} -O {target_path}")
            else:
                print(f"{filename} already exists - skipping.")

rule combine_rwi:
    """
    Combine all RWI CSVs into one standardized file.
    Removes "quadkey" column if it exists
    """
    input:
        raw_folder = "data/inputs/rwi/raw/"
    output:
        combined = "data/inputs/rwi/relative_wealth_index_combined.csv"
    run:
        import pandas as pd
        from pathlib import Path

        input_dir = Path(input.raw_folder)
        csv_files = sorted(input_dir.glob("*_relative_wealth_index.csv"))

        dfs = []

        for csv_file in csv_files:
            df = pd.read_csv(csv_file)

            # Drop 'quadkey' column if it exists
            if "quadkey" in df.columns:
                df = df.drop(columns=["quadkey"])
            
            dfs.append(df)
        
        combined_df = pd.concat(dfs, ignore_index=True)
        combined_df.to_csv(output.combined, index=False)

rule process_rwi_points:
    input:
        csv = "data/inputs/rwi/relative_wealth_index_combined.csv",
    output:
        gpkg = "data/inputs/rwi/relative_wealth_index_combined.gpkg",
    run:
        import pandas as pd
        import geopandas as gpd
        df = pd.read_csv(input.csv)
        gdf = gpd.GeoDataFrame(
            data=df[["rwi", "error"]],
            geometry=gpd.points_from_xy(df.longitude, df.latitude),
            crs="EPSG:4326"
        )
        # reproject to 3857 (point coordinates are provided in lat/lon but the
        # underlying regular grid seems to be based on Web Mercator via Quadkeys)
        gdf.to_crs("EPSG:3857", inplace=True)
        gdf.to_file(output.gpkg, engine="pyogrio")


rule process_rwi_grid:
    """Rasterise points to a grid data format

    -tr specifies x/y resolution, guessed from most common difference between point
    coordinates:

    def guess_resolution(coords):
        coords.sort()
        diffs = np.diff(coords)
        diffs.sort()
        # vals, counts = np.unique(diffs, return_counts=True)
        counts, vals = np.histogram(diffs)
        return vals[np.argmax(counts)]
    print("x res:", guess_resolution(gdf.geometry.x.unique())
    print("y res:", guess_resolution(gdf.geometry.y.unique())
    """
    input:
        gpkg = "data/inputs/rwi/relative_wealth_index_combined.gpkg",
    output:
        tiff = "data/inputs/rwi/rwi.tif",
    shell:
        """
        gdal_rasterize \
            -a rwi \
            -init -999 \
            -a_nodata -999 \
            -tr 2445.9786434 2445.96770335 \
            -ot Float64 \
            -of GTiff \
            {input.gpkg} \
            {output.tiff}
        """

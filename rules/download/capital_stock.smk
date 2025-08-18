"""
Download IMF capital stock data

References
---------
IMF capital stock: https://infrastructuregovern.imf.org/content/dam/PIMA/Knowledge-Hub/dataset/IMFInvestmentandCapitalStockDataset2021.xlsx
"""

rule download_IMF_capital_stock:
    output:
        capital_stock = "data/inputs/imf/IMFInvestmentandCapitalStockDataset2021.xlsx"
    shell:
        """
        mkdir -p $(dirname {output.capital_stock})
        wget -nc https://infrastructuregovern.imf.org/content/dam/PIMA/Knowledge-Hub/dataset/IMFInvestmentandCapitalStockDataset2021.xlsx -O {output.capital_stock}
        """
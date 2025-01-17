import os
import re
import sqlite3
import pandas as pd


DB_PATH = "resource/swgdrugdb.db"


def parseSearchTerm(searchTerm: str) -> str:

    conn   = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()


    if (not re.fullmatch(r"CL\d+", searchTerm)):
        
        cursor.execute("SELECT craftsLabEntry FROM idxTable WHERE SMILES = ?",
                       (searchTerm,))
        
        craftsLabsEntry = cursor.fetchone()[0]
    else:

        craftsLabsEntry = searchTerm

    return craftsLabsEntry


def download(searchTerm: str) -> None:

    allFragsDF    = searchAndFetch(searchTerm)
    downloadsPath = os.path.join(os.path.expanduser("~"), "Downloads")
    fullPath      = os.path.join(downloadsPath, f"{parseSearchTerm(searchTerm)}.csv")

    allFragsDF.to_csv(fullPath, index=False)


def searchAndFetch(searchTerm: str) -> pd.DataFrame:

    conn           = sqlite3.connect(DB_PATH)
    cursor         = conn.cursor()
    craftsLabEntry = parseSearchTerm(searchTerm)
    allFragsDF     = pd.read_sql(f"SELECT * FROM {craftsLabEntry}", conn)

    cursor.close()
    conn.close()

    return allFragsDF


def viewIdxTable() -> pd.DataFrame:

    conn         = sqlite3.connect(DB_PATH)
    idxTableAsDF = pd.read_sql(f"SELECT * FROM idxTable", conn)

    conn.close()

    return idxTableAsDF
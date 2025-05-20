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

        craftsLabEntry = cursor.fetchone()[0]

    else:

        craftsLabEntry = searchTerm

    cursor.close()
    conn.close()

    return craftsLabEntry


def searchAndFetch(searchTerm: str) -> pd.DataFrame:

    try:

        conn           = sqlite3.connect(DB_PATH)
        cursor         = conn.cursor()
        craftsLabEntry = parseSearchTerm(searchTerm)
        allFragsDF     = pd.read_sql(f"SELECT * FROM {craftsLabEntry}", conn)

        return allFragsDF

    except Exception as e:

        raise e

    finally:

        cursor.close()
        conn.close()


def viewIdxTable() -> pd.DataFrame:

    conn         = sqlite3.connect(DB_PATH)
    idxTableAsDF = pd.read_sql(f"SELECT * FROM idxTable", conn)

    conn.close()

    return idxTableAsDF
from fragmentor import Fragmentor
import os
import re
import sqlite3
import pandas as pd
from rdkit import Chem


DB_PATH = "resource/swgdrugdb.db"


def addEntryFromSmiles(smiles: str) -> None:

    conn   = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    cursor.execute("SELECT craftsLabEntry FROM idxTable WHERE SMILES = ?",
                   (smiles,))

    res = cursor.fetchone()

    if (res is not None):

        print(f"<*> {smiles} already exists with craftsLabEntry: {res[0]}")

        conn.close()

        return

    mol = Chem.MolFromSmiles(smiles)

    if (mol is None):

        print(f"<*> {smiles} is not a valid SMILES string")

        conn.close()

        return

    cursor.execute("INSERT INTO idxTable (SMILES) VALUES (?)", (smiles,))
    conn.commit()

    cursor.execute("SELECT craftsLabEntry FROM idxTable WHERE SMILES = ?",
                   (smiles,))

    craftsLabEntry = cursor.fetchone()[0]

    fragmentor     = Fragmentor()
    fragmentor.mol = mol
    allFragsDF     = fragmentor.fetchAllFragsData()

    allFragsDF.to_sql(craftsLabEntry, conn, index=False)

    print(f"{craftsLabEntry} successfully added to the MFD.")

    conn.commit()
    cursor.close()
    conn.close()


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
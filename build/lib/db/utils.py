from fragmentor import Fragmentor
import os
import re
import sqlite3
import pandas as pd
from rdkit import Chem


##########
# CONSTS #
##########

CURRENT_DIR     = os.path.dirname(os.path.abspath(__file__))
DB_PATH         = os.path.abspath(
    os.path.join(CURRENT_DIR, "..", "..", "resource", "swgdrugdb.db")
)

##########


def addEntryFromSmiles(smiles: str) -> None:

    if (not smiles):

        raise ValueError("<!> Error: SMILES string must not be empty.")

    mol = Chem.MolFromSmiles(smiles)

    if (mol is None):

        raise ValueError(f"<!> Error: {smiles} is not a valid SMILES string.")

    with sqlite3.connect(DB_PATH) as conn:

        cursor = conn.cursor()

        cursor.execute("SELECT craftsLabEntry FROM idxTable WHERE SMILES = ?",
                       (smiles,))

        res = cursor.fetchone()

        if (res is not None):

            print(f"<*> {smiles} already exists with craftsLabEntry: {res[0]}")

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

        print(f"<*> {craftsLabEntry} was successfully added to the MFD.")


def parseSearchTerm(searchTerm: str) -> str:

    with sqlite3.connect(DB_PATH) as conn:

        cursor = conn.cursor()

        if (not re.fullmatch(r"CL\d+", searchTerm)):
            
            cursor.execute("SELECT craftsLabEntry FROM idxTable WHERE SMILES = ?",
                           (searchTerm,))

            result = cursor.fetchone()

            if (not result):

                raise ValueError(f"<!> Error: {searchTerm} cannot be found.")

            craftsLabEntry = result[0]

            return craftsLabEntry

        craftsLabEntry = searchTerm

        return craftsLabEntry


def searchAndFetch(searchTerm: str) -> pd.DataFrame:

    with sqlite3.connect(DB_PATH) as conn:

        try:

            craftsLabEntry = parseSearchTerm(searchTerm)
            allFragsDF     = pd.read_sql(f"SELECT * FROM {craftsLabEntry}", conn)

            return allFragsDF

        except Exception:

            raise


def viewIdxTable() -> pd.DataFrame:

    with sqlite3.connect(DB_PATH) as conn:

        idxTableAsDF = pd.read_sql(f"SELECT * FROM idxTable", conn)

        return idxTableAsDF
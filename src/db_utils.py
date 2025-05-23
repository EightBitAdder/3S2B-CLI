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
    os.path.join(CURRENT_DIR, "..", "resource", "swgdrugdb.db")
)
##########


def addEntryFromSmiles(smiles: str) -> None:
    """
    Adds an entry to the Molecular Fragment Database (M.F.D.)
    by inserting it into the Index Table via a SMILES string, as well as
    the corresponding Fragment Table.


    params:
    -------

    smiles:
        The SMILES string of the molecule to be inserted.


    returns:
    --------

    None
    """

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

        print(f"{craftsLabEntry} successfully added to the MFD.")


def parseSearchTerm(searchTerm: str) -> str:
    """
    Parses the search term from either:

        i) A Crafts Lab Entry in the form --- CL# ---; or,

        ii) A SMILES string,

    to the corresponding Crafts Lab Entry.


    params:
    -------

    searchTerm: str
        The search term to be parsed.


    returns:
    --------

    craftsLabEntry: str
        The Crafts Lab Entry.
    """

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
    """
    Fetches the corresponding fragment table for the given search term.


    params:
    -------

    searchTerm: str
        The name of the corresponding fragment table for the given search term
        to be fetched.


    returns:
    --------
    
    allFragsDF: pd.DataFrame
        The dataframe of the fragment table.
    """

    with sqlite3.connect(DB_PATH) as conn:

        try:

            craftsLabEntry = parseSearchTerm(searchTerm)
            allFragsDF     = pd.read_sql(f"SELECT * FROM {craftsLabEntry}", conn)

            return allFragsDF

        except Exception:

            raise


def viewIdxTable() -> pd.DataFrame:
    """
    Fetches the Index Table.


    params:
    -------

    None


    returns:
    --------

    idxTableAsDF: pd.DataFrame
        The Index Table as a pandas dataframe object.
    """

    with sqlite3.connect(DB_PATH) as conn:

        idxTableAsDF = pd.read_sql(f"SELECT * FROM idxTable", conn)

        return idxTableAsDF
from src.fragmentor import Fragmentor
import os
import sqlite3
from rdkit import Chem
from rdkit import RDLogger
from tqdm import tqdm


RDLogger.DisableLog('rdApp.*')


CURRENT_DIR     = os.path.dirname(os.path.abspath(__file__))
DB_PATH         = os.path.abspath(os.path.join(CURRENT_DIR, "..", "..", "resource", "swgdrugdb.db"))
SWGDRUG313_PATH = os.path.abspath(os.path.join(CURRENT_DIR, "..", "..", "resource", "SWGDRUG313.SDF"))


def makeDB() -> None:

    conn   = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS idxTable (
            SMILES TEXT PRIMARY KEY,
            craftsLabEntry TEXT
        );
    """)

    cursor.execute("""
        CREATE TRIGGER IF NOT EXISTS autoIncrementCraftsLabEntry
        AFTER INSERT ON idxTable
        FOR EACH ROW
        WHEN NEW.craftsLabEntry IS NULL
        BEGIN
            UPDATE idxTable
            SET craftsLabEntry = 'CL' || NEW.rowid
            WHERE rowid = NEW.rowid;
        END;
    """)

    suppl = Chem.SDMolSupplier(SWGDRUG313_PATH)

    for mol in tqdm(suppl, desc="Processing SDF . . ."):

        if (not mol):

            continue
        
        
        smiles         = Chem.MolToSmiles(mol)
        fragmentor     = Fragmentor()
        fragmentor.mol = mol
        allFragsDF     = fragmentor.fetchAllFragsData()

        cursor.execute("SELECT craftsLabEntry FROM idxTable WHERE SMILES = ?", (smiles,))

        row = cursor.fetchone()

        if (row):

            craftsLabEntry = row[0]

        else:

            cursor.execute("INSERT INTO idxTable (SMILES) VALUES (?)", (smiles,))
            conn.commit()
            cursor.execute("SELECT craftsLabEntry FROM idxTable WHERE SMILES = ?", (smiles,))

            craftsLabEntry = cursor.fetchone()[0]
        
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name = ?", (craftsLabEntry,))
        
        if (cursor.fetchone() is None):

            if (not allFragsDF.empty):

                try:
            
                    allFragsDF.to_sql(craftsLabEntry, conn, index=False)

                    tqdm.write(f"<*> Created fragment table: {craftsLabEntry}.")

                except Exception as e:

                    tqdm.write(f"<*> Failed to create fragment table: {craftsLabEntry} >>> {e}")

            else:

                tqdm.write(f"<*> Skipped empty fragment table: {craftsLabEntry}.")

    conn.commit()
    cursor.close()
    conn.close()

if __name__ == "__main__":

    makeDB()
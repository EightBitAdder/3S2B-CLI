from src.fragmentor import Fragmentor
import sqlite3
from rdkit import Chem, RDLogger


DB_PATH = "resource/swgdrugdb.db"
SWGDRUG313_PATH = "resource/SWGDRUG313.SDF"


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

    for mol in suppl:

        if (not mol):

            continue
        
        
        smiles         = Chem.MolToSmiles(mol)
        fragmentor     = Fragmentor()
        fragmentor.mol = mol
        allFragsDF     = fragmentor.fetchAllFragsData()
        
        cursor.execute("INSERT OR IGNORE INTO idxTable (SMILES) VALUES (?)",
                       (smiles,))
        
        if (cursor.rowcount > 0):
            
            allFragsDF.to_sql(f"CL{cursor.lastrowid}", conn, index=False)

    conn.commit()
    cursor.close()
    conn.close()


RDLogger.DisableLog("rdApp.*") # rdkit warnings are really annoying!
makeDB()
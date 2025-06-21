from db.utils import (
    viewIdxTable,
    searchAndFetch
)
from main.comparator import (
    MassList,
    MSData,
    Comparator
)
import os
import pandas as pd
import numpy as np
from tqdm import tqdm


def fetchMassList(allFragsDF: pd.DataFrame) -> MassList:

    massList = []

    for row in allFragsDF.itertuples():

        Iso_Wts = list(map(float, row.Iso_Wts.split(", ")))
        masses  = [row.Exact_Mol_Wt, *Iso_Wts]

        massList.extend([mass, row.Mult] for mass in masses)

    df = pd.DataFrame(massList, columns=["Mass", "Multiplicity"])

    return MassList(df, weighted=True)


def searchAndFetchByMass(mz: float) -> pd.DataFrame:

    mz           = round(float(mz), 2)
    idxTableDF   = viewIdxTable()
    # filteredRows = []
    rows = []

    for craftsLabEntry in tqdm(idxTableDF.iloc[:, 0], desc="<*> Fetching fragments . . ."):

        for row in searchAndFetch(craftsLabEntry).itertuples():

            Iso_Wts = list(map(float, row.Iso_Wts.split(", ")))
            masses  = [row.Exact_Mol_Wt, *Iso_Wts]

            if (mz in np.round(masses, 2)):

                rows.append(row)

    return pd.DataFrame(rows)
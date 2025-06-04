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
    filteredRows = []

    for idx, entry in tqdm(idxTableDF.iloc[:, 1].items(), desc="<*> Fetching fragments . . ."):

        massList = fetchMassList(searchAndFetch(entry))

        if (mz in np.round(massList.masses, 2)):

            filteredRows.append(idx)
            tqdm.write(f"<*> Found a match!")

    return idxTableDF.loc[filteredRows].reset_index(drop=True)


def compare(msDataPath: str, massListPath: str, *, tol: float) -> float:
    
    msData                  = MSData.fromFile(msDataPath)
    massList                = MassList.fromFile(massListPath)
    comparator              = Comparator(msData, massList, tol)
    FPIEScore, plotMetaData = comparator.calculateFPIE()

    comparator.plotFPIE(plotMetaData, FPIEScore, os.path.splitext(os.path.basename(msDataPath))[0])

    return FPIEScore


def compareAll(msDataPath: str, *, tol: float) -> pd.DataFrame:

    msData          = MSData.fromFile(msDataPath)
    idxTableDF      = viewIdxTable()
    craftsLabEntrys = idxTableDF.iloc[:, 1]
    FPIEs           = []

    for entry in tqdm(craftsLabEntrys, desc=f"<*> Calculating FPIEs . . ."):

        massList     = fetchMassList(searchAndFetch(entry))
        comparator   = Comparator(msData, massList, tol)
        FPIEScore, _ = comparator.calculateFPIE()

        FPIEs.append(FPIEScore)

    df = pd.concat([idxTableDF.iloc[:, 0], pd.DataFrame({"FPIE": FPIEs})], axis=1)

    return df
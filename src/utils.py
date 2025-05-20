from db_utils import viewIdxTable, searchAndFetch
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm


def fetchMassList(allFragsDF: pd.DataFrame) -> np.ndarray:

    massList = []

    for row in allFragsDF.itertuples():

        Iso_Wts = list(map(float, row.Iso_Wts.split(", ")))

        massList.extend([row.Exact_Mol_Wt, *Iso_Wts, row.Mult])

    return np.array(massList)


def isHeader(row: pd.Series) -> bool:

    return row.apply(lambda x: isinstance(x, str)).all()


def checkFile(filePath: str, *, expectedCols: int) -> pd.DataFrame:

    if (not os.path.exists(filePath)):

        raise FileNotFoundError(f"<!> Error in reading {filePath}.")
    
    if (not (filePath.endswith(".csv") or filePath.endswith(".txt"))):

        raise ValueError(f"<!> Error in reading non-.csv or -.txt files.")
    
    # TODO:
    # Write helper method for programmatically determining sep
    tempDF = pd.read_csv(filePath, header=None, sep=" ")

    if (isHeader(tempDF.iloc[0])):
        
        df = tempDF[1:]

    else:

        df = tempDF

    if (df.shape[1] != expectedCols):

        print(df)

        raise ValueError(f"<!> Error in parsing {filePath}.")
    
    if (not np.issubdtype(df.values.dtype, np.number)):

        raise ValueError(f"<!> Error in parsing {filePath}.")
    
    return df


def compare(msDataPath: str, massListPath: str, *, tol: float) -> float:
    
    msData   = checkFile(msDataPath, expectedCols=2)
    massList = checkFile(massListPath, expectedCols=1)
    
    return FPIE(msData.to_numpy(), massList.to_numpy(), tol=tol, plotting=True)


def compareAll(msDataPath: str, *, tol: float) -> pd.DataFrame:

    msData          = checkFile(msDataPath, expectedCols=2)
    FPIEs           = []
    idxTableDF      = viewIdxTable()
    craftsLabEntrys = idxTableDF.iloc[:, 1]

    for entry in tqdm(craftsLabEntrys, desc=f"<*> Calculating FPIEs..."):

        massList = fetchMassList(searchAndFetch(entry))

        FPIEs.append(FPIE(msData.to_numpy(), massList, tol=tol, weighted=True))

    return pd.concat([idxTableDF.iloc[:, 0], pd.DataFrame({"FPIE": FPIEs})],
                     axis=1)


def FPIE(msData: np.ndarray, massList: np.ndarray, *, tol: float, weighted: bool=False, plotting: bool=False) -> float:
    """
        <*> Reminder:
            tol must not be 0 if msData is high-res...
    """

    if (massList.ndim == 1):

        massList = np.column_stack((massList, np.ones_like(massList)))

    elif (massList.ndim == 2):

        if (massList.shape[1] != 2):

            raise ValueError(f"<!> Error in parsing mass list with shape: {massList.shape}." \
                             f"Expected (n,) or (n, 2).")

    else:

        raise ValueError(f"<!> Error in parsing mass list with shape: {massList.shape}.")

    masses = massList[:, 0]
    mults  = massList[:, 1]


    massesMultDict = {}

    for mass, mult in zip(masses, mults):

        if ((mass not in massesMultDict) or (mult < massesMultDict[mass])):

            massesMultDict[mass] = mult

    masses   = np.array(list(massesMultDict.keys()))
    weights  = np.array([
        np.exp(-(massesMultDict[m] - 1)) if weighted else 1.0
        for m in masses
    ])


    msDataMasses = np.unique(msData[:, 0])

    if (tol == 0):

        masses = masses.astype(int)

    else:

        masses       = np.round(masses, 2)
        msDataMasses = np.round(msDataMasses, 2)

    commonMasses   = []
    matchedWeights = {}

    for msDataMass in msDataMasses:

        for i, m in enumerate(masses):

            if (m - tol <= msDataMass <= m + tol):

                commonMasses.append(msDataMass)
                matchedWeights[msDataMass] = weights[i]
                break

    explainedIntensities = []

    for mz, intensity in msData:

        if (mz in matchedWeights):

            explainedIntensities.append(intensity * matchedWeights[mz])

    msDataIntensities = msData[:, 1]
    FPIEScore         = float(sum(explainedIntensities) / np.sum(msDataIntensities))

    if (plotting):

        for i in range(len(msDataMasses)):

            if (msDataMasses[i] in commonMasses):
                
                plt.stem([msDataMasses[i]],
                         [msDataIntensities[i]],
                         linefmt="r-",
                         markerfmt=" ",
                         basefmt=" ")
            else:

                plt.stem([msDataMasses[i]],
                         [msDataIntensities[i]],
                         linefmt="b-",
                         markerfmt=" ",
                         basefmt=" ")

        plt.xlabel("m/z")
        plt.ylabel("Intensity")
        plt.show()

    return FPIEScore
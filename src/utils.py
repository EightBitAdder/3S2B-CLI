from db_utils import viewIdxTable, searchAndFetch
import os
import numpy as np
import pandas as pd
from matplotlib import cm
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from tqdm import tqdm


def fetchMassList(allFragsDF: pd.DataFrame) -> np.ndarray:

    massList = []

    for row in allFragsDF.itertuples():

        Iso_Wts = list(map(float, row.Iso_Wts.split(", ")))
        masses  = [row.Exact_Mol_Wt, *Iso_Wts]

        massList.extend([mass, row.Mult] for mass in masses)

    return np.array(massList)


def isHeader(row: pd.Series) -> bool:

    return row.apply(lambda x: isinstance(x, str)).all()


def checkFile(filePath: str, *, expectedCols: int | None=None, isMassList: bool=False) -> pd.DataFrame:

    if (not os.path.exists(filePath)):

        raise FileNotFoundError(f"<!> Error in reading {filePath}.")
    
    if (not (filePath.endswith(".csv") or filePath.endswith(".txt"))):

        raise ValueError(f"<!> Error in reading non-.csv or -.txt files.")
    
    # TODO:
    # Write helper method for programmatically determining sep
    tempDF = pd.read_csv(filePath, header=None, sep=" ")

    if (isHeader(tempDF.iloc[0])):
        
        df = tempDF[1:].copy()

    else:

        df = tempDF

    if (isMassList):

        if (df.shape[1] not in {1, 2}):

            raise ValueError(f"<!> Error parsing {filePath} as mass list.")

        df[0] = pd.to_numeric(df[0], errors="raise")

        if (df.shape[1] == 2):

            df[1] = pd.to_numeric(df[1], errors="raise", downcast="integer")

    elif (expectedCols is not None):

        if (df.shape[1] != expectedCols):

            print(df)

            raise ValueError(f"<!> Error in parsing {filePath}.")

    if (not np.issubdtype(df.values.dtype, np.number)):

        raise ValueError(f"<!> Error in parsing {filePath}.")
    
    return df


def compare(msDataPath: str, massListPath: str, *, tol: float) -> float:
    
    msData   = checkFile(msDataPath, expectedCols=2)
    massList = checkFile(massListPath, isMassList=True)
    
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


def FPIE(msData: np.ndarray,
         massList: np.ndarray,
         *,
         tol: float,
         weighted: bool=False,
         plotting: bool=False) -> float:
    """
        <*> Reminder:
            tol must not be 0 if msData is high-res...
    """

    if (massList.ndim == 1):

        massList = np.column_stack((massList, np.ones_like(massList)))

    masses       = massList[:, 0]
    mults        = massList[:, 1]
    msDataMasses = np.unique(msData[:, 0])

    if (tol == 0):

        masses = masses.astype(int)

    else:

        masses       = np.round(masses, 2)
        msDataMasses = np.round(msDataMasses, 2)

    massesMultDict = {}

    for mass, mult in zip(masses, mults):

        if ((mass not in massesMultDict) or (mult < massesMultDict[mass])):

            massesMultDict[mass] = mult

    masses         = np.array(list(massesMultDict.keys()))
    weights        = np.array([
        np.exp(-(massesMultDict[m] - 1)) if weighted else 1.0
        for m in masses
    ])
    matchedMults   = {}
    matchedWeights = {}

    for msDataMass in msDataMasses:

        for i, m in enumerate(masses):

            if (m - tol <= msDataMass <= m + tol):

                matchedMults[msDataMass]   = massesMultDict[m]
                matchedWeights[msDataMass] = weights[i]

                break

    explainedIntensities = []

    for mz, intensity in msData:

        if (mz in matchedWeights):

            explainedIntensities.append(intensity * matchedWeights[mz])

    msDataIntensities = msData[:, 1]
    FPIEScore         = float(sum(explainedIntensities) / np.sum(msDataIntensities))

    if (plotting):

        fig, ax      = plt.subplots()
        mults        = sorted(set(massesMultDict.values()))
        numColors    = len(mults)
        colorMap     = cm.get_cmap("tab20", numColors)
        multColorMap = {mult: colorMap(i) for i, mult in enumerate(mults)}

        for i, mz in enumerate(msDataMasses):

            intensity = msDataIntensities[i]

            if (mz in matchedWeights):

                mult  = matchedMults.get(mz, 1)
                color = multColorMap[mult]

            else:

                color = "lightgray"

            ax.stem([mz], [intensity],
                     linefmt=color,
                     markerfmt=" ",
                     basefmt=" ")

        ax.set_xlabel("m/z")
        ax.set_ylabel("Intensity")
        ax.set_title(f"FPIE Annotated Plot")

        legendHandles = [
            Patch(color=multColorMap[mult], label=f"mult = {int(mult)}")
            for mult in mults
        ]

        ax.legend(handles=legendHandles, title="Multiplicity")

        plt.show()

    return FPIEScore
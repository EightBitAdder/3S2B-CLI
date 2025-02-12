from db_utils import viewIdxTable, searchAndFetch
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def fetchMassList(allFragsDF: pd.DataFrame) -> np.ndarray:

    massList = []

    for row in allFragsDF.itertuples():

        Iso_Wts = list(map(float, row.Iso_Wts.split(", ")))

        massList.extend([row.Exact_Mol_Wt, *Iso_Wts])

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


def compare(msDataPath: str, massListPath: str, tol: float) -> float:
    
    msData = checkFile(msDataPath, expectedCols=2)
    massList = checkFile(massListPath, expectedCols=1)
    
    return FPIE(msData.to_numpy(), massList.to_numpy(), tol=tol, plotting=True)


def compareAll(msDataPath: str, tol: float) -> pd.DataFrame:

    msData          = checkFile(msDataPath, expectedCols=2)
    FPIEs           = []
    idxTableDF      = viewIdxTable()
    craftsLabEntrys = idxTableDF.iloc[:, 1]

    for entry in craftsLabEntrys:

        massList = fetchMassList(searchAndFetch(entry))

        FPIEs.append(FPIE(msData.to_numpy(), massList, tol=tol))

    return pd.concat([idxTableDF.iloc[:, 0], pd.DataFrame({"FPIE": FPIEs})],
                     axis=1)


def FPIE(msData: np.ndarray, massList: np.ndarray, *, tol: float, plotting: bool=False) -> float:

    msDataMasses = np.unique(msData[:, 0])

    if (tol == 0):

        massList = massList.astype(int)

    else:

        massList = np.round(massList, 5)
        msDataMasses = np.round(msDataMasses, 5)

    massList = np.unique(massList)
    commonMasses = []

    for msDataMass in msDataMasses:

        for m in massList:

            if (m - tol <= msDataMass <= m + tol):

                commonMasses.append(msDataMass)

    explainedIntensities = []

    for row in msData:

        if (row[0] in commonMasses):

            explainedIntensities.append(row[1])

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
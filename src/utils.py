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


def compare(msDataPath: str, massListPath: str) -> float:
    
    if (not os.path.exists(msDataPath) or not os.path.exists(massListPath)):

        raise FileNotFoundError(f"<!> Error in reading {msDataPath} and {massListPath}.")
    
    if (not (msDataPath.endswith(".csv") or msDataPath.endswith(".txt")) or not (massListPath.endswith(".csv") or massListPath.endswith(".txt"))):

        raise ValueError(f"<!> Error in reading non-.csv or -.txt files.")
    
    # TODO:
    # Write helper method for programmatically determining sep
    tempMSData = pd.read_csv(msDataPath, header=None, sep=" ")

    if (isHeader(tempMSData.iloc[0])):
        
        msData = tempMSData[1:]

    else:

        msData = tempMSData

    if (msData.shape[1] != 2):

        print(msData)

        raise ValueError(f"<!> Error in parsing {msDataPath}.")
    
    if (not np.issubdtype(msData.values.dtype, np.number)):

        raise ValueError(f"<!> Error in parsing {msDataPath}.")
    
    tempMassList = pd.read_csv(massListPath, header=None)

    if (isHeader(tempMassList.iloc[0])):

        massList = tempMassList[1:]

    else:

        massList = tempMassList

    if (massList.shape[1] != 1):

        raise ValueError(f"<!> Error in parsing {massListPath}.")
    
    if (not np.issubdtype(massList.values.dtype, np.number)):

        raise ValueError(f"<!> Error in parsing {massListPath}.")
    
    return FPIE(msData.to_numpy(), massList.to_numpy(), plotting=True)


def compareAll(msDataPath: str) -> str:

    if (not os.path.exists(msDataPath)):

        raise FileNotFoundError(f"<!> Error in reading {msDataPath}.")
    
    if (not (msDataPath.endswith(".csv") or msDataPath.endswith(".txt"))):

        raise ValueError(f"<!> Error in reading non-.csv or -.txt files.")
    
    # TODO:
    # Write helper method for programmatically determining sep
    tempMSData = pd.read_csv(msDataPath, header=None, sep=" ")

    if (isHeader(tempMSData.iloc[0])):

        msData = tempMSData[1:]

    else:

        msData = tempMSData

    if (msData.shape[1] != 2):

        print(msData)

        raise ValueError(f"<!> Error in parsing {msDataPath}.")
    
    if (not np.issubdtype(msData.values.dtype, np.number)):

        raise ValueError(f"<!> Error in parsing {msDataPath}.")
    
    FPIEs           = []
    idxTableDF      = viewIdxTable()
    craftsLabEntrys = idxTableDF.iloc[:, 1]

    for entry in craftsLabEntrys:

        massList = fetchMassList(searchAndFetch(entry))

        FPIEs.append(FPIE(msData.to_numpy(), massList))

    return pd.concat([idxTableDF.iloc[:, 0], pd.DataFrame({"FPIE": FPIEs})],
                     axis=1)


def FPIE(msData: np.ndarray, massList: np.ndarray, *, tol: float=0, plotting: bool=False) -> float:

    if (tol == 0):

        massList = massList.astype(int)

    msDataMasses = np.unique(msData[:, 0])
    massList = np.unique(massList)
    commonMasses = []

    for msDataMass in msDataMasses:

        for m in massList:

            if (m + tol < msDataMass):

                continue

            elif(m - tol >= msDataMass):

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
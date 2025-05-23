from db_utils import viewIdxTable, searchAndFetch
import os
from enum import Enum, auto
from dataclasses import dataclass
import csv
import numpy as np
import pandas as pd
from matplotlib import cm
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from tqdm import tqdm


class FileType(Enum):
    """
    """

    MS_DATA               = auto()
    MASS_LIST_DATA_SINGLE = auto()
    MASS_LIST_DATA_DOBULE = auto()


@dataclass
class TypedDataFrame():
    """
    """

    df  : pd.DataFrame
    type: FileType


def fetchMassList(allFragsDF: pd.DataFrame) -> np.ndarray:
    """
    Fetches the corresponding mass list from the given dataframe
    of the implied fragment table.


    params:
    -------

    allFragsDF: pd.DataFrame
        The dataframe of the implied fragment table.


    returns:
    --------
    
    massList: np.ndarray
        The fetched mass list.
    """

    massList = []

    for row in allFragsDF.itertuples():

        Iso_Wts = list(map(float, row.Iso_Wts.split(", ")))
        masses  = [row.Exact_Mol_Wt, *Iso_Wts]

        massList.extend([mass, row.Mult] for mass in masses)

    return TypedDataFrame(pd.DataFrame(massList, columns=["Mass", "Multiplicity"]), FileType.MASS_LIST_DATA_DOBULE)


def isHeader(row: pd.Series) -> bool:
    """
    Checks whether a row of a pandas dataframe object is a header row, or not.


    params:
    -------

    row: pd.Series
        The (header) row to check.


    returns:
    --------

    _: bool
        True or False.
    """

    return row.apply(lambda x: isinstance(x, str)).all()


def checkFile(filePath: str, *, isMassList: bool=False) -> TypedDataFrame:
    """
    Parses the file from either:

        i) a mass-spectral data file;
            a) of dimension (n, 2),
               where each row is composed of a pair of m/z and intensity values.
               Both m/z and intensity values are either all integers (low-res)
               or floats (high-res).
           or,

        ii) a mass list data file;
            a) of dimesion (n, 1),
               where each row is composed of a single mass value.
               The mass values are either integers (low-res)
               or floats (high-res); or,
            b) of dimension (n, 2),
               where each row is composed of a pair of mass and multiplicity values.
               The mass values are either integers (low-res)
               or float (high-res),
               and the multiplicity values are always integers.


    params:
    -------

    filePath: str
        The path to the file that is to be parsed.

    isMassList: bool
        A flag to mark if the file is to be parsed as a mass list data file.


    returns:
    --------

    df: TypedDataFrame
        The parsed file as a TypedDataFrame object;
        i.e., if the file is parsed as a mass-spectral data file,
        then it will be of type: TypedFile.MS_DATA.
        If the file is parsed as a mass list data file,
        then it will either be of type:

            i) TypedFile.MASS_LIST_DATA_SINGLE,
               whenever the mass list data file
               is not annotated with multiplicities; or,

            ii) TypedFile.MASS_LIST_DATA_DOUBLE,
                whenever the mass list data file
                is annotated with multiplicites.
    """

    if (not os.path.exists(filePath)):

        raise FileNotFoundError(f"<!> Error: {filePath} not found.")
    
    if (not (filePath.endswith(".csv") or filePath.endswith(".txt"))):

        raise ValueError(f"<!> Error: {filePath} is not of type .csv or -.txt.")

    tempDF = pd.read_csv(filePath, header=None, sep=r"\s+")
    df     = tempDF[1:].copy() if isHeader(tempDF.iloc[0]) else tempDF

    if (isMassList):

        if (df.shape[1] not in {1, 2}):

            raise ValueError(f"<!> Error: {filePath} must have 1 or 2 columns. " \
                             f"{df.shape[1]} were given.")

        df[0] = pd.to_numeric(df[0], errors="raise")

        if (df.shape[1] == 2):

            df[1] = pd.to_numeric(df[1], errors="raise", downcast="integer")

            return TypedDataFrame(df, FileType.MASS_LIST_DATA_DOBULE)

        else:

            return TypedDataFrame(df, FileType.MASS_LIST_DATA_SINGLE)

    elif (df.shape[1] != 2):

        raise ValueError(f"<!> Error: {filePath} must have 2 columns. " \
                         f"{df.shape[1]} was/were given.")

    if (not np.issubdtype(df.values.dtype, np.number)):

        raise ValueError(f"<!> Error: {filePath} contains non-numeric data.")
    
    return TypedDataFrame(df, FileType.MS_DATA)


def compare(msDataPath: str, massListPath: str, *, tol: float) -> float:
    """
    Compares the mass-spectral data associated with the given path
    and the mass list data associated with the given path,
    calculating the F.P.I.E. score between them
    and displaying an annotated stem plot of the spectrum.


    params:
    -------

    msDataPath: str
        The path to the mass-spectral data.

    massListPath: str
        The path to the mass list data.


    returns:
    --------

        The F.P.I.E. score.
    """
    
    msDataTyped             = checkFile(msDataPath)
    massListTyped           = checkFile(massListPath, isMassList=True)
    FPIEScore, plotMetaData = FPIE(
        msDataTyped.df.to_numpy(),
        massListTyped.df.to_numpy(),
        tol=tol,
        weighted=(massListTyped.type == FileType.MASS_LIST_DATA_DOBULE),
    )

    plotFPIE(plotMetaData, FPIEScore)

    return FPIEScore


def compareAll(msDataPath: str, *, tol: float) -> pd.DataFrame:
    """
    Compares the mass-spectral data associated with the given path
    to the mass list of each molecule in the M.F.D.,
    calculating the F.P.I.E. score between them.


    params:
    -------

    msDataPath: str
        The path to the mass-spectral data.


    returns:
    --------

    df: pd.DataFrame
        The table of F.P.I.E. scores.
    """

    msDataTyped     = checkFile(msDataPath)
    FPIEs           = []
    idxTableDF      = viewIdxTable()
    craftsLabEntrys = idxTableDF.iloc[:, 1]

    for entry in tqdm(craftsLabEntrys, desc=f"<*> Calculating FPIEs . . ."):

        massListTyped = fetchMassList(searchAndFetch(entry))
        FPIEScore, _  = FPIE(
            msDataTyped.df.to_numpy(),
            massListTyped.df.to_numpy(),
            tol=tol,
            weighted=(massListTyped.type == FileType.MASS_LIST_DATA_DOBULE)
        )
        
        FPIEs.append(FPIEScore)

    df = pd.concat([idxTableDF.iloc[:, 0], pd.DataFrame({"FPIE": FPIEs})],
                     axis=1)

    return df


def searchAndFetchByMass(mz: float) -> pd.DataFrame:
    """
    """

    mz           = round(float(mz), 2)
    idxTableDF   = viewIdxTable()
    filteredRows = []

    for idx, entry in tqdm(idxTableDF.iloc[:, 1].items(), desc=f"<*> Fetching fragments . . ."):

        massList = fetchMassList(searchAndFetch(entry)).df.iloc[:, 0].round(2)

        if (mz in massList.values):

            filteredRows.append(idx)
            tqdm.write(f"<*> Found a match!")

    return idxTableDF.loc[filteredRows].reset_index(drop=True)


def plotFPIE(plotMetaData: dict, FPIEScore: float) -> None:

        msDataMasses      = plotMetaData["msDataMasses"]
        msDataIntensities = plotMetaData["msDataIntensities"]
        matchedWeights    = plotMetaData["matchedWeights"]
        matchedMults      = plotMetaData["matchedMults"]
        massesMultDict    = plotMetaData["massesMultDict"]

        fig, ax = plt.subplots()

        fig.patch.set_facecolor("darkblue")
        ax.set_facecolor("lightblue")

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

                color = "black"

            _, stemline, _ = ax.stem([mz], [intensity],
                                     linefmt=color,
                                     markerfmt=" ",
                                     basefmt=" ")

            stemline.set_linewidth(10)

        ax.set_xlabel("m/z", color="yellow")
        ax.set_ylabel("Intensity", color="yellow")
        ax.set_title(f"FPIE Annotated Plot", color="yellow")

        legendHandles = [
            Patch(color=multColorMap[mult], label=f"{int(mult)}")
            for mult in mults
        ]

        ax.legend(handles=legendHandles, title="Multiplicity")
        ax.tick_params(axis="both", colors="yellow")

        plt.show()


def FPIE(msData: np.ndarray,
         massList: np.ndarray,
         *,
         tol: float,
         weighted: bool=False) -> float:
    """
    Calculates the F.P.I.E. score between the given mass-spectral data
    and the given mass list data.


    params:
    -------

    msData: np.ndarray
        The mass-spectral data as a np.ndarray object.

    massList: np.ndarray
        The mass list data as a np.ndarray object.

    tol: float
        The tolerance by which to compare masses in the mass-spectral data
        to masses on the mass list. If the mass-spectral data is high-res,
        then a tolerance of 0.005 is recommended.
        If the mass-spectral data is low-res,
        then a tolerance of 0 is recommended.

    weighted: bool
        A flag to mark if the mass list data is annotated with multiplicities;
        i.e., of tpye: FileType.MASS_LIST_DATA_DOUBLE.


    returns:
    --------

    FPIEScore: float
        The F.P.I.E. score.
    """
    if (massList.ndim == 1):

        massList = np.column_stack((massList, np.zeros_like(massList)))

    masses       = massList[:, 0]
    mults        = massList[:, 1]
    msDataMasses = np.unique(msData[:, 0])

    if (tol == 0):

        masses = masses.round().astype(int)

    else:

        masses       = np.round(masses, 2)
        msDataMasses = np.round(msDataMasses, 2)

    massesMultDict = {}

    for mass, mult in zip(masses, mults):

        if ((mass not in massesMultDict) or (mult < massesMultDict[mass])):

            massesMultDict[mass] = mult

    masses         = np.array(list(massesMultDict.keys()))
    weights        = np.array([
        np.exp(-(0.01 * massesMultDict[m])) if weighted else 1.0
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

    return FPIEScore, {
        "msDataMasses": msDataMasses,
        "msDataIntensities": msDataIntensities,
        "matchedWeights": matchedWeights,
        "matchedMults": matchedMults,
        "massesMultDict": massesMultDict,
        "weights": weights
    }
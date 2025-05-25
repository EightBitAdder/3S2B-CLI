from db.utils import (
    viewIdxTable,
    searchAndFetch
)
from file_readers.utils import (
    FileType,
    TypedDataFrame,
    FileReader,
    DelimitedFileReader,
    MassListReader,
    MSDataReader
)
import os
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from tqdm import tqdm


class MassList():


    def __init__(self, df: pd.DataFrame, weighted: bool):

        self.df       = df
        self.weighted = weighted
        self.masses   = df.iloc[:, 0].to_numpy()
        self.mults    = df.iloc[:, 1].to_numpy() if df.shape[1] > 1 else np.ones(len(df))


    @classmethod
    def fromFile(cls, filePath: str) -> "MassList":

        typedDF = MassListReader(filePath).read()
        weighted = (typedDF.type == FileType.MASS_LIST_DATA_DOBULE)

        return cls(typedDF.df, weighted)


class MSData():


    def __init__(self, df: pd.DataFrame):

        self.df          = df
        self.masses      = df.iloc[:, 0].to_numpy()
        self.intensities = df.iloc[:, 1].to_numpy()


    @classmethod
    def fromFile(cls, filePath: str) -> "MSData":

        typedDF = MSDataReader(filePath).read()

        return cls(typedDF.df)


def fetchMassList(allFragsDF: pd.DataFrame) -> MassList:

    massList = []

    for row in allFragsDF.itertuples():

        Iso_Wts = list(map(float, row.Iso_Wts.split(", ")))
        masses  = [row.Exact_Mol_Wt, *Iso_Wts]

        massList.extend([mass, row.Mult] for mass in masses)

    df = pd.DataFrame(massList, columns=["Mass", "Multiplicity"])

    return MassList(df, weighted=True)


def compare(msDataPath: str, massListPath: str, *, tol: float) -> float:
    
    msData                  = MSData.fromFile(msDataPath)
    massList                = MassList.fromFile(massListPath)
    FPIEScore, plotMetaData = FPIE(
        msData,
        massList,
        tol=tol,
        weighted=massList.weighted
    )

    plotFPIE(
        plotMetaData,
        FPIEScore,
        os.path.splitext(os.path.basename(msDataPath))[0]
    )

    return FPIEScore


def compareAll(msDataPath: str, *, tol: float) -> pd.DataFrame:

    msData          = MSData.fromFile(msDataPath)
    FPIEs           = []
    idxTableDF      = viewIdxTable()
    craftsLabEntrys = idxTableDF.iloc[:, 1]

    for entry in tqdm(craftsLabEntrys, desc=f"<*> Calculating FPIEs . . ."):

        massList      = fetchMassList(searchAndFetch(entry))
        FPIEScore, _  = FPIE(
            msData,
            massList,
            tol=tol,
            weighted=massList.weighted
        )
        
        FPIEs.append(FPIEScore)

    df = pd.concat([idxTableDF.iloc[:, 0], pd.DataFrame({"FPIE": FPIEs})],
                     axis=1)

    return df


def searchAndFetchByMass(mz: float) -> pd.DataFrame:

    mz           = round(float(mz), 2)
    idxTableDF   = viewIdxTable()
    filteredRows = []

    for idx, entry in tqdm(idxTableDF.iloc[:, 1].items(), desc=f"<*> Fetching fragments . . ."):

        massList = fetchMassList(searchAndFetch(entry))

        if (mz in np.round(massList.masses, 2)):

            filteredRows.append(idx)
            tqdm.write(f"<*> Found a match!")

    return idxTableDF.loc[filteredRows].reset_index(drop=True)


def plotFPIE(plotMetaData: dict, FPIEScore: float, title: str="") -> None:

        msDataMasses      = plotMetaData["msDataMasses"]
        msDataIntensities = plotMetaData["msDataIntensities"]
        matchedWeights    = plotMetaData["matchedWeights"]
        matchedMults      = plotMetaData["matchedMults"]
        massesMultDict    = plotMetaData["massesMultDict"]

        fig, ax = plt.subplots()

        fig.patch.set_facecolor("darkblue")
        ax.set_facecolor("white")

        mults        = sorted(set(massesMultDict.values()))
        numColors    = len(mults)
        colorMap     = cm.get_cmap("tab20b", numColors)
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

            stemline.set_linewidth(1)

            if (mz not in matchedWeights):

                ax.text(mz,
                        intensity + intensity * 0.05,
                        "*",
                        color="red",
                        ha="center",
                        fontsize=10,
                        clip_on=True)

        ax.set_xlabel("m/z", color="yellow")
        ax.set_ylabel("Intensity", color="yellow")
        ax.set_title(f"{title} FPIE Annotated Plot", color="yellow")
        ax.grid(True, linestyle='--', linewidth=1, color='black', alpha=0.5)

        legendHandles     = [
            Patch(color=multColorMap[mult], label=f"{int(mult)}")
            for mult in mults
        ]
        unexplainedHandle = Line2D([0],
                                   [0],
                                   marker="*",
                                   color="red",
                                   linestyle="None",
                                   markersize=10,
                                   label="Unexplained")

        legendHandles.append(unexplainedHandle)

        ax.legend(handles=legendHandles, title="Multiplicity")
        ax.tick_params(axis="both", colors="yellow")

        plt.show()


def FPIE(msData: MSData,
         massList: MassList,
         *,
         tol: float,
         weighted: bool=False) -> float:

    masses       = massList.masses
    mults        = massList.mults
    msDataMasses = np.unique(msData.masses)

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

    for mz, intensity in zip(msData.masses, msData.intensities):

        mzKey = round(mz) if tol == 0 else np.round(mz, 2)

        if (mzKey in matchedWeights):

            explainedIntensities.append(intensity * matchedWeights[mzKey])

    msDataIntensities = msData.intensities
    FPIEScore         = float(sum(explainedIntensities) / np.sum(msDataIntensities))

    return FPIEScore, {
        "msDataMasses": msDataMasses,
        "msDataIntensities": msDataIntensities,
        "matchedWeights": matchedWeights,
        "matchedMults": matchedMults,
        "massesMultDict": massesMultDict,
        "weights": weights
    }
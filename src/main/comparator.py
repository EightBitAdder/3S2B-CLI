from typing import Tuple, Dict
from file_readers.utils import (
    FileType,
    TypedDataFrame,
    FileReader,
    DelimitedFileReader,
    MassListReader,
    MSDataReader
)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


class MassList():


    def __init__(self, df: pd.DataFrame, weighted: bool):

        self.df       = df
        self.weighted = weighted
        self.masses   = df.iloc[:, 0].to_numpy()
        self.mults    = df.iloc[:, 1].to_numpy() if df.shape[1] > 1 else np.ones(len(df))


    @classmethod
    def fromFile(cls, filePath: str) -> "MassList":

        typedDF = MassListReader(filePath, delimiter=r"\s+").typedDF
        weighted = (typedDF.fileType == FileType.MASS_LIST_DATA_DOUBLE)

        return cls(typedDF.df, weighted)


class MSData():


    def __init__(self, df: pd.DataFrame):

        self.df          = df
        self.masses      = df.iloc[:, 0].to_numpy()
        self.intensities = df.iloc[:, 1].to_numpy()


    @classmethod
    def fromFile(cls, filePath: str) -> "MSData":

        return cls(MSDataReader(filePath, delimiter=r"\s+").typedDF.df)


class Comparator():

	
	def __init__(self, msData: MSData, massList: MassList, tol: float):

		self.msData   = msData
		self.massList = massList
		self.tol      = tol

		if (tol == 0):

			self.msDataMasses   = np.unique(np.round(msData.masses).astype(int))
			self.massListMasses = np.round(massList.masses).astype(int)

		else:

			self.msDataMasses   = np.unique(np.round(msData.masses, 2))
			self.massListMasses = np.round(massList.masses, 2)

		self.weights, self.massesMultsDict = self._fetchWeights()


	def _fetchWeights(self):

		massesMultsDict = {}

		for mass, mult in zip(self.massListMasses, self.massList.mults):

			if ((mass not in massesMultsDict) or (mult < massesMultsDict[mass])):

				massesMultsDict[mass] = mult

		masses  = np.array(list(massesMultsDict.keys()))
		weights = np.array([
			np.exp(-(0.01 * massesMultsDict[m])) if self.massList.weighted else 1.0
			for m in masses
		])

		return weights, massesMultsDict


	def calculateFPIE(self) -> Tuple[float, Dict]:

		matchedWeights = {}
		matchedMults   = {}
		masses         = np.array(list(self.massesMultsDict.keys()))

		for mz in self.msDataMasses:

			for i, m in enumerate(masses):

				if (m - self.tol <= mz <= m + self.tol):

					matchedMults[mz]   = self.massesMultsDict[m]
					matchedWeights[mz] = self.weights[i]

					break


		explainedIntensities = []

		for mz, intensity in zip(self.msDataMasses, self.msData.intensities):

			if (mz in matchedWeights):

				explainedIntensities.append(intensity * matchedWeights[mz])

		FPIEScore = float(sum(explainedIntensities) / np.sum(self.msData.intensities))

		return FPIEScore, {
			"msDataMasses": self.msDataMasses,
			"msDataIntensities": self.msData.intensities,
			"matchedMults": matchedMults,
			"massesMultsDict": self.massesMultsDict,
			"weights": self.weights,
			"matchedWeights": matchedWeights
		}


	def plotFPIE(self, plotMetaData: Dict, FPIEScore: float, title: str="") -> None:

		fig, ax = plt.subplots()

		fig.patch.set_facecolor("darkblue")
		ax.set_facecolor("white")

		mults         = sorted(set(plotMetaData["massesMultsDict"].values()))
		colorMap      = cm.get_cmap("tab20b", len(mults))
		multsColorMap = {m: colorMap(i) for i, m in enumerate(mults)}

		for i, mz in enumerate(plotMetaData["msDataMasses"]):

			intensity      = plotMetaData["msDataIntensities"][i]
			color          = multsColorMap.get(
				plotMetaData["matchedMults"].get(mz, None),
				"black"
			)
			_, stemline, _ = ax.stem(
				[mz],
				[intensity],
				linefmt=color,
				markerfmt=" ",
				basefmt=" "
			)

			if (mz not in plotMetaData["matchedWeights"]):

				ax.text(
					mz,
					intensity + intensity * 0.05,
					"*",
					color="red",
					ha="center",
					fontsize=10
				)

		ax.set_xlabel("m/z", color="yellow")
		ax.set_ylabel("Intensity", color="yellow")
		ax.set_title(f"{title} FPIE Annotated Plot", color="yellow")
		ax.grid(True, linestyle="--", linewidth=1, color="black", alpha=0.5)

		legendHandles = [
			Patch(
				color=multsColorMap[m],
				label=f"{int(m)}"
			)
			for m in mults
		]

		legendHandles.append(
			Line2D(
				[0],
				[0],
				marker="*",
				color="red",
				linestyle="None",
				markersize=10,
				label="Unexplained"
			)
		)

		ax.legend(handles=legendHandles, title="Multiplicity")
		ax.tick_params(axis="both", colors="yellow")

		plt.show()
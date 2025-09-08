from typing import Tuple, Dict
from dataclasses import dataclass
from .weight_functions import (
	ConstantWeight,
	ExponentialDecayWeight,
	fetchWeightFunction
)
from .mass_list import MassList
from .ms_data import MSData
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


@dataclass
class PlotMetaData:

	msDataMasses     : np.ndarray
	matchedMults     : Dict[float, int]
	massesMultsDict  : Dict[float, int]
	matchedWeights   : Dict[float, float]


class Comparator():

	
	def __init__(
		self,
		msData: MSData,
		massList: MassList,
		*,
		tol: float,
		weightFunction: str | None=None
	):

		self.msData   = msData
		self.massList = massList
		self.tol      = tol

		try:

			self.weightFunction = fetchWeightFunction(weightFunction)

		except:

			self.weightFunction = (
				ExponentialDecayWeight() if massList.weighted else ConstantWeight()
			)


	@property
	def tol(self) -> float:

		return self._tol


	@tol.setter
	def tol(self, val: float) -> None:

		if (val < 0):

			raise ValueError(f"<!> Error: {val} is not positive.")

		self._tol = val



	def _round(self) -> Tuple[np.ndarray, np.ndarray]:

		if (self.tol == 0):

			return (
				np.unique(np.round(self.msData.masses).astype(int)),
				np.round(self.massList.masses).astype(int)
			)
		else:

			return(
				np.unique(np.round(self.msData.masses, 2)),
				np.round(self.massList.masses, 2)
			)


	def calculateFPIE(self) -> Tuple[float, Dict]:

		msDataMasses, massListMasses = self._round()
		massListMults                = self.massList.mults

		massesMultsDict = {}

		for mass, mult in zip(massListMasses, massListMults):

			if (mass not in massesMultsDict or mult < massesMultsDict[mass]):

				massesMultsDict[mass] = mult

		weightsDict = {
			mass: self.weightFunction.compute(mass, massesMultsDict[mass])
			for mass in massesMultsDict
		}

		matchedWeights = {}
		matchedMults   = {}

		for mz in msDataMasses:

			for mass in massesMultsDict:

				if (mass - self.tol <= mz <= mass + self.tol):

					matchedMults[mz]   = massesMultsDict[mass]
					matchedWeights[mz] = weightsDict[mass]

					break

		explainedIntensities = [
			intensity * matchedWeights[mz]
			for mz, intensity in zip(msDataMasses, self.msData.intensities)
			if mz in matchedWeights
		]

		FPIEScore    = float(sum(explainedIntensities) / np.sum(self.msData.intensities))
		plotMetaData = PlotMetaData(
			msDataMasses=msDataMasses,
			matchedMults=matchedMults,
			massesMultsDict=massesMultsDict,
			matchedWeights=matchedWeights
		)

		return FPIEScore, plotMetaData


	def plotFPIE(self, plotMetaData: PlotMetaData, FPIEScore: float, title: str="") -> None:

		fig, ax = plt.subplots()

		fig.patch.set_facecolor("darkblue")
		ax.set_facecolor("white")

		mults         = sorted(set(plotMetaData.massesMultsDict.values()))
		colorMap      = cm.get_cmap("tab20b", len(mults))
		multsColorMap = {mult: colorMap(i) for i, mult in enumerate(mults)}

		for i, mz in enumerate(plotMetaData.msDataMasses):

			intensity      = self.msData.intensities[i]
			color          = multsColorMap.get(
				plotMetaData.matchedMults.get(mz, None),
				"black"
			)
			_, stemline, _ = ax.stem(
				[mz],
				[intensity],
				linefmt=color,
				markerfmt=" ",
				basefmt=" "
			)

			if (mz not in plotMetaData.matchedWeights):

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
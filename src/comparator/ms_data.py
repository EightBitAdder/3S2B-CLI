from file_readers import (
    fetchFileReader
)
import pandas as pd
import numpy as np


class MSData():


    def __init__(self, df: pd.DataFrame):

        self.df = self._mount(df)

        self.masses     : np.ndarray  = df.iloc[:, 0].to_numpy()
        self.intensities: np.ndarray  = df.iloc[:, 1].to_numpy()


    def _mount(self, df: pd.DataFrame) -> pd.DataFrame:

        if (df.shape[1] != 2):

           raise ValueError(f"<!> Error: {val}.shape is not (n, 2).")

        return df


    @classmethod
    def fromFile(cls, filePath: str, *, delimiter: str | None=None, fileReader: str="MS-DATA") -> "MSData":

        return cls(fetchFileReader(fileReader, filePath, delimiter=delimiter).typedDF.df)
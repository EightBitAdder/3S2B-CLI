from typing import Dict
from file_readers import (
    FileType,
    fetchFileReader
)
import pandas as pd
import numpy as np


class MassList():


    def __init__(self, df: pd.DataFrame, *, weighted: bool):

        self.df       = self._mount(df)
        self.weighted = weighted

        self.masses          : np.ndarray       = self.df.iloc[:, 0].to_numpy()
        self.mults           : np.ndarray       = (
            self.df.iloc[:, 1].to_numpy()
            if self.weighted else np.ones(len(self.df))
        )


    def _mount(self, df: pd.DataFrame) -> pd.DataFrame:

        if (df.shape[1] not in (1, 2)):

            raise ValueError(f"<!> Error: {val}.shape is not (n, 1) or (n, 2).")

        return df


    @classmethod
    def fromFile(cls, filePath: str, *, delimiter: str | None=None, fileReader: str="MASS-LIST") -> "MassList":

        reader  = fetchFileReader(fileReader, filePath, delimiter=delimiter)
        weighted = (reader.typedDF.fileType == FileType.MASS_LIST_DATA_DOUBLE)

        return cls(reader.typedDF.df, weighted=weighted)
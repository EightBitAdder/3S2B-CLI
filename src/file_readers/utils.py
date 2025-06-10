import os
from abc import ABC, abstractmethod
from enum import Enum, auto
from dataclasses import dataclass
import csv
import pandas as pd


class FileType(Enum):

    MASS_SPEC_DATA        = auto()
    MASS_LIST_DATA_SINGLE = auto()
    MASS_LIST_DATA_DOUBLE = auto()


@dataclass
class TypedDataFrame():

    df      : pd.DataFrame
    fileType: FileType


class FileReader(ABC):


    def __init__(self, filePath: str):

        self.filePath                       = filePath
        self._typedDF: TypedDataFrame | None = None


    @property
    def filePath(self) -> str:

        return self._filePath


    @filePath.setter
    def filePath(self, val: str) -> None:

        if (not os.path.exists(val)):

            raise FileNotFoundError()

        self._filePath = val


    @property
    def typedDF(self) -> TypedDataFrame:

        if (self._typedDF is None):

            self._typedDF = self._read()

        return self._typedDF


    @abstractmethod
    def _mount(self) -> pd.DataFrame:

        pass


    @abstractmethod
    def _scrub(self, df: pd.DataFrame) -> pd.DataFrame:

        pass


    @abstractmethod
    def _parse(self, df: pd.DataFrame) -> FileType:

        pass


    def _read(self) -> TypedDataFrame:

        df       = self._mount()
        fileType = self._parse(self._scrub(df))

        return TypedDataFrame(df, fileType)


class DelimitedFileReader(FileReader):


    def __init__(self, filePath: str, *, delimiter: str | None=None):

        self._delimiter = delimiter

        super().__init__(filePath)


    @property
    def delimiter(self) -> str:

        if (self._delimiter is None):

            self._delimiter = self._sniff()

        return self._delimiter


    def _sniff(self) -> str:

        with open(self.filePath, "r", encoding="utf-8") as fo:

            sample = fo.read(1024)

        try:

            dialect = csv.Sniffer().sniff(sample, delimiters=" \t,;|")

            return dialect.delimiter

        except:

            raise ValueError()


    def _isheader(self, row: pd.Series) -> bool:

        return row.apply(lambda char: isinstance(char, str)).all()


    def _mount(self) -> pd.DataFrame:

        df = pd.read_csv(self.filePath,
                         header=None,
                         sep=self.delimiter)
        df = df.iloc[1:] if self._isheader(df.iloc[0]) else df

        return df


class MassListReader(DelimitedFileReader):


    def __init__(self, filePath: str, *, delimiter: str):

        super().__init__(filePath, delimiter=delimiter)


    def _mount(self) -> pd.DataFrame:

        df = super()._mount()

        match(df.shape[1]):

            case 1:

                df.columns = ["mass"]

            case 2:

                df.columns = ["mass", "mult"]

            case _:

                raise ValueError()

        return df


    def _scrub(self, df: pd.DataFrame) -> pd.DataFrame:

        df["mass"] = pd.to_numeric(df["mass"], errors="raise")

        if ("mult" in df.columns):

            df["mult"] = pd.to_numeric(df["mult"], errors="raise", downcast="integer")

        return df


    def _parse(self, df: pd.DataFrame) -> FileType:

        return (
            FileType.MASS_LIST_DATA_DOUBLE
            if df.shape[1] == 2
            else FileType.MASS_LIST_DATA_SINGLE
        )


class MSDataReader(DelimitedFileReader):


    def __init__(self, filePath: str, *, delimiter: str):

        super().__init__(filePath, delimiter=delimiter)


    def _mount(self) -> pd.DataFrame:

        df = super()._mount()

        if (df.shape[1] != 2):

            raise ValueError()

        df.columns = ["mz", "intensity"]

        return df


    def _scrub(self, df: pd.DataFrame) -> None:

        df["mz"] = pd.to_numeric(df["mz"], errors="raise")
        df["intensity"] = pd.to_numeric(df["intensity"], errors="raise")

        return df


    def _parse(self, df: pd.DataFrame) -> FileType:

        return FileType.MASS_SPEC_DATA
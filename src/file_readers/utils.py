import os
from abc import ABC, abstractmethod
from enum import Enum, auto
from dataclasses import dataclass
import pandas as pd


class FileType(Enum):

    MS_DATA               = auto()
    MASS_LIST_DATA_SINGLE = auto()
    MASS_LIST_DATA_DOBULE = auto()


@dataclass
class TypedDataFrame():

    df  : pd.DataFrame
    type: FileType


class FileReader(ABC):


    def __init__(self, filePath: str):

        self.filePath = filePath


    @abstractmethod
    def read(self) -> TypedDataFrame:

        pass


class DelimitedFileReader(FileReader):


    def __init__(self, filePath: str, delimiter: str=r"\s+"):

        super().__init__(filePath)
        self.delimiter = delimiter

        self.df: pd.DataFrame | None = None


    def _isHeader(self, row: pd.Series) -> bool:

        return row.apply(lambda x: isinstance(x, str)).all()


    def load(self) -> pd.DataFrame:

        if (not os.path.exists(self.filePath)):

            raise FileNotFoundError(f"<!> Error: file not found.")

        if (not self.filePath.endswith((".csv", ".txt"))):

            raise ValueError(f"<!> Error: unsupported file extension.")

        df      = pd.read_csv(self.filePath, header=None, sep=self.delimiter)
        self.df = df[1:] if self._isHeader(df.iloc[0]) else df


    @abstractmethod
    def validate(self) -> None:

        pass


    @abstractmethod
    def convert(self) -> pd.DataFrame:

        pass


    @abstractmethod
    def fileType(self) -> FileType:

        pass


    def read(self) -> TypedDataFrame:

        self.load()
        self.validate()
        self.convert()

        return TypedDataFrame(self.df, self.fileType())


class MassListReader(DelimitedFileReader):


    def __init__(self, filePath: str, delimiter: str=r"\s+"):

        super().__init__(filePath)
        self.delimiter = delimiter


    def validate(self) -> None:

        if (self.df.shape[1] not in {1, 2}):

            raise ValueError(f"<!> Error: file must have 1 or 2 columns; but, " \
                             f"{self.df.shape[1]} were given.")


    def convert(self) -> pd.DataFrame:

        self.df[0] = pd.to_numeric(self.df[0], errors="raise")

        if (self.df.shape[1] == 2):

            self.df[1] = pd.to_numeric(self.df[1], errors="raise", downcast="integer")

        return self.df


    def fileType(self) -> FileType:

        return FileType.MASS_LIST_DATA_DOBULE if self.df.shape[1] == 2 else FileType.MASS_LIST_DATA_SINGLE


class MSDataReader(DelimitedFileReader):


    def __init__(self, filePath: str, delimiter: str=r"\s+"):

        super().__init__(filePath)
        self.delimiter = delimiter


    def validate(self) -> None:

        if (self.df.shape[1] != 2):

            raise ValueError(f"<!> Error: file must have 2 columns; but, " \
                             f"{df.shape[1]} was/were given.")


    def convert(self) -> None:

        self.df[0] = pd.to_numeric(self.df[0], errors="raise")
        self.df[1] = pd.to_numeric(self.df[1], errors="raise")


    def fileType(self) -> FileType:

        return FileType.MS_DATA
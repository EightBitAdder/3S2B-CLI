import os
from abc import ABC, abstractmethod
from enum import Enum, auto
from dataclasses import dataclass
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

        self.filePath                        = filePath
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
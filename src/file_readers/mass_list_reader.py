from .file_reader import FileType
from .delimited_file_reader import DelimitedFileReader
from .registry import registerFileReader
import pandas as pd


@registerFileReader("MASS-LIST")
class MassListReader(DelimitedFileReader):


    def __init__(self, filePath: str, *, delimiter: str | None=None):

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
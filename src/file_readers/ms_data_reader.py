from .file_reader import FileType
from .delimited_file_reader import DelimitedFileReader
from .registry import registerFileReader
import pandas as pd


@registerFileReader("MS-DATA")
class MSDataReader(DelimitedFileReader):


    def __init__(self, filePath: str, *, delimiter: str | None=None):

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
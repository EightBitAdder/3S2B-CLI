from .file_reader import FileReader
import csv
import pandas as pd


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

        except csv.Error:

            return r"\s+"


    def _isheader(self, row: pd.Series) -> bool:

        return row.apply(lambda char: isinstance(char, str)).all()


    def _mount(self) -> pd.DataFrame:

        df = pd.read_csv(self.filePath,
                         header=None,
                         sep=self.delimiter)
        df = df.iloc[1:] if self._isheader(df.iloc[0]) else df

        return df
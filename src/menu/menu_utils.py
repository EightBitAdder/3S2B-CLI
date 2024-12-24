import os
from utils import compare, compareAll
from db_utils import download, searchAndFetch, viewIdxTable
from menu.menu import Handler
import sqlite3
import pandas as pd
import textwrap
from prettytable import PrettyTable


DB_PATH = "resource/swgdrugdb.db"


def dfToPrettyTable(df: pd.DataFrame, header: bool=True, colWidth: int=50):

    table = PrettyTable()

    if (header):
        table.field_names = df.columns.tolist()
    else:
        table.header=False

    table.junction_char = "-"

    for row in df.itertuples(index=False):
        wrappedRow = [
            textwrap.fill(str(cell), width=colWidth) for cell in row
        ]

        table.add_row(wrappedRow)

    for col in table.field_names:
        table.align[col] = "l"

    table.hrules = True

    return table


def paginateDF(df: pd.DataFrame, preamble: str | PrettyTable=None, minRows: int=5, colWidth: int=50):

    numRows = len(df)
    numPage = 0

    def clearConsole():

        # Clearing PrettyTable objects from the console is a bloody mess!
        os.system("cls" if os.name == "nt" else "clear")

        print("\n" * 5)


    while True:

        start = numPage * minRows
        end = start + minRows

        if (start >= numRows):

            print("<*> End of Dataframe.")

            break

        currPage = df.iloc[start:end]
        table = dfToPrettyTable(currPage, colWidth=colWidth)

        clearConsole()

        if (preamble):

            if (isinstance(preamble, PrettyTable)):

                print(preamble, "\n")

            elif (isinstance(preamble, str)):

                print(dfToPrettyTable(pd.DataFrame({"": [preamble]}), header=False),
                      "\n")

        else:

            print(f"<!> Unsupported preamble type: {type(preamble)}\n")

        print(f"Page {numPage + 1} (Rows {start + 1}-{min(end, numRows)} of {numRows})")
        print(table)

        command = input("Press [n] for next page, [p] for previous page; or, [q] to quit: ").strip().lower()

        match command:

            case "n":

                if (end < numRows):

                    numPage += 1
                else:

                    print("<!> Page does not exist.")

            case "p":

                if (numPage > 0):

                    numPage -= 1

                else:

                    print("<!> Page does not exist.")  

            case "q":

                break

            case _:

                print("<!> Invalid choice. Please try again.")


class CompareHandler(Handler):

    def handle(self, *args) -> None:

        print(dfToPrettyTable(pd.DataFrame({"FPIE": [compare(*args)]})))


class CompareAllHandler(Handler):

    def handle(self, *args) -> None:

        result = compareAll(*args)
        maxFPIEs = result[result["FPIE"] == result["FPIE"].max()]

        paginateDF(maxFPIEs, "SWGDRUG Mols. w/ Max. FPIE:")
        paginateDF(result, "SWGDRUG Mols. w/ FPIE:")


class DownloadHandler(Handler):

    def handle(self, *args) -> None:

        try:

            print("Donwloading {}...".format(*args))

            result = download(*args)
        except:

            print("<!> Error in downloading {}.".format(*args))
        else:

            print("<*> Download complete.")


class SearchAndFetchHandler(Handler):

    def handle(self, *args) -> None:

        conn = sqlite3.connect(DB_PATH)
        cursor = conn.cursor()

        result = searchAndFetch(*args)
        titleDF = pd.read_sql(f"SELECT * FROM idxTable WHERE craftsLabEntry = ?",
                              conn,
                              params=(*args,))

        paginateDF(result, dfToPrettyTable(titleDF))

        cursor.close()
        conn.close()


class ViewIdxTableHandler(Handler):

    def handle(self, *args) -> None:

        result = viewIdxTable(*args)

        paginateDF(result, "MFD Index Table")
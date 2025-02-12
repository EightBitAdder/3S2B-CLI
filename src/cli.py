from utils import compare, compareAll
from db_utils import searchAndFetch, viewIdxTable
import os
import sqlite3
import click
from textual.app import App
from textual.containers import Vertical
from textual.widgets import DataTable, Footer


DB_PATH = "resource/swgdrugdb.db"


class ScrollableTable(App):

    BINDINGS = [("escape", "quit", "Quit"), ("d", "download", "Download")]


    def __init__(self, df, title=""):

        super().__init__()
        self.df = df
        self.title = title

    
    def compose(self):

        yield Vertical(id="table-container")
        yield Footer()

    
    def on_mount(self):

        container = self.query_one("#table-container")

        table = DataTable()

        table.add_columns(*self.df.columns.astype(str))

        for row in self.df.itertuples(index=False):
            
            table.add_row(*map(str, row))

        table.cursor_type       = "row"
        table.styles.height     = 50
        table.styles.background = "darkblue"
        table.styles.color      = "white"
        table.border_title      = self.title
        table.styles.border     = ("heavy", "yellow")

        container.mount(table)


    def action_quit(self):

        self.exit()


    def action_download(self):

        downloads_path = os.path.join(os.path.expanduser("~"), "Downloads")
        file_name      = "table.csv"
        full_path      = os.path.join(downloads_path, file_name)

        self.df.to_csv(full_path, index=False)

        self.notify(f"<*> Table saved to {full_path}.", severity="info")


@click.group()
def cli():

    pass


@click.command()
@click.argument("file_paths", nargs=-1, type=click.Path(exists=True))
@click.argument("tol", type=float)
def c(file_paths, tol):

    print(f"FPIE: {compare(*file_paths, tol)}")


@click.command()
@click.argument("ms_data_path", type=click.Path(exists=True))
@click.argument("tol", type=float)
def a(ms_data_path, tol):

    result = compareAll(ms_data_path, tol)
    result = result.sort_values(by="FPIE", ascending=False)

    ScrollableTable(result, "SWGDRUG SMILES w/ FPIEs").run()


@click.command()
@click.argument("search_term")
def f(search_term):

    conn    = sqlite3.connect(DB_PATH)
    
    ScrollableTable(searchAndFetch(search_term), f"{search_term} Fragment Table").run()

    conn.close()


@click.command()
@click.argument("args", nargs=-1)
def i(args):

    ScrollableTable(viewIdxTable(*args), "MFD Index Table").run()


cli.add_command(c)
cli.add_command(a)
cli.add_command(f)
cli.add_command(i)


if __name__ == "__main__":

    cli()
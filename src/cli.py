from fragmentor import Fragmentor
from utils import fetchMassList, compare, compareAll, searchAndFetchByMass
from db_utils import addEntryFromSmiles, searchAndFetch, viewIdxTable
import os
import pandas as pd
import sqlite3
import click
from textual.app import App
from textual import on
from textual.containers import Vertical
from textual.widgets import DataTable, Footer
from rdkit import Chem


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
DB_PATH     = os.path.abspath(os.path.join(CURRENT_DIR, "..", "resource", "swgdrugdb.db"))


class ScrollableTable(App):

    BINDINGS = [("escape", "go_back", "Go Back"),
                ("d", "download", "Download"),
                ("m", "fetch_mass_list", "Fetch Mass List")]


    def __init__(self, df, title="", *, parent_df=None, parent_title=""):

        super().__init__()
        self.df           = df
        self.title        = title
        self.parent_df    = parent_df
        self.parent_title = parent_title
        self.curr         = None

    
    def compose(self):

        yield Vertical(id="table-container")
        yield Footer()

    
    async def on_mount(self):

        await self.load_table(self.df, self.title)


    async def load_table(self, df, title=""):

        self.df    = df
        self.title = title

        container = self.query_one("#table-container")

        await container.remove_children()

        self.curr = DataTable()

        self.curr.add_columns(*self.df.columns.astype(str))

        for idx, row in enumerate(self.df.itertuples(index=False)):
            
            self.curr.add_row(*map(str, row), key=idx)

        self.curr.cursor_type       = "row"
        self.curr.styles.height     = 50
        self.curr.styles.background = "darkblue"
        self.curr.styles.color      = "white"
        self.curr.border_title      = self.title
        self.curr.styles.border     = ("heavy", "yellow")

        container.mount(self.curr)


    def action_download(self):

        downloads_path = os.path.join(os.path.expanduser("~"), "Downloads")
        file_name      = f"{self.title}.csv"
        full_path      = os.path.join(downloads_path, file_name)

        try:

            self.df.to_csv(full_path, index=False, sep=" ")
            self.notify(f"<*> Table saved to {full_path}.", severity="info")

        except:

            self.notify(f"<!> Error failed to save Table to {full_path}.")


    def action_fetch_mass_list(self):

        if ("Fragment Table" not in self.title):

            return

        mass_list      = fetchMassList(self.df).df
        mass_list_df   = pd.DataFrame(mass_list, columns=["Mass", "Multiplicity"])
        downloads_path = os.path.join(os.path.expanduser("~"), "Downloads")
        file_name      = f"{self.title.replace("Fragment Table", "")} Mass List.csv"
        full_path      = os.path.join(downloads_path, file_name)

        try:

            mass_list_df.to_csv(full_path, index=False, sep=" ")
            self.notify(f"<*> Mass List saved to {full_path}.", severity="info")

        except:

            self.notify(f"<!> Error failed to save Mass List to {full_path}.")


    async def action_go_back(self):

        if (self.parent_df is not None):

            await self.load_table(self.parent_df, self.parent_title)

            self.parent_df    = None
            self.parent_title = ""

        else:

            self.exit()


    async def _switch_table(self, new_df, new_title):

        await self.load_table(new_df, new_title)


    @on(DataTable.RowSelected)
    async def row_selected(self, event: DataTable.RowSelected):

        if (self.title == "MFD Index Table"):

            selected_idx = event.row_key.value
            selected_val = self.df.iloc[selected_idx, 1]
            new_df       = searchAndFetch(selected_val)

            self.parent_df    = self.df
            self.parent_title = self.title

            await self.load_table(new_df, f"{selected_val} Fragment Table")


@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx):

    if ctx.invoked_subcommand is None:

        click.echo("========================================================")
        click.echo("+ Welcome to the 3S2B-CLI v. 1.0.0                     +")
        click.echo("+                                                      +")
        click.echo("+ ~~~~                                                 +")
        click.echo("+                                                      +")
        click.echo("+ Jesse Fraser M.Sc. &                                 +")
        click.echo("+ Dr. Arun Moorthy Ph.D.                               +")
        click.echo("+                                                      +")
        click.echo("+ CRAFTS Lab | Trent University | 2025                 +")
        click.echo("+                                                      +")
        click.echo("+ <*>type 3s2b --help for a list of available commands +")
        click.echo("========================================================")

        return

    else:

        pass


@click.command()
@click.argument("file_paths", nargs=-1, type=click.Path(exists=True))
@click.argument("tol", type=float)
def c(file_paths, tol):

    ms_data_path, mass_list_path = file_paths

    # TODO:
    # Re-name file.
    ScrollableTable(pd.DataFrame([compare(*file_paths, tol=tol)], columns=["FPIE"]),
                    f"FPIE").run()


@click.command()
@click.argument("ms_data_path", type=click.Path(exists=True))
@click.argument("tol", type=float)
def a(ms_data_path, tol):

    result = compareAll(ms_data_path, tol=tol)
    result = result.sort_values(by="FPIE", ascending=False)

    ScrollableTable(result, "SWGDRUG SMILES w/ FPIEs").run()


@click.command()
@click.argument("search_term")
def f(search_term):
    
    try:

        ScrollableTable(searchAndFetch(search_term), f"{search_term} Fragment Table").run()

    except Exception as e:

        print(f"<*> {search_term} does not exist . . .")

@click.command()
@click.argument("mz")
def m(mz):

    ScrollableTable(searchAndFetchByMass(mz), f"Entries Containing: mz = {mz}").run()


@click.command()
@click.argument("args", nargs=-1)
def i(args):

    ScrollableTable(viewIdxTable(*args), "MFD Index Table").run()


@click.command()
@click.argument("smiles", type=str)
@click.argument("maxMult", type=int)
def fr(smiles, maxmult):

    # TODO:
    # Error handling
    fragmentor     = Fragmentor(maxMult=maxmult)
    fragmentor.mol = Chem.MolFromSmiles(smiles)

    return ScrollableTable(fragmentor.fetchAllFragsData(), f"{smiles} Fragment Table").run()


@click.command()
@click.argument("smiles")
def add(smiles):

    addEntryFromSmiles(smiles)


cli.add_command(i)
cli.add_command(f)
cli.add_command(m)
cli.add_command(c)
cli.add_command(a)
cli.add_command(fr)
cli.add_command(add)


if __name__ == "__main__":

    cli()
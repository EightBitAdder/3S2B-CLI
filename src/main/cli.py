from main.fragmentor import Fragmentor
from main.utils import fetchMassList, searchAndFetchByMass
from db.utils import addEntryFromSmiles, searchAndFetch, viewIdxTable
from main.comparator import (
    MassList,
    MSData,
    Comparator
)
import os
import re
import numpy as np
import pandas as pd
import sqlite3
import click
from textual.app import App
from textual import on
from textual.containers import Vertical
from textual.widgets import DataTable, Footer
from tqdm import tqdm
from rdkit import Chem


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
        file_name      = f"table.csv"
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

        if (self.title == "MFD Index Table" or self.title.startswith("Compare All >>> ")):

            selected_idx = event.row_key.value
            selected_val = self.df.iloc[event.row_key.value]["craftsLabEntry"]
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
@click.argument("file_paths", nargs=2, type=click.Path(exists=True))
@click.argument("tol", type=float)
@click.option("--sr", default="MS-DATA", help="<*> MS-Data file reader")
@click.option("--lr", default="MASS-LIST", help="<*> Mass-List file reader")
@click.option("--dr", default=None, help="<*> MS-Data file and Mass-List file delimiter")
@click.option("--wf", default=None, help="<*> Weight Function")
@click.option("--plot", is_flag=True, help="<*> Plot Annotated FPIE")
def c(file_paths, tol, sr, lr, dr, wf, plot):

    ms_data_path, mass_list_path = file_paths

    ms_data                 = MSData.fromFile(ms_data_path, delimiter=dr, fileReader=sr)
    mass_list               = MassList.fromFile(mass_list_path, delimiter=dr, fileReader=lr)
    comparator              = Comparator(ms_data, mass_list, tol, weightFunction=wf)
    FPIEScore, plotMetaData = comparator.calculateFPIE()

    if (plot):

        comparator.plotFPIE(
            plotMetaData,
            FPIEScore,
            os.path.splitext(os.path.basename(ms_data_path))[0]
        )

    ScrollableTable(pd.DataFrame([[FPIEScore]], columns=["FPIE"]), f"FPIE").run()


@click.command()
@click.argument("ms_data_path", type=click.Path(exists=True))
@click.argument("tol", type=float)
@click.option("--sr", default="MS-DATA", help="<*> MS-Data file reader")
@click.option("--dr", default=None, help="<*> MS-Data file and Mass-List file delimiter")
@click.option("--wf", default=None, help="<*> Weight Function")
def a(ms_data_path, tol, sr, dr, wf):

    ms_data           = MSData.fromFile(ms_data_path, delimiter=dr, fileReader=sr)
    idx_table_df      = viewIdxTable()
    crafts_lab_entrys = idx_table_df.iloc[:, 2]
    FPIEs             = []
    exact_mol_weights = []

    for entry in tqdm(crafts_lab_entrys, desc=f"<*> Calculating FPIEs . . ."):

        mass_list    = fetchMassList(searchAndFetch(entry))
        comparator   = Comparator(ms_data, mass_list, tol, weightFunction=wf)
        FPIEScore, _ = comparator.calculateFPIE()
        allFragsDF   = searchAndFetch(entry)

        FPIEs.append(FPIEScore)
        exact_mol_weights.append(allFragsDF.iloc[0, 2])

    df = pd.concat([idx_table_df.iloc[:, 0], pd.DataFrame({"FPIE": np.round(FPIEs, 2), "Exact_Mol_Wt": np.round(exact_mol_weights, 2)}), idx_table_df.iloc[:, 1]], axis=1)
    df = df.sort_values(by="FPIE", ascending=False)

    ScrollableTable(df, f"Compare All >>> {os.path.splitext(os.path.basename(ms_data_path))[0]}").run()


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

    ScrollableTable(searchAndFetchByMass(mz), f"MFD Index Table >>> By mz").run()


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
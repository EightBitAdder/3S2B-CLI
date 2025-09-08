from comparator import MassList
import pandas as pd


def fetchMassList(allFragsDF: pd.DataFrame) -> MassList:

    massList = []

    for row in allFragsDF.itertuples():

        Iso_Wts = list(map(float, row.Iso_Wts.split(", ")))
        masses  = [row.Exact_Mol_Wt, *Iso_Wts]

        massList.extend([mass, row.Mult] for mass in masses)

    df = pd.DataFrame(massList, columns=["Mass", "Multiplicity"])

    return MassList(df, weighted=True)
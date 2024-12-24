from typing import Tuple, Generator
from itertools import combinations, chain
import re
import pandas as pd
import IsoSpecPy as iso
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors


class Fragmentor():
    """
    A molecular mass fragment calculator.


    Attributes:
    -----------

    maxMult: int
        The maximum depth of recursive bond breaking.

    mol: Chem.rdchem.Mol
        The molecule to fragment.  
    """

    def __init__(self, *, maxMult: int=2) -> None:
        
        self.maxMult = maxMult


    @property
    def maxMult(self) -> int:

        return self._maxMult
    
    
    @maxMult.setter
    def maxMult(self, val: int) -> None:

        if (val < 0):
            
            raise ValueError(f"<!> {val} is not a positive integer.")
        
        self._maxMult = val


    @property
    def mol(self) -> Chem.rdchem.Mol:

        return self._mol
    

    @mol.setter
    def mol(self, val: Chem.rdchem.Mol) -> None:

        # TODO:
        # Error handling
        Chem.Kekulize(val, clearAromaticFlags=True)
        val = Chem.AddHs(val)

        self._mol = val
    

    #################
    # Fetch Methods #
    #################
    

    def _fetchBondIdxs(self) -> Generator[Tuple[int, ...], None, None]:
        """
        Fetches the bond indices to break.


        params:
        -------

        None


        returns:
        --------

        bondIdxs: Generator[Tuple[int, ...], None, None]
            The bond indices to break.
        """

        singleBondIdxs = [bond.GetIdx()
                          for bond in self.mol.GetBonds()
                          if (bond.GetBondType() == Chem.BondType.SINGLE)]
        bondIdxs = chain(*[combinations(singleBondIdxs, mult)
                          for mult in range(1, self.maxMult + 1)])

        return bondIdxs
    

    def fetchAllFrags(self) -> list[Tuple[Chem.rdchem.Mol, int]]:
        """
        Fetches all of the fragments, as well as their multiplicities.


        params:
        -------

        None


        returns:
        --------

        fragsWithMult: list[Tuple[Chem.rdchem.Mol, int]]
            The fragments, as well as their multiplicities.
        """

        fragsWithMult = [(self.mol, 0)]

        for bondIdx in self._fetchBondIdxs():

            frags = Chem.FragmentOnBonds(self.mol, [*bondIdx])
            fragMols = Chem.GetMolFrags(frags, asMols=True)
            fragsWithMult.extend([(fragMol, len(bondIdx))
                                  for fragMol in fragMols])

        return fragsWithMult
    

    def _cleanMolSmiles(self, mol: Chem.rdchem.Mol) -> str:
        """
        Cleans the SMILES output; nominally, after calls to:
        Chem.FragmentOnBonds() and Chem.GetMolFrags(),
        in fetchAllFrags().


        params:
        -------

        mol: Chem.rdchem.Mol
            The molecule from which the SMILES output is generated.

        
        returns:
        --------

        cleanMolSmiles: str
            The cleaned SMILES output.
        """

        cleanMolSmiles = re.sub(r"\(\[\d+\*\]\)|\[\d+\*\]|\*",
                                "",
                                Chem.MolToSmiles(mol))
        
        return cleanMolSmiles
    

    def _cleanMolFormula(self, mol: Chem.rdchem.Mol) -> str:
        """
        Cleans the formula output; nominally, after calls to:
        Chem.FragmentOnBonds() and Chem.GetMolFrags(),
        in fetchAllFrags().


        params:
        -------
        
        mol: Chem.rdchem.Mol
            The molecule from which the formula is generated.

        
        returns:
        --------

        cleanMolFormula: str
            The cleaned formula output.
        """

        cleanMolFormula = re.sub(r"\*\d+|\*|\+|\-",
                                 "",
                                 rdMolDescriptors.CalcMolFormula(mol))

        return cleanMolFormula
    

    def _fetchIsoMasses(self, mol: Chem.rdchem.Mol) -> list[float]:
        """
        Fetches the masses of all of the isotopes
        (up to a given cover-probability)
        of a given molecule.

        
        params:
        -------
        
        mol: Chem.rdchem.Mol
            The molecule from which the masses of all of its isotopes
            (up to a given cover-probability)
            will be calculated.

        
        returns:
        --------

        isoMasses: list[float]:
            The masses of all of the isotopes
            (up to a given cover-probability)
            of the given molecule.
        """

        isoMasses = []
        sp = iso.IsoTotalProb(formula=self._cleanMolFormula(mol),
                              prob_to_cover=0.9999)

        for mass, prob in sp:
            if (prob > 0.01):
                isoMasses.append(mass)

        return isoMasses
    

    def fetchAllFragsData(self) -> pd.DataFrame:
        """
        Fetches all of the relevant data
        for each of the fragments generated by:
        fetchAllFrags().


        params:
        -------

        None


        returns:
        --------

        allFragsDF: pd.DataFrame
            All of the relevant data
            for each of the fragments generated by:
            fetchAllFrags().
        """

        allFragsData = [
            *[
                [self._cleanMolSmiles(frag[0]),
                 self._cleanMolFormula(frag[0]),
                 Descriptors.ExactMolWt(frag[0]),
                 ", ".join(map(str, self._fetchIsoMasses(frag[0]))),
                 frag[1]]
                 for frag in self.fetchAllFrags()
            ]
        ]

        allFragsDF = pd.DataFrame(allFragsData,
                                  columns=["SMILES",
                                          "Formula",
                                          "Exact_Mol_Wt",
                                          "Iso_Wts",
                                          "Mult"])
        allFragsDF.drop_duplicates(inplace=True)

        return allFragsDF
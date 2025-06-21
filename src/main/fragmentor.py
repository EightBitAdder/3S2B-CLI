from typing import List, Tuple, Dict
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
        
        self.maxMult                      = maxMult
        self._smilesCache: Dict[str, str] = {}


    @property
    def maxMult(self) -> int:

        return self._maxMult
    
    
    @maxMult.setter
    def maxMult(self, val: int) -> None:

        if (val < 1):
            
            raise ValueError(f"<!> Error: {val} is not a positive integer.")
        
        self._maxMult = val


    @property
    def mol(self) -> Chem.rdchem.Mol:

        return self._mol
    

    @mol.setter
    def mol(self, val: Chem.rdchem.Mol) -> None:

        try:

            Chem.Kekulize(val, clearAromaticFlags=True)
        
        except Chem.rdchem.KekulizeException as e:

            raise ValueError(f"<!> Error: kekulization failed >>> {e}")

        val       = Chem.AddHs(val)
        self._mol = val
    

    #################
    # Fetch Methods #
    #################
    
    '''
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
        bondIdxs       = chain(*[combinations(singleBondIdxs, mult)
                                for mult in range(1, self.maxMult + 1)])

        return bondIdxs
    '''
    def _fetchBondIdxs(self) -> List[int]:

        return [bond.GetIdx() for bond in self._mol.GetBonds()
                if bond.GetBondType() == Chem.BondType.SINGLE]
    

    '''
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

            frags    = Chem.FragmentOnBonds(self.mol, [*bondIdx])
            fragMols = Chem.GetMolFrags(frags, asMols=True)

            fragsWithMult.extend([(fragMol, len(bondIdx))
                                  for fragMol in fragMols])

        return fragsWithMult
    '''
    def fetchAllFrags(self) -> Dict[str, Tuple[List[Chem.Mol], Tuple[int, ...]]]:

        seenResultants = {}

        try:

            seenResultants[self._cleanMolSmiles(self._mol)] = ([self._mol], tuple())

        except Exception as e:

            print(f"<*> Skipped fragment due to RDKit error >>> {e}")

        bondIdxs = self._fetchBondIdxs()

        for k in range(1, self.maxMult + 1):

            for combo in combinations(bondIdxs, k):

                try:

                    frags      = Chem.FragmentOnBonds(self._mol, list(combo))
                    fragMols   = Chem.GetMolFrags(frags, asMols=True)
                    fragSMILES = sorted(self._cleanMolSmiles(frag) for frag in fragMols if frag is not None)
                    key        = "::".join(fragSMILES)

                except:

                    continue

                if (not key or any(not smiles for smiles in fragSMILES)):

                    continue

                if (key not in seenResultants):

                    seenResultants[key] = (fragMols, combo)

        return seenResultants


    r'''
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
    '''
    def _cleanMolSmiles(self, mol: Chem.Mol) -> str:

        smiles = Chem.MolToSmiles(mol, canonical=True)

        if (smiles not in self._smilesCache):

            cleanMolSmiles            = re.sub(r"\(\[\d+\*\]\)|\[\d+\*\]|\*",
                                               "",
                                               smiles)
            self._smilesCache[smiles] = cleanMolSmiles

        return self._smilesCache[smiles]
    

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
        formula   = formula=self._cleanMolFormula(mol)

        if (not formula):

            return []

        try:

            sp = iso.IsoTotalProb(formula=formula, prob_to_cover=0.9999)

            for mass, prob in sp:

                if (prob > 0.01):

                    isoMasses.append(mass)

        except Exception:

            return []

        return isoMasses
    
    '''
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

        allFragsData = []

        for frag in self.fetchAllFrags():

            try:

                smiles  = self._cleanMolSmiles(frag[0])
                formula = self._cleanMolFormula(frag[0])

                if (not smiles or not formula):

                    continue

                exactWt   = Descriptors.ExactMolWt(frag[0])
                isoMasses = ", ".join(map(str, self._fetchIsoMasses(frag[0])))

                allFragsData.append([smiles, formula, exactWt, isoMasses, frag[1]])

            except Exception as e:

                print(f"<*> Skipped fragment due to RDKit error >>> {e}")

                continue

        allFragsDF = pd.DataFrame(allFragsData,
                                  columns=["SMILES",
                                           "Formula",
                                           "Exact_Mol_Wt",
                                           "Iso_Wts",
                                           "Mult"])
        
        allFragsDF.drop_duplicates(inplace=True)

        return allFragsDF
    '''
    def fetchAllFragsData(self) -> pd.DataFrame:

        allFragsData = []
        resultants   = self.fetchAllFrags()

        for fragMols, combo in resultants.values():

            for frag in fragMols:

                try:

                    smiles  = self._cleanMolSmiles(frag)
                    formula = self._cleanMolFormula(frag)

                    if (not smiles or not formula):

                        continue

                    exactWt   = Descriptors.ExactMolWt(frag)
                    isoMasses = ", ".join(map(str, self._fetchIsoMasses(frag)))

                    allFragsData.append([smiles, formula, exactWt, isoMasses, len(combo)])

                except Exception as e:

                    print(f"<*> Skipped fragment due to RDKit error >>> {e}")

                    continue

        return pd.DataFrame(
            allFragsData,
            columns=["SMILES",
                     "Formula",
                     "Exact_Mol_Wt",
                     "Iso_Wts",
                     "Mult"]
        ).drop_duplicates()
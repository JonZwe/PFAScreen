"""
Most functions from GraphMZ (written by Robert Young)
Structure-related calculation utilities.

This module provides functions for working with molecular structures,
including SMILES parsing, formula generation, and fragment analysis.
Further down: Additions from Jonathan Zweigle
"""

from typing import TypeAlias

import pandas as pd
from rdkit.Chem import Descriptors, inchi, rdmolops, RWMol, GetPeriodicTable
from rdkit.Chem.MolStandardize.rdMolStandardize import Uncharger
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles
from rdkit.Chem.rdmolops import GetFormalCharge, GetMolFrags
import pubchempy as pcp
import pyopenms as oms
from tqdm import tqdm

# Type aliases
FragmentList: TypeAlias = list[str]


def canonicalize_smiles(smiles: str) -> str:
    """
    Convert SMILES to canonical form using RDKit.

    Args:
        smiles: Input SMILES string to canonicalize

    Returns:
        Canonical SMILES string if successful, empty string if conversion fails

    Notes:
        This function acts as a validator for SMILES strings.
        If RDKit cannot parse the SMILES, an empty string is returned.
    """
    try:
        mol = MolFromSmiles(smiles)
        if mol is None:
            print(f"Cannot parse SMILES: {smiles}")
            return ""
        return MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception as e:
        print(f"Error canonicalizing SMILES '{smiles}': {str(e)}")
        return ""


def get_mol(smiles: str) -> Mol:
    """
    Convert SMILES string to RDKit molecule.

    Args:
        smiles: Input SMILES string

    Returns:
        RDKit Mol object

    Raises:
        ValueError: If SMILES string is invalid
    """
    mol = MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return mol


def get_inchi_from_smiles(smiles: str) -> str:
    """
    Return the standard InChI string for a SMILES string.

    Args:
        smiles: SMILES string to convert to InChI

    Returns:
        InChI string representation of the molecule

    Raises:
        ValueError: If SMILES string is invalid or cannot be converted
    """
    try:
        mol = MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        inchi_str = inchi.MolToInchi(mol)
        if not inchi_str:
            raise ValueError(f"Failed to generate InChI from SMILES: {smiles}")
        return inchi_str
    except Exception as e:
        print(f"Error converting SMILES to InChI '{smiles}': {str(e)}")
        raise ValueError(f"Failed to convert SMILES to InChI: {smiles}") from e


def get_inchikey_from_smiles(smiles: str) -> str:
    """
    Return the standard InChIKey for a SMILES string.

    Args:
        smiles: SMILES string to convert to InChIKey

    Returns:
        InChIKey string representation of the molecule

    Raises:
        ValueError: If SMILES string is invalid or cannot be converted
    """
    try:
        mol = MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        inchikey = inchi.MolToInchiKey(mol)
        if not inchikey or len(inchikey) != 27:  # Standard InChIKey length
            raise ValueError(f"Failed to generate valid InChIKey from SMILES: {smiles}")
        return inchikey
    except Exception as e:
        print(f"Error converting SMILES to InChIKey '{smiles}': {str(e)}")
        raise ValueError(f"Failed to convert SMILES to InChIKey: {smiles}") from e


def get_mol_from_inchi(
    inchi_str: str, sanitize: bool = True, removeHs: bool = True
) -> Mol:
    """
    Construct an RDKit Mol object from an InChI string.

    Args:
        inchi_str: InChI string to convert to RDKit molecule
        sanitize: Whether to sanitize the molecule during creation
        removeHs: Whether to remove hydrogens from the molecule

    Returns:
        RDKit Mol object

    Raises:
        ValueError: If InChI string is invalid or cannot be parsed
    """
    try:
        mol = inchi.MolFromInchi(inchi_str, sanitize=sanitize, removeHs=removeHs)
        if mol is None:
            raise ValueError(f"Unable to parse InChI string: {inchi_str}")
        return mol
    except Exception as e:
        print(f"Error creating molecule from InChI '{inchi_str}': {str(e)}")
        raise ValueError(f"Failed to create molecule from InChI: {inchi_str}") from e


def has_wildcards(smiles: str) -> bool:
    """
    Check if SMILES contains wildcards.

    Args:
        smiles: SMILES string to check

    Returns:
        Boolean indicating if wildcards are present
    """
    if pd.isna(smiles):  # Check for NaN/None
        return False
    return "*" in smiles


def neutralize_smiles(smiles: str) -> str:
    """
    Convert charged SMILES to neutral form.

    Args:
        smiles: SMILES string to neutralize

    Returns:
        Neutralized SMILES string

    Raises:
        ValueError: If SMILES is invalid or cannot be neutralized
    """
    try:
        mol = get_mol(smiles)
        uncharger = Uncharger()
        neutral_mol = uncharger.uncharge(mol)
        return MolToSmiles(neutral_mol, isomericSmiles=True)
    except ValueError as e:
        print(f"Error neutralizing {smiles}: {str(e)}")
        return smiles


def get_charge_smiles(smiles: str) -> int | None:
    """
    Get formal charge from SMILES.

    Args:
        smiles: SMILES string to calculate charge for

    Returns:
        Formal charge as integer, or None if calculation fails
    """
    try:
        mol = get_mol(smiles)
        return GetFormalCharge(mol)
    except ValueError as e:
        print(f"Error calculating charge for {smiles}: {str(e)}")
        return None


def smiles_to_formula(smiles: str) -> str:
    """
    Get molecular formula from SMILES.

    Args:
        smiles: SMILES string to convert

    Returns:
        Molecular formula as string, or empty string if conversion fails
    """
    try:
        mol = get_mol(smiles)
        return CalcMolFormula(mol, separateIsotopes=True, abbreviateHIsotopes=False)
    except ValueError as e:
        print(f"Error calculating formula for {smiles}: {str(e)}")
        return ""


def get_base_formula(smiles: str) -> str:
    """
    Get formula with isotopes converted to base elements.

    Args:
        smiles: SMILES string to convert

    Returns:
        Molecular formula with isotopes converted to base elements

    Raises:
        ValueError: If SMILES is invalid
    """
    mol = get_mol(smiles)
    return CalcMolFormula(mol, separateIsotopes=False)


def smiles_to_mass(smiles: str) -> float | None:
    """
    Get exact mass from SMILES.

    Args:
        smiles: SMILES string to calculate mass for

    Returns:
        Exact molecular mass as float, or None if calculation fails
    """
    try:
        mol = get_mol(smiles)
        return Descriptors.ExactMolWt(mol)  # type: ignore
    except ValueError as e:
        print(f"Error calculating mass for {smiles}: {str(e)}")
        return None


def get_unique_fluorine_fragments(smiles: str) -> list[str]:
    """
    Get F-containing fragments from SMILES.

    Args:
        smiles: SMILES string to analyze

    Returns:
        List of unique SMILES fragments containing fluorine atoms
    """
    try:
        mol = get_mol(smiles)
        fragments = GetMolFrags(mol, asMols=True)
        fluorine_fragments = [
            MolToSmiles(frag, isomericSmiles=True)
            for frag in fragments
            if any(atom.GetSymbol() == "F" for atom in frag.GetAtoms())
        ]
        unique = set(fluorine_fragments)
        return list(unique)
    except ValueError as e:
        print(f"Error processing fragments for {smiles}: {str(e)}")
        return []


def get_unique_carbon_fragments(smiles: str) -> list[str]:
    """
    Get C-containing fragments from SMILES and ignore elemental ions.

    Args:
        smiles: SMILES string to analyze

    Returns:
        List of unique SMILES fragments containing carbon atoms
    """
    try:
        mol = get_mol(smiles)
        fragments = GetMolFrags(mol, asMols=True)
        carbon_fragments = [
            MolToSmiles(frag, isomericSmiles=True)
            for frag in fragments
            if any(atom.GetSymbol() == "C" for atom in frag.GetAtoms())
        ]
        unique = set(carbon_fragments)
        return list(unique)
    except ValueError as e:
        print(f"Error processing fragments for {smiles}: {str(e)}")
        return []
    

# From here on my function come:

def is_charge_possible(smiles: str) -> bool:
    """
    Heuristic check if a molecule can be or already is charged.
    Rules:
      - If the SMILES encodes formal charges → return True
      - If all atoms are in (C, H, F, Cl, Br, I) → return False
      - Otherwise (heteroatoms, metals, etc.) → return True
    """
    # Define the safe set: atoms that form neutral molecules with carbon
    safe_set = {"C", "H", "F", "Cl", "Br", "I"}
    
    mol = MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    
    # Check for explicit charges
    if any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms()):
        return True
    
    # Collect unique element symbols
    elements = {atom.GetSymbol() for atom in mol.GetAtoms()}
    
    # If only safe atoms are present → cannot ionize
    if elements.issubset(safe_set):
        return False
    
    # Otherwise, charge is possible
    return True


def is_salt(smiles: str) -> bool:
    """
    Returns True if the SMILES contains multiple disconnected fragments (a salt).
    """
    mol = MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    
    # Get fragments as separate molecules
    frags = GetMolFrags(mol, asMols=True)
    
    # Salt if more than one fragment
    return len(frags) > 1


def has_fluorine(smiles: str) -> bool:
    """
    Check if SMILES contains fluorine atoms.

    Args:
        smiles: SMILES string to check

    Returns:
        Boolean indicating if fluorine atoms are present
    """
    mol = MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return any(atom.GetSymbol() == "F" for atom in mol.GetAtoms())


def has_element_formula(formula: str, element: str) -> bool:
    """
    Check if a given element is present in a molecular formula string.

    Args:
        formula: Molecular formula as a string (e.g., 'C7H8O4')
        element: Element symbol as a string (e.g., 'C', 'N', 'F')
    Returns:
        True if the element is present, False otherwise.
    """
    try:
        ef = oms.EmpiricalFormula(formula)
        # EmpiricalFormula.getElementalComposition() returns a dict with byte keys
        comp = ef.getElementalComposition()
        return element.encode() in comp
    except Exception as e:
        print(f"Error parsing formula '{formula}': {e}")
        return False
    

def has_element(smiles: str, element: str) -> bool:
    """
    Check if SMILES contains a specific element.

    Args:
        smiles: SMILES string to check
        element: Element symbol to look for (e.g., 'F', 'Cl', 'N')
    Returns:
        Boolean indicating if the specified element is present
    """
    mol = MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return any(atom.GetSymbol() == element for atom in mol.GetAtoms())


def has_formal_charge(formula: str) -> bool:
    """
    Check if a molecular formula indicates a charged species.
    """
    # Look for common charge indicators in the formula
    charge_indicators = ['+', '-']
    return any(indicator in formula for indicator in charge_indicators)


def get_charge_str(formula: str) -> str:
    """
    Get "positive", "negative", or "neutral" from formula string.
    """
    if '+' in formula and '-' in formula:
        return "neutral"  # Both charges present, likely neutral overall
    elif '+' in formula:
        return "positive"
    elif '-' in formula:
        return "negative"
    else:
        return "neutral"


def neutralize_smiles_fluorine(smiles: str) -> str:
    """
    NOTE: DEVELOPMENT STAGE, NOT FULLY TESTED YET!

    Neutralize a salt SMILES:
      - Keeps largest fluorine containing fragment
      - Adds hydrogens to negative charges, removes extra hydrogens from positives
      - Quaternary ammonium groups are not neutralized
      - Returns single neutral SMILES
    """
    # Pick largest fluorine-containing fragment
    fluorine_fragments = get_unique_fluorine_fragments(smiles)
    fluorine_fragments_mol = [get_mol(frag) for frag in fluorine_fragments]

    if len(fluorine_fragments_mol) > 1:
        print('Warning: SMILES contains multiple F fragments, neutralizing to largest fluorine-containing fragment.')
    if not fluorine_fragments_mol:
        # fallback: take largest fragment
        target_frag = max(fluorine_fragments_mol, key=lambda m: m.GetNumAtoms())
    else:
        target_frag = max(fluorine_fragments_mol, key=lambda m: m.GetNumAtoms())

    # Make a copy to modify
    mol_writable = RWMol(target_frag)
    for atom in mol_writable.GetAtoms():
        charge = atom.GetFormalCharge()
        if charge == 0:
            continue

        # Calculate total bonds
        total_bonds = atom.GetExplicitValence() + atom.GetImplicitValence()
        max_valence = GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())

        if charge > 0 and total_bonds >= max_valence:
            print(f"Warning: permanent positive charge on atom {atom.GetSymbol()} at index {atom.GetIdx()} kept.")
            continue

        if charge < 0 and total_bonds >= max_valence:
            print(f"Warning: permanent negative charge on atom {atom.GetSymbol()} at index {atom.GetIdx()} kept.")
            continue

        # Neutralize
        atom.SetFormalCharge(0)
        if charge > 0:
            atom.SetNumExplicitHs(atom.GetTotalNumHs() + charge)
        else:
            atom.SetNumExplicitHs(atom.GetTotalNumHs() - charge)


    # Sanitize and return SMILES
    rdmolops.SanitizeMol(mol_writable)
    return MolToSmiles(mol_writable, isomericSmiles=True)


def pubchem_names_to_df(names, properties=None) -> pd.DataFrame:
    """
    One row per input name (always).
    Flags multiple PubChem hits instead of expanding rows.
    Uses the first hit deterministically when multiple hits occur.
    """

    if properties is None:
        properties = [
                    "cid",
                    "molecular_formula",
                    "molecular_weight",
                    "exact_mass",
                    "inchi",
                    "inchikey",
                    "smiles",
                    "xlogp",
                    "tpsa",
                    "charge",
                    "rotatable_bond_count",
                    "h_bond_donor_count",
                    "heavy_atom_count"
                     ]

    rows = []

    for name in tqdm(names, desc="Running pubchempy queries"):
        try:
            cs = pcp.get_compounds(name, namespace="name")
        except Exception as e:
            rows.append({
                "query_name": name,
                "found": False,
                "n_hits": 0,
                "multiple_hits": False,
                "error": str(e)
            })
            continue

        # no hits
        if not cs:
            rows.append({
                "query_name": name,
                "found": False,
                "n_hits": 0,
                "multiple_hits": False
            })
            continue

        # one or more hits
        df_hits = pcp.compounds_to_frame(cs, properties=properties)

        rows.append({
            "query_name": name,
            "found": True,
            "n_hits": len(df_hits),
            "multiple_hits": len(df_hits) > 1,
            # take first hit deterministically
            **df_hits.iloc[0].to_dict()
        })

    return pd.DataFrame(rows)

from pyteomics.mass import calculate_mass, Composition
from itertools import product
from collections import defaultdict

def match_subformulas_to_peaks(formula_str, 
                               peaks, 
                               mass_tolerance=0.005, 
                               max_combinations=100_000):
    """
    Match subformulas of a molecular formula to MS2 peaks within a mass tolerance.
    Uses peak indices as keys instead of m/z values.

    Parameters:
        formula_str (str): e.g. "C8H10O3S"
        peaks (list): list of MS2 m/z values
        mass_tolerance (float): allowed error in Da (e.g., 0.005)
        max_combinations (int): maximum allowed subformulas before stopping

    Returns:
        dict: {peak_index: [ {formula, mass, delta}, ... ] }
    """
    # ensure that charges are removed!
    formula_str = formula_str.replace("+", "").replace("-", "")
    
    base_formula = Composition(formula_str)
    max_counts = {atom: base_formula[atom] for atom in base_formula}

    atoms = sorted(max_counts.keys())
    ranges = [range(max_counts[a] + 1) for a in atoms]

    # --- Estimate number of combinations
    total_combinations = 1
    for r in ranges:
        total_combinations *= len(r)

    if total_combinations > max_combinations:
        print(f"Stopped! Formula {formula_str} would create {total_combinations:,} combinations (limit = {max_combinations:,})")
        return {}

    # --- Store results using peak indices
    matches = defaultdict(list)

    for counts in product(*ranges):
        if sum(counts) == 0:
            continue
        sub = {atom: count for atom, count in zip(atoms, counts)}
        sub_mass = calculate_mass(sub)

        for idx, peak in enumerate(peaks):
            delta = abs(peak - sub_mass)
            if delta <= mass_tolerance:
                # Remove "1" from formula string (e.g., C6H5O1 → C6H5O)
                sub_formula_str = ''.join(
                    f"{k}" if v == 1 else f"{k}{v}" for k, v in sub.items() if v > 0
                )
                matches[idx].append({
                    'formula': sub_formula_str,
                    'mass': round(sub_mass, 5),
                    'delta': round(delta, 5)
                })

    return dict(matches)

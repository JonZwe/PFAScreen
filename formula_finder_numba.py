import numpy as np
from IsoSpecPy import IsoThreshold
from collections import defaultdict
import pubchempy as pcp
import re
import numba
from typing import List, Dict
from numba import njit, types
from numba.typed import Dict as NumbaDict
from numba.typed import List as NumbaList

allowed_elements = ['C', 'H', 'N', 'O', 'F', 'Cl', 'Br', 'I', 'P', 'S', 'Si']

# Element limits
MAX_ELEMENT_COUNTS: Dict[str, int] = {
    'C': 50, 'H': 100, 'N': 6, 'O': 6, 'S': 4, 'P': 1,
    'F': 36, 'Cl': 4, 'Br': 4, 'I': 4, 'Si': 1
}

# Monoisotopic masses
ELEMENT_MASSES: Dict[str, float] = {
    'C': 12.0000000,
    'H': 1.00782503223,
    'N': 14.00307400443,
    'O': 15.99491461957,
    'Si': 27.976928,
    'S': 31.9720711743,
    'P': 30.9737619984,
    'F': 18.9984031627,
    'Cl': 34.968852682,
    'Br': 78.9183376,
    'I': 126.904477
}
def parse_formula_regex(formula_str):
    """
    Parses a chemical formula string into a dictionary of element counts.
    E.g., "C12H14N2F12Cl" -> {'C': 12, 'H': 14, 'N': 2, 'F': 12, 'Cl': 1}
    """
    # Regex explanation: 
    # ([A-Z][a-z]?) : Captures element symbol (e.g., C, Cl, Mg). 
    # (\d*)         : Captures the following number (count).
    # g modifier ensures it finds all occurrences.
    matches = re.findall(r'([A-Z][a-z]?)(\d*)', formula_str)
    
    counts = defaultdict(int)
    for element, count_str in matches:
        # Convert the count string to an integer, defaulting to 1 if empty
        count = int(count_str) if count_str else 1
        counts[element] += count
        
    return counts

def get_fragment_match_count(fragment_list, main_formula_str):
    """
    Counts how many fragments in a list are a subset of the main formula using built-in modules
    """
    try:
        main_counts = parse_formula_regex(main_formula_str)
    except Exception as e:
        print(f"Error parsing main formula: {e}")
        return "0/0"

    try:
        total_fragments = len(fragment_list)
    except TypeError:
        return "0/0" 

    match_count = 0
    total_fragments = len(fragment_list)

    for sub_formula_str in fragment_list:
        try:
            sub_counts = parse_formula_regex(sub_formula_str)
            is_a_match = True

            for element, sub_count in sub_counts.items():
                main_count = main_counts.get(element, 0)
                if sub_count > main_count:
                    is_a_match = False
                    break # Stop checking this fragment

            if is_a_match:
                match_count += 1
                
        except Exception as e:
            # Handle cases where the regex fails (e.g., malformed input)
            print(f"Skipping invalid fragment '{sub_formula_str}': {e}")
            continue

    return f"{match_count}/{total_fragments}"


def bin_isotopes(formula, resolution, threshold=0.0001):
    sp = IsoThreshold(formula=formula, threshold=threshold)
    peaks = list(sp)
    peaks.sort(key=lambda x: x[0])

    binned = []
    current_mass, current_prob = peaks[0]

    for mass, prob in peaks[1:]:
        delta_m = current_mass / resolution
        if abs(mass - current_mass) <= delta_m:
            current_mass = (current_mass * current_prob + mass * prob) / (current_prob + prob)
            current_prob += prob
        else:
            binned.append((current_mass, current_prob))
            current_mass, current_prob = mass, prob
    binned.append((current_mass, current_prob))
    return binned

def isotope_score(intens_array_theoretical, intens_array_measured):

    # Adjust length of theoretic intensities (M+1, M+2, ...)
    target_len = len(intens_array_measured)
    if len(intens_array_theoretical) > target_len:
        intens_array_theoretical = intens_array_theoretical[:target_len]
    else:
        intens_array_theoretical += [0.0] * (target_len - len(intens_array_theoretical))

    # Normalize
    intens_array_measured_normalized = intens_array_measured / np.sum(intens_array_measured)
    intens_difference = abs(intens_array_theoretical - intens_array_measured_normalized)
    score = 1
    for n in range(len(intens_difference)):
        score *= (1 - intens_difference[n])
            
    return score

def match_spectrum_to_formulas_numba(mz, ppm, mode, isotope_score_lim, unsat_min, unsat_max, measured_intensities, dia_frags, resolution=1000, PubChem=True):

    if mode == '[M+H]+':
        target_mass = mz - 1.0072
    if mode == '[M-H]-':
        target_mass = mz + 1.0072

    print('Start formula finder...')

    formulas_with_details = find_formulas_with_info_numba(target_mass, ppm, unsat_min, unsat_max, allowed_elements, MAX_ELEMENT_COUNTS)
    scored_matches = []

    for entry in formulas_with_details:
        formula = entry['formula']
        unsat = entry['unsaturation_dbe']
        ppm_error = entry['ppm_error']

        binned = bin_isotopes(formula, resolution)
        theoretical_intensities = [p[1] for p in binned]

        score = isotope_score(theoretical_intensities, measured_intensities)
        penalty = mismatch_penalty(theoretical_intensities, measured_intensities, target_mass, formula)

        adjusted_score = score - penalty
        scored_matches.append({
            'formula': formula, 
            'score': adjusted_score, 
            'ppm_error': ppm_error, 
            'unsaturation': unsat
        })

    # Sort the full list by score
    scored_matches.sort(key=lambda x: x['score'], reverse=True)

    # Create a new list containing only matches above the threshold
    filtered_matches = [match for match in scored_matches if match['score'] > isotope_score_lim]
    
    print(f"Total formulas found matching mass/DBE: {len(scored_matches)}")
    print(f"Formulas matching score limit ({isotope_score_lim}): {len(filtered_matches)}")

    
    
    from pubchempy import get_cids
    formula_cache = {}
    
    # Iterate over the filtered list
    for match in filtered_matches:
        formula = match['formula']
        score = match['score']
        ppm_error = match['ppm_error']
        unsat = match['unsaturation']

        if PubChem==True:
            if formula in formula_cache:
                cids = formula_cache[formula]
            else:
                try:
                    cids = get_cids(formula, 'formula')
                    formula_cache[formula] = cids
                except Exception as e:
                    cids = []
                    formula_cache[formula] = []
        else:
            cids = []

        
        # Assuming get_fragment_match_count is defined
        possible_fragments = get_fragment_match_count(dia_frags, formula)
        
        print(f"Formula: {formula}, Score: {score:.5f}, ppm error: {ppm_error:.2f}, Unsaturation: {unsat}, Pubchem entries: {len(cids)}, Possible fragments: {possible_fragments}")
    print('Done!')
    # return filtered_matches

 

# Mismatch-Penalty: Checks if a theoretical peak has a relevant intensity (e.g., > 5%) but is missing in the measured spectrum (e.g., < 1%). Thresholds can be set in the function below
def mismatch_penalty(theoretical, measured, target_mass, formula, threshold_theoretical=0.05, threshold_measured=0.01):
    
    penalty = 0.0
    import re
    def count_atoms(formula, element):
        matches = re.findall(rf'{element}(\d*)', formula)
        count = 0
        for match in matches:
            count += int(match) if match else 1
        return count

    for i in range(len(theoretical)):
        theo = theoretical[i]
        measured_norm = [x / sum(measured) for x in measured]
        meas = measured_norm[i] if i < len(measured) else 0.0

        # Mismatch: strong theoretical peak, but no measured peak
        if theo > threshold_theoretical and meas < threshold_measured:
            penalty += theo

    # M+2-Peak-Check: Measured M+2 < 50% of theoretical M+2
    if len(theoretical) > 2:
        theo_m2 = theoretical[2]
        meas_m2 = measured_norm[2] if len(measured) > 2 else 0.0
        if theo_m2 > 0.05 and meas_m2 < 0.5 * theo_m2:
            penalty += theo_m2 # NOTE hoher Faktor so!! Aufpassen!

    #### Apply Rules from Fiehn & Kind

    # Penalty (N-rule) NOTE: Hier evtl noch anpassen für Massen über 500 Da weil die Rundung nicht stimmt?
    rounded_mass = round(target_mass)
    n_count = count_atoms(formula, 'N')
    if (rounded_mass % 2 == 0 and n_count % 2 != 0) or (rounded_mass % 2 != 0 and n_count % 2 == 0):
        penalty += 0.5 
    
    # Penalty (element ratios)
    c_count = count_atoms(formula, 'C')
    if c_count == 0:
        penalty += 1
    else:
        h_count = count_atoms(formula, 'H')
        f_count = count_atoms(formula, 'F')
        cl_count = count_atoms(formula, 'Cl')
        br_count = count_atoms(formula, 'Br')
        o_count = count_atoms(formula, 'O')
        p_count = count_atoms(formula, 'P')
        s_count = count_atoms(formula, 'S')
        if h_count / c_count > 3 or f_count/c_count > 6 or cl_count/c_count > 2 or br_count/c_count > 2 or n_count/c_count > 4 or o_count/c_count > 3 or p_count/c_count > 2 or s_count/c_count > 3: #removed:  or h_count / c_count < 0.125
            penalty += 0.5


    return penalty


# Convert Python Dict to Numba Dict
NUMBA_ELEMENT_MASSES = NumbaDict.empty(key_type=types.unicode_type, value_type=types.float64)
for k, v in ELEMENT_MASSES.items():
    NUMBA_ELEMENT_MASSES[k] = v

# ------------------------------------------------------------
# 1) Numba-Helper: DBE = (2C + 2 + N - H - X) / 2
#    - Indices c_idx, h_idx, n_idx (=-1 if element is missing)
#    - X = sum counts[j] for all halogens (halogen_mask[j] == True)
# ------------------------------------------------------------
@njit(cache=True)
def compute_dbe_counts_CHNX(counts: np.ndarray,
                            c_idx: int, h_idx: int, n_idx: int,
                            halogen_mask: np.ndarray) -> float:
    C = counts[c_idx] if c_idx >= 0 else 0
    H = counts[h_idx] if h_idx >= 0 else 0
    N = counts[n_idx] if n_idx >= 0 else 0

    X = 0
    for j in range(counts.shape[0]):
        if halogen_mask[j]:
            X += counts[j]

    # DBE-formula
    return (2.0 * C + 2.0 + N - H - X) / 2.0


# ------------------------------------------------------------
# 2) Numba-Helper: recursive search with pruning
# ------------------------------------------------------------
@njit(cache=True)
def recursive_search_numba(current_counts: np.ndarray,
                           current_mass: float,
                           element_index: int,
                           element_masses: np.ndarray,
                           element_max_counts: np.ndarray,
                           suffix_max_mass: np.ndarray,
                           lower_bound: float,
                           upper_bound: float,
                           target_mz: float,
                           unsat_min: float,
                           unsat_max: float,
                           # DBE-parameters:
                           c_idx: int, h_idx: int, n_idx: int,
                           halogen_mask: np.ndarray,
                           # results:
                           results_ppm,        # typed.List[float64]
                           results_dbe,        # typed.List[float64]
                           results_counts      # typed.List[np.ndarray[int64, 1d]]
                           ):
    n = element_masses.shape[0]

    # Pruning
    if current_mass + suffix_max_mass[element_index] < lower_bound:
        return
    if current_mass > upper_bound:
        return

    if element_index == n:
        if lower_bound <= current_mass <= upper_bound:
            ppm_err = 1e6 * (current_mass - target_mz) / target_mz
            dbe = compute_dbe_counts_CHNX(current_counts, c_idx, h_idx, n_idx, halogen_mask)
            if (dbe >= unsat_min) and (dbe <= unsat_max):
                results_ppm.append(ppm_err)
                results_dbe.append(dbe)
                results_counts.append(current_counts.copy())  # wichtig: Kopie!
        return

    mass_i = element_masses[element_index]
    max_k = element_max_counts[element_index]
    base_mass = current_mass

    for k in range(max_k + 1):
        new_mass = base_mass + k * mass_i

        if new_mass > upper_bound:
            break  # höhere k werden nur schwerer

        if new_mass + suffix_max_mass[element_index + 1] < lower_bound:
            # k ist noch zu klein; versuche größeres k
            current_counts[element_index] = k
            continue

        current_counts[element_index] = k
        recursive_search_numba(current_counts, new_mass, element_index + 1,
                               element_masses, element_max_counts, suffix_max_mass,
                               lower_bound, upper_bound, target_mz,
                               unsat_min, unsat_max,
                               c_idx, h_idx, n_idx, halogen_mask,
                               results_ppm, results_dbe, results_counts)

    current_counts[element_index] = 0  # Backtrack


# ------------------------------------------------------------
# 3) Public Wrapper-funktion (Python)
#    - builds indices/masks for C,H,N,Halogens
# ------------------------------------------------------------
def find_formulas_with_info_numba(
    target_mz: float,
    mass_tolerance_ppm: float,
    unsat_min: float,
    unsat_max: float,
    elements: List[str],
    max_counts: Dict[str, int],
) -> List[Dict]:
    """
    DBE = (2C + 2 + N - H - X)/2
    X ist die Summe monovalenter Halogene (F, Cl, Br, I).
    """

    # --- Prepare inputs ---
    element_data = sorted(
        [(e, NUMBA_ELEMENT_MASSES[e]) for e in elements],
        key=lambda x: x[1]
    )
    element_symbols = [e[0] for e in element_data]
    element_masses_arr = np.array([e[1] for e in element_data], dtype=np.float64)
    element_max_counts_arr = np.array([max_counts[e] for e in element_symbols], dtype=np.int64)

    # Indices for C, H, N (=-1 if not present)
    def idx_of(sym: str) -> int:
        try:
            return element_symbols.index(sym)
        except ValueError:
            return -1

    c_idx = idx_of("C")
    h_idx = idx_of("H")
    n_idx = idx_of("N")

    # Halogen-Mask (F, Cl, Br, I)
    halogens = set(["F", "Cl", "Br", "I"])
    halogen_mask = np.array([sym in halogens for sym in element_symbols], dtype=np.bool_)

    # --- Tolerance ---
    tolerance_da = target_mz * (mass_tolerance_ppm / 1e6)
    lower_bound = target_mz - tolerance_da
    upper_bound = target_mz + tolerance_da

    # --- Suffix-Max-Mass ---
    n = len(element_symbols)
    suffix_max_mass = np.zeros(n + 1, dtype=np.float64)
    running = 0.0
    for i in range(n - 1, -1, -1):
        running += element_max_counts_arr[i] * element_masses_arr[i]
        suffix_max_mass[i] = running

    # --- Initial-Counts ---
    current_counts0 = np.zeros(n, dtype=np.int64)

    # --- resultcontainer ---
    results_ppm = NumbaList.empty_list(types.float64)
    results_dbe = NumbaList.empty_list(types.float64)
    results_counts = NumbaList.empty_list(types.int64[:])  # 1D arrays of int64

    # --- JIT-Call ---
    recursive_search_numba(
        current_counts=current_counts0,
        current_mass=0.0,
        element_index=0,
        element_masses=element_masses_arr,
        element_max_counts=element_max_counts_arr,
        suffix_max_mass=suffix_max_mass,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        target_mz=target_mz,
        unsat_min=unsat_min,
        unsat_max=unsat_max,
        c_idx=c_idx, h_idx=h_idx, n_idx=n_idx,
        halogen_mask=halogen_mask,
        results_ppm=results_ppm,
        results_dbe=results_dbe,
        results_counts=results_counts
    )

    # --- Post-Processing: Formulas & Sorting ---
    final_results = []
    for i in range(len(results_ppm)):
        counts_arr = results_counts[i]
        parts = []
        for j, el in enumerate(element_symbols):
            c = int(counts_arr[j])
            if c > 0:
                parts.append(f"{el}{c if c > 1 else ''}")
        formula_str = "".join(parts)

        final_results.append({
            "formula": formula_str,
            "ppm_error": float(results_ppm[i]),
            "unsaturation_dbe": float(results_dbe[i]),
        })

    final_results.sort(key=lambda x: abs(x["ppm_error"]))
    return final_results
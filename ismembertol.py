"""
Function to find identicals within a tolerance
input: two arrays (vec_1 and vec_2) and absolute tolerance (tol)
returns booleans with indices of identical values 
e.g.: vec_1[bool_1] = values in vec_1 that are identical with vec_2 within tol
and boolean matrix diff_bool
"""
import numpy as np

def ismembertol(
        vec1, 
        vec2, 
        tol
        ):

    vec1 = np.array(vec1)
    vec2 = np.array(vec2)
    # Compute the absolute difference between vec1 and vec2
    diff = np.abs(np.subtract(vec1, vec2[:, np.newaxis]))

    # Create a boolean array indicating which elements are within the tolerance
    diff_bool = np.less_equal(diff, tol)

    # Check whether any elements in each row of bool1 are True
    bool1 = np.any(diff_bool, axis = 0)
    bool2 = np.any(diff_bool, axis = 1)

    return bool1, bool2, diff_bool
import re
import pyopenms as oms
from pyopenms.Constants import ELECTRON_MASS_U

# NOTE: Current error: [M+2H] does not work! FIX!
# TODO: Use and test second function!

def get_adduct_data(adduct):
    """
    Parse common adduct string (e.g., '[M+H]+', or '[M-Cl]-') and calculate the adduct mass

    Parameters
    ----------
    adduct : str
        The adduct string in the format of '[M+H]+' 

    Returns
    -------
    adduct_mass : float
        The adduct mass
    element : str
        The element symbol of the adduct (e.g. 'H', 'Na', 'K')
    polarity : str
        The polarity of the adduct (e.g. '+', '-')
    element_bool : bool
        Whether the element needs to be added or subtracted
    """
    electron_mass = ELECTRON_MASS_U

    polarity = adduct[-1]
    
    inner_data = adduct[adduct.index('[') + 1:adduct.index(']')]

    def check_plus_minus(s):
        if '+' in s:
            return '+'
        elif '-' in s:
            return '-'
        else:
            return None
    element_bool = check_plus_minus(inner_data)

    if not element_bool:
        element = None
    else:
        element = re.split(r'(\+|-)', inner_data)[-1]

        try:
            edb = oms.ElementDB()
            element_mass = edb.getElement(element).getMonoWeight()
        except AttributeError:
            element_mass = oms.EmpiricalFormula(element).getMonoWeight()

    if (polarity == '+') and (element_bool == '+'):
        adduct_mass = element_mass - electron_mass

    elif (polarity == '-') and (element_bool == '+'):
        adduct_mass = element_mass + electron_mass

    elif (polarity == '+') and (element_bool == '-'):
        adduct_mass = - element_mass - electron_mass

    elif (polarity == '-') and (element_bool == '-'):
        adduct_mass = - element_mass + electron_mass

    elif (polarity == '+') and (not element_bool):
        adduct_mass = - electron_mass

    elif (polarity == '-') and (not element_bool):
        adduct_mass = electron_mass

    return adduct_mass, element, polarity, element_bool




def get_adduct_data_new(adduct):
    """
    OMS-based version of get_adduct_data supporting multiplicities and multiple elements.

    Parameters
    ----------
    adduct : str
        Adduct string like '[M+H]+', '[M+2H]+', '[M+Na+H]+', or '[M-Cl]-'

    Returns
    -------
    adduct_mass : float
        Total adduct mass contribution
    element : str
        First element in adduct (for backward compatibility)
    polarity : str
        '+' or '-'
    element_bool : str
        '+' if added, '-' if removed, None if none
    """
    polarity = adduct[-1]
    inner_data = adduct[adduct.index('[')+1 : adduct.index(']')]
    inner_data = inner_data.replace('M', '')

    # Regex to extract sign, multiplicity, element
    pattern = r'([+-])(\d*)([A-Z][a-z]?)'
    matches = re.findall(pattern, inner_data)

    if not matches:
        # No elements, just [M]+ or [M]-
        if polarity == '+':
            return -ELECTRON_MASS_U, None, polarity, None
        else:
            return ELECTRON_MASS_U, None, polarity, None

    total_mass = 0.0
    first_element = None
    first_element_bool = None
    edb = oms.ElementDB()

    for i, (sign, mult, elem) in enumerate(matches):
        mult = int(mult) if mult else 1

        # Store first element for backward compatibility
        if i == 0:
            first_element = elem
            first_element_bool = sign

        try:
            element_mass = edb.getElement(elem).getMonoWeight()
        except AttributeError:
            # fallback to empirical formula
            element_mass = oms.EmpiricalFormula(elem).getMonoWeight()

        if sign == '+':
            total_mass += element_mass * mult
        else:
            total_mass -= element_mass * mult

    # Adjust for electron mass
    if polarity == '+':
        total_mass -= ELECTRON_MASS_U
    else:
        total_mass += ELECTRON_MASS_U

    return total_mass, first_element, polarity, first_element_bool

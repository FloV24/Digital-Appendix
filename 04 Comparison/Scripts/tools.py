import numpy as np
import re

# Define what we want to remove from SMILES strings
bond_symbols = ["=", "-", "~", "#", ":", "$", "/", "\\", "@"]

# Define our target elements
metals = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Hf', 'Ta', 'W', 'Re']
non_metals = ['C', 'N', 'O', 'P', 'S', 'Cl', 'H']
wanted_elements = metals + non_metals

# Complete list of ALL two-letter elements from periodic table
two_letter_elements = ['He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar', 'Ca', 'Sc', 'Ti', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

def convert_lowercase_to_uppercase(smiles):
    """Convert standalone c,n,s,o to C,N,S,O (avoid Cl->CL, Cr->CR, etc.)"""
    result = ""
    i = 0
    
    while i < len(smiles):
        char = smiles[i]
        
        if char in ['c', 'n', 's', 'o']:
            # Check if it's part of a multi-letter element
            is_part_of_larger_element = False
            
            # Look backward: is previous char + this char = any two-letter element?
            if i > 0 and smiles[i-1:i+1] in two_letter_elements:
                is_part_of_larger_element = True
            
            # If it's standalone, convert to uppercase
            if not is_part_of_larger_element:
                result += char.upper()
            else:
                result += char
        else:
            result += char
        i += 1
    
    return result

def extract_element_from_brackets(bracketed_text):
    """Extract element from [Cr+3] -> Cr, [NH3+] -> N, etc."""
    content = bracketed_text[1:-1]  # Remove [ and ]
    
    # Sort elements by length (longest first) to avoid Cl being matched as C
    sorted_elements = sorted(wanted_elements, key=len, reverse=True)
    
    for element in sorted_elements:
        if content.startswith(element):
            return element
    
    return bracketed_text  # If we can't identify it, keep original

def check_all_elements_wanted(key):
    """Check if ALL elements in the key are in wanted_elements"""
    # Sort wanted elements by length (longest first) to avoid partial matches
    sorted_elements = sorted(wanted_elements, key=len, reverse=True)
    
    i = 0
    while i < len(key):
        found_element = False
        
        # Try to match each wanted element at current position
        for element in sorted_elements:
            if key[i:i+len(element)] == element:
                # Check if this is a complete element (not part of a larger one)
                # Either at end of string or next char is not lowercase
                if (i + len(element) >= len(key) or 
                    not key[i + len(element)].islower() or
                    key[i + len(element)] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']):
                    i += len(element)
                    found_element = True
                    break
        
        if not found_element:
            # Check if it's a digit (allowed)
            if key[i].isdigit():
                i += 1
            else:
                # Found an unwanted character/element
                return False
    
    return True

def simplify_smiles(dictionary):
    """Main function: clean SMILES and keep only wanted elements"""
    cleaned_dict = {}
    
    for original_key in dictionary:
        key = original_key
        # Step 1: Convert lowercase letters
        key = convert_lowercase_to_uppercase(key)
        
        # Step 2: Remove bond symbols
        for bond in bond_symbols:
            key = key.replace(bond, "")
        
        # Step 3: Process bracketed elements [Cr+3] -> Cr
        brackets = re.findall(r'\[[^\]]+\]', key)
        for bracket in brackets:
            element = extract_element_from_brackets(bracket)
            key = key.replace(bracket, element)
        
        # Step 4: Skip if still has brackets (couldn't clean fully)
        if '[' in key or ']' in key:
            continue
        
        # Step 5: Check if ALL elements in key are wanted elements
        if not check_all_elements_wanted(key):
            continue
        
        # Step 6: Check if key contains at least one element (not just digits)
        has_any_element = any(element in key for element in wanted_elements)
        if not has_any_element:
            continue
        
        # Step 7: Add to cleaned dictionary
        if key not in cleaned_dict:
            cleaned_dict[key] = np.array([])
        cleaned_dict[key] = np.append(cleaned_dict[key], dictionary[original_key])
    
    return cleaned_dict

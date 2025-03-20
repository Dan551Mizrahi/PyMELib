from typing import Dict, List
from itertools import product


def generate_dictionaries_from_sets(data: Dict[str, List]) -> List:
    """
    Generate a list of dictionaries from a dictionary of lists of values.
    :param data: A dictionary where each key is a string and each value is a list of values.
    :return: A list of dictionaries where each dictionary has the same keys as the input dictionary and each value is a
    """
    # Extract keys and sets of values
    keys = list(data.keys())
    values = list(data.values())

    # Generate Cartesian product of the sets of values
    combinations = product(*values)

    # Create a list of dictionaries for each combination
    result = [dict(zip(keys, combination)) for combination in combinations]

    return result

def reduce_dict_by_function(data: Dict, property_function) -> set:
    """
    Reduce a dictionary by a function.
    :param data: A dictionary.
    :param property_function: A property function that takes a value and returns a boolean.
    :return: A set of keys for which the function returned True regarding their values.
    """
    return_set = set()
    for key, value in data.items():
        if property_function.__get__(value):
            return_set.add(key)
    return return_set
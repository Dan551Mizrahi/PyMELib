from abc import ABC, abstractmethod
from pprint import pformat

class Factor(ABC):
    """
    Abstract class for factors which are used in the preprocessing phase and then in the enumeration phase.
    """
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def get_value(self, key) -> bool:
        pass

    @abstractmethod
    def set_value(self, key, value):
        pass


class MemoTable(Factor):
    """Basic dictionary-based implementation of a factor."""
    def __init__(self):

        # Memoization table
        # key: tuple of assignments in order of node[order_for_factor]
        # value: bool
        self._hash_table = dict()

    def get_only_true_keys(self) -> list:
        return [key for key in self._hash_table.keys() if self._hash_table[key] == 1]

    def get_value(self, key) -> bool:
        return self._hash_table[key]

    def set_value(self, key, value):
        self._hash_table[key] = value

    def get_all_keys(self):
        return self._hash_table.keys()

    def __str__(self):
        return_str = pformat([{k: v.name for k, v in key.items()} for key in self.get_only_true_keys()], width=100)
        return return_str.replace('\n', '<br>')